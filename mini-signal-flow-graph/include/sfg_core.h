/**
 * sfg_core.h — Signal Flow Graph Core Data Structures
 *
 * Signal Flow Graphs (SFG), introduced by Claude Shannon (1942) and
 * formalized by Samuel J. Mason (1953), are directed graphs that represent
 * systems of linear algebraic equations. An SFG consists of:
 *   - Nodes: representing system variables (signals)
 *   - Branches: directed edges with associated transmittance (gain)
 *
 * This module defines the fundamental data types and basic graph
 * manipulation operations.
 *
 * References:
 *  - Mason, S.J. "Feedback Theory — Some Properties of Signal Flow Graphs"
 *    Proc. IRE, Vol. 41, No. 9, pp. 1144-1156, 1953.
 *  - Shannon, C.E. "The Theory and Design of Linear Differential Equation
 *    Machines" OSRD Report 411, 1942.
 *
 * Author: Knowledge-first construction
 * Date:   2026-06-20
 */

#ifndef SFG_CORE_H
#define SFG_CORE_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: Core Type Definitions
 * ================================================================ */

/**
 * sfg_node_type_t — Classification of SFG nodes
 *
 * L1 Definition: SFG nodes are classified by their role in the graph.
 * - Source (input) nodes have only outgoing branches
 * - Sink (output) nodes have only incoming branches
 * - Internal nodes have both incoming and outgoing branches
 * - Summing junction nodes combine multiple incoming signals
 * - Pickoff nodes branch a signal to multiple destinations
 */
typedef enum {
    SFG_NODE_SOURCE          = 0,  /* Input node — only outgoing branches */
    SFG_NODE_SINK            = 1,  /* Output node — only incoming branches */
    SFG_NODE_INTERNAL        = 2,  /* Internal variable node */
    SFG_NODE_SUMMING         = 3,  /* Summing junction (+/- weights) */
    SFG_NODE_PICKOFF         = 4,  /* Signal pickoff / branch point */
    SFG_NODE_STATE           = 5,  /* State variable (integrator output) */
    SFG_NODE_DUMMY           = 6   /* Auxiliary / dummy node */
} sfg_node_type_t;

/**
 * sfg_gain_t — Transmittance (branch gain) type
 *
 * Branches carry a transmittance (gain) that can be real or complex.
 * Complex gains represent frequency-dependent or phase-shifting elements.
 *
 * L1 Definition: A branch from node i to node j with transmittance g_ij
 * means: x_j = g_ij * x_i  (signal at j is g_ij times signal at i).
 *
 * L3 Math Structure: The gain field is the complex numbers C.
 */
typedef double complex sfg_gain_t;

/**
 * sfg_real_t — Real-valued gain for purely resistive/resistive-equivalent
 * systems (DC gain, static sensitivity, etc.)
 */
typedef double sfg_real_t;

/**
 * sfg_node_id_t — Unique node identifier (0-based index)
 */
typedef int32_t sfg_node_id_t;

#define SFG_INVALID_NODE_ID  ((sfg_node_id_t)(-1))
#define SFG_MAX_NODES        256
#define SFG_MAX_BRANCHES     1024

/**
 * sfg_node_t — A single node in the signal flow graph
 *
 * Each node represents a system variable. Nodes carry:
 * - A unique identifier
 * - A type classification
 * - An optional label (variable name)
 * - A flag indicating whether an input branch with gain 1 connects
 *   to this node (for mixed nodes in Mason's formula)
 */
typedef struct {
    sfg_node_id_t   id;              /* Unique node index */
    sfg_node_type_t type;            /* Node classification */
    char            label[32];       /* Variable name label */
    int             is_input_tapped; /* Has unity-gain input? (for mixed nodes) */
} sfg_node_t;

/**
 * sfg_branch_t — A directed branch (edge) in the signal flow graph
 *
 * L1 Definition: A branch from node `src` to node `dst` with
 * transmittance (gain) `gain`. The signal equation is:
 *   x_dst += gain * x_src
 *
 * Multiple branches may enter the same destination node (summation).
 */
typedef struct {
    sfg_node_id_t src;              /* Source node index */
    sfg_node_id_t dst;              /* Destination node index */
    sfg_gain_t    gain;             /* Transmittance (branch gain) */
    char          label[32];        /* Optional branch label */
} sfg_branch_t;

/**
 * sfg_path_t — A directed path (sequence of branches)
 *
 * L1 Definition: A forward path is a directed path from a source node
 * to a sink node that passes through each node at most once.
 * A loop is a closed path that starts and ends at the same node,
 * passing through other nodes at most once.
 */
typedef struct {
    sfg_node_id_t nodes[SFG_MAX_NODES]; /* Node sequence */
    int           length;               /* Number of nodes in the path */
    sfg_gain_t    gain;                 /* Product of branch gains along path */
    int           is_loop;              /* 1 if loop, 0 if forward path */
} sfg_path_t;

/**
 * sfg_graph_t — A complete signal flow graph
 *
 * The fundamental data structure. An SFG consists of a set of nodes
 * and a set of directed branches connecting them.
 */
typedef struct {
    sfg_node_t   nodes[SFG_MAX_NODES];
    int          num_nodes;

    sfg_branch_t branches[SFG_MAX_BRANCHES];
    int          num_branches;

    sfg_node_id_t input_node;       /* Primary source (input) node */
    sfg_node_id_t output_node;      /* Primary sink (output) node */

    char name[64];                  /* Graph name/label */
} sfg_graph_t;

/* ================================================================
 * L2: Core Graph Operations
 * ================================================================ */

/**
 * sfg_graph_init — Initialize an empty signal flow graph
 *
 * Concept: Graph construction begins with an empty container.
 *
 * Complexity: O(1)
 * Returns: pointer to initialized graph (same as `g`)
 */
sfg_graph_t* sfg_graph_init(sfg_graph_t* g);

/**
 * sfg_graph_clear — Reset graph to empty state (free all nodes/branches)
 *
 * Complexity: O(N + B) where N = nodes, B = branches
 */
void sfg_graph_clear(sfg_graph_t* g);

/**
 * sfg_add_node — Add a node to the graph
 *
 * L1: Creates a new system variable node in the SFG.
 *
 * Complexity: O(1)
 * Returns: node ID on success, SFG_INVALID_NODE_ID if graph is full
 */
sfg_node_id_t sfg_add_node(sfg_graph_t* g, sfg_node_type_t type,
                           const char* label);

/**
 * sfg_add_source — Convenience: add a source (input) node
 */
sfg_node_id_t sfg_add_source(sfg_graph_t* g, const char* label);

/**
 * sfg_add_sink — Convenience: add a sink (output) node
 */
sfg_node_id_t sfg_add_sink(sfg_graph_t* g, const char* label);

/**
 * sfg_add_state — Convenience: add a state-variable node (1/s integrator)
 */
sfg_node_id_t sfg_add_state(sfg_graph_t* g, const char* label);

/**
 * sfg_add_summing — Convenience: add a summing junction node
 */
sfg_node_id_t sfg_add_summing(sfg_graph_t* g, const char* label);

/**
 * sfg_remove_node — Remove a node and all incident branches
 *
 * Complexity: O(B) where B = number of branches
 * Returns: 0 on success, -1 if node not found
 */
int sfg_remove_node(sfg_graph_t* g, sfg_node_id_t node_id);

/**
 * sfg_add_branch — Add a directed branch between two nodes
 *
 * L1: Creates a signal transmission path from src to dst with given gain.
 * Multiple branches between the same nodes are allowed (parallel).
 *
 * Complexity: O(1)
 * Returns: branch index on success, -1 if invalid or graph is full
 */
int sfg_add_branch(sfg_graph_t* g, sfg_node_id_t src, sfg_node_id_t dst,
                   sfg_gain_t gain, const char* label);

/**
 * sfg_add_branch_real — Add a branch with real-valued gain
 */
int sfg_add_branch_real(sfg_graph_t* g, sfg_node_id_t src, sfg_node_id_t dst,
                        sfg_real_t gain, const char* label);

/**
 * sfg_remove_branch — Remove a specific branch
 *
 * Complexity: O(B) due to array compaction
 * Returns: 0 on success, -1 if not found
 */
int sfg_remove_branch(sfg_graph_t* g, sfg_node_id_t src, sfg_node_id_t dst);

/**
 * sfg_find_branch — Find a branch between src and dst
 *
 * Returns: pointer to branch, or NULL if not found
 * Complexity: O(B) linear scan
 */
const sfg_branch_t* sfg_find_branch(const sfg_graph_t* g,
                                     sfg_node_id_t src, sfg_node_id_t dst);

/**
 * sfg_set_input — Designate the primary source node
 */
void sfg_set_input(sfg_graph_t* g, sfg_node_id_t node_id);

/**
 * sfg_set_output — Designate the primary sink node
 */
void sfg_set_output(sfg_graph_t* g, sfg_node_id_t node_id);

/* ================================================================
 * L2: Graph Query Operations
 * ================================================================ */

/**
 * sfg_node_count — Return number of nodes in the graph
 */
int sfg_node_count(const sfg_graph_t* g);

/**
 * sfg_branch_count — Return number of branches in the graph
 */
int sfg_branch_count(const sfg_graph_t* g);

/**
 * sfg_get_node — Get node by ID
 *
 * Returns: pointer to node, or NULL if ID invalid
 */
const sfg_node_t* sfg_get_node(const sfg_graph_t* g, sfg_node_id_t id);

/**
 * sfg_get_branch — Get branch by index
 */
const sfg_branch_t* sfg_get_branch(const sfg_graph_t* g, int idx);

/**
 * sfg_node_in_degree — Count incoming branches to a node
 *
 * L2 Concept: In-degree reflects how many signals contribute to a variable.
 * For a sink node, in-degree = number of contributing forward paths.
 *
 * Complexity: O(B)
 */
int sfg_node_in_degree(const sfg_graph_t* g, sfg_node_id_t node_id);

/**
 * sfg_node_out_degree — Count outgoing branches from a node
 *
 * Complexity: O(B)
 */
int sfg_node_out_degree(const sfg_graph_t* g, sfg_node_id_t node_id);

/**
 * sfg_is_source — Check if node has only outgoing branches
 */
int sfg_is_source(const sfg_graph_t* g, sfg_node_id_t node_id);

/**
 * sfg_is_sink — Check if node has only incoming branches
 */
int sfg_is_sink(const sfg_graph_t* g, sfg_node_id_t node_id);

/**
 * sfg_node_exists — Check if a node ID is valid in the graph
 */
int sfg_node_exists(const sfg_graph_t* g, sfg_node_id_t id);

/**
 * sfg_branch_exists — Check if a branch exists between two nodes
 */
int sfg_branch_exists(const sfg_graph_t* g,
                      sfg_node_id_t src, sfg_node_id_t dst);

/* ================================================================
 * L2: Graph Manipulation
 * ================================================================ */

/**
 * sfg_graph_copy — Deep copy of a signal flow graph
 *
 * Complexity: O(N + B)
 * Returns: pointer to destination graph
 */
sfg_graph_t* sfg_graph_copy(sfg_graph_t* dst, const sfg_graph_t* src);

/**
 * sfg_graph_transpose — Create transpose (reverse all branch directions)
 *
 * L2 Concept: The transpose of an SFG corresponds to the adjoint system.
 * If original SFG represents equations A x = b, the transpose represents
 * A^T y = c. Used in sensitivity analysis and adjoint methods.
 *
 * Complexity: O(N + B)
 */
sfg_graph_t* sfg_graph_transpose(sfg_graph_t* dst, const sfg_graph_t* src);

/**
 * sfg_graph_reverse — In-place reversal of all branch directions
 *
 * Complexity: O(B)
 */
void sfg_graph_reverse(sfg_graph_t* g);

/**
 * sfg_merge_nodes — Merge two nodes into one (combine their incident branches)
 *
 * L5 Concept: Node merging is used in SFG reduction. Branches incident to
 * both nodes are redirected to the merged node.
 *
 * Returns: 0 on success, -1 on error
 */
int sfg_merge_nodes(sfg_graph_t* g, sfg_node_id_t keep,
                    sfg_node_id_t absorb);

/* ================================================================
 * L3: Adjacency / Connection Matrix
 * ================================================================ */

/**
 * sfg_connection_matrix — Build the connection (adjacency) matrix
 *
 * Concept: C[i][j] = gain of branch from node i to node j.
 * If no direct branch exists, C[i][j] = 0.
 *
 * The connection matrix maps directly to the system of equations:
 *   x_j = sum_i C[i][j] * x_i  for all j
 *
 * @param g    The signal flow graph
 * @param C    Pre-allocated N×N matrix (row-major, size N*N)
 * @param N    Number of nodes (must equal sfg_node_count(g))
 * @return     0 on success, -1 on dimension mismatch
 *
 * Complexity: O(N² + B)
 */
int sfg_connection_matrix(const sfg_graph_t* g, sfg_gain_t* C, int N);

/**
 * sfg_gain_matrix — Build the coefficient matrix for x = G x + B u form
 *
 * Concept: For state-space-like representation, G[i][j] is the gain
 * from state j to state i.
 *
 * Complexity: O(N² + B)
 */
int sfg_gain_matrix(const sfg_graph_t* g, sfg_gain_t* G, int N);

/**
 * sfg_node_index_map — Build array mapping node_id → sequential index
 *
 * Complexity: O(N)
 */
void sfg_node_index_map(const sfg_graph_t* g, int* map);

/* ================================================================
 * L3: Graph Display / Debugging
 * ================================================================ */

/**
 * sfg_graph_print — Print graph structure (node list + branch list)
 *
 * Complexity: O(N + B)
 */
void sfg_graph_print(const sfg_graph_t* g);

/**
 * sfg_graph_to_dot — Export graph in Graphviz DOT format
 *
 * Complexity: O(N + B)
 * @param filename  Output file path (NULL for stdout)
 */
void sfg_graph_to_dot(const sfg_graph_t* g, const char* filename);

/**
 * sfg_graph_validate — Check graph integrity
 *
 * Verifies: all branch endpoints reference existing nodes,
 * no duplicate node IDs, node count consistency.
 *
 * Returns: 0 if valid, negative error code otherwise
 */
int sfg_graph_validate(const sfg_graph_t* g);

#ifdef __cplusplus
}
#endif

#endif /* SFG_CORE_H */
