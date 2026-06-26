/**
 * @file blockdiagram.h
 * @brief Block Diagram — core data structures for control system modeling
 *
 * Knowledge Coverage (L1-L9):
 *   L1: Node types (input, output, block, summer, takeoff, gain),
 *       TransferFunction struct (ratio of polynomials), edge/path types
 *   L2: Construction/lifecycle (bd_create/destroy), property assignment
 *   L3: Graph representation (nodes as vertices, edges as directed connections)
 *   L4: Topological queries for Mason's formula — forward paths via DFS,
 *       loop detection with backtracking, non-touching analysis
 *   L5: Clone (deep copy), validation (structural consistency),
 *       printing (debug/education), degree queries
 *
 * Block diagrams are the foundational graphical language of control
 * engineering. They model signal flow through interconnected dynamic
 * elements: transfer function blocks process signals, summing junctions
 * combine multiple signals, and take-off points branch one signal to
 * multiple destinations.
 *
 * Ref: K. Ogata, "Modern Control Engineering", 5th ed., Ch. 3, Prentice Hall
 * Ref: R.C. Dorf & R.H. Bishop, "Modern Control Systems", 13th ed., Ch. 2, Pearson
 * Ref: G.F. Franklin et al., "Feedback Control of Dynamic Systems", 8th ed., Ch. 3
 * MIT 6.302 Feedback Systems / Stanford ENGR105 Feedback Control
 * Berkeley ME132 Dynamic Systems / ETH 151-0591 Control Systems I
 */
#ifndef BLOCKDIAGRAM_H
#define BLOCKDIAGRAM_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================
 * L1: Core Type Definitions — the fundamental building blocks
 * ======================================================================== */

/** The six types of nodes that can appear in a block diagram.
 *  Each type has specific degree constraints:
 *  INPUT:  in=0, out>=1   — only source of external signals
 *  OUTPUT: in>=1, out=0   — only sink for measured signals
 *  BLOCK:  in=1, out=1    — transfer function element G(s)
 *  SUMMER: in>=1, out=1   — addition/subtraction of signals
 *  TAKEOFF:in=1, out>=1   — signal branching point
 *  GAIN:   in=1, out=1    — constant multiplier K (special case of BLOCK) */
typedef enum { BD_NODE_INPUT=0, BD_NODE_OUTPUT=1, BD_NODE_BLOCK=2,
               BD_NODE_SUMMER=3, BD_NODE_TAKEOFF=4, BD_NODE_GAIN=5 } bd_node_type_t;

/**
 * Transfer Function: G(s) = N(s)/D(s) where
 *   N(s) = num[0] + num[1]*s + ... + num[m]*s^m  (ascending powers)
 *   D(s) = den[0] + den[1]*s + ... + den[n]*s^n  (ascending powers)
 *
 * Proper TF condition: den_order >= num_order (n >= m).
 * The delay field models transport lag: G(s) * exp(-delay*s).
 * is_discrete flag distinguishes s-domain (continuous) from z-domain.
 */
typedef struct { double *num,*den; int num_order,den_order; double delay;
                 int is_discrete; double sample_time; } TransferFunction;

/**
 * A node in the block diagram graph. Fields:
 *   id: unique identifier within the diagram
 *   type: determines the node's behavior and degree constraints
 *   tf: transfer function for BLOCK nodes (NULL for non-block types)
 *   gain: multiplicative factor for GAIN nodes, default 1.0
 *   input_signs: +/-1 array for SUMMER node input ports
 *   input_count/output_count: actual number of connected edges
 *   n_inputs/n_outputs: allocated capacity for edge connections
 *   initial_cond: initial value for state (time-domain simulation)
 */
typedef struct bd_node_t { int id; bd_node_type_t type; char label[64];
    TransferFunction *tf; double gain; int n_inputs,n_outputs;
    int *input_signs; int input_count,output_count; double initial_cond; } bd_node_t;

/** A directed edge in the block diagram carrying signal from src to dst.
 *  weight: multiplicative factor on the edge (default 1.0).
 *  dst_port: which input port at the destination node this edge connects to. */
typedef struct { int id,src_node,dst_node,dst_port; double weight; } bd_edge_t;

/** A path through the block diagram — sequence of nodes and edges.
 *  Forward paths go from input to output (is_loop=0).
 *  Loops are closed paths where first==last node (is_loop=1).
 *  gain: product of all edge weights along the path.
 *  Used as the primary data structure for Mason's gain formula. */
typedef struct { int *node_ids,*edge_ids,length,is_loop; double gain; } bd_path_t;

/** Top-level block diagram: a directed graph of signal processing elements.
 *  nodes/edges: dynamically resized arrays of pointers.
 *  input_node_id/output_node_id: designated I/O terminals.
 *  next_node_id/next_edge_id: auto-increment counters for unique IDs. */
typedef struct { int n_nodes,max_nodes; bd_node_t **nodes; int n_edges,max_edges;
    bd_edge_t **edges; int input_node_id,output_node_id,next_node_id,next_edge_id;
    char name[128]; } BlockDiagram;

/* ========================================================================
 * L2: Construction and Lifecycle
 * ======================================================================== */

/** Create an empty block diagram with pre-allocated capacity.
 *  initial_nodes/edges set the starting array sizes (auto-grow on demand).
 *  Returns NULL on memory allocation failure. */
BlockDiagram* bd_create(const char *name, int ini_n, int ini_e);

/** Free all memory associated with a block diagram.
 *  Recursively frees nodes (and their transfer functions), edges, and arrays. */
void bd_destroy(BlockDiagram *bd);

/** Add a node of the given type. Returns the new node ID, or -1 on failure.
 *  Automatically sets node as input/output if it's the first of that type. */
int bd_add_node(BlockDiagram *bd, bd_node_type_t t, const char *label);

/** Remove a node and all its incident edges from the diagram.
 *  Complexity: O(N+E) due to edge scan and array compaction. */
int bd_remove_node(BlockDiagram *bd, int nid);

/** Add a directed edge between two existing nodes with given weight.
 *  Returns edge ID, or -1 if either node does not exist. */
int bd_add_edge(BlockDiagram *bd, int src, int dst, double w);

/** Remove an edge from the diagram. Updates node degree counts. */
int bd_remove_edge(BlockDiagram *bd, int eid);

/* ========================================================================
 * L3: Node Property Operations
 * ======================================================================== */

/** Assign a transfer function to a BLOCK node (deep-copied internally). */
int bd_set_transfer_function(BlockDiagram *bd, int nid, const TransferFunction *tf);

/** Set the constant gain of a GAIN or BLOCK node. */
int bd_set_gain(BlockDiagram *bd, int nid, double g);

/** Configure the +/- signs for a SUMMER node's inputs.
 *  signs[i] >= 0 means addition (+), < 0 means subtraction (-). */
int bd_set_summer_signs(BlockDiagram *bd, int nid, const double *s, int ns);

/** Set the initial condition value for time-domain simulations. */
int bd_set_initial_condition(BlockDiagram *bd, int nid, double ic);

/* ========================================================================
 * L4: Topological Queries — Foundation for Mason's Formula
 * ======================================================================== */

/** Look up a node by its unique ID. Returns NULL if not found. O(N). */
bd_node_t* bd_get_node(const BlockDiagram *bd, int nid);

/** Find all forward paths from input to output using DFS with backtracking.
 *  Each path is a simple walk (no repeated nodes). Returns count found. */
int bd_find_forward_paths(const BlockDiagram *bd, int maxp, bd_path_t *out);

/** Find all individual loops (closed paths where start==end, no other repeats).
 *  Uses DFS from each non-I/O node. Returns count found. */
int bd_find_loops(const BlockDiagram *bd, int maxl, bd_path_t *out);

/** Check if two loops share no nodes (are non-touching).
 *  Essential for computing higher-order terms in Mason's determinant:
 *  Delta = 1 - Sum(L_i) + Sum(L_i*L_j) - Sum(L_i*L_j*L_k) + ... */
int bd_loops_are_non_touching(const bd_path_t *a, const bd_path_t *b);

/** Compute the product of gains for a set of loops.
 *  For a set of k non-touching loops, the contribution is (-1)^k * product. */
double bd_loop_gain_product(const bd_path_t *loops, int n);

/* ========================================================================
 * L5: Utility Operations
 * ======================================================================== */

/** Compute the total gain along a path by multiplying all edge weights. */
double bd_path_gain(const BlockDiagram *bd, const bd_path_t *p);

/** Free a dynamically allocated path structure (caller must not reuse). */
void bd_path_free(bd_path_t *p);

/** Create a deep copy of the entire block diagram including all TFs. */
BlockDiagram* bd_clone(const BlockDiagram *bd);

/** Print a human-readable representation of the diagram.
 *  detailed=0: summary only. detailed=1: full node/edge listing. */
void bd_print(const BlockDiagram *bd, int det);

/** Count the number of nodes of each type. counts[BD_NODE_INPUT..GAIN]. */
void bd_count_node_types(const BlockDiagram *bd, int c[6]);

/** Validate structural consistency: no dangling edges, I/O nodes have
 *  correct edge patterns. Returns 0 if valid, negative error code otherwise. */
int bd_validate(const BlockDiagram *bd);

/** Return the in-degree (incoming edge count) of a node. */
int bd_node_degree_in(const BlockDiagram *bd, int nid);

/** Return the out-degree (outgoing edge count) of a node. */
int bd_node_degree_out(const BlockDiagram *bd, int nid);

/** Get the string name of a node type (for display/debugging).
 *  Returns "IN","OUT","BLK","SUM","TKO","GAIN" or "???" for invalid type. */
const char* bd_node_type_name(bd_node_type_t t);

/** Find a node by its label string. Returns node ID or -1 if not found.
 *  Useful for constructing diagrams by name rather than ID. */
int bd_find_node_by_label(const BlockDiagram *bd, const char *label);

/** Check whether a path exists from src to dst (simple BFS reachability). */
int bd_path_exists(const BlockDiagram *bd, int src_id, int dst_id);

#ifdef __cplusplus
}
#endif
#endif
