/**
 * sfg_path.h — Signal Flow Graph Path Enumeration
 *
 * Path and loop enumeration is the computational core of Mason's Gain
 * Formula. This module implements depth-first search (DFS) based
 * algorithms for finding all forward paths, all loops (cycles), and
 * identifying non-touching loop sets.
 *
 * L5 Engineering Method: Johnson's algorithm (1975) for cycle enumeration
 * in directed graphs adapted for SFG loop finding. Forward path
 * enumeration uses recursive DFS with visited-node tracking.
 *
 * References:
 *  - Johnson, D.B. "Finding All the Elementary Circuits of a Directed
 *    Graph" SIAM J. Comput., Vol. 4, No. 1, pp. 77-84, 1975.
 *  - Mason, S.J. "Feedback Theory — Further Properties of Signal Flow
 *    Graphs" Proc. IRE, Vol. 44, No. 7, pp. 920-926, 1956.
 *  - Lorens, C.S. "Flowgraphs for the Modeling and Analysis of Linear
 *    Systems" McGraw-Hill, 1964.
 */

#ifndef SFG_PATH_H
#define SFG_PATH_H

#include "sfg_core.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SFG_MAX_PATHS      256
#define SFG_MAX_LOOPS      256
#define SFG_MAX_NT_GROUPS  128

/* ================================================================
 * L1: Path and Loop Container Types
 * ================================================================ */

/**
 * sfg_path_list_t — Collection of forward paths
 *
 * L1 Definition: A set of directed paths from a source to a sink,
 * where each node appears at most once in each path.
 */
typedef struct {
    sfg_path_t paths[SFG_MAX_PATHS];
    int        count;
} sfg_path_list_t;

/**
 * sfg_loop_list_t — Collection of loops (feedback cycles)
 *
 * L1 Definition: A loop is a closed directed path that starts and ends
 * at the same node, passing through other nodes at most once.
 * The loop gain is the product of all branch gains around the loop.
 */
typedef struct {
    sfg_path_t loops[SFG_MAX_LOOPS];
    int        count;
} sfg_loop_list_t;

/**
 * sfg_loop_group_t — A group of mutually non-touching loops
 *
 * L1 Definition: Two loops are non-touching if they share no nodes.
 * A non-touching loop group is a set of loops where every pair is
 * non-touching. The group gain is the product of individual loop gains.
 */
typedef struct {
    int        loop_indices[SFG_MAX_LOOPS]; /* Indices into loop list */
    int        num_loops;                   /* Number of loops in group */
    sfg_gain_t group_gain;                  /* Product of loop gains */
} sfg_loop_group_t;

/**
 * sfg_nt_groups_t — All non-touching loop groups, organized by size
 */
typedef struct {
    sfg_loop_group_t groups[SFG_MAX_NT_GROUPS];
    int              count;
    int              sizes[SFG_MAX_LOOPS + 1]; /* size → count mapping */
} sfg_nt_groups_t;

/* ================================================================
 * L5: Path Enumeration (DFS-based)
 * ================================================================ */

/**
 * sfg_find_forward_paths — Enumerate all forward paths from source to sink
 *
 * L5 Algorithm: Modified DFS that tracks visited nodes. A forward path
 * is a simple directed path from the source to the sink node.
 * Branches are traversed in order; nodes already in the current path
 * are skipped to prevent cycles within the forward path.
 *
 * L4 Connection: The set of all forward paths forms the basis for the
 * numerator of Mason's transfer function formula.
 *
 * Complexity: O(N + B + P·L) where P = number of paths, L = avg path length.
 *             Worst-case exponential in N for dense graphs.
 *
 * @param g         The signal flow graph
 * @param source    Source node ID
 * @param sink      Sink node ID
 * @param paths     Output: filled with all forward paths
 * @return          Number of forward paths found, -1 on error
 */
int sfg_find_forward_paths(const sfg_graph_t* g,
                           sfg_node_id_t source, sfg_node_id_t sink,
                           sfg_path_list_t* paths);

/**
 * sfg_find_all_loops — Enumerate all distinct loops (cycles) in the graph
 *
 * L5 Algorithm: Johnson's algorithm for elementary circuit enumeration.
 * Starting from each node, perform DFS to find all cycles that begin
 * and end at that node. Strongly connected component (SCC) based
 * pruning is used to avoid redundant exploration.
 *
 * A loop's gain is the product of all branch gains around the loop.
 *
 * Complexity: O((N + B)(C + 1)) where C = number of cycles.
 *             Worst-case exponential for complete digraphs.
 */
int sfg_find_all_loops(const sfg_graph_t* g, sfg_loop_list_t* loops);

/**
 * sfg_find_self_loops — Find all self-loops (branch from node to itself)
 *
 * L1 Definition: A self-loop is a branch whose source and destination
 * are the same node. Self-loop gain contributes directly to the
 * denominator of the transfer function.
 *
 * Complexity: O(B)
 */
int sfg_find_self_loops(const sfg_graph_t* g, sfg_loop_list_t* loops);

/* ================================================================
 * L5: Non-Touching Loop Analysis
 * ================================================================ */

/**
 * sfg_loops_touch — Determine if two loops share any nodes
 *
 * L2 Concept: Two loops are "touching" if they have at least one node
 * in common. If they share no nodes, they are "non-touching."
 * Non-touching loops contribute to higher-order terms in Mason's formula.
 *
 * Complexity: O(min(L1_len, L2_len))
 * Returns: 1 if loops touch (share nodes), 0 if non-touching
 */
int sfg_loops_touch(const sfg_path_t* loop1, const sfg_path_t* loop2);

/**
 * sfg_find_non_touching_groups — Enumerate all non-touching loop groups
 *
 * L5 Algorithm: For each loop, find all other loops that don't touch it.
 * Extend recursively to find groups of size k (k = 2, 3, ...).
 * Stop when no new non-touching groups can be formed.
 *
 * This is the most computationally intensive part of Mason's formula
 * for graphs with many feedback loops.
 *
 * Complexity: O(L² · 2^L) worst case, but typically much less due to
 *             sparsity of non-touching relationships.
 */
int sfg_find_non_touching_groups(const sfg_graph_t* g,
                                  const sfg_loop_list_t* all_loops,
                                  sfg_nt_groups_t* groups);

/**
 * sfg_path_loops_touch — Check if a forward path touches a loop
 *
 * Concept: A forward path touches a loop if they share at least one node.
 * If they do not touch, the loop is excluded from the path's cofactor Δₖ.
 *
 * Complexity: O(P_len + L_len)
 * Returns: 1 if touching, 0 if non-touching
 */
int sfg_path_loops_touch(const sfg_path_t* path, const sfg_path_t* loop);

/* ================================================================
 * L5: Path/Loop Utility Functions
 * ================================================================ */

/**
 * sfg_path_gain_compute — Compute the gain of a path (product of branch gains)
 *
 * Complexity: O(L) where L = path length
 */
sfg_gain_t sfg_path_gain_compute(const sfg_graph_t* g,
                                  const sfg_path_t* path);

/**
 * sfg_path_contains_node — Check if a node is on a given path
 *
 * Complexity: O(L)
 */
int sfg_path_contains_node(const sfg_path_t* path, sfg_node_id_t node);

/**
 * sfg_paths_touch — Check if two paths share any node
 *
 * Complexity: O(L1 + L2)
 * Returns: 1 if touching, 0 if non-touching
 */
int sfg_paths_touch(const sfg_path_t* p1, const sfg_path_t* p2);

/**
 * sfg_path_print — Print a path (node sequence + gain)
 */
void sfg_path_print(const sfg_path_t* path);

/**
 * sfg_loop_print — Print a loop with its gain
 */
void sfg_loop_print(const sfg_path_t* loop);

/**
 * sfg_path_list_sort_by_gain — Sort paths by gain magnitude (ascending)
 *
 * Complexity: O(P log P)
 */
void sfg_path_list_sort_by_gain(sfg_path_list_t* paths);

/**
 * sfg_loop_list_sort_by_gain — Sort loops by gain magnitude (ascending)
 */
void sfg_loop_list_sort_by_gain(sfg_loop_list_t* loops);

/**
 * sfg_path_list_clear — Reset a path list to empty
 */
void sfg_path_list_clear(sfg_path_list_t* pl);

/**
 * sfg_loop_list_clear — Reset a loop list to empty
 */
void sfg_loop_list_clear(sfg_loop_list_t* ll);

/**
 * sfg_nt_groups_clear — Reset non-touching groups to empty
 */
void sfg_nt_groups_clear(sfg_nt_groups_t* nt);

#ifdef __cplusplus
}
#endif

#endif /* SFG_PATH_H */
