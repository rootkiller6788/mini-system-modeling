/**
 * sfg_reduction.h — Signal Flow Graph Topological Reduction
 *
 * SFG reduction is a systematic method for simplifying signal flow
 * graphs to a single equivalent transmittance between input and output.
 * The reduction rules preserve the overall transfer function while
 * eliminating internal nodes and combining parallel paths.
 *
 * The reduction rules correspond to algebraic manipulations of the
 * underlying linear equations and provide an alternative to Mason's
 * formula for computing transfer functions.
 *
 * L5 Engineering Method: SFG reduction rules provide a visual/graphical
 * approach to system simplification, often preferred in control
 * engineering practice.
 *
 * Reduction Rules (Lorenz 1964, D'Azzo & Houpis 1995):
 *   R1. Series (cascade) branches:  a→b→c  →  a→c (gain = g1·g2)
 *   R2. Parallel branches:          a⇉b    →  a→b (gain = g1+g2)
 *   R3. Self-loop elimination:      a⟲a    →  divide incoming by (1 - g_self)
 *   R4. Node absorption:            eliminate intermediate node
 *   R5. Node splitting:             duplicate node to untangle paths
 *
 * References:
 *  - Lorens, C.S. "Flowgraphs for the Modeling and Analysis of Linear
 *    Systems" McGraw-Hill, 1964. Chapters 2-3.
 *  - D'Azzo, J.J., Houpis, C.H. "Linear Control System Analysis and
 *    Design" McGraw-Hill, 4th Ed., 1995. Section 3.6.
 */

#ifndef SFG_REDUCTION_H
#define SFG_REDUCTION_H

#include "sfg_core.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: Reduction Result Types
 * ================================================================ */

/**
 * sfg_reduction_step_t — Record of a single reduction step
 *
 * Records the type of reduction and what was done, for audit trail
 * and educational display.
 */
typedef enum {
    SFG_REDUCE_SERIES           = 0,
    SFG_REDUCE_PARALLEL         = 1,
    SFG_REDUCE_SELF_LOOP        = 2,
    SFG_REDUCE_NODE_ABSORB      = 3,
    SFG_REDUCE_NODE_SPLIT       = 4,
    SFG_REDUCE_SOURCE_SHIFT     = 5,
    SFG_REDUCE_COMPLETE         = 6
} sfg_reduce_op_t;

typedef struct {
    sfg_reduce_op_t operation;      /* Type of reduction performed */
    sfg_node_id_t   node_a;         /* First affected node */
    sfg_node_id_t   node_b;         /* Second affected node (if applicable) */
    sfg_gain_t      old_gain;       /* Gain before reduction */
    sfg_gain_t      new_gain;       /* Gain after reduction */
    char            description[128]; /* Human-readable description */
} sfg_reduction_step_t;

#define SFG_MAX_REDUCTION_STEPS 256

typedef struct {
    sfg_graph_t           reduced_graph;  /* Final reduced graph */
    sfg_gain_t            overall_gain;   /* Equivalent transmittance */
    sfg_reduction_step_t  steps[SFG_MAX_REDUCTION_STEPS];
    int                   num_steps;
    int                   success;        /* 1 if reduction completed */
} sfg_reduction_result_t;

/* ================================================================
 * L5: Reduction Operations (Rule R1-R5)
 * ================================================================ */

/**
 * sfg_reduce_series — Eliminate series (cascade) branches
 *
 * R1: If node B has exactly one incoming branch (A→B) and one outgoing
 * branch (B→C), eliminate B and create direct branch A→C with
 * gain = g_AB * g_BC.
 *
 * This corresponds to algebraic substitution: x_B = g_AB·x_A,
 * x_C = g_BC·x_B = g_BC·g_AB·x_A.
 *
 * Complexity: O(B)
 * Returns: number of series eliminations performed
 */
int sfg_reduce_series(sfg_graph_t* g);

/**
 * sfg_reduce_parallel — Combine parallel branches
 *
 * R2: If multiple branches exist between the same pair of nodes
 * (A→B with gains g1, g2, ..., gk), replace with single branch
 * A→B with gain = g1 + g2 + ... + gk.
 *
 * This corresponds to: x_B = g1·x_A + g2·x_A = (g1+g2)·x_A.
 *
 * Complexity: O(B²) for pairwise comparison
 * Returns: number of parallel combinations performed
 */
int sfg_reduce_parallel(sfg_graph_t* g);

/**
 * sfg_reduce_self_loop — Eliminate self-loops
 *
 * R3: If node A has a self-loop with gain g_self, remove the self-loop
 * and divide all incoming branch gains by (1 - g_self). Outgoing
 * branches are unaffected.
 *
 * Rationale: x_A = Σ(g_in·x_in) + g_self·x_A
 *           → x_A(1 - g_self) = Σ(g_in·x_in)
 *           → x_A = Σ(g_in/(1 - g_self))·x_in
 *
 * Complexity: O(B)
 * Returns: number of self-loops eliminated
 */
int sfg_reduce_self_loop(sfg_graph_t* g);

/**
 * sfg_reduce_node_absorb — Absorb (eliminate) an internal node
 *
 * R4: Eliminate an internal node by redirecting all paths through it.
 * For each incoming branch (X→N) with gain g_in and each outgoing
 * branch (N→Y) with gain g_out, add a new branch X→Y with gain
 * g_in·g_out. Then remove node N and all its incident branches.
 *
 * Self-loops on N are eliminated first (R3 must precede R4).
 *
 * Complexity: O(d_in · d_out) where d_in, d_out are in/out degrees
 * Returns: 0 on success, -1 if node is a source or sink
 */
int sfg_reduce_node_absorb(sfg_graph_t* g, sfg_node_id_t node_id);

/**
 * sfg_reduce_node_split — Split a node to isolate non-touching paths
 *
 * R5: Duplicate a node so that paths that were touching through this
 * node become non-touching. All incoming branches go to the original
 * node; select outgoing branches go to the duplicate.
 *
 * This is the most subtle reduction rule and is rarely needed for
 * standard control system SFGs.
 *
 * Complexity: O(B)
 * Returns: 0 on success, -1 on error
 */
int sfg_reduce_node_split(sfg_graph_t* g, sfg_node_id_t node_id);

/* ================================================================
 * L5: Automated Reduction
 * ================================================================ */

/**
 * sfg_reduce_to_single — Reduce an SFG to a single branch input→output
 *
 * Full automated reduction: repeatedly apply R1-R4 until only the
 * source and sink nodes remain with a single branch between them.
 *
 * The final branch gain equals the system transfer function T.
 * This provides an alternative method to Mason's formula for computing
 * transfer functions, and serves as a verification.
 *
 * Algorithm:
 *   while (nodes > 2 || branches > 1):
 *     try: reduce_parallel()
 *     try: reduce_series()
 *     try: reduce_self_loop()
 *     if no progress: reduce_node_absorb(most_connected_internal_node)
 *
 * Complexity: O(N²·B) worst case
 * Returns: 0 on full reduction, -1 if stuck
 */
int sfg_reduce_to_single(sfg_graph_t* g, sfg_node_id_t source,
                          sfg_node_id_t sink,
                          sfg_reduction_result_t* result);

/**
 * sfg_reduce_interactive — Perform one reduction step, return what was done
 *
 * For educational and debugging purposes. Finds the best next reduction
 * step and applies it.
 *
 * Returns: 1 if a step was performed, 0 if fully reduced, -1 if stuck
 */
int sfg_reduce_interactive(sfg_graph_t* g, sfg_reduction_step_t* step);

/**
 * sfg_reduction_result_init — Initialize reduction result structure
 */
void sfg_reduction_result_init(sfg_reduction_result_t* result);

/**
 * sfg_reduction_result_print — Print the sequence of reduction steps
 */
void sfg_reduction_result_print(const sfg_reduction_result_t* result);

/**
 * sfg_is_reducible — Check if the graph can be reduced to single branch
 *
 * A graph is reducible iff it is connected and has no unreachable
 * internal nodes.
 *
 * Returns: 1 if reducible, 0 otherwise
 */
int sfg_is_reducible(const sfg_graph_t* g, sfg_node_id_t source,
                      sfg_node_id_t sink);

/**
 * sfg_reduce_mason_compare — Compare reduction result with Mason's formula
 *
 * Verifies that sfg_reduce_to_single produces the same transfer function
 * as sfg_mason_compute.
 *
 * Returns: relative error between the two methods (should be < 1e-10)
 */
double sfg_reduce_mason_compare(const sfg_graph_t* g,
                                 sfg_node_id_t source,
                                 sfg_node_id_t sink);

#ifdef __cplusplus
}
#endif

#endif /* SFG_REDUCTION_H */
