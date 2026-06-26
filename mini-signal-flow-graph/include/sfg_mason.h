/**
 * sfg_mason.h — Mason's Gain Formula Implementation
 *
 * Mason's Gain Formula (MGF) is the central theorem of signal flow
 * graph theory. It provides a closed-form expression for the transfer
 * function (gain) from any source to any sink in a linear system
 * represented as an SFG.
 *
 * Mason's Formula (L4 Theorem):
 *   T = (1/Δ) * Σ_k P_k * Δ_k
 *
 * where:
 *   T     = overall transfer function (system gain)
 *   P_k   = gain of the k-th forward path
 *   Δ     = graph determinant
 *         = 1 - ΣL_i + ΣL_i·L_j - ΣL_i·L_j·L_k + ...
 *         (sum of loop gains, minus sum of non-touching pair gains,
 *          plus sum of non-touching triple gains, etc.)
 *   Δ_k   = cofactor of the k-th forward path
 *         = value of Δ for the subgraph obtained by removing all
 *           nodes on the k-th forward path
 *         = 1 - (sum of loop gains not touching path k) + ...
 *
 * References:
 *  - Mason, S.J. "Feedback Theory — Some Properties of Signal Flow Graphs"
 *    Proc. IRE, Vol. 41, No. 9, pp. 1144-1156, September 1953.
 *  - Mason, S.J. "Feedback Theory — Further Properties of Signal Flow
 *    Graphs" Proc. IRE, Vol. 44, No. 7, pp. 920-926, July 1956.
 *  - D'Azzo, J.J., Houpis, C.H. "Linear Control System Analysis and
 *    Design" McGraw-Hill, 4th Ed., 1995. Chapter 3.
 */

#ifndef SFG_MASON_H
#define SFG_MASON_H

#include "sfg_core.h"
#include "sfg_path.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: Mason's Formula Result Types
 * ================================================================ */

/**
 * sfg_mason_result_t — Complete Mason's formula computation result
 *
 * Stores all intermediate quantities for transparency and verification.
 */
typedef struct {
    sfg_gain_t transfer_function;     /* T = overall gain from input to output */

    /* Numerator components */
    sfg_path_list_t forward_paths;    /* All forward paths P_k */
    sfg_gain_t      path_gains[SFG_MAX_PATHS]; /* Gain of each forward path */
    sfg_gain_t      cofactors[SFG_MAX_PATHS];  /* Δ_k for each forward path */
    int             num_forward_paths;

    /* Determinant components */
    sfg_loop_list_t all_loops;        /* All individual loops */
    sfg_nt_groups_t nt_groups;        /* Non-touching loop groups */
    sfg_gain_t      determinant;      /* Δ = graph determinant */

    /* Stratified determinant terms */
    sfg_gain_t sum_L1;               /* Σ L_i (sum of individual loop gains) */
    sfg_gain_t sum_L2;               /* Σ L_i L_j (non-touching pairs) */
    sfg_gain_t sum_L3;               /* Σ L_i L_j L_k (non-touching triples) */
    sfg_gain_t sum_L4;               /* Σ L_i L_j L_k L_m (non-touching quads) */
    sfg_gain_t sum_L5;               /* Sum of quintets and higher (combined) */

    /* Status */
    int    success;                   /* 1 if computation completed */
    char   error_msg[256];            /* Error description if failed */
} sfg_mason_result_t;

/* ================================================================
 * L5: Mason's Formula Computation
 * ================================================================ */

/**
 * sfg_mason_compute — Compute transfer function using Mason's Gain Formula
 *
 * L5 Algorithm: Full implementation of Mason's formula.
 *
 * Steps:
 *   1. Find all forward paths from input to output (DFS)
 *   2. Find all loops in the graph (Johnson's algorithm)
 *   3. Compute graph determinant Δ:
 *      Δ = 1 - ΣL₁ + ΣL₂ - ΣL₃ + ΣL₄ - ... (alternating signs)
 *   4. For each forward path P_k:
 *      a. Identify all loops that do NOT touch P_k
 *      b. Compute Δ_k from these non-touching loops
 *   5. Sum over paths: T = (Σ P_k · Δ_k) / Δ
 *
 * L4 Theorem: This formula is an algebraic consequence of Cramer's rule
 * applied to the system of linear equations represented by the SFG.
 * It is exact for any linear, time-invariant system.
 *
 * Complexity: O(P·2^L) worst case where P = forward paths, L = loops.
 *             For typical control systems (≤10 loops), very fast.
 *
 * @param g       The signal flow graph
 * @param source  Source node ID
 * @param sink    Sink node ID
 * @param result  Output: complete Mason's formula result
 * @return        0 on success, -1 on error
 */
int sfg_mason_compute(const sfg_graph_t* g,
                      sfg_node_id_t source, sfg_node_id_t sink,
                      sfg_mason_result_t* result);

/**
 * sfg_mason_compute_simple — Simplified: return just the transfer function
 *
 * Convenience wrapper around sfg_mason_compute.
 *
 * Complexity: Same as sfg_mason_compute.
 * Returns: transfer function gain, or 0+0i on error
 */
sfg_gain_t sfg_mason_compute_simple(const sfg_graph_t* g,
                                     sfg_node_id_t source,
                                     sfg_node_id_t sink);

/**
 * sfg_mason_result_init — Initialize a result structure
 */
void sfg_mason_result_init(sfg_mason_result_t* result);

/**
 * sfg_mason_result_free — Free internal allocations in result
 */
void sfg_mason_result_free(sfg_mason_result_t* result);

/* ================================================================
 * L5: Determinant Decomposition
 * ================================================================ */

/**
 * sfg_determinant_compute — Compute graph determinant Δ only
 *
 * L4 Conservation: The determinant Δ captures the complete feedback
 * structure of the system. Δ = 0 indicates the presence of a pole at
 * the origin or an algebraic loop with infinite gain.
 *
 * Δ = 0 → system is not well-posed (singular).
 *
 * Complexity: O(L·2^L) worst case, O(L²) typical
 */
sfg_gain_t sfg_determinant_compute(const sfg_graph_t* g,
                                    const sfg_loop_list_t* loops,
                                    const sfg_nt_groups_t* groups);

/**
 * sfg_cofactor_compute — Compute cofactor Δ_k for a specific forward path
 *
 * Δ_k = determinant of the subgraph formed by removing all nodes
 * on the given forward path. Only loops not touching the path contribute.
 *
 * Complexity: O(L·2^L) worst case
 */
sfg_gain_t sfg_cofactor_compute(const sfg_path_t* path,
                                 const sfg_loop_list_t* all_loops,
                                 const sfg_nt_groups_t* all_groups);

/* ================================================================
 * L5: Stratified Determinant Terms
 * ================================================================ */

/**
 * sfg_sum_L1 — Sum of all individual loop gains: Σ L_i
 *
 * This is the first correction term in the determinant expansion.
 * For systems with only single loops (no non-touching pairs),
 * Δ = 1 - ΣL₁.
 */
sfg_gain_t sfg_sum_L1(const sfg_loop_list_t* loops);

/**
 * sfg_sum_L2 — Sum of products of non-touching loop pairs: Σ L_i · L_j
 *
 * Only non-touching pairs contribute (touching pairs are excluded).
 * Sign is positive in the determinant expansion.
 */
sfg_gain_t sfg_sum_L2(const sfg_nt_groups_t* groups);

/**
 * sfg_sum_L3 — Sum of products of non-touching loop triples: Σ L_i·L_j·L_k
 *
 * Sign is negative in the determinant expansion.
 */
sfg_gain_t sfg_sum_L3(const sfg_nt_groups_t* groups);

/**
 * sfg_sum_L4 — Sum of products of non-touching loop quadruplets
 */
sfg_gain_t sfg_sum_L4(const sfg_nt_groups_t* groups);

/**
 * sfg_sum_L5_and_higher — Combined sum of quintets and above
 */
sfg_gain_t sfg_sum_L5_and_higher(const sfg_nt_groups_t* groups);

/* ================================================================
 * L5: Multi-Input Multi-Output (MIMO) Mason
 * ================================================================ */

/**
 * sfg_mason_mimo — Compute MIMO transfer function matrix via Mason
 *
 * L8 Advanced: For a system with M inputs and N outputs,
 * compute the complete N×M transfer function matrix T(s)
 * where T[i][j] = gain from input j to output i.
 *
 * Each entry is computed independently using Mason's formula.
 *
 * @param g          Signal flow graph
 * @param inputs     Array of M source node IDs
 * @param M          Number of inputs
 * @param outputs    Array of N sink node IDs
 * @param N          Number of outputs
 * @param T          Output: N×M matrix (row-major, T[i*M + j])
 * @return           0 on success, -1 on error
 *
 * Complexity: O(M·N·P·2^L)
 */
int sfg_mason_mimo(const sfg_graph_t* g,
                   const sfg_node_id_t* inputs, int M,
                   const sfg_node_id_t* outputs, int N,
                   sfg_gain_t* T);

/**
 * sfg_mason_verify — Verify Mason's formula against direct matrix inversion
 *
 * L4 Validation: For a given SFG, construct the connection matrix,
 * solve (I - C)x = b directly via LU decomposition, and compare
 * against Mason's formula result.
 *
 * Returns: relative error between Mason and direct method
 */
double sfg_mason_verify(const sfg_graph_t* g,
                        sfg_node_id_t source, sfg_node_id_t sink);

/* ================================================================
 * L5: Mason's Formula Print / Debug
 * ================================================================ */

/**
 * sfg_mason_result_print — Print complete Mason's formula derivation
 *
 * Displays: all forward paths with gains, all loops, determinant terms,
 * cofactors, and final transfer function. Useful for educational
 * purposes and debugging.
 */
void sfg_mason_result_print(const sfg_mason_result_t* result);

/**
 * sfg_mason_summary_print — Print concise summary (T and Δ only)
 */
void sfg_mason_summary_print(const sfg_mason_result_t* result);

#ifdef __cplusplus
}
#endif

#endif /* SFG_MASON_H */
