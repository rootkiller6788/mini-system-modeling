/**
 * sfg_matrix.h — Matrix Methods in Signal Flow Graph Analysis
 *
 * Declares matrix-algebraic functions that connect SFG theory to
 * linear algebra: Cramer's rule verification, resolvent computation,
 * eigenvalue bounds, and condition number estimation.
 */

#ifndef SFG_MATRIX_H
#define SFG_MATRIX_H

#include "sfg_core.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * sfg_cramer_verify — Verify Mason's formula against Cramer's rule
 *
 * L4 Theorem: Mason's formula is equivalent to Cramer's rule for
 * the linear system represented by the SFG.
 *
 * Returns: relative error (should be < 1e-10 for well-conditioned systems)
 */
double sfg_cramer_verify(const sfg_graph_t* g,
                          sfg_node_id_t source, sfg_node_id_t sink);

/**
 * sfg_resolvent_entry — Compute entry of resolvent matrix (sI-A)^(-1)
 */
double complex sfg_resolvent_entry(const sfg_graph_t* g,
                                    sfg_node_id_t state_src,
                                    sfg_node_id_t state_dst);

/**
 * sfg_eigenvalue_bounds — Gershgorin circle bounds from connection matrix
 */
int sfg_eigenvalue_bounds(const sfg_graph_t* g,
                           double* centers, double* radii, int N);

/**
 * sfg_condition_number — Estimate condition number κ(I-G) from SFG
 */
double sfg_condition_number(const sfg_graph_t* g);

/**
 * sfg_matrix_print_stats — Print matrix statistics from SFG
 */
void sfg_matrix_print_stats(const sfg_graph_t* g);

#ifdef __cplusplus
}
#endif

#endif /* SFG_MATRIX_H */
