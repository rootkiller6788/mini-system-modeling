/**
 * sfg_matrix.c — Matrix Methods in Signal Flow Graph Analysis
 *
 * Explores the deep connection between SFGs and linear algebra.
 * The signal flow graph is fundamentally a graphical representation
 * of the matrix equation (I - G)x = b, where G is the gain matrix.
 *
 * Topics:
 *   - Connection matrix and its relationship to SFG structure
 *   - Determinant computation via SFG (Cramer's rule connection)
 *   - Resolvent matrix (sI - A)^(-1) from SFG
 *   - Eigenvalue extraction from SFG structure
 *   - Gaussian elimination interpreted as SFG reduction
 *
 * Knowledge coverage:
 *   L4: Cramer's rule ↔ Mason's formula equivalence
 *   L3: Matrix representations of SFG topology
 *   L8: Resolvent matrix and eigenvalue analysis via SFG
 */

#include "sfg_matrix.h"
#include "sfg_core.h"
#include "sfg_complex.h"
#include "sfg_mason.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * L4: Cramer's Rule Verification via SFG
 * ================================================================ */

/**
 * sfg_cramer_verify — Verify that Mason's formula gives the same
 * result as Cramer's rule for the linear system represented by SFG.
 *
 * L4 Theorem: For an SFG representing (I - G)x = Bu,
 * the transfer function from input u_j to output y_i computed by
 * Mason's formula equals the ratio of determinants from Cramer's rule:
 *   T = det(M_j) / det(I - G)
 * where M_j is I - G with the j-th column replaced by the input vector.
 *
 * Returns: relative error between the two methods.
 */
double sfg_cramer_verify(const sfg_graph_t* g,
                          sfg_node_id_t source, sfg_node_id_t sink) {
    if (!g) return -1.0;

    int N = g->num_nodes;
    if (N <= 0 || N > 32) return -1.0; /* Limit for determinant computation */

    /* Compute via Mason */
    sfg_gain_t T_mason = sfg_mason_compute_simple(g, source, sink);

    /* Build matrices for Cramer's rule: (I - G)x = e_src, T = x[sink] */

    /* Connection matrix C: C[src*N + dst] = gain from src to dst */
    double complex* M = (double complex*)calloc((size_t)(N * N),
                                                  sizeof(double complex));
    if (!M) return -1.0;

    /* M = I - C^T (where C[src*N + dst] = gain) */
    for (int i = 0; i < N; i++) {
        M[i * N + i] = 1.0 + 0.0 * I;
    }
    for (int b = 0; b < g->num_branches; b++) {
        int src = (int)g->branches[b].src;
        int dst = (int)g->branches[b].dst;
        /* For (I - G)x = b form, G has entry at (dst, src) */
        M[dst * N + src] -= g->branches[b].gain;
    }

    /* Compute determinant of M = det(I - G) */
    /* Using recursive cofactor expansion for small matrices,
     * or LU decomposition for larger ones. */

    /* Copy M for LU */
    double complex* LU = (double complex*)calloc((size_t)(N * N),
                                                   sizeof(double complex));
    if (!LU) { free(M); return -1.0; }
    memcpy(LU, M, (size_t)(N * N) * sizeof(double complex));

    /* LU decomposition with partial pivoting */
    int* pivot = (int*)calloc((size_t)N, sizeof(int));
    if (!pivot) { free(M); free(LU); return -1.0; }
    for (int i = 0; i < N; i++) pivot[i] = i;

    int sign = 1;
    for (int col = 0; col < N; col++) {
        /* Find pivot */
        int max_row = col;
        double max_val = cabs(LU[col * N + col]);
        for (int row = col + 1; row < N; row++) {
            double val = cabs(LU[row * N + col]);
            if (val > max_val) {
                max_val = val;
                max_row = row;
            }
        }

        if (max_val < 1e-15) {
            free(M); free(LU); free(pivot);
            return -1.0; /* Singular */
        }

        if (max_row != col) {
            sign = -sign;
            for (int j = 0; j < N; j++) {
                double complex tmp = LU[col * N + j];
                LU[col * N + j] = LU[max_row * N + j];
                LU[max_row * N + j] = tmp;
            }
        }

        for (int row = col + 1; row < N; row++) {
            double complex factor = LU[row * N + col] / LU[col * N + col];
            LU[row * N + col] = factor;
            for (int j = col + 1; j < N; j++) {
                LU[row * N + j] -= factor * LU[col * N + j];
            }
        }
    }

    /* det = product of diagonal elements × sign */
    double complex det_M = (double complex)sign;
    for (int i = 0; i < N; i++) {
        det_M *= LU[i * N + i];
    }

    /* For Cramer's rule: T = x[sink] where (I-G)x = b
     * The solution vector x is computed directly. */
    /* Build augmented system (I-G)x = e_src and solve directly */

    double complex* b_cramer = (double complex*)calloc((size_t)N,
                                                         sizeof(double complex));
    if (!b_cramer) { free(M); free(LU); free(pivot); return -1.0; }
    b_cramer[(int)source] = 1.0 + 0.0 * I;

    /* Apply forward substitution using LU */
    double complex* y = (double complex*)calloc((size_t)N,
                                                  sizeof(double complex));
    if (!y) { free(M); free(LU); free(pivot); free(b_cramer); return -1.0; }

    /* Forward: Ly = Pb */
    for (int i = 0; i < N; i++) {
        y[i] = b_cramer[i];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i * N + j] * y[j];
        }
    }

    /* Backward: Ux = y */
    double complex* x_cramer = (double complex*)calloc((size_t)N,
                                                         sizeof(double complex));
    if (!x_cramer) {
        free(M); free(LU); free(pivot); free(b_cramer); free(y);
        return -1.0;
    }

    for (int i = N - 1; i >= 0; i--) {
        x_cramer[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x_cramer[i] -= LU[i * N + j] * x_cramer[j];
        }
        x_cramer[i] /= LU[i * N + i];
    }

    double complex T_cramer = x_cramer[(int)sink];

    double error = cabs(T_mason - T_cramer) / (cabs(T_cramer) + 1e-15);

    free(M); free(LU); free(pivot);
    free(b_cramer); free(y); free(x_cramer);

    return error;
}

/* ================================================================
 * L3: Resolvent Matrix via SFG
 * ================================================================ */

/**
 * sfg_resolvent_entry — Compute a single entry of the resolvent
 * matrix (sI - A)^(-1) using SFG methods.
 *
 * In the SFG representation of a state-space system, the transfer
 * function from state j to state i equals the (i,j) entry of the
 * resolvent matrix: [(sI - A)^(-1)]_{ij}.
 *
 * This provides an alternative to matrix inversion for computing
 * the resolvent when the system has a sparse SFG representation.
 */
double complex sfg_resolvent_entry(const sfg_graph_t* g,
                                    sfg_node_id_t state_src,
                                    sfg_node_id_t state_dst) {
    if (!g) return 0.0;
    return sfg_mason_compute_simple(g, state_src, state_dst);
}

/**
 * sfg_eigenvalue_bound — Compute Gershgorin bounds on eigenvalues
 * from SFG connection matrix.
 *
 * L3: Gershgorin circle theorem: Every eigenvalue λ of matrix A
 * satisfies: |λ - a_{ii}| ≤ Σ_{j≠i} |a_{ij}| for some i.
 *
 * From the SFG connection matrix, we can bound the eigenvalues
 * without computing them explicitly.
 */
int sfg_eigenvalue_bounds(const sfg_graph_t* g,
                           double* centers, double* radii, int N) {
    if (!g || !centers || !radii) return -1;
    if (N != g->num_nodes) return -1;

    /* Build connection matrix and compute Gershgorin disks */
    for (int i = 0; i < N; i++) {
        centers[i] = 0.0;
        radii[i]   = 0.0;
    }

    for (int b = 0; b < g->num_branches; b++) {
        int src = (int)g->branches[b].src;
        int dst = (int)g->branches[b].dst;

        if (src == dst) {
            /* Self-loop: diagonal element of connection matrix */
            centers[src] = creal(g->branches[b].gain);
        } else {
            /* Off-diagonal: contributes to radius of src row */
            radii[src] += cabs(g->branches[b].gain);
        }
    }

    return 0;
}

/**
 * sfg_condition_number — Estimate condition number of (I-G) from SFG
 *
 * L4: The condition number κ(I-G) = ||I-G|| · ||(I-G)^(-1)||
 * determines numerical sensitivity. A large condition number
 * indicates an ill-conditioned system.
 *
 * Using the SFG, we estimate the condition number from the
 * determinant Δ and the largest loop gain.
 */
double sfg_condition_number(const sfg_graph_t* g) {
    if (!g) return -1.0;

    /* Find all loops to estimate condition number */
    sfg_loop_list_t loops;
    sfg_loop_list_clear(&loops);
    sfg_find_all_loops(g, &loops);

    /* Estimate: κ ≈ max(1, Σ|L_i|) / Δ */
    double sum_loop_mag = 0.0;
    double max_loop_mag = 0.0;
    for (int i = 0; i < loops.count; i++) {
        double mag = cabs(loops.loops[i].gain);
        sum_loop_mag += mag;
        if (mag > max_loop_mag) max_loop_mag = mag;
    }

    /* Estimate determinant magnitude */
    sfg_gain_t sum_L1 = sfg_sum_L1(&loops);
    double delta_est = fabs(1.0 - cabs(sum_L1));
    /* More accurate: use the full determinant if available */
    /* For coarse estimate: Δ ≈ 1 - ΣL₁ */

    if (delta_est < 1e-15) return INFINITY;

    double norm_I_minus_G_est = 1.0 + sum_loop_mag;
    double norm_inv_est = 1.0 / delta_est;

    double kappa = norm_I_minus_G_est * norm_inv_est;
    return kappa;
}

/**
 * sfg_matrix_print_stats — Print matrix statistics derived from SFG
 */
void sfg_matrix_print_stats(const sfg_graph_t* g) {
    if (!g) return;
    printf("=== SFG Matrix Statistics ===\n");
    printf("  Nodes: %d, Branches: %d\n", g->num_nodes, g->num_branches);

    int N = g->num_nodes;
    if (N <= 0) return;

    /* Connection density */
    double density = (double)g->num_branches / (double)(N * N);
    printf("  Connection density: %.4f\n", density);

    /* Gershgorin bounds */
    double* centers = (double*)calloc((size_t)N, sizeof(double));
    double* radii   = (double*)calloc((size_t)N, sizeof(double));
    if (centers && radii) {
        sfg_eigenvalue_bounds(g, centers, radii, N);
        double max_radius = 0.0;
        for (int i = 0; i < N; i++) {
            if (radii[i] > max_radius) max_radius = radii[i];
        }
        printf("  Max Gershgorin radius: %.6f\n", max_radius);
        printf("  Spectral radius bound: %.6f\n", max_radius);
        free(centers);
        free(radii);
    }

    /* Condition number estimate */
    double kappa = sfg_condition_number(g);
    printf("  Condition number estimate κ: %.2e\n", kappa);

    /* Sparsity analysis */
    int n_self_loops = 0;
    int n_parallel  = 0;
    for (int i = 0; i < g->num_branches; i++) {
        if (g->branches[i].src == g->branches[i].dst) n_self_loops++;
    }
    /* Count parallel branches */
    for (int i = 0; i < g->num_nodes; i++) {
        for (int j = 0; j < g->num_nodes; j++) {
            int count = 0;
            for (int b = 0; b < g->num_branches; b++) {
                if ((int)g->branches[b].src == i
                    && (int)g->branches[b].dst == j) count++;
            }
            if (count > 1) n_parallel += (count - 1);
        }
    }
    printf("  Self-loops: %d\n", n_self_loops);
    printf("  Parallel branch groups: %d\n", n_parallel);

    printf("==============================\n");
}
