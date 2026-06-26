/**
 * sfg_mason.c — Mason's Gain Formula Implementation
 *
 * Full implementation of Mason's Gain Formula (MGF), the central
 * theorem of signal flow graph theory. Computes the exact transfer
 * function from any source to any sink in a linear system represented
 * as an SFG.
 *
 * Mason's Formula:
 *   T = (Σ P_k · Δ_k) / Δ
 *
 * where:
 *   Δ = 1 - ΣL₁ + ΣL₂ - ΣL₃ + ΣL₄ - ...  (graph determinant)
 *   P_k = gain of k-th forward path
 *   Δ_k = cofactor = determinant of subgraph not touching path k
 *
 * Knowledge coverage:
 *   L4: Mason's theorem, Cramer's rule connection, determinant structure
 *   L5: Full Mason evaluation algorithm
 *   L6: Transfer function computation for control systems
 *   L8: MIMO Mason extension
 */

#include "sfg_mason.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * L1: Result Structure Management
 * ================================================================ */

void sfg_mason_result_init(sfg_mason_result_t* result) {
    if (!result) return;
    memset(result, 0, sizeof(sfg_mason_result_t));
}

void sfg_mason_result_free(sfg_mason_result_t* result) {
    if (!result) return;
    sfg_path_list_clear(&result->forward_paths);
    sfg_loop_list_clear(&result->all_loops);
    sfg_nt_groups_clear(&result->nt_groups);
    memset(result, 0, sizeof(sfg_mason_result_t));
}

/* ================================================================
 * L5: Determinant Component Sums
 * ================================================================ */

sfg_gain_t sfg_sum_L1(const sfg_loop_list_t* loops) {
    sfg_gain_t sum = 0.0;
    if (!loops) return sum;
    for (int i = 0; i < loops->count; i++) {
        sum += loops->loops[i].gain;
    }
    return sum;
}

sfg_gain_t sfg_sum_L2(const sfg_nt_groups_t* groups) {
    sfg_gain_t sum = 0.0;
    if (!groups) return sum;
    for (int i = 0; i < groups->count; i++) {
        if (groups->groups[i].num_loops == 2) {
            sum += groups->groups[i].group_gain;
        }
    }
    return sum;
}

sfg_gain_t sfg_sum_L3(const sfg_nt_groups_t* groups) {
    sfg_gain_t sum = 0.0;
    if (!groups) return sum;
    for (int i = 0; i < groups->count; i++) {
        if (groups->groups[i].num_loops == 3) {
            sum += groups->groups[i].group_gain;
        }
    }
    return sum;
}

sfg_gain_t sfg_sum_L4(const sfg_nt_groups_t* groups) {
    sfg_gain_t sum = 0.0;
    if (!groups) return sum;
    for (int i = 0; i < groups->count; i++) {
        if (groups->groups[i].num_loops == 4) {
            sum += groups->groups[i].group_gain;
        }
    }
    return sum;
}

sfg_gain_t sfg_sum_L5_and_higher(const sfg_nt_groups_t* groups) {
    sfg_gain_t sum = 0.0;
    if (!groups) return sum;
    for (int i = 0; i < groups->count; i++) {
        if (groups->groups[i].num_loops >= 5) {
            sum += groups->groups[i].group_gain;
        }
    }
    return sum;
}

/* ================================================================
 * L4: Determinant and Cofactor Computation
 * ================================================================ */

sfg_gain_t sfg_determinant_compute(const sfg_graph_t* g,
                                    const sfg_loop_list_t* loops,
                                    const sfg_nt_groups_t* groups) {
    (void)g; /* g is reserved for future use (e.g., verification) */

    /*
     * Δ = 1 - ΣL₁ + ΣL₂ - ΣL₃ + ΣL₄ - ΣL₅_and_higher
     *
     * This is the alternating sum of non-touching loop group gains.
     * The formula is a direct consequence of Cramer's rule applied to
     * the system of linear equations (I - G)x = b, where G is the
     * gain matrix extracted from the SFG.
     *
     * L4 Theorem: For any linear SFG with at most finite number of loops,
     * the determinant series terminates (there are only finitely many
     * non-touching loop groups). The resulting Δ is the denominator
     * of the transfer function.
     */
    sfg_gain_t delta = 1.0 + 0.0 * I;  /* The "1" term */

    if (loops && loops->count > 0) {
        delta -= sfg_sum_L1(loops);             /* - ΣL₁ */
    }

    if (groups) {
        delta += sfg_sum_L2(groups);            /* + ΣL₂ */
        delta -= sfg_sum_L3(groups);            /* - ΣL₃ */
        delta += sfg_sum_L4(groups);            /* + ΣL₄ */
        delta -= sfg_sum_L5_and_higher(groups);  /* - ΣL₅₊ */
    }

    return delta;
}

sfg_gain_t sfg_cofactor_compute(const sfg_path_t* path,
                                 const sfg_loop_list_t* all_loops,
                                 const sfg_nt_groups_t* all_groups) {
    if (!path || !all_loops) return 1.0 + 0.0 * I;

    /*
     * Δ_k = cofactor of path k
     *     = value of Δ for subgraph obtained by removing all nodes on path k
     *     = 1 - (sum of loop gains not touching path k)
     *       + (sum of non-touching loop pairs not touching path k) - ...
     *
     * We compute this by filtering loops and non-touching groups
     * to those that don't touch the forward path.
     */

    /* Collect loops not touching the path */
    sfg_loop_list_t nt_loops;
    sfg_loop_list_clear(&nt_loops);

    for (int i = 0; i < all_loops->count; i++) {
        if (!sfg_path_loops_touch(path, &all_loops->loops[i])) {
            if (nt_loops.count < SFG_MAX_LOOPS) {
                nt_loops.loops[nt_loops.count++] = all_loops->loops[i];
            }
        }
    }

    /* Filter groups to those not touching the path.
     * A group doesn't touch the path if none of its loops touch it. */
    sfg_nt_groups_t nt_groups_filtered;
    sfg_nt_groups_clear(&nt_groups_filtered);

    if (all_groups) {
        for (int i = 0; i < all_groups->count; i++) {
            int touches = 0;
            for (int j = 0; j < all_groups->groups[i].num_loops; j++) {
                int loop_idx = all_groups->groups[i].loop_indices[j];
                if (loop_idx < all_loops->count
                    && sfg_path_loops_touch(path,
                        &all_loops->loops[loop_idx])) {
                    touches = 1;
                    break;
                }
            }
            if (!touches && nt_groups_filtered.count < SFG_MAX_NT_GROUPS) {
                nt_groups_filtered.groups[nt_groups_filtered.count++]
                    = all_groups->groups[i];
            }
        }
    }

    /* Compute determinant of the subgraph */
    sfg_gain_t cofactor = 1.0 + 0.0 * I;
    if (nt_loops.count > 0) {
        cofactor -= sfg_sum_L1(&nt_loops);
    }
    cofactor += sfg_sum_L2(&nt_groups_filtered);
    cofactor -= sfg_sum_L3(&nt_groups_filtered);
    cofactor += sfg_sum_L4(&nt_groups_filtered);
    cofactor -= sfg_sum_L5_and_higher(&nt_groups_filtered);

    return cofactor;
}

/* ================================================================
 * L5: Full Mason's Formula Computation
 * ================================================================ */

int sfg_mason_compute(const sfg_graph_t* g,
                      sfg_node_id_t source, sfg_node_id_t sink,
                      sfg_mason_result_t* result) {
    if (!g || !result) return -1;
    sfg_mason_result_init(result);

    if (!sfg_node_exists(g, source)) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "Source node %d does not exist", (int)source);
        return -1;
    }
    if (!sfg_node_exists(g, sink)) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "Sink node %d does not exist", (int)sink);
        return -1;
    }
    if (source == sink) {
        /* Trivial: signal passes through unchanged */
        result->transfer_function = 1.0 + 0.0 * I;
        result->determinant       = 1.0 + 0.0 * I;
        result->num_forward_paths = 1;
        result->path_gains[0]     = 1.0 + 0.0 * I;
        result->cofactors[0]      = 1.0 + 0.0 * I;
        result->success           = 1;
        return 0;
    }

    /* Step 1: Find all forward paths */
    sfg_path_list_clear(&result->forward_paths);
    int n_paths = sfg_find_forward_paths(g, source, sink,
                                          &result->forward_paths);
    if (n_paths < 0) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "Path enumeration failed");
        return -1;
    }
    result->num_forward_paths = n_paths;

    if (n_paths == 0) {
        /* No path from source to sink: transfer function is zero */
        result->transfer_function = 0.0;
        result->determinant       = 1.0 + 0.0 * I;
        result->success           = 1;
        return 0;
    }

    /* Extract path gains */
    for (int i = 0; i < n_paths; i++) {
        result->path_gains[i] = result->forward_paths.paths[i].gain;
    }

    /* Step 2: Find all loops */
    sfg_loop_list_clear(&result->all_loops);
    sfg_find_all_loops(g, &result->all_loops);

    /* Step 3: Find non-touching loop groups */
    sfg_nt_groups_clear(&result->nt_groups);
    sfg_find_non_touching_groups(g, &result->all_loops, &result->nt_groups);

    /* Step 4: Compute graph determinant Δ */
    result->sum_L1 = sfg_sum_L1(&result->all_loops);
    result->sum_L2 = sfg_sum_L2(&result->nt_groups);
    result->sum_L3 = sfg_sum_L3(&result->nt_groups);
    result->sum_L4 = sfg_sum_L4(&result->nt_groups);
    result->sum_L5 = sfg_sum_L5_and_higher(&result->nt_groups);

    result->determinant = sfg_determinant_compute(g, &result->all_loops,
                                                   &result->nt_groups);

    /* Check for singular graph (Δ = 0) */
    if (sfg_cmplx_is_zero(result->determinant, 1e-15)) {
        snprintf(result->error_msg, sizeof(result->error_msg),
                 "Graph determinant is zero (singular system)");
        result->success = 0;
        return -1;
    }

    /* Step 5: Compute cofactors and numerator */
    sfg_gain_t numerator = 0.0;
    for (int i = 0; i < n_paths; i++) {
        result->cofactors[i] = sfg_cofactor_compute(
            &result->forward_paths.paths[i],
            &result->all_loops, &result->nt_groups);
        numerator += result->path_gains[i] * result->cofactors[i];
    }

    /* Step 6: T = (Σ P_k Δ_k) / Δ */
    result->transfer_function = numerator / result->determinant;
    result->success = 1;

    return 0;
}

sfg_gain_t sfg_mason_compute_simple(const sfg_graph_t* g,
                                     sfg_node_id_t source,
                                     sfg_node_id_t sink) {
    sfg_mason_result_t result;
    sfg_mason_result_init(&result);
    int ret = sfg_mason_compute(g, source, sink, &result);
    sfg_gain_t tf = ret == 0 ? result.transfer_function : 0.0;
    sfg_mason_result_free(&result);
    return tf;
}

/* ================================================================
 * L8: MIMO Mason Extension
 * ================================================================ */

int sfg_mason_mimo(const sfg_graph_t* g,
                   const sfg_node_id_t* inputs, int M,
                   const sfg_node_id_t* outputs, int N,
                   sfg_gain_t* T) {
    if (!g || !inputs || !outputs || !T) return -1;
    if (M <= 0 || N <= 0) return -1;

    /* T[i][j] = transfer function from input j to output i
     * Stored row-major: T[i*M + j] */

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            sfg_mason_result_t result;
            sfg_mason_result_init(&result);

            int ret = sfg_mason_compute(g, inputs[j], outputs[i], &result);
            if (ret == 0) {
                T[i * M + j] = result.transfer_function;
            } else {
                T[i * M + j] = 0.0;  /* No path or singular */
            }

            sfg_mason_result_free(&result);
        }
    }

    return 0;
}

/* ================================================================
 * L4: Mason's Formula Verification
 * ================================================================ */

double sfg_mason_verify(const sfg_graph_t* g,
                        sfg_node_id_t source, sfg_node_id_t sink) {
    if (!g) return -1.0;
    if (!sfg_node_exists(g, source) || !sfg_node_exists(g, sink)) return -1.0;

    /* Compute via Mason's formula */
    sfg_mason_result_t mason_result;
    sfg_mason_result_init(&mason_result);
    int ret = sfg_mason_compute(g, source, sink, &mason_result);

    if (ret != 0) {
        sfg_mason_result_free(&mason_result);
        return -1.0;
    }

    sfg_gain_t T_mason = mason_result.transfer_function;

    /* Compute via direct matrix method:
     * (I - G) x = e_src, T = x[sink]
     * where G is the gain matrix (dst = G * src convention) */
    int N = g->num_nodes;
    if (N <= 0) {
        sfg_mason_result_free(&mason_result);
        return -1.0;
    }

    /* Build matrix (I - G) where G[j][i] = gain from i to j */
    double complex* M = (double complex*)calloc((size_t)(N * N),
                                                  sizeof(double complex));
    double complex* b = (double complex*)calloc((size_t)N,
                                                  sizeof(double complex));
    if (!M || !b) {
        free(M); free(b);
        sfg_mason_result_free(&mason_result);
        return -1.0;
    }

    /* M = I - G */
    for (int i = 0; i < N; i++) {
        M[i * N + i] = 1.0 + 0.0 * I;  /* Identity */
    }
    for (int k = 0; k < g->num_branches; k++) {
        int src = (int)g->branches[k].src;
        int dst = (int)g->branches[k].dst;
        /* G has entry at (dst, src) with value = gain */
        M[dst * N + src] -= g->branches[k].gain;
    }

    /* b = e_source (unit vector) */
    b[(int)source] = 1.0 + 0.0 * I;

    /* Solve M x = b using Gaussian elimination with partial pivoting */
    /* Augmented matrix [M | b] */
    double complex* aug = (double complex*)calloc((size_t)(N * (N + 1)),
                                                    sizeof(double complex));
    if (!aug) {
        free(M); free(b);
        sfg_mason_result_free(&mason_result);
        return -1.0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            aug[i * (N + 1) + j] = M[i * N + j];
        }
        aug[i * (N + 1) + N] = b[i];
    }

    /* Gaussian elimination */
    for (int col = 0; col < N; col++) {
        /* Find pivot */
        int pivot = col;
        double max_val = cabs(aug[col * (N + 1) + col]);
        for (int row = col + 1; row < N; row++) {
            double val = cabs(aug[row * (N + 1) + col]);
            if (val > max_val) {
                max_val = val;
                pivot = row;
            }
        }

        if (max_val < 1e-14) {
            /* Singular matrix */
            free(M); free(b); free(aug);
            sfg_mason_result_free(&mason_result);
            return -1.0;
        }

        /* Swap rows */
        if (pivot != col) {
            for (int j = 0; j <= N; j++) {
                double complex tmp = aug[col * (N + 1) + j];
                aug[col * (N + 1) + j] = aug[pivot * (N + 1) + j];
                aug[pivot * (N + 1) + j] = tmp;
            }
        }

        /* Eliminate below */
        for (int row = col + 1; row < N; row++) {
            double complex factor = aug[row * (N + 1) + col]
                                   / aug[col * (N + 1) + col];
            for (int j = col; j <= N; j++) {
                aug[row * (N + 1) + j] -= factor * aug[col * (N + 1) + j];
            }
        }
    }

    /* Back substitution */
    double complex* x = (double complex*)calloc((size_t)N,
                                                  sizeof(double complex));
    if (!x) {
        free(M); free(b); free(aug);
        sfg_mason_result_free(&mason_result);
        return -1.0;
    }

    for (int i = N - 1; i >= 0; i--) {
        double complex sum = aug[i * (N + 1) + N];
        for (int j = i + 1; j < N; j++) {
            sum -= aug[i * (N + 1) + j] * x[j];
        }
        x[i] = sum / aug[i * (N + 1) + i];
    }

    double complex T_direct = x[(int)sink];

    /* Compute relative error */
    double error = cabs(T_mason - T_direct) / (cabs(T_direct) + 1e-15);

    free(M); free(b); free(aug); free(x);
    sfg_mason_result_free(&mason_result);

    return error;
}

/* ================================================================
 * L5: Print / Display Functions
 * ================================================================ */

void sfg_mason_result_print(const sfg_mason_result_t* result) {
    if (!result) return;

    printf("\n==========================================\n");
    printf("  Mason's Gain Formula — Complete Result\n");
    printf("==========================================\n\n");

    if (!result->success) {
        printf("  ERROR: %s\n", result->error_msg);
        return;
    }

    /* Forward paths */
    printf("--- Forward Paths (P_k) ---\n");
    printf("  Number of forward paths: %d\n\n", result->num_forward_paths);
    for (int i = 0; i < result->num_forward_paths; i++) {
        printf("  P_%d: ", i + 1);
        sfg_path_print(&result->forward_paths.paths[i]);
    }

    /* Loops */
    printf("\n--- Loops (L_i) ---\n");
    printf("  Number of loops: %d\n\n", result->all_loops.count);
    for (int i = 0; i < result->all_loops.count; i++) {
        printf("  L_%d: ", i + 1);
        sfg_loop_print(&result->all_loops.loops[i]);
    }

    /* Non-touching groups */
    printf("\n--- Non-Touching Loop Groups ---\n");
    printf("  Number of groups: %d\n", result->nt_groups.count);
    for (int sz = 2; sz <= 5; sz++) {
        if (result->nt_groups.sizes[sz] > 0) {
            printf("  Size %d groups: %d\n", sz,
                   result->nt_groups.sizes[sz]);
        }
    }
    printf("\n");

    /* Determinant */
    printf("--- Determinant Δ ---\n");
    printf("  Δ = 1");
    printf(" - (");
    sfg_cmplx_print(result->sum_L1, NULL);
    printf(")");
    printf(" + (");
    sfg_cmplx_print(result->sum_L2, NULL);
    printf(")");
    printf(" - (");
    sfg_cmplx_print(result->sum_L3, NULL);
    printf(")");
    printf(" + (");
    sfg_cmplx_print(result->sum_L4, NULL);
    printf(")");
    printf(" - (");
    sfg_cmplx_print(result->sum_L5, NULL);
    printf(")\n");

    printf("  Δ = ");
    sfg_cmplx_print(result->determinant, NULL);
    printf("\n");

    /* Cofactors */
    printf("\n--- Cofactors (Δ_k) ---\n");
    for (int i = 0; i < result->num_forward_paths; i++) {
        printf("  Δ_%d = ", i + 1);
        sfg_cmplx_print(result->cofactors[i], NULL);
    }

    /* Transfer function */
    printf("\n--- Transfer Function ---\n");
    printf("  T = (");
    for (int i = 0; i < result->num_forward_paths; i++) {
        if (i > 0) printf(" + ");
        printf("P_%d·Δ_%d", i + 1, i + 1);
    }
    printf(") / Δ\n");

    printf("  T = ");
    sfg_cmplx_print(result->transfer_function, NULL);

    printf("  |T| = %.6g,  ∠T = %.2f°\n",
           cabs(result->transfer_function),
           carg(result->transfer_function) * 180.0 / 3.14159265358979323846);

    printf("\n==========================================\n\n");
}

void sfg_mason_summary_print(const sfg_mason_result_t* result) {
    if (!result) return;
    printf("Mason: T = ");
    sfg_cmplx_print(result->transfer_function, NULL);
    printf("  (Δ = ");
    sfg_cmplx_print(result->determinant, NULL);
    printf(", %d paths, %d loops)\n",
           result->num_forward_paths, result->all_loops.count);
}
