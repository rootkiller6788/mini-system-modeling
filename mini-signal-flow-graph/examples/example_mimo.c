/**
 * example_mimo.c — MIMO System via Multi-Input Multi-Output SFG
 *
 * L8 Advanced: Multi-input multi-output system analysis using SFG.
 * Demonstrates the MIMO extension of Mason's Gain Formula for
 * computing the complete transfer function matrix.
 *
 * System: 2-input, 2-output coupled system
 *   u1 → [+a11] → Σ → x1 → [+c11] → y1
 *   u2 → [+b12] → Σ → x2 → [+c22] → y2
 *   Cross-coupling: x1 ⇄ x2
 */

#include "sfg_core.h"
#include "sfg_mason.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <math.h>

int main(void) {
    printf("==========================================\n");
    printf("  Example 5: MIMO System via SFG\n");
    printf("==========================================\n\n");

    /*
     * Build a 2×2 MIMO SFG:
     *
     * u1 ──→[a11]→ Σ₁ → x1 ──→[c11]→ y1
     *              ↑    ↕  [a12/a21]
     * u2 ──→[b22]→ Σ₂ → x2 ──→[c22]→ y2
     *
     * With values: a11=2, a22=3, a12=0.5, a21=-0.3
     *              b11=1, b22=1
     *              c11=1, c22=1
     */
    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t u1 = sfg_add_source(&g, "u1");
    sfg_node_id_t u2 = sfg_add_source(&g, "u2");
    sfg_node_id_t s1 = sfg_add_summing(&g, "Σ1");
    sfg_node_id_t s2 = sfg_add_summing(&g, "Σ2");
    sfg_node_id_t x1 = sfg_add_state(&g, "x1");
    sfg_node_id_t x2 = sfg_add_state(&g, "x2");
    sfg_node_id_t y1 = sfg_add_sink(&g, "y1");
    sfg_node_id_t y2 = sfg_add_sink(&g, "y2");

    /* Input to summing junctions */
    sfg_add_branch_real(&g, u1, s1, 1.0, NULL);
    sfg_add_branch_real(&g, u2, s2, 1.0, NULL);

    /* Self-gains (a11, a22) */
    sfg_add_branch_real(&g, x1, s1, 2.0, NULL);
    sfg_add_branch_real(&g, x2, s2, 3.0, NULL);

    /* Cross-coupling (a12, a21) */
    sfg_add_branch_real(&g, x2, s1, 0.5, NULL);
    sfg_add_branch_real(&g, x1, s2, -0.3, NULL);

    /* State integrators */
    sfg_add_branch_real(&g, s1, x1, 1.0, NULL);
    sfg_add_branch_real(&g, s2, x2, 1.0, NULL);

    /* Output gains */
    sfg_add_branch_real(&g, x1, y1, 1.0, NULL);
    sfg_add_branch_real(&g, x2, y2, 1.0, NULL);

    printf("MIMO Signal Flow Graph:\n");
    sfg_graph_print(&g);

    /* Compute MIMO transfer function matrix */
    sfg_node_id_t inputs[2]  = {u1, u2};
    sfg_node_id_t outputs[2] = {y1, y2};
    sfg_gain_t T[4]; /* 2×2 matrix, row-major */

    int ret = sfg_mason_mimo(&g, inputs, 2, outputs, 2, T);

    if (ret == 0) {
        printf("\n--- MIMO Transfer Function Matrix ---\n");
        printf("T = [ T11  T12 ]\n");
        printf("    [ T21  T22 ]\n\n");

        for (int i = 0; i < 2; i++) {
            printf("  ");
            for (int j = 0; j < 2; j++) {
                sfg_cmplx_print(T[i * 2 + j], NULL);
                printf("  ");
            }
            printf("\n");
        }

        /* Individual TF computation for verification */
        printf("\n--- Individual TF via Mason ---\n");
        sfg_gain_t t11 = sfg_mason_compute_simple(&g, u1, y1);
        printf("T11 = u1→y1: ");
        sfg_cmplx_print(t11, NULL);

        sfg_gain_t t12 = sfg_mason_compute_simple(&g, u2, y1);
        printf("T12 = u2→y1: ");
        sfg_cmplx_print(t12, NULL);

        sfg_gain_t t21 = sfg_mason_compute_simple(&g, u1, y2);
        printf("T21 = u1→y2: ");
        sfg_cmplx_print(t21, NULL);

        sfg_gain_t t22 = sfg_mason_compute_simple(&g, u2, y2);
        printf("T22 = u2→y2: ");
        sfg_cmplx_print(t22, NULL);
    }

    printf("\nDone.\n");
    return 0;
}
