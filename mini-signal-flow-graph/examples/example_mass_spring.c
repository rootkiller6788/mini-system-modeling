/**
 * example_mass_spring.c — Mass-Spring-Damper System via SFG
 *
 * L6 Canonical Problem: Model a mechanical mass-spring-damper system
 * and compute the transfer function from force to displacement using
 * Mason's Gain Formula on the signal flow graph.
 *
 * System: F → [1/m] → a → [1/s] → v → [1/s] → x
 *         Feedback: x → [-k] → Σ; v → [-b] → Σ
 *
 * Transfer function: X(s)/F(s) = 1/(ms² + bs + k)
 *
 * This demonstrates the classic mechanical system modeling pattern.
 */

#include "sfg_core.h"
#include "sfg_mason.h"
#include "sfg_state_space.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("==========================================\n");
    printf("  Example 1: Mass-Spring-Damper via SFG\n");
    printf("==========================================\n\n");

    double m = 1.0;  /* Mass (kg) */
    double b = 0.5;  /* Damping coefficient (N·s/m) */
    double k = 2.0;  /* Spring constant (N/m) */

    printf("Parameters: m=%.1f kg, b=%.1f N·s/m, k=%.1f N/m\n\n", m, b, k);

    /* Build SFG:
     * Node 0: F (source)
     * Node 1: Σ (summing junction for forces)
     * Node 2: a = ẋ_dot = acceleration state (state derivative)
     * Node 3: v = ẋ (state — velocity)
     * Node 4: x = state (position)
     * Node 5: x_out (sink)
     */
    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t F_in   = sfg_add_source(&g, "F");
    sfg_node_id_t sum_F  = sfg_add_summing(&g, "ΣF");
    sfg_node_id_t a      = sfg_add_node(&g, SFG_NODE_INTERNAL, "a");
    sfg_node_id_t v      = sfg_add_state(&g, "v");
    sfg_node_id_t x_pos  = sfg_add_state(&g, "x");
    sfg_node_id_t x_out  = sfg_add_sink(&g, "x_out");

    /* Forward path: F → Σ → a → v → x → x_out */
    sfg_add_branch_real(&g, F_in, sum_F, 1.0, NULL);
    sfg_add_branch_real(&g, sum_F, a, 1.0 / m, NULL);  /* a = F/m */
    sfg_add_branch_real(&g, a, v, 1.0, NULL);            /* v = ∫a dt */
    sfg_add_branch_real(&g, v, x_pos, 1.0, NULL);        /* x = ∫v dt */
    sfg_add_branch_real(&g, x_pos, x_out, 1.0, NULL);

    /* Feedback paths (negative feedback) */
    sfg_add_branch_real(&g, x_pos, sum_F, -k, NULL);      /* Spring: -k·x */
    sfg_add_branch_real(&g, v, sum_F, -b, NULL);          /* Damper: -b·v */

    sfg_graph_print(&g);

    /* Compute transfer function via Mason's formula */
    sfg_mason_result_t result;
    sfg_mason_result_init(&result);

    int ret = sfg_mason_compute(&g, F_in, x_out, &result);
    if (ret == 0) {
        printf("\n--- Mason's Formula Result ---\n");
        printf("Number of forward paths: %d\n", result.num_forward_paths);
        printf("Number of loops: %d\n", result.all_loops.count);
        printf("Transfer function (DC gain): T = ");
        sfg_cmplx_print(result.transfer_function, NULL);

        /* Expected DC gain: 1/k = 1/2 = 0.5 */
        double dc_gain = creal(result.transfer_function);
        printf("\nExpected DC gain: 1/k = %.4f\n", 1.0 / k);
        printf("Computed DC gain: %.4f\n", dc_gain);
        printf("Match: %s\n", fabs(dc_gain - 1.0/k) < 1e-6 ? "YES" : "NO");
    }

    sfg_mason_result_free(&result);

    printf("\n--- State-Space Representation ---\n");
    sfg_ss_system_t ss;
    sfg_ss_init(&ss);
    sfg_ss_alloc(&ss, 2, 1, 1);
    /* ẋ = v, v̇ = (-k/m)x + (-b/m)v + (1/m)F */
    sfg_ss_set_A(&ss, 0, 0, 0.0);
    sfg_ss_set_A(&ss, 0, 1, 1.0);
    sfg_ss_set_A(&ss, 1, 0, -k/m);
    sfg_ss_set_A(&ss, 1, 1, -b/m);
    sfg_ss_set_B(&ss, 0, 0, 0.0);
    sfg_ss_set_B(&ss, 1, 0, 1.0/m);
    sfg_ss_set_C(&ss, 0, 0, 1.0);
    sfg_ss_set_C(&ss, 0, 1, 0.0);
    sfg_ss_print(&ss);

    /* Characteristic polynomial */
    double poly[3];
    sfg_ss_characteristic_polynomial(&ss, poly, 2);
    printf("\nCharacteristic polynomial: %.3f + %.3f s + %.3f s²\n",
           poly[0], poly[1], poly[2]);
    printf("Expected: k + b s + m s² = %.3f + %.3f s + %.3f s²\n",
           k, b, m);

    sfg_ss_free(&ss);

    printf("\nDone.\n");
    return 0;
}
