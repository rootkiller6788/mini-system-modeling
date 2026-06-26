/**
 * example_aircraft.c — Aircraft Pitch Control via SFG
 *
 * L7 Application: Aircraft pitch dynamics modeled as a signal flow
 * graph. The short-period approximation of longitudinal dynamics
 * is analyzed using Mason's Gain Formula.
 *
 * Boeing 747 reference (approximate cruise conditions):
 *   - Pitch rate to elevator transfer function
 *   - Angle of attack dynamics
 *
 * Longitudinal short-period equations:
 *   α̇ = Z_α·α + q + Z_δe·δe
 *   q̇ = M_α·α + M_q·q + M_δe·δe
 *
 * L7 Application keywords: Boeing, 747, aircraft, GNC
 */

#include "sfg_core.h"
#include "sfg_mason.h"
#include "sfg_state_space.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <math.h>

int main(void) {
    printf("==========================================\n");
    printf("  Example 4: Aircraft Pitch Control via SFG\n");
    printf("==========================================\n\n");

    /* Boeing 747 approximate short-period parameters (cruise) */
    double Z_alpha = -0.5;     /* Z-force due to α */
    double Z_delta = -0.05;    /* Z-force due to elevator */
    double M_alpha = -2.0;     /* Pitch moment due to α */
    double M_q     = -0.6;     /* Pitch damping */
    double M_delta = -10.0;    /* Control effectiveness */

    printf("Aircraft Short-Period Parameters:\n");
    printf("  Z_α=%.2f, Z_δe=%.2f, M_α=%.2f, M_q=%.2f, M_δe=%.2f\n\n",
           Z_alpha, Z_delta, M_alpha, M_q, M_delta);

    /*
     * SFG for short-period dynamics:
     *
     * δe → [Z_δe] → Σ → [1] → α
     *                 ↑
     * δe → [M_δe] → Σ → [1] → q
     *           ↑      ↑     ↓
     *           +--[M_q]--+  [1]→ α
     *           +--[M_α]--α
     */

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t de_in  = sfg_add_source(&g, "δe");
    sfg_node_id_t a_sum  = sfg_add_summing(&g, "Σα");
    sfg_node_id_t alpha  = sfg_add_state(&g, "α");
    sfg_node_id_t q_sum  = sfg_add_summing(&g, "Σq");
    sfg_node_id_t q_out  = sfg_add_state(&g, "q");
    sfg_node_id_t q_sink = sfg_add_sink(&g, "q_out");

    /* α equation: α̇ = Z_α·α + q + Z_δe·δe */
    sfg_add_branch_real(&g, de_in, a_sum, Z_delta, NULL);
    sfg_add_branch_real(&g, alpha, a_sum, Z_alpha, NULL);
    sfg_add_branch_real(&g, q_out, a_sum, 1.0, NULL);
    sfg_add_branch_real(&g, a_sum, alpha, 1.0, NULL);

    /* q equation: q̇ = M_α·α + M_q·q + M_δe·δe */
    sfg_add_branch_real(&g, de_in, q_sum, M_delta, NULL);
    sfg_add_branch_real(&g, alpha, q_sum, M_alpha, NULL);
    sfg_add_branch_real(&g, q_out, q_sum, M_q, NULL);
    sfg_add_branch_real(&g, q_sum, q_out, 1.0, NULL);

    /* Output */
    sfg_add_branch_real(&g, q_out, q_sink, 1.0, NULL);

    sfg_graph_print(&g);

    /* Compute transfer function q/δe */
    sfg_mason_result_t result;
    sfg_mason_result_init(&result);

    int ret = sfg_mason_compute(&g, de_in, q_sink, &result);
    if (ret == 0) {
        printf("\n--- Pitch Rate Transfer Function q/δe ---\n");
        sfg_mason_summary_print(&result);
        printf("\nDC gain (pitch rate per elevator degree): %.6f rad/s\n",
               creal(result.transfer_function));
    } else {
        printf("Computation error: %s\n", result.error_msg);
    }

    sfg_mason_result_free(&result);

    /* State-space verification */
    printf("\n--- State-Space Check ---\n");
    sfg_ss_system_t ss;
    sfg_ss_init(&ss);
    sfg_ss_alloc(&ss, 2, 1, 1);
    sfg_ss_set_A(&ss, 0, 0, Z_alpha);
    sfg_ss_set_A(&ss, 0, 1, 1.0);
    sfg_ss_set_A(&ss, 1, 0, M_alpha);
    sfg_ss_set_A(&ss, 1, 1, M_q);
    sfg_ss_set_B(&ss, 0, 0, Z_delta);
    sfg_ss_set_B(&ss, 1, 0, M_delta);
    sfg_ss_set_C(&ss, 0, 0, 0.0);
    sfg_ss_set_C(&ss, 0, 1, 1.0);

    double poly[3];
    sfg_ss_characteristic_polynomial(&ss, poly, 2);
    printf("Characteristic polynomial: %.3f + %.3f s + s²\n", poly[0], poly[1]);
    printf("Damping ratio ζ = %.3f\n", poly[1] / (2.0 * sqrt(poly[0])));
    printf("Natural frequency ω_n = %.3f rad/s\n", sqrt(poly[0]));

    sfg_ss_free(&ss);

    printf("\nDone.\n");
    return 0;
}
