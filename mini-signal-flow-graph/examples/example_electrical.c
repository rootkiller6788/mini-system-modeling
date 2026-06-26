/**
 * example_electrical.c — RLC Circuit Analysis via SFG
 *
 * L6 Canonical Problem: Analyze a series RLC circuit using signal
 * flow graphs and Mason's Gain Formula.
 *
 * Circuit: Voltage source V_in → R → L → C → ground
 * Output: Voltage across capacitor V_c
 *
 * State equations:
 *   L·di/dt = V_in - R·i - V_c
 *   C·dV_c/dt = i
 *
 * Transfer function: V_c(s)/V_in(s) = 1/(LCs² + RCs + 1)
 *
 * L7 Application: This SFG approach generalizes to complex circuit
 * analysis (Kirchhoff's laws → SFG), used in circuit simulators
 * and electronic design automation (EDA) tools.
 */

#include "sfg_core.h"
#include "sfg_mason.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("==========================================\n");
    printf("  Example 2: RLC Circuit via SFG\n");
    printf("==========================================\n\n");

    double R = 10.0;    /* Resistance (Ω) */
    double L = 0.001;   /* Inductance (H) = 1 mH */
    double C = 1e-6;    /* Capacitance (F) = 1 μF */

    printf("Parameters: R=%.0fΩ, L=%.0e H, C=%.0e F\n\n", R, L, C);

    /*
     * SFG for RLC circuit:
     *
     * V_in → [1] → Σ → [1/L] → i_dot → [1] → i → [1/C] → Vc_dot → [1] → Vc
     *            ↑              ↑
     *            |   [-R]       |   [-1]
     *            +--- i --------+--- Vc
     *
     * Where:
     *   i_dot = di/dt = (V_in - R·i - V_c) / L
     *   Vc_dot = dVc/dt = i / C
     */

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t V_in     = sfg_add_source(&g, "V_in");
    sfg_node_id_t sum_volt = sfg_add_summing(&g, "ΣV");
    sfg_node_id_t i_dot    = sfg_add_node(&g, SFG_NODE_INTERNAL, "i_dot");
    sfg_node_id_t i_l      = sfg_add_state(&g, "i_L");
    sfg_node_id_t Vc_dot   = sfg_add_node(&g, SFG_NODE_INTERNAL, "Vc_dot");
    sfg_node_id_t V_c      = sfg_add_state(&g, "V_c");
    sfg_node_id_t Vc_out   = sfg_add_sink(&g, "Vc_out");

    /* Forward path */
    sfg_add_branch_real(&g, V_in, sum_volt, 1.0, NULL);
    sfg_add_branch_real(&g, sum_volt, i_dot, 1.0 / L, NULL);
    sfg_add_branch_real(&g, i_dot, i_l, 1.0, NULL);
    sfg_add_branch_real(&g, i_l, Vc_dot, 1.0 / C, NULL);
    sfg_add_branch_real(&g, Vc_dot, V_c, 1.0, NULL);
    sfg_add_branch_real(&g, V_c, Vc_out, 1.0, NULL);

    /* Feedback branches (negative) */
    sfg_add_branch_real(&g, i_l, sum_volt, -R, NULL);   /* -R·i */
    sfg_add_branch_real(&g, V_c, sum_volt, -1.0, NULL); /* -V_c */

    printf("Signal Flow Graph:\n");
    sfg_graph_print(&g);

    /* Mason's formula */
    sfg_mason_result_t result;
    sfg_mason_result_init(&result);

    int ret = sfg_mason_compute(&g, V_in, Vc_out, &result);
    if (ret == 0) {
        printf("\n--- Mason's Formula Result ---\n");
        sfg_mason_result_print(&result);

        /* Expected DC gain: 1 (capacitor charges to Vin) */
        double dc_gain = creal(result.transfer_function);
        printf("DC gain = %.6f (expected 1.0)\n", dc_gain);

        /* Check determinant */
        printf("Determinant Δ = ");
        sfg_cmplx_print(result.determinant, NULL);
    }

    sfg_mason_result_free(&result);

    /* Numerical verification: frequency at resonance ω₀ = 1/√(LC) */
    double omega_0 = 1.0 / sqrt(L * C);
    printf("\nResonant frequency: ω₀ = %.2e rad/s (%.2f kHz)\n",
           omega_0, omega_0 / (2.0 * 3.14159265358979323846 * 1000.0));
    printf("Quality factor: Q = 1/R · √(L/C) = %.3f\n",
           1.0/R * sqrt(L/C));

    printf("\nDone.\n");
    return 0;
}
