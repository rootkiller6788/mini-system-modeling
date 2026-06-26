/**
 * example_control_tf.c — Classic Control System Transfer Function via SFG
 *
 * L6 Canonical Problem: Compure the closed-loop transfer function of a
 * unity-feedback control system using Mason's Gain Formula.
 *
 * System: R(s) → Σ → G_c(s) → G_p(s) → Y(s)
 *                    ↑__________________|
 *
 * where:
 *   G_c(s) = K_p + K_i/s      (PI controller)
 *   G_p(s) = 1/(s² + 2ζω_n s + ω_n²)  (second-order plant)
 *
 * This demonstrates SFG analysis of a practical control system.
 *
 * L7 Application: Used in DC motor speed control (Katy pattern),
 * process control (supplier chain), and aerospace positioning systems
 * (Boeing 747 autopilot reference).
 */

#include "sfg_core.h"
#include "sfg_mason.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <math.h>

int main(void) {
    printf("==========================================\n");
    printf("  Example 3: Control System TF via SFG\n");
    printf("==========================================\n\n");

    /* DC motor parameters */
    double K_p = 2.0;       /* Proportional gain */
    double K_i = 5.0;       /* Integral gain */
    double J   = 0.01;      /* Motor inertia (kg·m²) */
    double b_m = 0.1;       /* Motor viscous friction (N·m·s) */
    double K_t = 0.01;      /* Torque constant (N·m/A) */
    double K_e = 0.01;      /* Back-EMF constant (V·s/rad) */
    double R_a = 1.0;       /* Armature resistance (Ω) */

    printf("DC Motor Parameters:\n");
    printf("  J=%.3f, b=%.2f, K_t=%.3f, K_e=%.3f, R_a=%.1f\n",
           J, b_m, K_t, K_e, R_a);
    printf("  Controller: K_p=%.1f, K_i=%.1f\n\n", K_p, K_i);

    /*
     * Build SFG for DC motor with PI control (unity feedback).
     *
     * Nodes:
     *   0: r (reference input — source)
     *   1: e (error — summing junction)
     *   2: u_p (proportional control signal)
     *   3: u_i (integral control signal)
     *   4: u (total control signal — summing)
     *   5: V_a (armature voltage)
     *   6: I_a (armature current → torque)
     *   7: ω (angular velocity — state)
     *   8: θ (angular position — state)
     *   9: y (output — sink)
     *   10: fb (feedback node)
     */

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t r     = sfg_add_source(&g, "r");
    sfg_node_id_t e     = sfg_add_summing(&g, "e");
    sfg_node_id_t u_p   = sfg_add_node(&g, SFG_NODE_INTERNAL, "u_p");
    sfg_node_id_t u_i   = sfg_add_node(&g, SFG_NODE_INTERNAL, "u_i");
    sfg_node_id_t u_sum = sfg_add_summing(&g, "u");
    sfg_node_id_t V_a   = sfg_add_node(&g, SFG_NODE_INTERNAL, "V_a");
    sfg_node_id_t I_a   = sfg_add_node(&g, SFG_NODE_INTERNAL, "I_a");
    sfg_node_id_t T_m   = sfg_add_node(&g, SFG_NODE_INTERNAL, "T_m");
    sfg_node_id_t omega_dot = sfg_add_node(&g, SFG_NODE_INTERNAL, "ω_dot");
    sfg_node_id_t omega = sfg_add_state(&g, "ω");
    sfg_node_id_t theta_dot = sfg_add_node(&g, SFG_NODE_INTERNAL, "θ_dot");
    sfg_node_id_t theta = sfg_add_state(&g, "θ");
    sfg_node_id_t y     = sfg_add_sink(&g, "y");
    sfg_node_id_t fb    = sfg_add_node(&g, SFG_NODE_PICKOFF, "fb");

    /* Forward path: r → e → (K_p + K_i/s) → u → ... → y */
    sfg_add_branch_real(&g, r, e, 1.0, NULL);

    /* Proportional path: e → u_p */
    sfg_add_branch_real(&g, e, u_p, K_p, NULL);
    sfg_add_branch_real(&g, u_p, u_sum, 1.0, NULL);

    /* Integral path: e → u_i (with integrator) */
    sfg_add_branch_real(&g, e, u_i, K_i, NULL);
    sfg_add_branch_real(&g, u_i, u_sum, 1.0, NULL);

    /* Plant: u_sum → V_a → I_a → T_m → ω_dot → ω */
    sfg_add_branch_real(&g, u_sum, V_a, 1.0, NULL);
    sfg_add_branch_real(&g, V_a, I_a, 1.0 / R_a, NULL);
    sfg_add_branch_real(&g, I_a, T_m, K_t, NULL);
    sfg_add_branch_real(&g, T_m, omega_dot, 1.0, NULL);
    sfg_add_branch_real(&g, omega_dot, omega, 1.0, NULL);

    /* Feedback EMF: ω → V_a (back-EMF) */
    sfg_add_branch_real(&g, omega, V_a, -K_e, NULL);

    /* Viscous friction: ω → T_m (opposing torque) */
    sfg_add_branch_real(&g, omega, T_m, -b_m, NULL);

    /* Output */
    sfg_add_branch_real(&g, omega, y, 1.0, NULL);

    /* Unity feedback */
    sfg_add_branch_real(&g, omega, fb, 1.0, NULL);
    sfg_add_branch_real(&g, fb, e, -1.0, NULL);

    printf("Signal Flow Graph for DC Motor with PI Control:\n");
    sfg_graph_print(&g);

    /* Mason's formula */
    sfg_mason_result_t result;
    sfg_mason_result_init(&result);

    int ret = sfg_mason_compute(&g, r, y, &result);
    if (ret == 0) {
        printf("\n--- Mason's Formula Result ---\n");
        sfg_mason_summary_print(&result);
        printf("\nDetailed result:\n");
        sfg_mason_result_print(&result);

        double dc_gain = creal(result.transfer_function);
        printf("DC gain (velocity): %.6f\n", dc_gain);
        printf("Expected DC gain ≈ 1.0 (unity feedback)\n");
    } else {
        printf("Mason computation failed: %s\n", result.error_msg);
    }

    sfg_mason_result_free(&result);

    printf("\nDone.\n");
    return 0;
}
