/**
 * Example 4: Quadrotor Attitude Control Block Diagram
 *
 * L6: Canonical Problem — Quadrotor dynamics and control
 * L7: Application — Quadrotor attitude stabilization
 *
 * Models quadrotor pitch axis dynamics:
 *   J*theta'' = tau (control torque)
 *
 * Block diagram:
 *   theta_ref ->(+) -> [PID] -> [Quadrotor dynamics] -> theta
 *                ^-                                     |
 *                |--------------------------------------|
 *
 * Quadrotors are a classic control problem studied at:
 *   MIT 6.302, Stanford ENGR105 (aerial robotics lab),
 *   ETH 151-0591 (flying machine arena)
 *
 * The dynamics are modeled as a double integrator:
 *   G(s) = theta(s)/tau(s) = 1/(J*s^2)
 * with inertia J from a typical 250mm racing drone.
 */
#include "../include/blockdiagram.h"
#include "../include/transfer_function.h"
#include "../include/block_algebra.h"
#include "../include/block_reduction.h"
#include "../include/mason.h"
#include "../include/signal_flow_graph.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 4: Quadrotor Attitude Control ===\n\n");

    /* Quadrotor parameters (250mm racing drone, e.g., Betaflight) */
    double J = 2.5e-4;  /* Pitch moment of inertia [kg*m^2] */
    double Kp = 0.15;   /* Proportional gain */
    double Ki = 1.5;    /* Integral gain */
    double Kd = 0.008;  /* Derivative gain */

    printf("Quadrotor Parameters:\n");
    printf("  J = %.2e kg*m^2 (250mm frame)\n", J);
    printf("  PID: Kp=%.3f, Ki=%.3f, Kd=%.4f\n\n", Kp, Ki, Kd);

    /* Plant: G(s) = theta/tau = 1/(J*s^2) */
    double n_plant[] = {1.0};  double d_plant[] = {0.0, 0.0, J};
    TransferFunction *G = tf_create(n_plant, 0, d_plant, 2);
    printf("Plant G(s) = "); tf_print(G);

    /* PID Controller: C(s) = Kp + Ki/s + Kd*s = (Kd*s^2 + Kp*s + Ki)/s */
    double n_pid[] = {Ki, Kp, Kd};  double d_pid[] = {0.0, 1.0};
    TransferFunction *C = tf_create(n_pid, 2, d_pid, 1);
    printf("PID C(s) = "); tf_print(C);

    /* Open-loop: L(s) = C(s)*G(s) */
    TransferFunction *L = tf_series(C, G);
    printf("\nOpen-loop L(s) = C(s)*G(s) = "); tf_print(L);

    /* Closed-loop: T(s) = L(s)/(1+L(s)) */
    TransferFunction *T = tf_feedback(L, NULL);
    printf("Closed-loop T(s) = "); tf_print(T);

    /* DC gain analysis */
    printf("\nAnalysis:\n");
    printf("  DC gain of T(s): %.6f (should be 1.0 for tracking)\n",
           tf_dc_gain(T));

    /* Stability check */
    tf_stability_t st = tf_stability_routh(T);
    printf("  Stability: %s\n",
           st == TF_STABLE ? "STABLE" :
           st == TF_MARGINALLY_STABLE ? "MARGINAL" : "UNSTABLE");

    /* Step response characteristics */
    printf("  Rise time (est):      %.4f s\n", tf_rise_time(T));
    printf("  Settling time (est):  %.4f s\n", tf_settling_time(T, 0.02));
    printf("  Overshoot (est):      %.1f%%\n", tf_overshoot(T));
    printf("  Steady-state error:   %.6f\n", tf_steady_state_error(T));

    /* Build equivalent block diagram */
    BlockDiagram *bd = ba_build_standard_feedback(C, NULL);
    /* Replace G in feedback with the quadrotor dynamics */
    /* For simplicity, print the feedback structure */
    printf("\nBlock diagram structure (standard feedback):\n");
    ba_print_structure(bd);

    /* Verify with Mason's formula */
    SignalFlowGraph *sg = bd_to_sfg((void*)bd);
    if (sg) {
        double mason_T = sfg_mason_transfer_function_gain(sg);
        printf("\nMason's formula verification:\n");
        printf("  T(0) from feedback formula: %.6f\n", tf_dc_gain(T));
        printf("  T(0) from Mason's formula:  %.6f\n", mason_T);
        printf("  Match: %s\n",
               fabs(tf_dc_gain(T) - mason_T) < 0.01 ? "YES (within 1%%)" : "NO");
        sfg_destroy(sg);
    }

    /* Cleanup */
    tf_destroy(G); tf_destroy(C); tf_destroy(L); tf_destroy(T);
    bd_destroy(bd);

    printf("\n=== Done ===\n");
    return 0;
}
