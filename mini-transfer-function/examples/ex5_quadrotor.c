/**
 * Example 5: Quadrotor Attitude Control
 *
 * Demonstrates: L6 Canonical System, L7 Application
 *
 * Models a quadrotor UAV's single-axis attitude dynamics
 * (double integrator) and designs a PD controller.
 * Relevant to drone control and SpaceX landing stability.
 */

#include "transfer_function.h"
#include "tf_analysis.h"
#include "tf_interconnections.h"
#include "tf_advanced.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Quadrotor Attitude Control Analysis ===\n\n");

    /* Quadrotor parameters (small UAV):
     * I=0.01 kg*m^2 (moment of inertia about pitch/roll axis)
     */
    double I=0.01;
    TransferFunction *attitude = tf_quadrotor_attitude(I);

    printf("Quadrotor Attitude TF: Theta(s)/Torque(s) = 1/(I*s^2)\n");
    printf("  Moment of inertia I = %.4f kg*m^2\n", I);
    tf_print(attitude);
    printf("\n");

    /* Properties */
    printf("System Properties:\n");
    printf("  System Type: %d (double integrator)\n", tf_system_type(attitude));
    printf("  Relative Degree: %d\n", tf_relative_degree(attitude));
    printf("  Stability: MARGINALLY STABLE (poles at s=0)\n");
    printf("  Phase Margin: %.2f deg (infinite for pure integrator)\n",
           tf_phase_margin(attitude));
    printf("\n");

    /* PD Controller Design */
    printf("--- PD Controller Design ---\n");
    double Kp = 4.0, Kd = 2.0;
    TransferFunction *pd = tf_create_pid(Kp, 0.0, Kd);
    printf("PD Controller: Kp=%.1f, Kd=%.1f\n", Kp, Kd);
    tf_print(pd);
    printf("\n");

    /* Closed-loop analysis */
    TransferFunction *loop = tf_series(pd, attitude);
    TransferFunction *cl = tf_unity_feedback(loop);

    printf("Closed-Loop TF (PD + Quadrotor):\n");
    tf_print(cl);
    printf("\n");

    /* Closed-loop properties:
     * CL denominator: I*s^2 + Kd*s + Kp
     * wn = sqrt(Kp/I), zeta = Kd/(2*sqrt(Kp*I)) */
    double wn_cl = sqrt(Kp/I);
    double zeta_cl = Kd/(2.0*sqrt(Kp*I));
    printf("Closed-Loop Properties:\n");
    printf("  Natural Frequency: wn = %.2f rad/s (%.2f Hz)\n",
           wn_cl, wn_cl/(2.0*M_PI));
    printf("  Damping Ratio: zeta = %.4f\n", zeta_cl);
    printf("  Stability: %s\n",
           tf_stability_routh(cl) == TF_STABLE ? "STABLE" : "UNSTABLE");

    /* Step response */
    printf("\nStep Response (unit step torque command):\n");
    double po = tf_overshoot(cl);
    double tr = tf_rise_time(cl);
    double ts = tf_settling_time(cl, 0.02);
    printf("  Overshoot: %.2f %%\n", po);
    printf("  Rise Time: %.4f s\n", tr);
    printf("  Settling Time (2%%): %.4f s\n", ts);
    printf("  Steady-State Error: %.6f\n", tf_steady_state_error(cl));

    /* Frequency response */
    printf("\nFrequency Response around resonance:\n");
    printf("  %-12s %-12s %-12s\n", "w (rad/s)", "|G|", "Phase (deg)");
    double test_w[] = {0.1*wn_cl, 0.5*wn_cl, wn_cl, 2.0*wn_cl, 5.0*wn_cl};
    int i;
    for (i = 0; i < 5; i++) {
        double w = test_w[i];
        double mag = tf_magnitude_at(cl, w);
        double phase = tf_phase_at(cl, w) * 180.0 / M_PI;
        printf("  %-12.2f %-12.4f %-12.2f\n", w, mag, phase);
    }

    tf_destroy(attitude); tf_destroy(pd);
    tf_destroy(loop); tf_destroy(cl);

    printf("\n=== Quadrotor Analysis Complete ===\n");
    return 0;
}
