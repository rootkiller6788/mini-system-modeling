/**
 * Example 1: DC Motor Speed Control Transfer Function
 *
 * Demonstrates: L6 Canonical System, L7 Application
 *
 * Models an armature-controlled DC motor and analyzes its
 * transfer function, stability, and step response.
 *
 * Realistic parameters for a small industrial DC motor
 * (similar to those used in Detroit auto plants / Toyota robotics).
 */

#include "transfer_function.h"
#include "tf_analysis.h"
#include "tf_interconnections.h"
#include "tf_advanced.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== DC Motor Speed Control Analysis ===\n\n");

    /* DC motor parameters:
     * J=0.01 kg*m^2 (rotor inertia)
     * b=0.1 N*m*s (viscous friction)
     * Kt=0.01 N*m/A (torque constant)
     * Ke=0.01 V*s/rad (back EMF, =Kt in SI)
     * Ra=1.0 Ohm (armature resistance)
     * La=0.5 H (armature inductance)
     */
    double J=0.01, b=0.1, Kt=0.01, Ke=0.01, Ra=1.0, La=0.5;

    TransferFunction *motor = tf_dc_motor_speed(J, b, Kt, Ke, Ra, La);

    printf("DC Motor Speed TF: G(s) = Omega(s)/V_a(s)\n");
    tf_print(motor);
    printf("\n");

    /* Properties */
    printf("Properties:\n");
    printf("  DC Gain (steady-state speed/volt): %.4f rad/(V*s)\n", tf_dc_gain(motor));
    printf("  System Type: %d\n", tf_system_type(motor));
    printf("  Relative Degree: %d\n", tf_relative_degree(motor));
    printf("  Properness: %s\n",
           tf_properness(motor) == TF_STRICTLY_PROPER ? "strictly proper" :
           tf_properness(motor) == TF_PROPER ? "proper" : "improper");
    printf("\n");

    /* Stability */
    tf_stability_t st = tf_stability_routh(motor);
    printf("Stability (Routh-Hurwitz): %s\n",
           st == TF_STABLE ? "STABLE" :
           st == TF_MARGINALLY_STABLE ? "MARGINALLY STABLE" : "UNSTABLE");

    /* Time domain characteristics */
    double po = tf_overshoot(motor);
    double tr = tf_rise_time(motor);
    double ts = tf_settling_time(motor, 0.02);
    printf("\nStep Response (2nd-order approx):\n");
    printf("  Overshoot: %.2f %%\n", po);
    printf("  Rise Time: %.4f s\n", tr);
    printf("  Settling Time (2%%): %.4f s\n", ts);

    /* Frequency response sample */
    printf("\nFrequency Response (Bode magnitude sample):\n");
    double freqs[] = {0.1, 1.0, 10.0, 100.0};
    int i;
    for (i = 0; i < 4; i++) {
        double mag = tf_magnitude_at(motor, freqs[i]);
        double mag_db = 20.0 * log10(mag);
        printf("  w=%.1f rad/s: %.4f (%.2f dB)\n", freqs[i], mag, mag_db);
    }

    tf_destroy(motor);

    /* Also demonstrate PID control for this motor */
    printf("\n--- PID Controller for DC Motor ---\n");
    TransferFunction *pid = tf_create_pid(10.0, 5.0, 1.0);
    printf("PID Controller: ");
    tf_print(pid);

    /* Closed-loop with unity feedback */
    TransferFunction *loop = tf_multiply(pid,
        tf_dc_motor_speed(J, b, Kt, Ke, Ra, La));
    TransferFunction *cl = tf_unity_feedback(loop);
    printf("\nClosed-loop TF (PID + Motor, unity FB):\n");
    tf_print(cl);

    tf_stability_t cl_st = tf_stability_routh(cl);
    printf("Closed-loop stability: %s\n",
           cl_st == TF_STABLE ? "STABLE" : "UNSTABLE");

    tf_destroy(pid); tf_destroy(loop); tf_destroy(cl);

    printf("\n=== DC Motor Analysis Complete ===\n");
    return 0;
}
