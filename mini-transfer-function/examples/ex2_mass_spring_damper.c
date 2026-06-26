/**
 * Example 2: Mass-Spring-Damper System
 *
 * Demonstrates: L6 Canonical System, L7 Application
 *
 * Models a mechanical mass-spring-damper and analyzes
 * natural frequency, damping ratio, and oscillation.
 */

#include "transfer_function.h"
#include "tf_analysis.h"
#include "tf_advanced.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Mass-Spring-Damper System Analysis ===\n\n");

    /* Suspension-like parameters:
     * m=250 kg (quarter-car sprung mass)
     * b=1000 N*s/m (damper)
     * k=20000 N/m (spring)
     */
    double m=250.0, bv=1000.0, k=20000.0;
    TransferFunction *msd = tf_mass_spring_damper(m, bv, k);

    printf("Mass-Spring-Damper TF: G(s)=X(s)/F(s)\n");
    printf("  m=%.1f kg, b=%.1f N*s/m, k=%.1f N/m\n", m, bv, k);
    tf_print(msd);
    printf("\n");

    /* Analysis */
    double wn = sqrt(k/m);
    double zeta = bv/(2.0*sqrt(m*k));
    printf("Natural Frequency: wn = %.4f rad/s (%.4f Hz)\n", wn, wn/(2*M_PI));
    printf("Damping Ratio: zeta = %.4f\n", zeta);
    printf("DC Gain (static deflection): %.6f m/N\n", tf_dc_gain(msd));

    if (zeta < 1.0) {
        printf("System is UNDERDAMPED - will oscillate\n");
        double wd = wn * sqrt(1.0 - zeta*zeta);
        printf("Damped natural frequency: wd = %.4f rad/s\n", wd);
    } else if (zeta == 1.0) {
        printf("System is CRITICALLY DAMPED\n");
    } else {
        printf("System is OVERDAMPED - no oscillation\n");
    }

    /* Step response */
    printf("\nStep Response Characteristics:\n");
    printf("  Rise Time: %.4f s\n", tf_rise_time(msd));
    printf("  Overshoot: %.2f %%\n", tf_overshoot(msd));
    printf("  Settling Time (2%%): %.4f s\n", tf_settling_time(msd, 0.02));
    printf("  Peak Time: %.4f s\n", tf_peak_time(msd));

    /* Frequency response */
    printf("\nFrequency Response:\n");
    printf("  DC Gain: %.6f\n", tf_dc_gain(msd));
    printf("  Bandwidth (rad/s): %.4f\n", tf_bandwidth(msd));

    double rp = tf_resonant_peak(msd);
    double rf = tf_resonant_frequency(msd);
    if (rp > 1.0) {
        printf("  Resonant Peak: %.4f at %.4f rad/s\n", rp, rf);
    }

    /* Stability check */
    tf_stability_t st = tf_stability_routh(msd);
    printf("\nStability: %s (expected stable for positive b, m, k)\n",
           st == TF_STABLE ? "STABLE" : "UNSTABLE");

    /* Compare with quarter-car suspension */
    printf("\n--- Comparison: Quarter-Car Suspension ---\n");
    TransferFunction *susp = tf_quarter_car_suspension(m, 40.0, k, bv, 180000.0);
    printf("Quarter-car TF (4th order):\n");
    tf_print(susp);
    printf("Stability: %s\n",
           tf_stability_routh(susp) == TF_STABLE ? "STABLE" : "UNSTABLE");

    tf_destroy(msd); tf_destroy(susp);

    printf("\n=== MSD Analysis Complete ===\n");
    return 0;
}
