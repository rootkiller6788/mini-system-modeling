/**
 * @example example_sdof_vibration.c
 * @brief End-to-end SDOF vibration analysis --- free and forced response
 *
 * Demonstrates: Natural frequency, damping ratio, step response,
 * frequency response, transmissibility calculation.
 * Reference: Thomson & Dahleh (1998), Ch.2-3.
 */

#include <stdio.h>
#include <math.h>
#include "mechanical_elements.h"
#include "mechanical_transfer_function.h"

int main(void)
{
    printf("\n=== SDOF Vibration Analysis ===\n\n");

    /* System parameters: mass=100kg, k=100000 N/m, b=2000 N*s/m */
    double m = 100.0;
    double k = 100000.0;
    double b = 2000.0;

    /* Natural frequency and damping */
    double wn = natural_frequency_undamped(k, m);
    double fn = natural_frequency_hz(k, m);
    double zeta = damping_ratio(b, k, m);
    double b_cr = critical_damping_coefficient(k, m);
    double wd = damped_natural_frequency(k, m, b);

    printf("System Parameters:\n");
    printf("  Mass m      = %.2f kg\n", m);
    printf("  Stiffness k = %.0f N/m\n", k);
    printf("  Damping b   = %.0f N*s/m\n\n", b);

    printf("Dynamic Characteristics:\n");
    printf("  Natural freq fn = %.3f Hz (omega_n = %.3f rad/s)\n", fn, wn);
    printf("  Damping ratio zeta = %.4f\n", zeta);
    printf("  Critical damping = %.1f N*s/m\n", b_cr);

    if (zeta < 1.0) {
        printf("  Damped natural freq = %.3f rad/s\n", wd);
        printf("  System is UNDERDAMPED\n");
    } else if (fabs(zeta - 1.0) < 1e-9) {
        printf("  System is CRITICALLY DAMPED\n");
    } else {
        printf("  System is OVERDAMPED\n");
    }

    /* Static deflection */
    double delta_st = static_deflection(m, k, 9.81);
    printf("  Static deflection = %.4f m (%.1f mm)\n\n", delta_st, delta_st*1000.0);

    /* Time-domain response */
    double OS = percent_overshoot(zeta);
    double tp = peak_time_underdamped(zeta, wn);
    double ts = settling_time_2percent(zeta, wn);
    double tr = rise_time_underdamped(zeta, wn);

    printf("Step Response Specifications:\n");
    printf("  Percent overshoot = %.1f %%\n", OS);
    printf("  Peak time tp = %.4f s\n", tp);
    printf("  Settling time ts (2%%) = %.4f s\n", ts);
    printf("  Rise time tr = %.4f s\n\n", tr);

    /* Step response values */
    printf("Step Response Values x(t) [m]:\n");
    double times[] = {0.0, 0.5*tp, tp, 2.0*tp, ts, 2.0*ts};
    const char *labels[] = {"t=0", "t=0.5*tp", "t=tp (peak)", "t=2*tp", "t=ts", "t=2*ts"};
    for (int i = 0; i < 6; i++) {
        double x = sdof_step_response_displacement(delta_st, wn, zeta, times[i]);
        printf("  %-15s = %.6f m\n", labels[i], x);
    }

    /* Frequency response */
    printf("\nFrequency Response |H(j*w)| [m/N]:\n");
    printf("  Freq [rad/s]  |H|           |H|_dB       Phase [deg]\n");
    double freqs[] = {0.1*wn, 0.5*wn, wn, 1.5*wn, 5.0*wn};
    MechanicalTF tf = tf_sdof_displacement(m, b, k);
    for (int i = 0; i < 5; i++) {
        double H = tf_mag(&tf, freqs[i]);
        double Hdb = tf_mag_db(&tf, freqs[i]);
        double phase = tf_phase(&tf, freqs[i]) * 180.0 / M_PI;
        printf("  %10.2f    %10.6f    %8.2f dB    %8.1f deg\n", freqs[i], H, Hdb, phase);
    }

    /* Transmissibility */
    printf("\nVibration Isolation Analysis:\n");
    printf("  r     TR     Isolation %%\n");
    for (int i = 1; i <= 5; i++) {
        double r = 0.5 * i;
        double TR = force_transmissibility(r, zeta);
        double eff = isolation_efficiency(r, zeta);
        printf("  %.1f   %.4f   %.1f %%\n", r, TR, eff);
    }

    /* Quality factor */
    double Q = quality_factor_from_damping(zeta);
    double bw = half_power_bandwidth(zeta, wn);
    printf("\n  Quality factor Q = %.3f\n", Q);
    printf("  Half-power bandwidth = %.3f rad/s\n", bw);

    printf("\n=== Analysis Complete ===\n");
    return 0;
}
