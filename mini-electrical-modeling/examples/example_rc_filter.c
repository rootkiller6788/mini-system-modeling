/**
 * @file example_rc_filter.c
 * @brief End-to-end example: RC low-pass filter analysis and design
 *
 * Demonstrates: transfer function creation, Bode analysis, step response,
 * cutoff frequency search, filter design with standard components.
 *
 * Problem: Design a first-order RC low-pass filter with fc = 159 Hz,
 * analyze its frequency and time-domain response.
 *
 * Reference: Dorf & Svoboda (2018) Chapter 14 — Frequency Response
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "transfer_function_electrical.h"
#include "electrical_elements.h"

int main(void)
{
    printf("=== RC Low-Pass Filter Design Example ===\n\n");

    /* Problem: Design filter for fc = 159 Hz */
    double target_fc = 159.0;
    printf("Target cutoff frequency: %.1f Hz\n", target_fc);

    /* Design using standard component values */
    FilterDesign fd = design_rc_lowpass_to_cutoff(target_fc);
    printf("Selected R = %.0f ohm, C = %.2f nF\n",
           fd.r_standard, fd.c_standard * 1e9);
    printf("Actual fc = %.1f Hz (error: %.1f%%)\n\n",
           fd.actual_cutoff, fd.error_percent);

    /* Create transfer function */
    TransferFunction tf = tf_rc_lowpass(fd.r_standard, fd.c_standard);

    /* DC gain */
    printf("DC Gain: %.3f\n\n", tf_dc_gain(&tf));

    /* Frequency response at key points */
    printf("=== Frequency Response ===\n");
    printf("  Freq [Hz]    |H| [dB]    Phase [deg]\n");
    printf("  ---------    --------    -----------\n");
    double freqs[] = {1.0, 10.0, 50.0, 100.0, 159.0, 500.0, 1000.0, 10000.0};
    for (int i = 0; i < 8; i++) {
        double w = 2.0 * M_PI * freqs[i];
        double mag_db = tf_magnitude_db(&tf, w);
        double phase_deg = tf_phase_deg(&tf, w);
        printf("  %8.0f    %8.2f    %8.1f\n", freqs[i], mag_db, phase_deg);
    }

    /* Find precise -3dB cutoff */
    double wc = tf_find_cutoff_freq(&tf, 1.0, 10000.0);
    printf("\nMeasured -3dB cutoff: %.2f Hz (%.2f rad/s)\n",
           wc/(2.0*M_PI), wc);

    /* Step response */
    printf("\n=== Step Response ===\n");
    printf("  Time [ms]    v_out(t)/V_in\n");
    printf("  ---------    ------------\n");
    double tau = time_constant_rc(fd.r_standard, fd.c_standard);
    double times[] = {0.0, 0.1*tau, 0.5*tau, 1.0*tau, 2.0*tau, 3.0*tau, 5.0*tau};
    for (int i = 0; i < 7; i++) {
        double y = tf_step_response(&tf, times[i]);
        printf("  %8.3f    %8.3f\n", times[i]*1e3, y);
    }
    printf("  (tau = %.3f ms)\n", tau * 1e3);

    /* Group delay at low frequency */
    double tau_g = tf_group_delay(&tf, 10.0, 1.0);
    printf("\nGroup delay at 10 rad/s: %.3f ms\n", tau_g * 1e3);

    /* Bode asymptotic */
    printf("\n=== Asymptotic Bode Check ===\n");
    printf("  Freq [rad/s]  Exact [dB]  Asympt [dB]\n");
    printf("  ------------  ----------  -----------\n");
    double afreqs[] = {0.1/tau, 1.0/tau, 10.0/tau, 100.0/tau};
    for (int i = 0; i < 4; i++) {
        double exact = tf_magnitude_db(&tf, afreqs[i]);
        double asympt = tf_bode_magnitude_asymptotic(&tf, afreqs[i]);
        printf("  %12.1f  %10.2f  %10.2f\n", afreqs[i], exact, asympt);
    }

    printf("\n=== Design Complete ===\n");
    return 0;
}
