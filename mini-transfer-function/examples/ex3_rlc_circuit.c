/**
 * Example 3: RLC Circuit Analysis
 *
 * Demonstrates: L6 Canonical System, L7 Application
 *
 * Models a series RLC circuit as a low-pass filter and analyzes
 * its frequency response and step response characteristics.
 * Relevant to Tesla inverter EMI filter design.
 */

#include "transfer_function.h"
#include "tf_analysis.h"
#include "tf_advanced.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== RLC Circuit Transfer Function Analysis ===\n\n");

    /* Tesla inverter EMI filter parameters:
     * R=0.1 Ohm (parasitic resistance)
     * L=100 uH = 0.0001 H (filter inductor)
     * C=10 uF = 0.00001 F (filter capacitor)
     */
    double R=0.1, L=0.0001, C=0.00001;
    TransferFunction *rlc = tf_rlc_lowpass(R, L, C);

    printf("RLC Low-Pass Filter: Vc(s)/Vin(s)\n");
    printf("  R=%.4f Ohm, L=%.6f H, C=%.8f F\n", R, L, C);
    tf_print(rlc);
    printf("\n");

    /* Natural frequency and damping */
    double wn = 1.0/sqrt(L*C);
    double zeta = R/2.0 * sqrt(C/L);
    printf("Natural Frequency: wn = %.2f rad/s (%.2f kHz)\n",
           wn, wn/(2000.0*M_PI));
    printf("Damping Ratio: zeta = %.6f\n", zeta);
    printf("Quality Factor: Q = 1/(2*zeta) = %.2f\n",
           zeta > 0.0 ? 1.0/(2.0*zeta) : INFINITY);
    printf("\n");

    /* Frequency response sweep */
    printf("Frequency Response (Bode):\n");
    printf("  %-12s %-12s %-12s %-12s\n", "w (rad/s)", "|G|", "|G| dB", "Phase (deg)");
    double freqs[] = {0.1*wn, 0.5*wn, 1.0*wn, 2.0*wn, 10.0*wn};
    int i;
    for (i = 0; i < 5; i++) {
        double w = freqs[i];
        double mag = tf_magnitude_at(rlc, w);
        double mag_db = 20.0 * log10(mag);
        double phase = tf_phase_at(rlc, w) * 180.0 / M_PI;
        printf("  %-12.1f %-12.6f %-12.2f %-12.2f\n",
               w, mag, mag_db, phase);
    }

    /* Bandwidth and resonance */
    printf("\nBandwidth (-3 dB): %.2f rad/s\n", tf_bandwidth(rlc));
    double rp = tf_resonant_peak(rlc);
    if (rp > 1.0)
        printf("Resonant Peak: %.4f at %.2f rad/s\n",
               rp, tf_resonant_frequency(rlc));

    /* Also demonstrate bandpass configuration */
    printf("\n--- RLC Bandpass (across resistor) ---\n");
    TransferFunction *bp = tf_rlc_bandpass(R, L, C);
    printf("Bandpass TF: ");
    tf_print(bp);
    printf("Bandpass center frequency: %.2f rad/s (resonance)\n", wn);

    tf_destroy(rlc); tf_destroy(bp);

    printf("\n=== RLC Analysis Complete ===\n");
    return 0;
}
