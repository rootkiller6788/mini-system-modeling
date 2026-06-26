/**
 * Example 4: Quarter-Car Suspension Analysis
 *
 * Demonstrates: L6 Canonical System, L7 Application
 *
 * Models a vehicle suspension system (quarter-car model) used in
 * Tesla, Toyota, and other automotive designs (ISO 2631).
 */

#include "transfer_function.h"
#include "tf_analysis.h"
#include "tf_advanced.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Quarter-Car Suspension Analysis ===\n\n");

    /* Tesla Model S-like parameters:
     * Ms=450 kg (sprung mass per corner)
     * Mus=50 kg (unsprung mass: wheel, tire, brake)
     * Ks=25000 N/m (suspension spring)
     * Cs=2000 N*s/m (damper)
     * Kt=200000 N/m (tire vertical stiffness)
     */
    double Ms=450.0, Mus=50.0, Ks=25000.0, Cs=2000.0, Kt=200000.0;

    TransferFunction *susp = tf_quarter_car_suspension(Ms, Mus, Ks, Cs, Kt);

    printf("Suspension TF: Xs(s)/Xr(s) - sprung mass to road\n");
    printf("  Tesla Model S parameters:\n");
    printf("  Ms=%.0f kg, Mus=%.0f kg, Ks=%.0f N/m, Cs=%.0f N*s/m, Kt=%.0f N/m\n",
           Ms, Mus, Ks, Cs, Kt);
    printf("\nTransfer Function:\n");
    tf_print(susp);
    printf("\n");

    /* Properties */
    printf("System Properties:\n");
    printf("  Order: %d poles\n", tf_num_poles(susp));
    printf("  Properness: %s\n",
           tf_properness(susp) == TF_STRICTLY_PROPER ? "strictly proper" :
           tf_properness(susp) == TF_PROPER ? "proper" : "improper");
    printf("  DC Gain (static): %.6f (long wavelength road)\n", tf_dc_gain(susp));
    printf("\n");

    /* Stability */
    tf_stability_t st = tf_stability_routh(susp);
    printf("Stability: %s\n",
           st == TF_STABLE ? "STABLE (passive suspension)" : "UNSTABLE");

    /* Frequency response (ride comfort analysis) */
    printf("\nFrequency Response (ride comfort):\n");
    printf("  %-12s %-12s %-12s\n", "f (Hz)", "Gain", "Gain (dB)");

    double test_freqs[] = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0};
    int i;
    for (i = 0; i < 8; i++) {
        double w = test_freqs[i] * 2.0 * M_PI;
        double mag = tf_magnitude_at(susp, w);
        double mag_db = 20.0 * log10(mag);
        printf("  %-12.2f %-12.4f %-12.2f\n", test_freqs[i], mag, mag_db);
    }

    /* Body resonance (sprung mass) and wheel hop (unsprung mass) */
    double w_body = sqrt(Ks/Ms);  /* approx */
    double w_wheel = sqrt((Ks+Kt)/Mus);  /* approx */
    printf("\nApproximate Natural Frequencies:\n");
    printf("  Body (sprung mass) resonance: %.2f rad/s (%.2f Hz)\n",
           w_body, w_body/(2.0*M_PI));
    printf("  Wheel hop (unsprung mass): %.2f rad/s (%.2f Hz)\n",
           w_wheel, w_wheel/(2.0*M_PI));

    /* ISO 2631 comfort: human sensitivity peaks at 4-8 Hz */
    double comfort_w = 5.0 * 2.0 * M_PI;
    double comfort_gain = tf_magnitude_at(susp, comfort_w);
    printf("\nISO 2631 Ride Comfort:\n");
    printf("  Gain at 5 Hz (human sensitivity peak): %.4f (%.2f dB)\n",
           comfort_gain, 20.0*log10(comfort_gain));

    tf_destroy(susp);

    printf("\n=== Suspension Analysis Complete ===\n");
    return 0;
}
