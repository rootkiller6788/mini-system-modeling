/**
 * @example example_quarter_car.c
 * @brief Quarter-car suspension design and analysis
 *
 * Demonstrates: Ride comfort vs handling tradeoff, suspension tuning,
 * natural frequency target selection, transmissibility at wheel hop.
 * Reference: Gillespie, Fundamentals of Vehicle Dynamics (1992).
 * Applications: Detroit/Tesla automotive suspension engineering.
 */

#include <stdio.h>
#include <math.h>
#include "mechanical_applications.h"
#include "mechanical_elements.h"

int main(void)
{
    printf("\n=== Quarter-Car Suspension Design ===\n\n");

    /* Tesla Model 3-like parameters (rear corner, estimated) */
    double sprung_mass = 380.0;      /* kg per rear corner */
    double unsprung_mass = 45.0;     /* kg (wheel + tire + brake + half suspension) */
    double tire_stiffness = 220000.0; /* N/m (radial) */
    double tire_damping = 80.0;       /* N*s/m */

    printf("Vehicle Parameters:\n");
    printf("  Sprung mass (per corner)  = %.1f kg\n", sprung_mass);
    printf("  Unsprung mass             = %.1f kg\n", unsprung_mass);
    printf("  Tire radial stiffness     = %.0f N/m\n", tire_stiffness);
    printf("\n");

    /* Design targets */
    printf("Design Targets:\n");
    printf("  Body bounce frequency target: 1.2 Hz (ride comfort)\n");
    printf("  Damping ratio target: 0.3 (balanced ride/handling)\n");
    printf("\n");

    /* Compute suspension parameters */
    double ks = suspension_stiffness_for_ride(sprung_mass, 1.2);
    double bs = suspension_damping_for_ratio(sprung_mass, ks, 0.3);

    printf("Designed Suspension Parameters:\n");
    printf("  Spring stiffness ks = %.0f N/m (%.1f N/mm)\n", ks, ks/1000.0);
    printf("  Damping coefficient bs = %.0f N*s/m\n\n", bs);

    /* Analyze the design */
    QuarterCarModel qcm = {sprung_mass, unsprung_mass, ks, bs, tire_stiffness, tire_damping};

    double body_hz, wheel_hz;
    quarter_car_natural_freqs(&qcm, &body_hz, &wheel_hz);

    printf("System Natural Frequencies:\n");
    printf("  Body bounce  = %.2f Hz\n", body_hz);
    printf("  Wheel hop    = %.2f Hz\n", wheel_hz);

    double defl = quarter_car_static_deflection(&qcm, 9.81);
    printf("  Static deflection = %.1f mm\n\n", defl*1000.0);

    double zeta_actual = quarter_car_damping_ratio(&qcm);
    printf("  Actual damping ratio = %.3f\n", zeta_actual);

    /* Ride comfort analysis */
    printf("\nRide Comfort Analysis (ISO 8608 road):\n");
    double speeds[] = {30.0, 60.0, 90.0, 120.0}; /* km/h */
    for (int i = 0; i < 4; i++) {
        double v_ms = speeds[i] / 3.6;
        double accel = quarter_car_body_accel_rms(&qcm, 5e-6, v_ms);
        double travel = quarter_car_suspension_travel_rms(&qcm, 5e-6, v_ms) * 1000.0;
        double dload = quarter_car_dynamic_load_rms(&qcm, 5e-6, v_ms);
        double comfort = quarter_car_ride_comfort(&qcm, v_ms, 3600.0);

        printf("  %.0f km/h: accel=%.3f m/s^2, travel=%.1f mm, DFL=%.0f N, comfort=%.3f\n",
               speeds[i], accel, travel, dload, comfort);
    }

    /* Skyhook comparison */
    double c_sky = skyhook_damping_ideal(&qcm);
    printf("\n  Skyhook ideal damping = %.0f N*s/m (actual = %.0f)\n", c_sky, bs);

    /* Isolation at wheel hop */
    double wn_body = body_hz * 2.0 * M_PI;
    double wn_wheel = wheel_hz * 2.0 * M_PI;
    double r_wheel = wn_wheel / wn_body;
    double TR_wheel = force_transmissibility(r_wheel, zeta_actual);
    printf("  Transmissibility at wheel hop (r=%.2f) = %.4f\n", r_wheel, TR_wheel);

    /* Half-car pitch estimate */
    double pitch_inertia = sprung_mass * 2.0 * 2.5 * 2.5; /* approximate: m * L^2 */
    double pitch_freq = half_car_pitch_freq(sprung_mass*2, pitch_inertia, 2.8, ks);
    printf("\nHalf-Car Pitch Estimate:\n");
    printf("  Pitch natural frequency = %.2f Hz\n", pitch_freq);

    /* Alternative designs for comparison */
    printf("\nAlternative Designs Comparison:\n");
    printf("  Target fn [Hz]  ks [N/m]  bs [N*s/m]  defl [mm]\n");
    double target_fns[] = {0.8, 1.0, 1.2, 1.5, 2.0};
    for (int i = 0; i < 5; i++) {
        double k_alt = suspension_stiffness_for_ride(sprung_mass, target_fns[i]);
        double b_alt = suspension_damping_for_ratio(sprung_mass, k_alt, 0.3);
        defl = sprung_mass * 9.81 / k_alt;
        printf("  %.1f           %.0f     %.0f         %.1f\n",
               target_fns[i], k_alt, b_alt, defl*1000.0);
    }

    printf("\n=== Suspension Design Complete ===\n");
    printf("\nNote: Typical sedan targets: body 1.0-1.5 Hz, damping 0.2-0.4,\n");
    printf("       deflection 100-200 mm, wheel hop 10-15 Hz.\n");
    printf("       Tesla and Detroit OEM applications.\n");
    return 0;
}
