/**
 * @file example_ball_screw.c
 * @brief L6 Canonical Problem — Ball Screw Linear Axis Design
 *
 * Designs a ball screw linear axis for a CNC/precision positioning stage:
 *   1. Define lead screw, motor, and linear load
 *   2. Convert linear requirements to rotary motor specs
 *   3. Compute required torque considering gravity and friction
 *   4. Build state-space model
 *   5. Analyze eigenvalues and stability
 *
 * Build: make examples && ./build/example_ball_screw
 */

#include "mechatronic_definitions.h"
#include "electromechanical_coupling.h"
#include "mechatronic_state_space.h"
#include "actuator_modeling.h"
#include "sensor_modeling.h"
#include "mechatronic_system_design.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

int main(void)
{
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║ Ball Screw Linear Axis — CNC Stage Design    ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    /* Step 1: Define the screw axis ------------------------------------ */
    printf("━━━ Step 1: Ball Screw Specifications ━━━\n");

    /* NSK compact ball screw: 16mm dia, 10mm lead */
    LeadScrew screw = lead_screw_create(
        0.010,    /* pitch = 10 mm/rev */
        0.92,     /* efficiency = 92% (ball screw) */
        0.016,    /* diameter = 16 mm */
        2.5e-5,   /* J_screw = screw inertia [kg·m²] */
        0.003     /* μ = 0.003 rolling friction coefficient */
    );

    printf("Screw: pitch=%.0f mm, η=%.0f%%, D=%.0f mm\n",
           screw.pitch * 1000.0, screw.efficiency * 100.0,
           screw.diameter * 1000.0);
    printf("       J_screw=%.2e kg·m², μ=%.3f\n",
           screw.J_screw, screw.mu_v);

    /* Linear stage: 50 kg payload on horizontal guides */
    double M_payload = 50.0;       /* [kg] */
    double B_guide = 50.0;         /* Guide damping [N·s/m] */
    double accel_target = 5.0;     /* Acceleration [m/s²] (≈0.5G) */
    double v_max = 0.5;            /* Max velocity [m/s] */
    double F_process = 500.0;      /* Cutting force [N] (worst case) */

    printf("Payload: M=%.0f kg, a_max=%.1f m/s², v_max=%.2f m/s\n",
           M_payload, accel_target, v_max);
    printf("         F_cutting=%.0f N\n", F_process);

    /* Step 2: Linear-to-rotary conversion ------------------------------ */
    printf("\n━━━ Step 2: Linear-to-Rotary Conversion ━━━\n");

    double omega_motor, T_motor_linear;
    linear_to_rotary_conversion(v_max, F_process, screw.pitch,
                                 screw.efficiency, &omega_motor, &T_motor_linear);

    double RPM = omega_motor * 30.0 / M_PI;
    printf("Required motor speed: %.1f rad/s (%.0f RPM)\n", omega_motor, RPM);
    printf("Process force at motor: %.3f N·m\n", T_motor_linear);

    /* Step 3: Compute required torque ---------------------------------- */
    printf("\n━━━ Step 3: Torque Requirements ━━━\n");

    /* Horizontal axis (θ_tilt = 0) */
    double T_req = ball_screw_required_torque(
        M_payload, screw.pitch, screw.efficiency,
        screw.mu_v,
        0.0,              /* θ_tilt = 0 (horizontal) */
        accel_target,     /* acceleration */
        screw.J_screw,
        F_process);

    /* Inclined axis (θ_tilt = 30°) */
    double T_req_inclined = ball_screw_required_torque(
        M_payload, screw.pitch, screw.efficiency,
        screw.mu_v,
        M_PI / 6.0,       /* θ_tilt = 30° */
        accel_target,
        screw.J_screw,
        F_process);

    printf("Horizontal axis → T_req = %.3f N·m\n", T_req);
    printf("  30° inclined   → T_req = %.3f N·m\n", T_req_inclined);

    /* Gravity component breakdown */
    double g = 9.80665;
    double r = screw.pitch / (2.0 * M_PI);
    double F_gravity_30 = M_payload * g * sin(M_PI / 6.0);
    double T_gravity_30 = F_gravity_30 * r / screw.efficiency;
    printf("  Gravity component (30°): %.3f N·m\n", T_gravity_30);

    /* Step 4: Select motor and build state-space ----------------------- */
    printf("\n━━━ Step 4: Motor Selection & State-Space ━━━\n");

    /* Yaskawa SGM7J-02A — 200W AC servo */
    DCMotorParams motor = dc_motor_create(
        2.5,       /* R_a [Ω] */
        0.004,     /* L_a [H] */
        0.0637,    /* K_e [V·s/rad] — 200V/3000RPM */
        0.26e-4,   /* J_m [kg·m²] */
        5e-5,      /* B_m [N·m·s/rad] */
        200.0,     /* V_max [V] */
        2.0        /* I_max [A] */
    );

    printf("Motor: %.0f W servo, K_e=%.4f V·s/rad, K_t=%.4f N·m/A\n",
           motor.P_rated, motor.K_e, motor.K_t);
    printf("       Rated torque: %.3f N·m, Max speed: %.0f RPM\n",
           motor.max_torque, motor.max_speed * 30.0 / M_PI);

    /* Torque margin check */
    double margin = (motor.max_torque - T_req) / motor.max_torque * 100.0;
    printf("       Torque margin: %.1f%%  →  %s\n",
           margin, margin > 20.0 ? "GOOD" : "MARGINAL");

    /* Build state-space model */
    LinearStateSpace ss = ss_ball_screw_axis(
        &motor, &screw, M_payload, B_guide);

    printf("\nBall screw axis state-space: %d states × %d inputs × %d outputs\n",
           ss.n_states, ss.n_inputs, ss.n_outputs);

    /* Step 5: Eigenvalue analysis -------------------------------------- */
    printf("\n━━━ Step 5: Dynamic Analysis ━━━\n");

    double eval_re[SS_MAX_DIM], eval_im[SS_MAX_DIM];
    int ret = ss_eigenvalues(ss.A, ss.n_states, eval_re, eval_im);
    if (ret == 0) {
        printf("Eigenvalues:\n");
        for (int i = 0; i < ss.n_states; i++) {
            printf("  λ_%d = %.2f", i + 1, eval_re[i]);
            if (fabs(eval_im[i]) > 1e-6)
                printf(" ± j%.2f", fabs(eval_im[i]));
            double wn, z;
            ss_mode_properties(eval_re[i], eval_im[i], &wn, &z);
            printf("  →  ω_n=%.1f rad/s (%.1f Hz), ζ=%.3f",
                   wn, wn / (2.0 * M_PI), z);
            printf("\n");
        }

        int stable = ss_is_stable(eval_re, eval_im, ss.n_states);
        printf("System is %s\n", stable ? "STABLE ✓" : "UNSTABLE ✗");
    }

    /* Step 6: Controllability & Observability -------------------------- */
    printf("\n━━━ Step 6: Structural Properties ━━━\n");

    int ctrl = ss_is_controllable(&ss);
    int obs  = ss_is_observable(&ss);
    printf("Controllable: %s\n", ctrl ? "YES ✓" : "NO ✗");
    printf("Observable:   %s\n", obs ? "YES ✓" : "NO ✗");

    /* Step 7: Cascade controller for this axis ------------------------- */
    printf("\n━━━ Step 7: Cascade Controller ━━━\n");

    RotarySensor enc = rotary_sensor_create(
        SENSOR_ENCODER_INCREMENTAL, 2048, 100e3, 1e-6);

    ServoAxis axis;
    /* Direct drive (no gear) */
    GearTrain direct = gear_train_create(1.0, 1.0, 0.0, 0.0, 0.0, 1e6);
    MechatronicLoad mech_load = mechatronic_load_create(
        0.0, M_payload, B_guide, 5.0, 7.0);

    design_servo_axis(&motor, &direct, &mech_load, &enc,
                       20.0,  /* 20 Hz position bandwidth */
                       5.0, 10.0, &axis);

    printf("Position loop: %.1f Hz (Kp=%.2f 1/s)\n",
           axis.bandwidth_pos, axis.Kp_position);
    printf("Velocity loop: %.1f Hz (Kp=%.4f, Ki=%.3f)\n",
           axis.bandwidth_vel, axis.Kp_velocity, axis.Ki_velocity);
    printf("Current loop:  %.1f Hz (Kp=%.3f, Ki=%.1f)\n",
           axis.bandwidth_cur, axis.Kp_current, axis.Ki_current);

    /* Step 8: Position resolution -------------------------------------- */
    printf("\n━━━ Step 8: Positioning Accuracy ━━━\n");

    double linear_res = enc.resolution * screw.pitch / (2.0 * M_PI);
    printf("Encoder: %d PPR, quadrature -> %d counts/rev\n",
           enc.ppr, enc.ppr * 4);
    printf("Angular resolution: %.2e rad\n", enc.resolution);
    printf("Linear resolution:  %.2f um\n", linear_res * 1e6);

    printf("\n╔══════════════════════════════════════════════╗\n");
    printf("║     Ball Screw Axis Design Complete           ║\n");
    printf("╚══════════════════════════════════════════════╝\n");

    return 0;
}
