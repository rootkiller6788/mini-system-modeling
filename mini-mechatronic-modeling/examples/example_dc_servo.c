/**
 * @file example_dc_servo.c
 * @brief L6 Canonical Problem — DC Servo Motor Design
 *
 * Designs a complete DC servo axis for an industrial positioning system:
 *   1. Define motor, gear, load parameters
 *   2. Reflect dynamics to motor shaft
 *   3. Design cascaded PI-P-P controller
 *   4. Analyze step response performance
 *   5. Simulate a motion profile
 *
 * This demonstrates the end-to-end mechatronic design workflow.
 *
 * Build: make examples && ./build/example_dc_servo
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
    printf("║   DC Servo Motor Design — Positioning Axis   ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    /* Step 1: Define components ---------------------------------------- */
    printf("━━━ Step 1: Component Parameters ━━━\n");

    /* Maxon RE 30 (60W) — representative small servo motor */
    DCMotorParams motor = dc_motor_create(
        0.85,      /* R_a  [Ω]  — armature resistance */
        0.00034,   /* L_a  [H]  — armature inductance */
        0.0255,    /* K_e  [V·s/rad] — back-EMF constant */
        3.38e-6,   /* J_m  [kg·m²] — rotor inertia */
        1.5e-6,    /* B_m  [N·m·s/rad] — viscous damping */
        24.0,      /* V_max [V] */
        3.0        /* I_max [A] */
    );

    printf("Motor: R_a=%.3f Ω, L_a=%.2f mH, K_e=%.4f V·s/rad\n",
           motor.R_a, motor.L_a * 1000.0, motor.K_e);
    printf("       J_m=%.2e kg·m², B_m=%.2e N·m·s/rad\n",
           motor.J_m, motor.B_m);
    printf("       Rated: %.0f W, max %.0f RPM\n",
           motor.P_rated, motor.max_speed * 30.0 / M_PI);

    printf("       τ_e = %.4f ms (electrical)\n",
           motor.tau_e * 1000.0);
    printf("       τ_m = %.2f s (mechanical, motor only)\n",
           motor.tau_m);

    /* Harmonic Drive gearhead, 50:1 */
    GearTrain gear = gear_train_create(
        50.0,   /* ratio N */
        0.85,   /* efficiency η */
        1e-6,   /* J_input [kg·m²] */
        8e-6,   /* J_output [kg·m²] */
        0.0005, /* backlash [rad] ≈ 1.7 arcmin */
        2e4     /* K_stiffness [N·m/rad] */
    );

    printf("\nGear: N=%0.f:1, η=%.0f%%, backlash=%.6f rad\n",
           gear.ratio, gear.efficiency * 100.0, gear.backlash);

    /* Load: a small rotary table */
    MechatronicLoad load = mechatronic_load_create(
        0.005,  /* J_load [kg·m²] */
        2.0,    /* M_load [kg] */
        0.01,   /* B_load [N·m·s/rad] — viscous */
        0.1,    /* T_coulomb [N·m] */
        0.15    /* T_static [N·m] */
    );

    printf("Load: J=%.4f kg·m², T_coulomb=%.3f N·m, T_static=%.3f N·m\n",
           load.J_load, load.T_coulomb, load.T_static);

    /* Encoder: 2000 PPR incremental with 4x quadrature */
    RotarySensor sensor = rotary_sensor_create(
        SENSOR_ENCODER_INCREMENTAL, 2000, 100e3, 1e-6);

    printf("Sensor: %d PPR encoder, res=%.2e rad, BW=%.0f kHz\n",
           sensor.ppr, sensor.resolution, sensor.bandwidth / 1000.0);

    /* Step 2: Reflect dynamics ----------------------------------------- */
    printf("\n━━━ Step 2: Reflected Dynamics ━━━\n");

    double J_total, B_total;
    compute_reflected_dynamics(&motor, &gear,
                                load.J_load, load.B_load,
                                &J_total, &B_total);

    printf("Motor shaft: J_total=%.6f kg·m², B_total=%.6f N·m·s/rad\n",
           J_total, B_total);

    /* Inertia ratio check */
    double J_load_reflected = load.J_load / (gear.ratio * gear.ratio);
    double lambda = inertia_ratio(J_load_reflected, motor.J_m);
    printf("Inertia ratio λ = J_load_reflected/J_motor = %.2f\n", lambda);
    printf("Classification: %s\n", inertia_ratio_classification(lambda));

    /* Electromechanical time constant */
    double tau_em = dc_motor_tau_electromechanical(&motor, J_total, B_total);
    printf("Electromechanical τ_em = %.4f s\n", tau_em);

    /* Step 3: Cascade controller design -------------------------------- */
    printf("\n━━━ Step 3: Cascade Controller Design ━━━\n");

    ServoAxis axis;
    design_servo_axis(&motor, &gear, &load, &sensor,
                       15.0,  /* ω_pos = 15 Hz position bandwidth */
                       5.0,   /* ω_vel/ω_pos = 5 */
                       10.0,  /* ω_cur/ω_vel = 10 */
                       &axis);

    printf("Cascaded bandwidths:\n");
    printf("  Position loop: %.1f Hz (Kp=%.3f 1/s)\n",
           axis.bandwidth_pos, axis.Kp_position);
    printf("  Velocity loop: %.1f Hz (Kp=%.4f, Ki=%.4f)\n",
           axis.bandwidth_vel, axis.Kp_velocity, axis.Ki_velocity);
    printf("  Current loop:  %.1f Hz (Kp=%.3f, Ki=%.1f)\n",
           axis.bandwidth_cur, axis.Kp_current, axis.Ki_current);

    /* Step 4: Step response performance -------------------------------- */
    printf("\n━━━ Step 4: Step Response Analysis ━━━\n");

    double rise, settle, overshoot, ss_err;
    servo_axis_step_analysis(&axis, 0.5, 0.001, 1.0,
                              &rise, &settle, &overshoot, &ss_err);
    printf("Predicted step response:\n");
    printf("  Rise time (10-90%%): %.2f ms\n", rise * 1000.0);
    printf("  Settling time (2%%): %.2f ms\n", settle * 1000.0);
    printf("  Overshoot: %.1f%%\n", overshoot);
    printf("  Steady-state error: %.4f rad\n", ss_err);

    /* Step 5: Motion profile simulation -------------------------------- */
    printf("\n━━━ Step 5: Motion Profile ━━━\n");

    DutyCycle cycle;
    memset(&cycle, 0, sizeof(cycle));
    cycle.n_segments = 5;
    cycle.cycle_time = 2.0;

    /* Segment 1: Accelerate */
    cycle.segments[0].torque = 0.3;
    cycle.segments[0].speed = 50.0;
    cycle.segments[0].duration = 0.2;

    /* Segment 2: Constant speed */
    cycle.segments[1].torque = 0.12;
    cycle.segments[1].speed = 100.0;
    cycle.segments[1].duration = 0.5;

    /* Segment 3: Decelerate */
    cycle.segments[2].torque = -0.15;
    cycle.segments[2].speed = 50.0;
    cycle.segments[2].duration = 0.2;

    /* Segment 4: Dwell (positioning) */
    cycle.segments[3].torque = 0.0;
    cycle.segments[3].speed = 0.0;
    cycle.segments[3].duration = 0.3;

    /* Segment 5: Return (accelerate reverse) */
    cycle.segments[4].torque = -0.25;
    cycle.segments[4].speed = -80.0;
    cycle.segments[4].duration = 0.8;

    double T_rms = rms_torque(&cycle);
    double T_peak = peak_torque_requirement(&cycle);
    double omega_avg = average_speed(&cycle);

    printf("RMS torque: %.4f N·m\n", T_rms);
    printf("Peak torque: %.4f N·m\n", T_peak);
    printf("Average speed: %.1f rad/s (%.0f RPM)\n",
           omega_avg, omega_avg * 30.0 / M_PI);

    /* Check against motor capability */
    printf("\nFeasibility check:\n");
    printf("  T_rms (%.3f) vs T_rated (%.3f) → %s\n",
           T_rms, motor.max_torque,
           T_rms <= motor.max_torque ? "OK" : "FAIL");
    printf("  T_peak (%.3f) vs T_stall (%.3f) → %s\n",
           T_peak, motor.max_torque,
           T_peak <= motor.max_torque ? "OK" : "FAIL");

    /* Energy analysis */
    double energy, peak_p, avg_p;
    motion_profile_energy(&cycle, &motor, &gear, J_total,
                           &energy, &peak_p, &avg_p);
    printf("\nEnergy per cycle: %.2f J\n", energy);
    printf("Peak electrical power: %.2f W\n", peak_p);
    printf("Average power: %.2f W\n", avg_p);

    printf("\n╔══════════════════════════════════════════════╗\n");
    printf("║         Servo Axis Design Complete           ║\n");
    printf("╚══════════════════════════════════════════════╝\n");

    return 0;
}
