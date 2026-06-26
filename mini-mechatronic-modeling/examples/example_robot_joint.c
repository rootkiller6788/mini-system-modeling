/**
 * @file example_robot_joint.c
 * @brief L6 Canonical Problem — Industrial Robot Joint Design
 *
 * Designs the shoulder joint of a 6-DOF industrial robot (similar to
 * KUKA KR 6 or ABB IRB 120):
 *   1. Compute link inertia and gravity loads
 *   2. Select motor + harmonic drive gear
 *   3. Verify torque and speed requirements
 *   4. Check encoder resolution for endpoint accuracy
 *   5. Analyze flexible gear dynamics
 *
 * Build: make examples && ./build/example_robot_joint
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
    printf("║  Industrial Robot Joint — Shoulder Design    ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    /* Step 1: Link parameters ------------------------------------------ */
    printf("━━━ Step 1: Link & Kinematics ━━━\n");

    /* Link parameters (representative of a 6 kg payload robot) */
    double M_link = 12.0;       /* Link + forearm + wrist + payload [kg] */
    double L_cm  = 0.45;        /* Center of mass distance from joint [m] */
    double J_link = M_link * L_cm * L_cm / 3.0;  /* Rod inertia ≈ ML²/3 */

    printf("Link mass: %.1f kg, CM distance: %.2f m\n", M_link, L_cm);
    printf("Link inertia about joint: J_link = %.4f kg·m²\n", J_link);

    /* Target motion profile */
    double accel_max = 15.0;    /* Max angular acceleration [rad/s²] */
    double omega_max = 3.0;     /* Max angular velocity [rad/s] (≈0.5 rev/s) */
    double reach = 0.8;         /* Robot reach [m] */

    printf("Motion: α_max=%.1f rad/s², ω_max=%.1f rad/s\n",
           accel_max, omega_max);
    printf("Reach: %.1f m, payload: 6 kg\n", reach);

    /* Step 2: Gear selection ------------------------------------------- */
    printf("\n━━━ Step 2: Transmission Selection ━━━\n");

    /* Harmonic Drive CSD-25-100: ratio 100:1, high precision */
    GearTrain harmonic_drive = gear_train_create(
        100.0,    /* ratio N = 100:1 */
        0.80,     /* efficiency (typical for harmonic drive at rated torque) */
        6.7e-6,   /* J_input [kg·m²] (wave generator) */
        1.3e-4,   /* J_output [kg·m²] (flexspline + output) */
        0.0002,   /* backlash [rad] ≈ 0.7 arcmin (near-zero for harmonic) */
        1.5e4     /* K_stiffness [N·m/rad] */
    );

    printf("Harmonic Drive 100:1, η=%.0f%%, backlash=%.4f rad\n",
           harmonic_drive.efficiency * 100.0, harmonic_drive.backlash);

    /* Inertia matching */
    double J_load_reflected = J_link
                              / (harmonic_drive.ratio * harmonic_drive.ratio);
    double inertia_match = inertia_ratio(J_load_reflected, 6.7e-6);
    printf("Reflected link inertia: %.2e kg·m²\n", J_load_reflected);
    printf("Inertia ratio: %.1f → %s\n",
           inertia_match, inertia_ratio_classification(inertia_match));

    /* Optimal gear ratio check */
    double N_opt = optimal_gear_ratio_inertia_match(6.7e-6, J_link);
    printf("Optimal N from inertia matching: %.0f:1\n", N_opt);

    /* Step 3: Motor selection ------------------------------------------ */
    printf("\n━━━ Step 3: Motor Selection ━━━\n");

    /* Kollmorgen AKM24C — 0.8 N·m rated, 3000 RPM */
    DCMotorParams motor = dc_motor_create(
        8.5,       /* R_a [Ω] */
        0.012,     /* L_a [H] */
        0.114,     /* K_e=K_t [V·s/rad, N·m/A] — 170V/1500RPM */
        2.4e-5,    /* J_m [kg·m²] */
        1e-4,      /* B_m [N·m·s/rad] */
        320.0,     /* V_dc [V] */
        8.0        /* I_max [A] (peak) */
    );
    /* Note: this motor is higher power — for shoulder joint */
    motor.max_torque = motor.K_t * motor.max_current;  /* peak torque */
    motor.max_speed = motor.max_voltage / motor.K_e;
    motor.P_rated = motor.max_voltage * motor.max_current * 0.7;

    printf("Motor: K_t=%.3f N·m/A, R_a=%.1f Ω\n",
           motor.K_t, motor.R_a);
    printf("       Rated torque: %.2f N·m (cont), Peak: %.2f N·m\n",
           motor.K_t * 4.0, motor.max_torque);
    printf("       No-load speed: %.0f RPM\n",
           motor.max_speed * 30.0 / M_PI);

    /* Step 4: Torque requirements -------------------------------------- */
    printf("\n━━━ Step 4: Joint Torque Analysis ━━━\n");

    double T_friction_joint = 2.0;  /* N·m at joint (seal + bearing friction) */
    double T_req = robot_joint_max_torque(
        M_link, L_cm, J_link, harmonic_drive.ratio,
        accel_max, motor.J_m, T_friction_joint);

    printf("Required peak motor torque: %.3f N·m\n", T_req);

    /* Breakdown */
    double T_inertial = (motor.J_m + J_link / 10000.0) * 100.0 * accel_max;
    double T_gravity = M_link * 9.80665 * L_cm / 100.0;
    printf("  Inertial component:   %.3f N·m\n", T_inertial);
    printf("  Gravity component:    %.3f N·m (worst case: horizontal)\n",
           T_gravity);
    printf("  Friction component:   %.3f N·m\n",
           T_friction_joint / harmonic_drive.ratio);

    /* Safety factor */
    double SF = motor.max_torque / T_req;
    printf("  Safety factor: %.2f → %s\n",
           SF, SF >= 1.5 ? "ACCEPTABLE" : "INSUFFICIENT");

    /* Step 5: Flexible joint dynamics ---------------------------------- */
    printf("\n━━━ Step 5: Flexible Joint Dynamics ━━━\n");

    LinearStateSpace ss_joint = ss_motor_flexible_gear(
        &motor, motor.J_m, motor.B_m,
        &harmonic_drive,
        J_link, 0.1);  /* J_L, B_L */

    printf("Flexible joint state-space: %d states\n", ss_joint.n_states);
    printf("States: [I_a, ω_motor, ω_link, Δθ_windup]^T\n");

    double eval_re[SS_MAX_DIM], eval_im[SS_MAX_DIM];
    ss_eigenvalues(ss_joint.A, ss_joint.n_states, eval_re, eval_im);

    printf("Eigenvalues:\n");
    for (int i = 0; i < ss_joint.n_states; i++) {
        printf("  λ_%d = %.2f", i + 1, eval_re[i]);
        if (fabs(eval_im[i]) > 1e-6) {
            printf(" ± j%.2f", fabs(eval_im[i]));
            double wn, z;
            ss_mode_properties(eval_re[i], eval_im[i], &wn, &z);
            printf("  →  f_n=%.1f Hz, ζ=%.4f", wn / (2.0 * M_PI), z);

            /* Identify flexible mode (typically high freq, low damping) */
            if (wn / (2.0 * M_PI) > 50.0 && z < 0.1)
                printf("  ← FLEXIBLE MODE");
        }
        printf("\n");
    }

    /* Step 6: Endpoint accuracy ---------------------------------------- */
    printf("\n━━━ Step 6: Endpoint Accuracy ━━━\n");

    /* Motor-side encoder: 20-bit absolute (1,048,576 counts/rev) */
    RotarySensor motor_encoder = rotary_sensor_create(
        SENSOR_ENCODER_ABSOLUTE, 1 << 20, 1000.0, 1e-6);

    double motor_res = motor_encoder.resolution; /* rad at motor */
    double joint_res = motor_res / harmonic_drive.ratio; /* rad at joint */
    double endpoint_res = joint_res * reach; /* linear at endpoint */

    printf("Motor encoder: %d counts/rev → %.2e rad resolution\n",
           motor_encoder.ppr, motor_res);
    printf("Joint resolution (after 100:1 gear): %.2e rad\n", joint_res);
    printf("Endpoint resolution at %.1f m reach: %.1f μm\n",
           reach, endpoint_res * 1e6);

    /* Step 7: Thermal considerations ----------------------------------- */
    printf("\n━━━ Step 7: Thermal Analysis ━━━\n");

    double T_amb = 25.0;          /* Ambient [°C] */
    double R_th = 1.2;            /* Thermal resistance [K/W] */
    double I_cont = 4.0;          /* Continuous current for RMS torque */
    double T_ss = dc_motor_steady_state_temp(I_cont, motor.R_a, T_amb, R_th);

    printf("Continuous current: %.1f A → I²R = %.1f W copper loss\n",
           I_cont, I_cont * I_cont * motor.R_a);
    printf("Steady-state winding temperature: %.1f °C\n", T_ss);

    /* Resistance increase at operating temp */
    double R_hot = resistance_temp_correction(motor.R_a, 25.0, T_ss, 0.00393);
    printf("Resistance at %.1f°C: %.3f Ω (cold: %.3f Ω, +%.1f%%)\n",
           T_ss, R_hot, motor.R_a, (R_hot/motor.R_a - 1.0)*100.0);

    printf("\n╔══════════════════════════════════════════════╗\n");
    printf("║      Robot Joint Design Complete             ║\n");
    printf("╚══════════════════════════════════════════════╝\n");

    return 0;
}
