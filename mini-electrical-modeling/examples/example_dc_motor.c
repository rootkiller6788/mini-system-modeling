/**
 * @file example_dc_motor.c
 * @brief End-to-end example: DC motor modeling and analysis
 *
 * Demonstrates: DC motor steady-state analysis, speed-torque curve,
 * transfer function derivation, state-space simulation, PWM analysis.
 *
 * Problem: Model a small 12V DC motor (Maxon RE 25 style) and analyze
 * its behavior under various operating conditions.
 *
 * Reference: MIT 6.302 — DC motor as canonical control plant
 */

#include <stdio.h>
#include <math.h>
#include "electromechanical_systems.h"
#include "electrical_elements.h"
#include "electrical_state_space.h"

int main(void)
{
    printf("=== DC Motor Modeling Example ===\n\n");

    /* Define motor parameters (typical small brushed DC motor) */
    DCMotorParams motor;
    motor.armature_resistance  = 2.5;     /* ohm */
    motor.armature_inductance  = 0.0015;  /* H */
    motor.back_emf_constant    = 0.023;   /* V/(rad/s) = 0.023 V*s/rad */
    motor.torque_constant      = 0.023;   /* Nm/A (same as Kb in SI) */
    motor.rotor_inertia        = 1.2e-6;  /* kg*m^2 */
    motor.viscous_friction     = 5.0e-7;  /* Nm/(rad/s) */
    motor.rated_voltage        = 12.0;
    motor.rated_current        = 1.5;
    motor.rated_speed          = 500.0;    /* rad/s (~4800 rpm) */
    motor.rated_torque          = 0.03;   /* Nm */

    printf("Motor Parameters:\n");
    printf("  Ra = %.2f ohm, La = %.2f mH\n",
           motor.armature_resistance, motor.armature_inductance*1e3);
    printf("  Kt = Kb = %.3f Nm/A = V/(rad/s)\n", motor.torque_constant);
    printf("  J = %.2e kg*m^2, B = %.2e Nm/(rad/s)\n",
           motor.rotor_inertia, motor.viscous_friction);

    /* Steady-state analysis */
    printf("\n=== Steady-State Analysis ===\n");

    double V_nom = 12.0;
    double T_l = 0.01; /* 10 mNm load */

    double w_ss = dc_motor_steady_state_speed(&motor, V_nom, T_l);
    double i_ss = dc_motor_steady_state_current(&motor, V_nom, w_ss);
    double T_ss = dc_motor_steady_state_torque(&motor, V_nom, w_ss);

    printf("At Va=%.0fV, TL=%.3f Nm:\n", V_nom, T_l);
    printf("  Speed:  %.1f rad/s (%.0f rpm)\n", w_ss, w_ss*60/(2*M_PI));
    printf("  Current: %.3f A\n", i_ss);
    printf("  Torque:  %.4f Nm\n", T_ss);

    /* Operating point with efficiency */
    DCMotorOperatingPoint op = dc_motor_operating_point(&motor, V_nom, T_l);
    printf("  Pin: %.2f W, Pout: %.2f W, Eff: %.1f%%\n",
           op.power_in, op.power_out, op.efficiency*100.0);

    /* Speed-torque characteristic */
    printf("\n=== Speed-Torque Curve ===\n");
    printf("  Torque [Nm]   Speed [rad/s]   Speed [rpm]\n");
    printf("  -----------   -------------   -----------\n");
    double T_stall = dc_motor_stall_torque(&motor, V_nom);
    double w_nl = dc_motor_no_load_speed(&motor, V_nom);

    for (int i = 0; i <= 10; i++) {
        double tq = T_stall * (double)i / 10.0;
        double sp = dc_motor_steady_state_speed(&motor, V_nom, tq);
        printf("  %11.4f   %13.1f   %11.0f\n", tq, sp, sp*60/(2*M_PI));
    }

    printf("  Stall torque: %.4f Nm\n", T_stall);
    printf("  No-load speed: %.1f rad/s (%.0f rpm)\n",
           w_nl, w_nl*60/(2*M_PI));

    /* Time constants */
    printf("\n=== Time Constants ===\n");
    double tau_e = dc_motor_electrical_time_constant(&motor);
    double tau_m = dc_motor_mechanical_time_constant(&motor);
    printf("  Electrical: %.3f ms\n", tau_e * 1e3);
    printf("  Mechanical: %.1f ms\n", tau_m * 1e3);
    printf("  Dominant: mechanical (tau_m/tau_e = %.0f:1)\n", tau_m/tau_e);

    /* Transfer function */
    printf("\n=== Transfer Function omega(s)/Va(s) ===\n");
    double num[3], den[3];
    int order = dc_motor_transfer_function_voltage_to_speed(&motor, num, den);
    if (order == 2) {
        printf("  H(s) = %.4g / (%.3e*s^2 + %.3e*s + %.4g)\n",
               num[0], den[2], den[1], den[0]);
        double wn = sqrt(den[0]/den[2]);
        double zeta = den[1]/(2*den[2]*wn);
        printf("  wn = %.1f rad/s, zeta = %.3f\n", wn, zeta);
    } else {
        printf("  H(s) = %.4g / (%.3e*s + %.4g) [1st-order approx]\n",
               num[0], den[1], den[0]);
    }

    /* State-space simulation: startup transient */
    printf("\n=== Startup Transient (RK4 Simulation) ===\n");
    printf("  Time [ms]   Current [A]   Speed [rad/s]\n");
    printf("  ---------   -----------   -------------\n");

    StateSpace ss;
    ss.n = 2; ss.m = 1; ss.p = 1;
    double A_data[4] = {
        -motor.armature_resistance/motor.armature_inductance,
        -motor.back_emf_constant/motor.armature_inductance,
        motor.torque_constant/motor.rotor_inertia,
        -motor.viscous_friction/motor.rotor_inertia
    };
    double B_data[2] = {1.0/motor.armature_inductance, 0.0};
    double C_data[2] = {0.0, 1.0};
    double D_data[1] = {0.0};

    for (int i = 0; i < 4; i++) ss.A[i] = A_data[i];
    ss.B[0] = B_data[0]; ss.B[1] = B_data[1];
    ss.C[0] = C_data[0]; ss.C[1] = C_data[1];
    ss.D[0] = D_data[0];

    double x[2] = {0.0, 0.0}; /* initial: zero current, zero speed */
    double u[1] = {12.0};
    double dt_ms = 1.0;
    double t_total = 0.030; /* 30 ms */
    int steps_per = 10;
    int total_steps = (int)(t_total / (dt_ms*1e-3));

    for (int k = 0; k <= total_steps; k += steps_per) {
        printf("  %8.1f   %11.3f   %13.1f\n",
               k*dt_ms, x[0], x[1]);
        ss_simulate_rk4(&ss, x, u, k*dt_ms*1e-3, (k+steps_per)*dt_ms*1e-3, steps_per*10);
    }

    /* PWM analysis */
    printf("\n=== PWM Drive Analysis ===\n");
    double Vdc = 24.0;
    double duty = 0.5;
    double f_pwm = 20000.0;
    double V_avg = pwm_average_voltage(Vdc, duty);
    double I_ripple = pwm_current_ripple(Vdc, duty, f_pwm, motor.armature_inductance);
    double f_min = pwm_minimum_frequency(Vdc, duty, 0.1, motor.armature_inductance);
    printf("  PWM freq: %.0f Hz, Duty: %.0f%%\n", f_pwm, duty*100);
    printf("  Avg voltage: %.1f V\n", V_avg);
    printf("  Current ripple: %.3f A\n", I_ripple);
    printf("  Min freq for 0.1A ripple: %.0f Hz\n", f_min);

    printf("\n=== DC Motor Analysis Complete ===\n");
    return 0;
}
