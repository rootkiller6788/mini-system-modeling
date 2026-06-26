/**
 * @file electromechanical_systems.c
 * @brief DC motor, generator, and electromechanical actuator models
 *
 * Knowledge points — each function = one independent concept:
 *   L1: DC motor/generator parameter structures and operating points
 *   L2: Steady-state analysis — speed, current, torque, power, efficiency
 *   L3: Transfer functions — voltage-to-speed, voltage-to-current
 *   L4: Time constants — electrical (La/Ra), mechanical (Ra*J/(Kt*Kb))
 *   L5: Motor parameter identification from measurements
 *   L6: PWM drive modeling — average voltage, current ripple
 *   L7: Motor sizing and gear ratio optimization
 *
 * Reference:
 *   MIT 6.302 Feedback Systems — DC motor as canonical plant
 *   Krause, Wasynczuk, Sudhoff, Analysis of Electric Machinery (2002)
 *   Fitzgerald, Kingsley, Umans, Electric Machinery (6th ed.)
 */

#include "electromechanical_systems.h"
#include <math.h>
#include <float.h>
#include <string.h>

/* ===== L2: Steady-State Analysis ===== */

double dc_motor_steady_state_speed(const DCMotorParams *params,
                                    double voltage, double load_torque)
{
    /* ω = (Kt*Va - Ra*TL) / (Kt*Kb + Ra*B) */
    if (!params) return 0.0;
    double num = params->torque_constant * voltage
                - params->armature_resistance * load_torque;
    double den = params->torque_constant * params->back_emf_constant
                + params->armature_resistance * params->viscous_friction;
    if (fabs(den) < 1e-30) return 0.0;
    return num / den;
}

double dc_motor_steady_state_current(const DCMotorParams *params,
                                      double voltage, double speed)
{
    /* ia = (Va - Kb*ω)/Ra */
    if (!params || fabs(params->armature_resistance) < 1e-30) return 0.0;
    return (voltage - params->back_emf_constant * speed)
           / params->armature_resistance;
}

double dc_motor_steady_state_torque(const DCMotorParams *params,
                                     double voltage, double speed)
{
    /* Te = Kt*ia = Kt*(Va - Kb*ω)/Ra */
    double ia = dc_motor_steady_state_current(params, voltage, speed);
    return params->torque_constant * ia;
}

/* ===== L2: Speed-Torque Characteristic ===== */

double dc_motor_no_load_speed(const DCMotorParams *params, double voltage)
{
    /* ω_nl = Va/Kb (when load torque and friction are zero) */
    if (!params || fabs(params->back_emf_constant) < 1e-30) return DBL_MAX;
    return voltage / params->back_emf_constant;
}

double dc_motor_stall_torque(const DCMotorParams *params, double voltage)
{
    /* T_stall = Kt*Va/Ra (at ω=0) */
    if (!params || fabs(params->armature_resistance) < 1e-30) return DBL_MAX;
    return params->torque_constant * voltage / params->armature_resistance;
}

DCMotorOperatingPoint dc_motor_operating_point(const DCMotorParams *params,
                                                 double voltage, double load_torque)
{
    DCMotorOperatingPoint op;
    if (!params) {
        memset(&op, 0, sizeof(op));
        return op;
    }

    op.speed  = dc_motor_steady_state_speed(params, voltage, load_torque);
    op.current = dc_motor_steady_state_current(params, voltage, op.speed);
    op.torque = params->torque_constant * op.current;
    op.voltage = voltage;

    op.power_in  = voltage * op.current;
    op.power_out = op.torque * op.speed;

    if (op.power_in > 1e-30)
        op.efficiency = op.power_out / op.power_in;
    else
        op.efficiency = 0.0;
    if (op.efficiency > 1.0) op.efficiency = 1.0;
    if (op.efficiency < 0.0) op.efficiency = 0.0;

    return op;
}

/* ===== L3: Transfer Functions ===== */

int dc_motor_transfer_function_voltage_to_speed(
    const DCMotorParams *params,
    double *num, double *den)
{
    /* H(s) = Kt / ((La*s + Ra)(J*s + B) + Kt*Kb)
     *      = Kt / (La*J*s^2 + (La*B + Ra*J)*s + (Ra*B + Kt*Kb))
     */
    if (!params || !num || !den) return -1;

    double a0 = params->armature_resistance * params->viscous_friction
              + params->torque_constant * params->back_emf_constant;
    double a1 = params->armature_inductance * params->viscous_friction
              + params->armature_resistance * params->rotor_inertia;
    double a2 = params->armature_inductance * params->rotor_inertia;

    num[0] = params->torque_constant; /* b0 = Kt */
    num[1] = 0.0;                     /* b1 = 0 */

    if (fabs(a2) > 1e-30) {
        /* Full second-order model */
        den[0] = a0;
        den[1] = a1;
        den[2] = a2;
        return 2;
    } else {
        /* Simplified first-order model (neglect La) */
        den[0] = a0;
        den[1] = params->armature_resistance * params->rotor_inertia;
        return 1;
    }
}

int dc_motor_transfer_function_voltage_to_current(
    const DCMotorParams *params,
    double *num, double *den)
{
    /* H(s) = (J*s + B) / (La*J*s^2 + (La*B+Ra*J)*s + (Ra*B+Kt*Kb)) */
    if (!params || !num || !den) return -1;

    double a0 = params->armature_resistance * params->viscous_friction
              + params->torque_constant * params->back_emf_constant;
    double a1 = params->armature_inductance * params->viscous_friction
              + params->armature_resistance * params->rotor_inertia;
    double a2 = params->armature_inductance * params->rotor_inertia;

    num[0] = params->viscous_friction; /* b0 = B */
    num[1] = params->rotor_inertia;    /* b1 = J */

    if (fabs(a2) > 1e-30) {
        den[0] = a0; den[1] = a1; den[2] = a2;
        return 2;
    } else {
        den[0] = a0;
        den[1] = params->armature_resistance * params->rotor_inertia;
        return 1;
    }
}

/* ===== L4: Time Constants ===== */

double dc_motor_electrical_time_constant(const DCMotorParams *params)
{
    /* tau_e = La/Ra */
    if (!params || fabs(params->armature_resistance) < 1e-30) return DBL_MAX;
    return params->armature_inductance / params->armature_resistance;
}

double dc_motor_mechanical_time_constant(const DCMotorParams *params)
{
    /* tau_m = Ra*J / (Kt*Kb) (dominant time constant when La << Ra) */
    if (!params) return DBL_MAX;
    double ktkb = params->torque_constant * params->back_emf_constant;
    if (fabs(ktkb) < 1e-30) return DBL_MAX;
    return params->armature_resistance * params->rotor_inertia / ktkb;
}

/* ===== L5: Parameter Identification ===== */

DCMotorParams dc_motor_identify_parameters(const DCMotorMeasurements *meas)
{
    DCMotorParams params;
    memset(&params, 0, sizeof(params));

    if (!meas || fabs(meas->test_voltage) < 1e-30) return params;

    double Va = meas->test_voltage;
    double Ra = meas->measured_resistance;
    double w_nl = meas->no_load_speed_at_voltage;
    double I_nl = meas->no_load_current_at_voltage;
    double I_stall = meas->stall_current_at_voltage;

    params.armature_resistance = Ra;

    /* Back-EMF constant: Kb = (Va - I_nl*Ra) / w_nl */
    if (fabs(w_nl) > 1e-30)
        params.back_emf_constant = (Va - I_nl * Ra) / w_nl;

    /* Torque constant = back-EMF constant (in SI units: Nm/A = V/(rad/s)) */
    params.torque_constant = params.back_emf_constant;

    /* Viscous friction: B = Kt * I_nl / w_nl */
    if (fabs(w_nl) > 1e-30)
        params.viscous_friction = params.torque_constant * I_nl / w_nl;

    /* Rated values from stall test */
    params.rated_voltage = Va;
    params.rated_current = I_stall;
    params.rated_torque = params.torque_constant * I_stall;
    params.rated_speed = w_nl;

    /* Rotor inertia estimated from mechanical time constant */
    /* Default estimate: tau_m ~ 0.1s for small motors */
    params.rotor_inertia = 1e-6;

    /* Armature inductance estimated from electrical time constant */
    params.armature_inductance = 1e-4;

    return params;
}

/* ===== L6: PWM Drive Modeling ===== */

double pwm_average_voltage(double supply_voltage, double duty_cycle)
{
    /* V_avg = D * V_supply */
    if (duty_cycle < 0.0) duty_cycle = 0.0;
    if (duty_cycle > 1.0) duty_cycle = 1.0;
    return supply_voltage * duty_cycle;
}

double pwm_current_ripple(double supply_voltage, double duty_cycle,
                           double pwm_freq, double inductance)
{
    /* Delta_I = V_supply * D * (1-D) / (f_pwm * L) */
    if (duty_cycle < 0.0) duty_cycle = 0.0;
    if (duty_cycle > 1.0) duty_cycle = 1.0;
    if (fabs(pwm_freq) < 1e-30 || fabs(inductance) < 1e-30) return DBL_MAX;
    return supply_voltage * duty_cycle * (1.0 - duty_cycle)
           / (pwm_freq * inductance);
}

double pwm_minimum_frequency(double supply_voltage, double duty_cycle,
                              double max_ripple, double inductance)
{
    /* f_min = V_supply * D * (1-D) / (Delta_I_max * L) */
    if (duty_cycle < 0.0) duty_cycle = 0.0;
    if (duty_cycle > 1.0) duty_cycle = 1.0;
    if (fabs(max_ripple) < 1e-30 || fabs(inductance) < 1e-30) return DBL_MAX;
    return supply_voltage * duty_cycle * (1.0 - duty_cycle)
           / (max_ripple * inductance);
}

/* ===== L2: DC Generator Model ===== */

double dc_generator_terminal_voltage(const DCGeneratorParams *params,
                                      double speed, double load_current)
{
    /* Vt = Kb*w - Ia*Ra */
    if (!params) return 0.0;
    return params->back_emf_constant * speed
           - load_current * params->armature_resistance;
}

double dc_generator_output_power(const DCGeneratorParams *params,
                                  double speed, double load_current)
{
    /* P_out = Vt * Ia */
    double vt = dc_generator_terminal_voltage(params, speed, load_current);
    return vt * load_current;
}

/* ===== L7: Motor Sizing ===== */

MotorSizingResult dc_motor_sizing(const MotorSizingRequirements *req)
{
    MotorSizingResult result;
    memset(&result, 0, sizeof(result));

    if (!req) return result;

    double omega_max = req->max_speed_rad_s;
    double torque_max = req->max_torque_nm;
    double Vdc = req->supply_voltage;

    if (omega_max < 1e-30 || fabs(Vdc) < 1e-30) return result;

    /* Peak mechanical power */
    double P_mech = torque_max * omega_max;
    result.estimated_power_w = P_mech;

    /* Recommended Kt = Vdc/omega_max (voltage-matched at max speed) */
    result.recommended_kt = Vdc / omega_max;

    /* Required current */
    if (result.recommended_kt > 1e-30)
        result.estimated_current_a = torque_max / result.recommended_kt;
    else
        result.estimated_current_a = torque_max;

    /* Recommended Ra for acceptable efficiency (~85% at max power) */
    double target_eff = 0.85;
    double P_loss = P_mech * (1.0 - target_eff) / target_eff;
    if (result.estimated_current_a > 1e-30)
        result.recommended_ra = P_loss
            / (result.estimated_current_a * result.estimated_current_a);
    else
        result.recommended_ra = 1.0;

    return result;
}

/* ===== L7: Optimal Gear Ratio ===== */

double optimal_gear_ratio(double motor_inertia, double load_inertia)
{
    /* N_opt = sqrt(J_load / J_motor) for minimum acceleration time */
    if (fabs(motor_inertia) < 1e-30) return DBL_MAX;
    if (load_inertia < 0.0) return 1.0;
    return sqrt(load_inertia / motor_inertia);
}
