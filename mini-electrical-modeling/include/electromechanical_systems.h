/**
 * @file electromechanical_systems.h
 * @brief Electromechanical system models DC motors, generators, actuators
 *
 * Knowledge Coverage:
 *   L1 - Definitions: back-EMF constant, torque constant, armature,
 *        field winding, mechanical time constant, electrical time constant
 *   L2 - Core Concepts: Lorentz force, Faraday induction in rotating machines,
 *        energy conversion efficiency, electromechanical coupling
 *   L3 - Mathematical Structures: coupled ODEs (electrical + mechanical),
 *        transfer functions with mechanical load, state-space models
 *   L4 - Fundamental Laws: Faraday law (back-EMF), Lorentz force (torque),
 *        Newton second law (mechanical), energy conservation (power balance)
 *   L5 - Engineering Methods: motor parameter identification,
 *        steady-state speed-torque curves, PWM drive modeling
 *
 * Reference:
 *   MIT 6.302 Feedback Systems — DC motor as canonical plant
 *   Stanford ENGR105 — Servomotor transfer functions
 *   Krause, Wasynczuk, and Sudhoff, Analysis of Electric Machinery (2002)
 *   Fitzgerald, Kingsley, and Umans, Electric Machinery (6th ed.)
 */

#ifndef ELECTROMECHANICAL_SYSTEMS_H
#define ELECTROMECHANICAL_SYSTEMS_H

#include <stdlib.h>
#include <math.h>

/* L1: DC motor parameters */
typedef struct {
    double armature_resistance;     /* Ra [ohm] */
    double armature_inductance;     /* La [H] */
    double back_emf_constant;       /* Kb [V/(rad/s)] = [V·s/rad] */
    double torque_constant;         /* Kt [N·m/A] */
    double rotor_inertia;           /* J [kg·m^2] */
    double viscous_friction;        /* B [N·m/(rad/s)] */
    double rated_voltage;           /* [V] */
    double rated_current;           /* [A] */
    double rated_speed;             /* [rad/s] */
    double rated_torque;            /* [N·m] */
} DCMotorParams;

/* L1: DC motor operating point */
typedef struct {
    double voltage;       /* Terminal voltage [V] */
    double current;       /* Armature current [A] */
    double speed;         /* Angular velocity [rad/s] */
    double torque;        /* Electromagnetic torque [N·m] */
    double power_in;      /* Electrical input power [W] */
    double power_out;     /* Mechanical output power [W] */
    double efficiency;    /* Power efficiency [0,1] */
} DCMotorOperatingPoint;

/**
 * L2: DC motor electrical dynamics:
 *   Va = Ra·ia + La·dia/dt + eb
 *   eb = Kb·ω (back-EMF)
 *
 * L4: DC motor mechanical dynamics:
 *   J·dω/dt = Te - B·ω - TL
 *   Te = Kt·ia (electromagnetic torque)
 *
 * L2: Steady-state analysis
 */

/**
 * @brief Steady-state speed from applied voltage and load torque
 *
 * At steady state (dia/dt = 0, dω/dt = 0):
 *   ω = (Kt·Va - Ra·TL) / (Kt·Kb + Ra·B)
 *
 * Theorem: Derived from coupled electrical-mechanical equilibrium.
 *   No-load speed: ω_nl = Va/Kb (when TL = 0, B ≈ 0)
 *   Stall torque: T_stall = Kt·Va/Ra (when ω = 0)
 *   Speed-torque curve is linear for DC motors.
 *
 * Complexity: O(1)
 */
double dc_motor_steady_state_speed(const DCMotorParams *params,
                                    double voltage, double load_torque);

/**
 * @brief Steady-state current
 *
 * ia = (Va - Kb·ω) / Ra
 * Complexity: O(1)
 */
double dc_motor_steady_state_current(const DCMotorParams *params,
                                      double voltage, double speed);

/**
 * @brief Steady-state torque
 *
 * Te = Kt·ia = Kt·(Va - Kb·ω)/Ra
 * At ω = 0 (stall): T_stall = Kt·Va/Ra
 * Complexity: O(1)
 */
double dc_motor_steady_state_torque(const DCMotorParams *params,
                                     double voltage, double speed);

/* L2: Speed-torque characteristic */

/**
 * @brief No-load speed: ω_nl = Va/Kb
 *
 * The speed at which back-EMF equals applied voltage, zero net torque.
 * Complexity: O(1)
 */
double dc_motor_no_load_speed(const DCMotorParams *params, double voltage);

/**
 * @brief Stall torque: T_stall = Kt·Va/Ra
 *
 * Maximum torque at zero speed (stalled condition).
 * Complexity: O(1)
 */
double dc_motor_stall_torque(const DCMotorParams *params, double voltage);

/**
 * @brief Full operating point computation
 *
 * Solves coupled steady-state equations iteratively.
 * Complexity: O(1) for linear case
 */
DCMotorOperatingPoint dc_motor_operating_point(const DCMotorParams *params,
                                                 double voltage, double load_torque);

/* L3: Transfer function from voltage to speed */

/**
 * @brief DC motor transfer function: H(s) = ω(s)/Va(s)
 *
 * H(s) = Kt / ((La·s + Ra)(J·s + B) + Kt·Kb)
 *
 * Simplified (neglecting La when La << Ra):
 *   H(s) ≈ Kt / (Ra·J·s + Ra·B + Kt·Kb)
 *        = K / (τ·s + 1)  where  K = Kt/(Ra·B + Kt·Kb), τ = Ra·J/(Ra·B + Kt·Kb)
 *
 * Full second-order:
 *   H(s) = K / (s^2/ω_n^2 + 2ζ·s/ω_n + 1)
 *   with ω_n = sqrt((Ra·B + Kt·Kb)/(La·J))
 *        ζ = (Ra·J + La·B) / (2·sqrt(La·J·(Ra·B + Kt·Kb)))
 *
 * @param num Output numerator coefficients [b0, b1, ...] (caller allocates 3 doubles)
 * @param den Output denominator coefficients [a0, a1, a2] (caller allocates 3 doubles)
 * @return Order of transfer function (2 for full, 1 for simplified if La=0)
 * Complexity: O(1)
 */
int dc_motor_transfer_function_voltage_to_speed(
    const DCMotorParams *params,
    double *num, double *den);

/**
 * @brief Transfer function from voltage to current: H(s) = Ia(s)/Va(s)
 *
 * H(s) = (J·s + B) / ((La·s + Ra)(J·s + B) + Kt·Kb)
 * Complexity: O(1)
 */
int dc_motor_transfer_function_voltage_to_current(
    const DCMotorParams *params,
    double *num, double *den);

/* L3: Mechanical time constant and electrical time constant */

/**
 * @brief Electrical time constant: τ_e = La/Ra [s]
 *
 * Determines how fast current responds to voltage changes.
 * Complexity: O(1)
 */
double dc_motor_electrical_time_constant(const DCMotorParams *params);

/**
 * @brief Mechanical time constant: τ_m = Ra·J/(Kt·Kb) [s]
 *
 * Determines how fast speed responds to voltage changes.
 * Typically τ_m >> τ_e for most motors.
 * Complexity: O(1)
 */
double dc_motor_mechanical_time_constant(const DCMotorParams *params);

/* L5: Motor parameter identification */

/**
 * @brief Estimate motor parameters from measurements
 *
 * Given measured: Ra (resistance), ω_nl at Va (no-load speed),
 * I_nl at Va (no-load current), I_stall at Va (stall current).
 *
 * Computes: Kb = (Va - I_nl·Ra)/ω_nl, Kt = Kb (in SI units),
 *           B = Kt·I_nl/ω_nl
 *
 * Complexity: O(1)
 */
typedef struct {
    double measured_resistance;
    double no_load_speed_at_voltage;
    double no_load_current_at_voltage;
    double stall_current_at_voltage;
    double test_voltage;
} DCMotorMeasurements;

DCMotorParams dc_motor_identify_parameters(const DCMotorMeasurements *meas);

/* L5: PWM drive modeling */

/**
 * @brief Average voltage from PWM duty cycle
 *
 * V_avg = D · V_supply where D = t_on/(t_on + t_off) ∈ [0, 1]
 *
 * Theorem: For PWM frequency >> electrical time constant (f_pwm > 10/τ_e),
 *   the motor acts as a low-pass filter and responds to average voltage.
 *   Ripple current: Δi = V_supply·D·(1-D)/(f_pwm·La)
 *
 * Complexity: O(1)
 */
double pwm_average_voltage(double supply_voltage, double duty_cycle);

/**
 * @brief PWM current ripple amplitude
 *
 * ΔI = V_supply·D·(1-D) / (f_pwm·La)
 *
 * Maximum ripple at D = 0.5: ΔI_max = V_supply/(4·f_pwm·La)
 *
 * Complexity: O(1)
 */
double pwm_current_ripple(double supply_voltage, double duty_cycle,
                           double pwm_freq, double inductance);

/**
 * @brief Minimum PWM frequency for given ripple specification
 *
 * f_min = V_supply·D·(1-D) / (ΔI_spec·La)
 *
 * Complexity: O(1)
 */
double pwm_minimum_frequency(double supply_voltage, double duty_cycle,
                              double max_ripple, double inductance);

/* L2: DC generator model */

typedef struct {
    double armature_resistance;    /* Ra [ohm] */
    double armature_inductance;    /* La [H] */
    double back_emf_constant;      /* Kb [V/(rad/s)] */
    double field_resistance;       /* Rf [ohm] — for separately excited */
    double field_inductance;       /* Lf [H] */
    double rated_power;            /* [W] */
    double rated_speed;            /* [rad/s] */
    double rated_voltage;          /* [V] */
} DCGeneratorParams;

/**
 * @brief Generator terminal voltage: Vt = Ea - Ia·Ra = Kb·ω - Ia·Ra
 *
 * Theorem: In generator mode, mechanical power is converted to electrical.
 *   Terminal voltage equals generated EMF minus armature resistance drop.
 *   Ea = Kb·ω where ω is the prime mover speed.
 *
 * Complexity: O(1)
 */
double dc_generator_terminal_voltage(const DCGeneratorParams *params,
                                      double speed, double load_current);

/**
 * @brief Generator output power: P_out = Vt·Ia
 *
 * Efficiency: η = P_out/(P_out + P_losses)
 *   P_losses = Ia^2·Ra (copper) + P_mech (friction, windage) + P_core (iron losses)
 *
 * Complexity: O(1)
 */
double dc_generator_output_power(const DCGeneratorParams *params,
                                  double speed, double load_current);

/* L7: Application-oriented motor sizing */

/**
 * @brief Compute required motor parameters for a given load profile
 *
 * Load requirements: max speed, max torque, desired mechanical time constant.
 * Returns recommended motor parameters.
 *
 * Reference: Maxon/FAULHABER motor selection guides
 * Complexity: O(1)
 */
typedef struct {
    double max_speed_rad_s;        /* Maximum required speed */
    double max_torque_nm;          /* Maximum required torque */
    double desired_bandwidth_hz;   /* Desired control bandwidth */
    double supply_voltage;         /* Available DC bus voltage */
} MotorSizingRequirements;

typedef struct {
    double recommended_kt;         /* Recommended torque constant [N·m/A] */
    double recommended_ra;         /* Recommended armature resistance [ohm] */
    double estimated_power_w;      /* Estimated mechanical power [W] */
    double estimated_current_a;    /* Estimated current draw [A] */
} MotorSizingResult;

MotorSizingResult dc_motor_sizing(const MotorSizingRequirements *req);

/**
 * @brief Gear ratio selection for optimal power transfer
 *
 * Theorem: For a motor driving an inertial load through a gearbox,
 *   optimal gear ratio N = sqrt(J_load/J_motor) minimizes
 *   acceleration time (matches load inertia reflected to motor shaft).
 *
 *   Reflected inertia: J_ref = J_load/N^2
 *   Reflected torque: T_ref = T_load/N
 *
 * @return Optimal gear ratio
 * Complexity: O(1)
 */
double optimal_gear_ratio(double motor_inertia, double load_inertia);

#endif /* ELECTROMECHANICAL_SYSTEMS_H */
