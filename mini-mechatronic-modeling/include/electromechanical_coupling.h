#ifndef ELECTROMECHANICAL_COUPLING_H
#define ELECTROMECHANICAL_COUPLING_H

/**
 * @file electromechanical_coupling.h
 * @brief L2-L4 Electromechanical Coupling Physics
 *
 * Implements the fundamental laws linking electrical and mechanical domains:
 * Lorentz force, back-EMF, energy conversion, and DC motor dynamic equations.
 *
 * Key References:
 *  - Krause et al. "Analysis of Electric Machinery and Drive Systems" (2013)
 *  - Krishnan, R. "Electric Motor Drives" (2001)
 *  - Fitzgerald, Kingsley, Umans "Electric Machinery" (2014)
 *
 * 9-Layer Coverage:
 *   L2 — Electromechanical analogy (F ↔ V, v ↔ I, etc.)
 *   L3 — Coupled ODEs, complex power, Laplace domain
 *   L4 — Lorentz/Faraday/Ampere, conservation of energy
 */

#include "mechatronic_definitions.h"

/*===========================================================================*/
/* L2: Electromechanical Analogy                                             */
/*===========================================================================*/

/**
 * @brief The force-voltage (mobility) analogy mapping.
 *
 * Mechanical      ↔  Electrical
 *   Force   (F)   ↔  Voltage (V)
 *   Velocity (v)  ↔  Current (I)
 *   Mass    (m)   ↔  Capacitance (C)
 *   Spring  (k)   ↔  Inductance inverse (1/L)
 *   Damper  (b)   ↔  Resistance inverse (1/R)
 */

/** Convert mechanical impedance to electrical impedance via F-V analogy
 * @param Z_mech  Mechanical impedance (F/v) [N·s/m]
 * @return Electrical impedance [Ω] */
double mech_impedance_to_electrical_fv(double Z_mech);

/** Convert mechanical impedance to electrical impedance via F-I analogy
 *  (force-current / velocity-voltage analogy)
 * @param Z_mech  Mechanical impedance (F/v) [N·s/m]
 * @return Electrical impedance [1/Ω] */
double mech_impedance_to_electrical_fi(double Z_mech);

/** Mobility analogy: mass → capacitance
 *  @param m  Mass [kg]
 *  @return Equivalent capacitance [F] */
double mass_to_capacitance_fv(double m);

/** Mobility analogy: spring → inverse inductance
 *  @param k  Spring stiffness [N/m]
 *  @return Equivalent inductance [H] */
double spring_to_inductance_fv(double k);

/** Mobility analogy: damper → inverse resistance
 *  @param b  Damping coefficient [N·s/m]
 *  @return Equivalent resistance [Ω] */
double damper_to_resistance_fv(double b);

/*===========================================================================*/
/* L4.1: Lorentz Force Law                                                   */
/*===========================================================================*/

/** Lorentz force on a current-carrying conductor in magnetic field
 *  F = I · L × B (vector form, simplified to scalar here)
 *  @param I      Current [A]
 *  @param B      Magnetic flux density [T]
 *  @param L_len   Conductor length perpendicular to B [m]
 *  @param theta  Angle between current and B-field [rad]
 *  @return Force [N]
 *
 *  Ref: Griffiths "Introduction to Electrodynamics" (2017), Ch.5 */
double lorentz_force(double I, double B, double L_len, double theta);

/** Lorentz force density (force per unit volume)
 *  @param J      Current density [A/m²]
 *  @param B      Magnetic flux density [T]
 *  @return Force density [N/m³] */
double lorentz_force_density(double J, double B);

/*===========================================================================*/
/* L4.2: Faraday's Law of Induction (Back-EMF)                               */
/*===========================================================================*/

/** Motional EMF in a conductor moving through magnetic field
 *  V_emf = B · L · v  (for perpendicular motion)
 *  @param B    Magnetic flux density [T]
 *  @param L    Conductor length [m]
 *  @param v    Velocity perpendicular to B [m/s]
 *  @return Induced EMF [V]
 *
 *  Ref: Hayt & Buck "Engineering Electromagnetics" (2018), Ch.9 */
double motional_emf(double B, double L, double v);

/** Rotational back-EMF in a DC motor
 *  V_bemf = K_e · ω
 *  where K_e = N · φ / (2π) with N=conductors, φ=flux per pole
 *  @param K_e  Back-EMF constant [V·s/rad]
 *  @param omega Angular velocity [rad/s]
 *  @return Back-EMF voltage [V] */
double dc_motor_back_emf(double K_e, double omega);

/** Equivalence of torque constant and back-EMF constant in SI units
 *  K_t [N·m/A] == K_e [V·s/rad]  (derives from conservation of energy)
 *  @param K_e  Back-EMF constant [V·s/rad]
 *  @param K_t  Torque constant [N·m/A]
 *  @return Absolute relative error |K_t - K_e| / max(K_t, K_e) */
double kt_ke_equivalence_check(double K_e, double K_t);

/*===========================================================================*/
/* L4.3: Electromechanical Torque Generation                                 */
/*===========================================================================*/

/** DC motor electromagnetic torque
 *  T_em = K_t · I_a
 *  @param K_t  Torque constant [N·m/A]
 *  @param I_a  Armature current [A]
 *  @return Electromagnetic torque [N·m] */
double dc_motor_electromagnetic_torque(double K_t, double I_a);

/** BLDC motor torque (simplified, sinusoidal commutation)
 *  T_em = (3/2) · K_t · I_q   where I_q is q-axis current
 *  @param K_t  Torque constant [N·m/A]
 *  @param I_q  Quadrature-axis current [A]
 *  @return Torque [N·m] */
double bldc_motor_torque(double K_t, double I_q);

/*===========================================================================*/
/* L4.4: DC Motor Dynamic Equations                                          */
/*===========================================================================*/

/** DC motor electrical dynamics:
 *  V_a = R_a·I_a + L_a·dI_a/dt + K_e·ω
 *
 *  Compute dI_a/dt = (V_a - R_a·I_a - K_e·ω) / L_a
 *  @param V_a    Applied armature voltage [V]
 *  @param R_a    Armature resistance [Ω]
 *  @param I_a    Armature current [A]
 *  @param L_a    Armature inductance [H]
 *  @param K_e    Back-EMF constant [V·s/rad]
 *  @param omega  Rotor angular velocity [rad/s]
 *  @return dI_a/dt [A/s] */
double dc_motor_electrical_dynamics(double V_a, double R_a, double I_a,
                                    double L_a, double K_e, double omega);

/** DC motor mechanical dynamics:
 *  J·dω/dt = T_em - B·ω - T_L - T_friction
 *
 *  Compute dω/dt = (K_t·I_a - B·ω - T_load) / J
 *  @param K_t     Torque constant [N·m/A]
 *  @param I_a     Armature current [A]
 *  @param B       Viscous friction coefficient [N·m·s/rad]
 *  @param omega   Angular velocity [rad/s]
 *  @param T_load  Load torque [N·m]
 *  @param J       Total moment of inertia [kg·m²]
 *  @return dω/dt [rad/s²] */
double dc_motor_mechanical_dynamics(double K_t, double I_a,
                                    double B, double omega,
                                    double T_load, double J);

/** Combined state-space DC motor model [I_a, omega]^T
 *
 *  d/dt [I_a ] = [-R_a/L_a    -K_e/L_a  ] [I_a ] + [1/L_a   0   ] [V_a ]
 *       [omega]   [ K_t/J      -B/J     ] [omega]   [ 0    -1/J ] [T_L ]
 *
 *  @param motor  DC motor parameters
 *  @param V_a    Applied voltage [V]
 *  @param T_L    Load torque [N·m]
 *  @param I_a    Current state [A] (input/output)
 *  @param omega  Speed state [rad/s] (input/output)
 *  @param dI_dt  Output: current derivative [A/s]
 *  @param domega_dt Output: speed derivative [rad/s²]
 *
 *  Ref: Ogata "Modern Control Engineering" (2010), Ch.3 */
void dc_motor_state_derivatives(const DCMotorParams *motor,
                                double V_a, double T_L,
                                double I_a, double omega,
                                double *dI_dt, double *domega_dt);

/** Compute the steady-state operating point of a DC motor
 *  I_ss = (V_a - K_e·ω) / R_a
 *  ω_ss = (K_t·V_a - R_a·T_L) / (K_t·K_e + R_a·B)
 *
 *  @param motor   Motor parameters
 *  @param V_a     Applied voltage [V]
 *  @param T_L     Load torque [N·m]
 *  @param I_ss    Output: steady-state current [A]
 *  @param omega_ss Output: steady-state speed [rad/s]
 *  @param T_em_ss Output: steady-state torque [N·m] */
void dc_motor_steady_state(const DCMotorParams *motor,
                           double V_a, double T_L,
                           double *I_ss, double *omega_ss, double *T_em_ss);

/*===========================================================================*/
/* L4.5: DC Motor Electrical and Mechanical Time Constants                   */
/*===========================================================================*/

/** Electrical time constant: τ_e = L_a / R_a
 *  Characterizes the armature current response
 *  @param motor  Motor parameters (reads R_a, L_a)
 *  @return Electrical time constant [s] */
double dc_motor_tau_electrical(const DCMotorParams *motor);

/** Mechanical time constant: τ_m = J_total / B_total
 *  Characterizes the speed response with fixed voltage
 *  @param J_total  Total inertia [kg·m²]
 *  @param D_total  Total damping [N·m·s/rad]
 *  @return Mechanical time constant [s] */
double dc_motor_tau_mechanical(double J_total, double D_total);

/** Electromechanical time constant (coupled)
 *  τ_em accounts for both electrical and mechanical dynamics
 *  @param motor  Motor parameters
 *  @param J_total Total reflected inertia [kg·m²]
 *  @param B_total Total damping [N·m·s/rad]
 *  @return Electromechanical time constant [s] */
double dc_motor_tau_electromechanical(const DCMotorParams *motor,
                                       double J_total, double B_total);

/*===========================================================================*/
/* L4.6: Motor Speed-Torque Characteristic                                   */
/*===========================================================================*/

/** Motor speed-torque curve: ω = (V_a/K_e) - (R_a/(K_t·K_e)) · T_em
 *  @param motor  Motor parameters
 *  @param V_a    Armature voltage [V]
 *  @param T_em   Electromagnetic torque [N·m]
 *  @return Speed [rad/s] */
double dc_motor_speed_from_torque(const DCMotorParams *motor,
                                  double V_a, double T_em);

/** No-load speed: ω_nl = V_a / K_e
 *  @param motor  Motor parameters
 *  @param V_a    Armature voltage [V]
 *  @return No-load speed [rad/s] */
double dc_motor_no_load_speed(const DCMotorParams *motor, double V_a);

/** Stall torque: T_stall = K_t · V_a / R_a
 *  @param motor  Motor parameters
 *  @param V_a    Armature voltage [V]
 *  @return Stall torque [N·m] */
double dc_motor_stall_torque(const DCMotorParams *motor, double V_a);

/** Motor output mechanical power: P_mech = T_em · ω
 *  @param T_em  Electromagnetic torque [N·m]
 *  @param omega Angular velocity [rad/s]
 *  @return Mechanical output power [W] */
double dc_motor_mechanical_power(double T_em, double omega);

/** Motor electrical input power: P_elec = V_a · I_a
 *  @param V_a   Terminal voltage [V]
 *  @param I_a   Armature current [A]
 *  @return Electrical input power [W] */
double dc_motor_electrical_power(double V_a, double I_a);

/** Motor copper loss: P_cu = I_a² · R_a
 *  @param I_a   Armature current [A]
 *  @param R_a   Armature resistance [Ω]
 *  @return Copper loss [W] */
double dc_motor_copper_loss(double I_a, double R_a);

/** Motor efficiency: η = P_mech / P_elec
 *  @param T_em   Torque [N·m]
 *  @param omega  Speed [rad/s]
 *  @param V_a    Voltage [V]
 *  @param I_a    Current [A]
 *  @return Efficiency (0-1) */
double dc_motor_efficiency(double T_em, double omega,
                           double V_a, double I_a);

/** Maximum mechanical power point of DC motor
 *  At half stall torque, half no-load speed:
 *  P_max = (V_a²) / (4 · R_a)
 *  @param V_a    Armature voltage [V]
 *  @param R_a    Armature resistance [Ω]
 *  @return Maximum mechanical power [W] */
double dc_motor_max_mechanical_power(double V_a, double R_a);

#endif /* ELECTROMECHANICAL_COUPLING_H */
