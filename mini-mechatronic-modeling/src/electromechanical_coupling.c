/**
 * @file electromechanical_coupling.c
 * @brief L2-L4 Electromechanical coupling physics implementations.
 *
 * Covers the fundamental physical laws linking electrical and mechanical
 * domains in mechatronic actuators: Lorentz force, Faraday's law
 * (back-EMF), DC motor dynamic equations, energy conversion, and
 * speed-torque characteristics.
 *
 * 9-Layer Coverage:
 *   L2 — Electromechanical analogies
 *   L3 — Coupled ODE systems
 *   L4 — Conservation laws (energy, Lorentz, Faraday)
 *
 * References:
 *   Griffiths (2017), Krause et al. (2013), Fitzgerald et al. (2014)
 *   Ogata "Modern Control Engineering" (2010)
 */

#include "electromechanical_coupling.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================*/
/* L2: Electromechanical Analogy (Force-Voltage / Mobility)                  */
/*===========================================================================*/

double mech_impedance_to_electrical_fv(double Z_mech)
{
    /*
     * Force-Voltage analogy:
     *   F ↔ V, v ↔ I
     *   Z_elec = V/I = F/v = Z_mech
     *
     * Mechanical impedance Z_mech = F/v has same numerical value
     * as electrical impedance in this analogy.
     *
     * Ref: Firestone "A New Analogy Between Mechanical and
     *      Electrical Systems" (1933), JASA
     */
    return Z_mech;
}

double mech_impedance_to_electrical_fi(double Z_mech)
{
    /*
     * Force-Current analogy (dual of F-V):
     *   F ↔ I, v ↔ V
     *   Y_elec = I/V = F/v = Z_mech
     *
     * Here Z_mech maps to electrical admittance, so
     * electrical impedance = 1/Z_mech.
     */
    if (fabs(Z_mech) < 1e-30) return 1e30;
    return 1.0 / Z_mech;
}

double mass_to_capacitance_fv(double m)
{
    /*
     * In F-V analogy:
     *   F = m · dv/dt  ↔  I_eq = C_eq · dV/dt  where I = v, V = F
     * So: v = C_eq · dF/dt, but F = m·dv/dt
     * Solving: C_eq = m (mass maps to capacitance directly)
     *
     * Ref: Shearer et al. "Introduction to System Dynamics" (1967)
     */
    return m > 0.0 ? m : 0.0;
}

double spring_to_inductance_fv(double k)
{
    /*
     * In F-V analogy:
     *   Mechanical spring: F = k · ∫v dt
     *   Electrical: V = (1/C) · ∫I dt
     * So k ↔ 1/C. But we want equivalent inductance.
     *
     * Spring equation: v = (1/k) · dF/dt
     * Inductor equation: V = L · dI/dt
     * With F↔V, v↔I: F = L · dv/dt? No — careful:
     *
     * F↔V, v↔I:
     *   Spring: F = k ∫ v dt → V = k ∫ I dt → capacitor: V=(1/C)∫I dt
     *   So k ↔ 1/C
     *
     * For inductor equivalent we use F-I analogy or mechanical
     * inductance analogy: L_equiv = 1/k (units: H ↔ s²/kg)
     */
    if (k < 1e-30) return 1e30;
    return 1.0 / k;
}

double damper_to_resistance_fv(double b)
{
    /*
     * In F-V analogy:
     *   F = b · v  ↔  V = R · I
     * So: b ↔ R (damping maps to resistance)
     *
     * Electrical resistance = mechanical damping directly.
     */
    return b;
}

/*===========================================================================*/
/* L4.1: Lorentz Force                                                       */
/*===========================================================================*/

double lorentz_force(double I, double B, double L_len, double theta)
{
    /*
     * Lorentz force on a straight current-carrying wire:
     *   F = I · (L × B)
     *
     * Magnitude: F = I · L · B · sin(θ)
     * where θ is the angle between the conductor direction and B field.
     *
     * In motors, conductors are perpendicular to B (θ = π/2),
     * so sin(θ) = 1, giving F = I·B·L.
     *
     * Ref: Griffiths "Introduction to Electrodynamics" (2017) §5.1.2
     *      Eq. 5.12: F_mag = ∫ I(dl × B)
     */
    return I * B * L_len * sin(theta);
}

double lorentz_force_density(double J, double B)
{
    /*
     * Lorentz force per unit volume:
     *   f = J × B  [N/m³]
     *
     * For perpendicular J and B: f = J · B
     *
     * Ref: Griffiths (2017) §5.1.3
     */
    return J * B;
}

/*===========================================================================*/
/* L4.2: Faraday's Law — Back-EMF                                           */
/*===========================================================================*/

double motional_emf(double B, double L, double v)
{
    /*
     * Motional EMF: conductor moving through magnetic field:
     *   EMF = ∮ (v × B) · dl
     *
     * For straight wire moving perpendicular to B:
     *   EMF = v · B · L
     *
     * This is the microscopic origin of back-EMF in motors.
     *
     * Ref: Hayt & Buck (2018) §9.3, Faraday's Law
     */
    return B * L * v;
}

double dc_motor_back_emf(double K_e, double omega)
{
    /*
     * DC motor back-EMF:
     *   V_bemf = K_e · ω
     *
     * where K_e = N·φ/(2π) is the back-EMF constant,
     * N = number of armature conductors,
     * φ = flux per pole [Wb].
     *
     * This is derived from Faraday's law for a rotating machine.
     *
     * Ref: Krause et al. (2013) §2.3
     *
     * Complexity: O(1)
     */
    return K_e * omega;
}

double kt_ke_equivalence_check(double K_e, double K_t)
{
    /*
     * In SI units, the torque constant K_t and back-EMF constant
     * K_e are numerically equal:
     *
     *   K_t [N·m/A] = K_e [V·s/rad]
     *
     * This derives from power conservation:
     *   P_elec = V_bemf · I = (K_e·ω)·I
     *   P_mech = T_em · ω = (K_t·I)·ω
     *   P_elec = P_mech → K_e = K_t
     *
     * Ref: Fitzgerald et al. (2014) §7.2
     */
    double denom = fmax(fabs(K_e), fabs(K_t));
    if (denom < 1e-30) return 0.0;
    return fabs(K_e - K_t) / denom;
}

/*===========================================================================*/
/* L4.3: Electromechanical Torque                                            */
/*===========================================================================*/

double dc_motor_electromagnetic_torque(double K_t, double I_a)
{
    /*
     * DC motor torque:
     *   T_em = K_t · I_a
     *
     * This is the Lorentz force integrated over the armature winding
     * structure. For a machine with N conductors, pole flux φ:
     *   T_em = (N·φ/(2π)) · I_a = K_t · I_a
     *
     * Ref: Krause et al. (2013) §2.2
     */
    return K_t * I_a;
}

double bldc_motor_torque(double K_t, double I_q)
{
    /*
     * BLDC/PMSM torque (sinusoidal FOC):
     *   T_em = (3/2) · P_p · (λ_pm·I_q + (L_d - L_q)·I_d·I_q)
     *
     * For surface-mount PM (L_d ≈ L_q = L_s), simplified:
     *   T_em = (3/2) · P_p · λ_pm · I_q
     *        = (3/2) · K_t_per_phase · I_q
     *
     * For six-step (trapezoidal) BLDC: T_em = 2 · K_t_ph · I_dc
     *
     * This implementation uses the per-phase K_t convention where
     * K_t is the phase torque constant.
     *
     * Ref: Krishnan "Electric Motor Drives" (2001) §9.3
     */
    return 1.5 * K_t * I_q;
}

/*===========================================================================*/
/* L4.4: DC Motor Dynamic Equations                                          */
/*===========================================================================*/

double dc_motor_electrical_dynamics(double V_a, double R_a, double I_a,
                                    double L_a, double K_e, double omega)
{
    /*
     * DC motor armature circuit:
     *   V_a = R_a·I_a + L_a·(dI_a/dt) + K_e·ω
     *
     * Therefore:
     *   dI_a/dt = (V_a - R_a·I_a - K_e·ω) / L_a
     *
     * The three terms on RHS represent:
     *   1. V_a: applied voltage (control input)
     *   2. -R_a·I_a: resistive (ohmic) voltage drop
     *   3. -K_e·ω: back-EMF (speed-dependent, couples elec↔mech)
     *
     * Ref: Ogata (2010) Example 3.6
     */
    if (L_a < 1e-30) return 0.0;  /* pure resistive case */
    return (V_a - R_a * I_a - K_e * omega) / L_a;
}

double dc_motor_mechanical_dynamics(double K_t, double I_a,
                                    double B, double omega,
                                    double T_load, double J)
{
    /*
     * DC motor mechanical dynamics (Newton's 2nd law for rotation):
     *   J · (dω/dt) = T_em - B·ω - T_load
     *
     * Therefore:
     *   dω/dt = (K_t·I_a - B·ω - T_load) / J
     *
     * Terms:
     *   1. K_t·I_a: electromagnetic driving torque
     *   2. -B·ω: viscous friction torque (opposes motion)
     *   3. -T_load: external load torque
     *
     * Ref: Cetinkunt (2007) §6.2
     */
    if (J < 1e-30) return 0.0;
    return (K_t * I_a - B * omega - T_load) / J;
}

void dc_motor_state_derivatives(const DCMotorParams *motor,
                                double V_a, double T_L,
                                double I_a, double omega,
                                double *dI_dt, double *domega_dt)
{
    /*
     * Combined DC motor state-space model (2nd order):
     *
     *   d/dt [I_a ] = [ -R_a/L_a    -K_e/L_a ] [I_a ] + [ 1/L_a    0  ] [V_a]
     *        [omega]   [  K_t/J      -B/J    ] [omega]   [   0    -1/J ] [T_L]
     *
     * This is the fundamental model for DC motor control.
     *
     * Ref: Ogata (2010) Example 3.6, Eq. 3-28
     */
    if (!motor) return;

    double L_a = motor->L_a;
    double J_m = motor->J_m;

    if (dI_dt) {
        *dI_dt = (L_a > 1e-30)
                 ? (V_a - motor->R_a * I_a - motor->K_e * omega) / L_a
                 : 0.0;
    }
    if (domega_dt) {
        *domega_dt = (J_m > 1e-30)
                     ? (motor->K_t * I_a - motor->B_m * omega - T_L) / J_m
                     : 0.0;
    }
}

void dc_motor_steady_state(const DCMotorParams *motor,
                           double V_a, double T_L,
                           double *I_ss, double *omega_ss, double *T_em_ss)
{
    /*
     * Steady-state DC motor operating point. Set derivatives to zero:
     *
     * 0 = V_a - R_a·I_ss - K_e·ω_ss          (1)
     * 0 = K_t·I_ss - B·ω_ss - T_L          (2)
     *
     * From (2): I_ss = (B·ω_ss + T_L) / K_t
     * Substitute into (1):
     *   V_a = R_a·(B·ω_ss + T_L)/K_t + K_e·ω_ss
     *
     * Since K_t = K_e in SI:
     *   ω_ss = (K_t·V_a - R_a·T_L) / (K_t² + R_a·B)
     *
     * Then: I_ss = (V_a - K_e·ω_ss) / R_a
     *       T_em_ss = K_t · I_ss
     *
     * Ref: Fitzgerald et al. (2014) §7.4
     */
    if (!motor || !I_ss || !omega_ss || !T_em_ss) return;

    double R = motor->R_a;
    double K = motor->K_e;  /* = K_t in SI */
    double B = motor->B_m;

    double denom = K * K + R * B;
    if (denom < 1e-30) {
        *I_ss = 0.0;
        *omega_ss = 0.0;
        *T_em_ss = 0.0;
        return;
    }

    *omega_ss = (K * V_a - R * T_L) / denom;
    *I_ss = (V_a - K * (*omega_ss)) / R;
    *T_em_ss = K * (*I_ss);
}

/*===========================================================================*/
/* L4.5: DC Motor Time Constants                                             */
/*===========================================================================*/

double dc_motor_tau_electrical(const DCMotorParams *motor)
{
    /*
     * Electrical time constant τ_e = L_a / R_a
     *
     * Governs how fast armature current responds to voltage changes.
     * Typical values: 0.1–10 ms for small servo motors.
     *
     * Ref: Ogata (2010) §3.4
     */
    if (!motor || motor->R_a < 1e-30) return 0.0;
    return motor->L_a / motor->R_a;
}

double dc_motor_tau_mechanical(double J_total, double B_total)
{
    /*
     * Mechanical time constant τ_m = J_total / B_total
     *
     * Governs speed response under constant torque.
     * Much slower than τ_e typically (10-100x).
     *
     * Ref: Ogata (2010) §3.4
     */
    if (B_total < 1e-30) return 1e30;
    return J_total / B_total;
}

double dc_motor_tau_electromechanical(const DCMotorParams *motor,
                                       double J_total, double B_total)
{
    (void)B_total;  /* not needed for electromechanical τ */
    /*
     * The electromechanical time constant accounts for the coupling
     * between electrical and mechanical dynamics.
     *
     * For a DC motor where τ_e << τ_m, the effective response from
     * voltage to speed is approximately first-order with:
     *
     *   τ_em ≈ J_total · R_a / (K_t · K_e)
     *
     * This is the dominant pole location in the voltage-to-speed
     * transfer function when inductance is negligible.
     *
     * Ref: Cetinkunt (2007) §6.2.3
     */
    if (!motor) return 0.0;
    double denom = motor->K_t * motor->K_e;
    if (denom < 1e-30) return 1e30;
    return J_total * motor->R_a / denom;
}

/*===========================================================================*/
/* L4.6: Speed-Torque Characteristic and Power                               */
/*===========================================================================*/

double dc_motor_speed_from_torque(const DCMotorParams *motor,
                                  double V_a, double T_em)
{
    /*
     * DC motor speed-torque line:
     *   ω = V_a/K_e - (R_a/(K_t·K_e)) · T_em
     *
     * At T_em = 0: ω_nl = V_a/K_e (no-load speed)
     * At ω = 0: T_stall = K_t · V_a / R_a (stall torque)
     *
     * This is a straight line in the (T, ω) plane.
     *
     * Ref: Fitzgerald et al. (2014) §7.5
     */
    if (!motor || motor->K_e < 1e-30) return 0.0;
    return V_a / motor->K_e
           - (motor->R_a / (motor->K_t * motor->K_e)) * T_em;
}

double dc_motor_no_load_speed(const DCMotorParams *motor, double V_a)
{
    /*
     * No-load speed: ω_nl = V_a / K_e
     * At no-load: I_a ≈ 0, T_em ≈ 0
     */
    if (!motor || motor->K_e < 1e-30) return 0.0;
    return V_a / motor->K_e;
}

double dc_motor_stall_torque(const DCMotorParams *motor, double V_a)
{
    /*
     * Stall torque: T_stall = K_t · V_a / R_a
     * At stall (ω=0): I_a = V_a/R_a, T_em = K_t·I_a
     */
    if (!motor || motor->R_a < 1e-30) return 0.0;
    return motor->K_t * V_a / motor->R_a;
}

double dc_motor_mechanical_power(double T_em, double omega)
{
    /*
     * Mechanical output power: P_mech = T_em · ω [W]
     *
     * This is the power delivered to the load shaft.
     */
    return T_em * omega;
}

double dc_motor_electrical_power(double V_a, double I_a)
{
    /*
     * Electrical input power: P_elec = V_a · I_a [W]
     */
    return V_a * I_a;
}

double dc_motor_copper_loss(double I_a, double R_a)
{
    /*
     * Copper (I²R) losses in armature winding:
     *   P_cu = I_a² · R_a [W]
     *
     * This is the dominant loss mechanism in most motors.
     * At stall, all electrical power becomes copper loss.
     */
    return I_a * I_a * R_a;
}

double dc_motor_efficiency(double T_em, double omega,
                           double V_a, double I_a)
{
    /*
     * Motor efficiency: η = P_mech / P_elec
     *
     * Maximum efficiency occurs near rated load, typically 70-90%.
     * Efficiency drops to zero at no-load (P_mech=0) and at stall (ω=0).
     *
     * Ref: Cetinkunt (2007) §6.4
     */
    double P_elec = V_a * I_a;
    if (fabs(P_elec) < 1e-30) return 0.0;
    double P_mech = T_em * omega;
    return P_mech / P_elec;
}

double dc_motor_max_mechanical_power(double V_a, double R_a)
{
    /*
     * Maximum mechanical power occurs at half the no-load speed
     * and half the stall torque:
     *
     *   P_max = V_a² / (4 · R_a)
     *
     * At this point: ω = ω_nl/2, T = T_stall/2
     * Efficiency at P_max = 50% (half power lost as copper loss).
     *
     * Ref: Fitzgerald et al. (2014) §7.5
     */
    if (R_a < 1e-30) return 0.0;
    return (V_a * V_a) / (4.0 * R_a);
}
