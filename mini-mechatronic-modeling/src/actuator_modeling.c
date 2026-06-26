/**
 * @file actuator_modeling.c
 * @brief L5-L7 Actuator selection, sizing, and modeling.
 *
 * Implements motor selection (RMS torque, peak torque, inertia matching),
 * thermal analysis, stepper motor torque curves, voice coil and piezo
 * actuator models.
 *
 * 9-Layer Coverage:
 *   L5 — Selection and sizing algorithms
 *   L6 — Canonical actuator problems
 *   L7 — Application constraints (robot joint, CNC, antenna)
 *   L8 — Piezoelectric actuators (advanced)
 *
 * References:
 *   Bolton (2015), Cetinkunt (2007), Morar (2003), Uchino (2010)
 *   Armstrong-Hélouvry (1994)
 */

#include "actuator_modeling.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================*/
/* L5: Duty Cycle Analysis                                                   */
/*===========================================================================*/

double rms_torque(const DutyCycle *cycle)
{
    /*
     * Root Mean Square (RMS) torque:
     *
     *   T_rms = √( (1/T_cycle) · Σ T_i² · Δt_i )
     *
     * The RMS torque determines the thermal (continuous) torque
     * requirement for the motor. The motor's rated continuous
     * torque must exceed T_rms plus a safety margin (typically 20%).
     *
     * Derivation: Copper losses ∝ I² ∝ T² (since T = K_t·I).
     * RMS torque represents the equivalent constant torque that
     * produces the same average heating.
     *
     * Ref: Cetinkunt (2007) §7.3
     *
     * Complexity: O(N) where N = number of segments
     */
    if (!cycle || cycle->n_segments <= 0) return 0.0;

    double sum_sq = 0.0;
    double total_time = 0.0;

    for (int i = 0; i < cycle->n_segments; i++) {
        double T = cycle->segments[i].torque;
        double dt = cycle->segments[i].duration;
        sum_sq += T * T * dt;
        total_time += dt;
    }

    if (total_time < 1e-30) return 0.0;
    if (cycle->cycle_time > total_time)
        total_time = cycle->cycle_time;  /* include idle time */

    return sqrt(sum_sq / total_time);
}

double peak_torque_requirement(const DutyCycle *cycle)
{
    /*
     * The peak torque is the maximum absolute torque across all
     * segments. This must be less than the motor's stall torque.
     *
     * T_peak = max(|T_i|)
     */
    if (!cycle || cycle->n_segments <= 0) return 0.0;

    double peak = 0.0;
    for (int i = 0; i < cycle->n_segments; i++) {
        double abs_t = fabs(cycle->segments[i].torque);
        if (abs_t > peak) peak = abs_t;
    }
    return peak;
}

double average_speed(const DutyCycle *cycle)
{
    /*
     * Time-weighted average speed:
     *
     *   ω_avg = (1/T_cycle) · Σ ω_i · Δt_i
     */
    if (!cycle || cycle->n_segments <= 0) return 0.0;

    double sum = 0.0;
    double total_time = 0.0;

    for (int i = 0; i < cycle->n_segments; i++) {
        sum += cycle->segments[i].speed * cycle->segments[i].duration;
        total_time += cycle->segments[i].duration;
    }

    if (total_time < 1e-30) return 0.0;
    return sum / total_time;
}

/*===========================================================================*/
/* L5: Gear Reflection                                                       */
/*===========================================================================*/

void reflect_load_through_gear(double J_load, double T_load,
                               const GearTrain *gear, int direction,
                               double *J_out, double *T_out)
{
    /*
     * Load reflection through a gear train:
     *
     *   J_reflected = J_load / N²
     *
     *   Motor driving load (direction=0):
     *     T_reflected = T_load / (N · η)
     *
     *   Load back-driving motor (direction=1):
     *     T_reflected = T_load · η / N
     *
     * Derivation: Power conservation + gear kinematics.
     *   P_out = η · P_in
     *   T_load · ω_load = η · T_motor · ω_motor
     *   ω_motor = N · ω_load
     *   → T_motor = T_load / (N · η) (driving)
     *
     * Ref: Cetinkunt (2007) §7.2
     */
    if (!gear) {
        if (J_out) *J_out = J_load;
        if (T_out) *T_out = T_load;
        return;
    }

    double N = gear->ratio;
    double N2 = N * N;

    if (J_out) *J_out = J_load / N2;

    if (T_out) {
        if (direction == 0) {
            /* Motor drives load */
            *T_out = T_load / (N * gear->efficiency);
        } else {
            /* Load back-drives motor (motoring as generator) */
            *T_out = T_load * gear->efficiency / N;
        }
    }
}

double optimal_gear_ratio_inertia_match(double J_m, double J_L)
{
    /*
     * Optimal gear ratio for inertia matching:
     *
     *   N_opt = √(J_L / J_m)
     *
     * The torque required to accelerate load at rate α is:
     *   T_req = (J_m + J_L/N²) · N · α + T_load/N
     *
     * Minimizing w.r.t. N for pure inertial load (T_load=0):
     *   d/dN [α · (J_m·N + J_L/N)] = 0
     *   → J_m - J_L/N² = 0
     *   → N_opt = √(J_L / J_m)
     *
     * This maximizes the load acceleration for a given motor torque.
     *
     * Ref: Cetinkunt (2007) §7.3.1
     */
    if (J_m < 1e-30) return 100.0;  /* avoid division by zero */
    return sqrt(J_L / J_m);
}

/*===========================================================================*/
/* L5: Motor Selection                                                       */
/*===========================================================================*/

MotorSelection select_dc_motor(const DCMotorParams *motor,
                               const GearTrain *gear,
                               const DutyCycle *cycle)
{
    /*
     * Motor selection algorithm:
     *
     * 1. Reflect load dynamics through gear
     * 2. Compute RMS torque (thermal constraint)
     * 3. Compute peak torque (stall constraint)
     * 4. Check speed constraint
     * 5. Evaluate margins
     *
     * Feasibility criteria:
     *   T_rms_reflected ≤ motor.max_torque (continuous)
     *   T_peak_reflected ≤ motor.max_torque (at stall)
     *   ω_max_reflected ≤ motor.max_speed
     *
     * Ref: Bolton (2015) §3.5
     */
    MotorSelection result;
    memset(&result, 0, sizeof(result));

    if (!motor || !cycle) {
        result.feasible = 0;
        return result;
    }

    result.motor = *motor;

    /* Collect torque/speed from cycle, reflect to motor */
    double T_rms = rms_torque(cycle);
    double T_peak = peak_torque_requirement(cycle);

    /* Reflect through gear */
    double N = gear ? gear->ratio : 1.0;
    double eta = gear ? gear->efficiency : 1.0;

    double J_load = 0.001;  /* conservatively estimated load inertia */
    result.reflected_J = J_load / (N * N);

    double T_rms_ref = T_rms / (N * eta);
    double T_peak_ref = T_peak / (N * eta);

    result.root_mean_sq_torque = T_rms_ref;
    result.peak_torque = T_peak_ref;
    result.continuous_speed = average_speed(cycle) * N;

    /* Feasibility checks */
    double torque_margin_cont = (motor->max_torque - T_rms_ref)
                                / motor->max_torque;
    double torque_margin_peak = (motor->max_torque - T_peak_ref)
                                / motor->max_torque;
    double speed_margin = (motor->max_speed - result.continuous_speed)
                          / motor->max_speed;

    result.torque_margin = torque_margin_cont;
    result.thermal_margin = torque_margin_cont;  /* simplified */

    result.feasible = (torque_margin_cont >= 0.0)
                      && (torque_margin_peak >= 0.0)
                      && (speed_margin >= 0.0)
                      && (result.continuous_speed <= motor->max_speed);

    return result;
}

/*===========================================================================*/
/* L5: Motor Thermal Model                                                   */
/*===========================================================================*/

double dc_motor_temperature(double I_rms, double R_a,
                            double T_ambient, double R_th, double C_th,
                            double t, double T_initial)
{
    /*
     * First-order thermal model:
     *
     *   C_th · dT/dt = P_loss - (T - T_amb) / R_th
     *
     * where P_loss = I_rms² · R_a (copper loss, dominant).
     *
     * Solution: T(t) = T_amb + P_loss·R_th
     *                  + [T_initial - T_amb - P_loss·R_th] · exp(-t/τ_th)
     *
     * where τ_th = R_th · C_th is the thermal time constant.
     *
     * Ref: Cetinkunt (2007) §6.4.1
     */
    double P_loss = I_rms * I_rms * R_a;
    double T_ss = T_ambient + P_loss * R_th;
    double tau_th = R_th * C_th;

    if (tau_th < 1e-30) return T_ss;

    return T_ss + (T_initial - T_ss) * exp(-t / tau_th);
}

double dc_motor_steady_state_temp(double I_rms, double R_a,
                                   double T_ambient, double R_th)
{
    /*
     * Steady-state temperature:
     *   T_ss = T_amb + I²·R_a · R_th
     */
    return T_ambient + I_rms * I_rms * R_a * R_th;
}

double resistance_temp_correction(double R0, double T0, double T,
                                   double alpha)
{
    /*
     * Temperature-dependent resistance for copper:
     *
     *   R(T) = R₀ · [1 + α · (T - T₀)]
     *
     * Copper: α ≈ 0.00393 K⁻¹ at T₀ = 20°C
     *
     * This is important because winding resistance increases
     * with temperature, reducing available torque.
     *
     * Ref: IEC 60034-1 (Rotating electrical machines — Rating and performance)
     */
    return R0 * (1.0 + alpha * (T - T0));
}

/*===========================================================================*/
/* L5: Stepper Motor Analysis                                                */
/*===========================================================================*/

double stepper_torque_at_rate(const StepperMotorParams *motor,
                               double step_rate, double V_supply)
{
    /*
     * Stepper motor pull-out torque curve approximation:
     *
     *   T(ω) = T_hold · ω_c / (ω_c + ω)
     *
     * where ω_c is the corner speed determined by L/R time constant:
     *   ω_c ∝ V_supply / (L_coil · steps_per_rev)
     *
     * At low speeds: T ≈ T_hold
     * At high speeds: T ∝ 1/ω (constant power region)
     *
     * Ref: Kenjo & Sugawara "Stepping Motors" (1994)
     */
    if (!motor) return 0.0;

    double L = motor->L_coil;
    double R = motor->R_coil;
    if (R < 1e-30) R = 1.0;
    double tau_elec = L / R;

    /* Corner frequency in steps/s */
    double f_c = V_supply / (tau_elec * motor->steps_per_rev * R);

    return motor->holding_torque * f_c / (f_c + step_rate);
}

double stepper_max_start_rate(const StepperMotorParams *motor,
                               double T_load, double J_load)
{
    /*
     * Maximum start-stop (no-ramp) step rate.
     *
     * Must satisfy: T_motor ≥ T_load + J_total · α_max
     *
     * For a single step, the motor must overcome load torque
     * and accelerate inertia within one step interval.
     *
     * Approximate:
     *   T_hold - T_detent - T_load = J_total · (step_angle · f²) / (2π)
     *
     * → f_start = √((T_hold - T_detent - T_load) · 2π / (J_total · step_angle))
     *
     * Ref: Morar "Stepper Motor Control" (2003)
     */
    if (!motor) return 0.0;

    double J_total = motor->J_r + J_load;
    if (J_total < 1e-30) return motor->max_slew_rate;

    double T_available = motor->holding_torque - motor->detent_torque - T_load;
    if (T_available <= 0.0) return 0.0;

    double f_start_sq = T_available * 2.0 * M_PI
                        / (J_total * motor->step_angle);
    if (f_start_sq < 0.0) return 0.0;

    return sqrt(f_start_sq);
}

double stepper_max_slew_rate(const StepperMotorParams *motor,
                               double T_load)
{
    /*
     * Maximum slew rate (with acceleration ramp).
     * Determined by torque-speed curve crossing the load line.
     *
     * Simplified: intersection of pull-out curve and T_load.
     *
     * T_load = T_hold · ω_c / (ω_c + ω_slew)
     * → ω_slew = ω_c · (T_hold/T_load - 1)
     *
     * Convert to steps/s: f_slew = ω_slew / step_angle
     */
    if (!motor || T_load <= 0.0) return motor->max_slew_rate;

    double L = motor->L_coil;
    double R = motor->R_coil;
    if (R < 1e-30) R = 1.0;
    double tau_elec = L / R;

    /* Corner speed in rad/s */
    double omega_c = 1.0 / (tau_elec * motor->steps_per_rev);

    double ratio = motor->holding_torque / T_load;
    if (ratio <= 1.0) return 0.0;

    double omega_slew = omega_c * (ratio - 1.0);
    double f_slew = omega_slew * motor->steps_per_rev / (2.0 * M_PI);

    return f_slew < motor->max_slew_rate ? f_slew : motor->max_slew_rate;
}

int stepper_ramp_time(const StepperMotorParams *motor,
                       double target_rate, double T_load, double J_total,
                       double *t_ramp)
{
    /*
     * Required linear acceleration ramp time.
     *
     * α = (T_available - T_load) / J_total   [rad/s²]
     * ω_target = target_rate · step_angle   [rad/s]
     * t_ramp = ω_target / α
     *
     * Returns 0 if feasible (T_available > T_load),
     *         -1 if target torque exceeds motor capability.
     */
    if (!motor || !t_ramp) return -1;

    double T_avail = motor->holding_torque - motor->detent_torque - T_load;
    if (T_avail <= 0.0) { *t_ramp = 1e30; return -1; }

    double alpha = T_avail / J_total;  /* rad/s² */
    double omega_target = target_rate * motor->step_angle;  /* rad/s */

    *t_ramp = omega_target / alpha;
    return 0;
}

/*===========================================================================*/
/* L6: Voice Coil Actuator                                                   */
/*===========================================================================*/

double vca_force(const VoiceCoilActuator *vca, double I)
{
    /*
     * Voice coil actuator force (Lorentz force principle):
     *
     *   F = K_f · I
     *
     * where K_f = B·L is the force constant.
     *
     * VCAs are used in HDD head positioning, fast steering mirrors,
     * and precision micro-positioning.
     *
     * Ref: Cetinkunt (2007) §6.2
     */
    if (!vca) return 0.0;
    return vca->K_f * I;
}

double vca_current_derivative(const VoiceCoilActuator *vca,
                               double I, double v, double V)
{
    /*
     * VCA electrical dynamics:
     *
     *   V = R·I + L·dI/dt + K_f·v   (back-EMF from motion)
     *
     *   dI/dt = (V - R·I - K_f·v) / L
     *
     * Note: In VCA, K_f is both force constant AND back-EMF constant
     * (same electromechanical reciprocity as K_t = K_e for DC motors).
     */
    if (!vca || vca->L_coil < 1e-30) return 0.0;
    return (V - vca->R_coil * I - vca->K_f * v) / vca->L_coil;
}

double vca_bandwidth(const VoiceCoilActuator *vca)
{
    /*
     * VCA bandwidth estimation:
     *
     * Electrical BW: ω_elec = R/L
     * Mechanical BW: ω_mech = K_f / √(M·K_spring)  (if spring-return)
     *                or ω_mech = K_f / (M·ω)  (moving coil, no spring)
     *
     * Overall BW ≈ min(ω_elec, ω_mech)
     *
     * For typical HDD VCA: BW ~500 Hz–2 kHz
     */
    if (!vca) return 0.0;

    double we = (vca->L_coil > 1e-30) ? vca->R_coil / vca->L_coil : 1e6;
    double wm;

    if (vca->K_spring > 0.0 && vca->M_moving > 0.0) {
        wm = vca->K_f / sqrt(vca->M_moving * vca->K_spring);
    } else if (vca->M_moving > 0.0) {
        wm = vca->K_f / vca->M_moving;  /* rough estimate */
    } else {
        wm = 1e6;
    }

    return we < wm ? we : wm;
}

/*===========================================================================*/
/* L8: Piezoelectric Actuator                                                */
/*===========================================================================*/

double piezo_free_stroke(const PiezoActuator *piezo, double V)
{
    /*
     * Piezoelectric free displacement (no external load):
     *
     *   ΔL = d₃₃ · n_layers · V
     *
     * where d₃₃ is the piezoelectric strain coefficient [m/V]
     * in the 3-direction (polarization axis).
     *
     * For PZT-5H: d₃₃ ≈ 593 pm/V
     * For PZT-8:  d₃₃ ≈ 225 pm/V
     *
     * Ref: Uchino "Piezoelectric Actuators and Ultrasonic Motors" (2010)
     */
    if (!piezo) return 0.0;
    return piezo->d33 * V;
}

double piezo_blocked_force(const PiezoActuator *piezo, double V)
{
    /*
     * Blocked force (zero displacement):
     *
     *   F_block = K_piezo · ΔL_free
     *           = K_piezo · d₃₃ · V
     *
     * where K_piezo is the actuator stiffness.
     *
     * Typical PZT block force: 100–10000 N
     */
    if (!piezo) return 0.0;
    return piezo->K_stiffness * piezo->d33 * V;
}

double piezo_stroke_under_load(const PiezoActuator *piezo,
                                double V, double K_load, double F_ext)
{
    /*
     * Piezo stroke under elastic load:
     *
     *   ΔL_eff = (K_p · ΔL_free - F_ext) / (K_p + K_load)
     *
     * where K_p = piezo stiffness, K_load = load stiffness.
     *
     * This is a series-spring model: the total displacement is shared
     * between compressing the load spring and the actuator's own
     * compliance.
     *
     * Derivation: Force balance
     *   K_p · (ΔL_free - ΔL_eff) = K_load · ΔL_eff + F_ext
     *   → ΔL_eff = (K_p · ΔL_free - F_ext) / (K_p + K_load)
     */
    if (!piezo) return 0.0;

    double K_p = piezo->K_stiffness;
    double dL_free = piezo_free_stroke(piezo, V);

    return (K_p * dL_free - F_ext) / (K_p + K_load);
}

double piezo_electrical_bandwidth(const PiezoActuator *piezo,
                                   double R_drive)
{
    /*
     * Piezo electrical bandwidth (RC limited):
     *
     *   f_bw = 1 / (2π · R_drive · C_piezo)
     *
     * Piezo actuators have significant capacitance (100 nF — 10 μF),
     * which limits the slew rate when driven through a finite
     * output resistance.
     *
     * High-speed piezo requires high-current drivers to charge/discharge
     * this capacitance rapidly.
     *
     * Ref: Uchino (2010) Ch.4
     */
    if (!piezo || R_drive < 1e-30 || piezo->C_cap < 1e-30)
        return 0.0;
    return 1.0 / (2.0 * M_PI * R_drive * piezo->C_cap);
}

double piezo_power_dissipation(const PiezoActuator *piezo,
                                double f, double V_pp, double tan_delta)
{
    /*
     * Piezoelectric power dissipation:
     *
     *   P = π · f · C · V_pp² · tan(δ)
     *
     * where tan(δ) is the dielectric loss tangent.
     *
     * This is the dielectric heating, which becomes significant
     * at high frequencies and large voltage swings.
     *
     * Typical PZT: tan(δ) ≈ 0.02 at 1 kHz
     *
     * Ref: Uchino (2010) §4.3
     */
    if (!piezo) return 0.0;
    return M_PI * f * piezo->C_cap * V_pp * V_pp * tan_delta;
}
