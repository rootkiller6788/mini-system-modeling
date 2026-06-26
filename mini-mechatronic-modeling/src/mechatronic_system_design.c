/**
 * @file mechatronic_system_design.c
 * @brief L5-L7 Mechatronic system design, integration, and analysis.
 *
 * Implements end-to-end design workflows: motor-gear-load matching,
 * cascade controller tuning, bandwidth allocation, friction compensation,
 * backlash modeling, and application-specific analyses (HDD, pick-and-place,
 * robot joint, antenna scanner, printer carriage).
 *
 * 9-Layer Coverage:
 *   L5 — Design algorithms (sizing, tuning, anti-windup)
 *   L6 — Canonical problems (servo axis, ball screw, robot joint)
 *   L7 — Real-world applications
 *
 * References:
 *   Shetty & Kolk (2011), Isermann (2005), Altintas (2012)
 *   Åström & Hägglund (2006), Armstrong-Hélouvry (1994)
 */

#include "mechatronic_system_design.h"
#include "mechatronic_state_space.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================*/
/* L5: Motor-Gear-Load Dynamics                                              */
/*===========================================================================*/

void compute_reflected_dynamics(const DCMotorParams *motor,
                                const GearTrain *gear,
                                double J_load, double B_load,
                                double *J_total, double *B_total)
{
    /*
     * Reflect all inertias and damping to motor shaft:
     *
     *   J_total = J_motor + J_gear_in + (J_gear_out + J_load)/N²
     *   B_total = B_motor + B_gear_in + (B_gear_out + B_load)/N²
     *
     * Derivation: Kinetic energy conservation
     *   KE = ½ J_m·ω_m² + ½ J_g1·ω_m² + ½ (J_g2 + J_L)·(ω_m/N)²
     *      = ½ [J_m + J_g1 + (J_g2 + J_L)/N²] · ω_m²
     *
     * Ref: Cetinkunt (2007) §7.2
     */
    if (!motor) return;

    double N = gear ? gear->ratio : 1.0;
    double N2 = N * N;

    double Jm = motor->J_m;
    double Bm = motor->B_m;
    double Jg1 = gear ? gear->J_input : 0.0;
    double Jg2 = gear ? gear->J_output : 0.0;
    double Bg = gear ? gear->B_damping : 0.0;

    if (J_total)
        *J_total = Jm + Jg1 + (Jg2 + J_load) / N2;
    if (B_total)
        *B_total = Bm + Bg + (Bg + B_load) / N2;
}

double inertia_ratio(double J_load_reflected, double J_motor)
{
    /*
     * Inertia ratio:
     *
     *   λ = J_load_reflected / J_motor
     *
     * Industry guidelines (from servo motor application notes):
     *   λ ≤ 1:1  — Direct drive, highest performance
     *   1:1 < λ ≤ 3:1 — Excellent control
     *   3:1 < λ ≤ 5:1 — Good, standard servo apps
     *   5:1 < λ ≤ 10:1 — Acceptable with careful tuning
     *   λ > 10:1 — Avoid or use advanced control (resonance issues)
     *
     * Ref: Yaskawa "Servo Motor Selection Guide" (2015)
     */
    if (J_motor < 1e-30) return 100.0;
    return J_load_reflected / J_motor;
}

const char *inertia_ratio_classification(double ratio)
{
    if (ratio <= 1.0)   return "Direct drive — Excellent";
    if (ratio <= 3.0)   return "Low inertia — Good";
    if (ratio <= 5.0)   return "Medium — Standard servo";
    if (ratio <= 10.0)  return "High — Needs tuning";
    return "Very high — Avoid or use advanced control";
}

/*===========================================================================*/
/* L5: Cascade Bandwidth Allocation                                          */
/*===========================================================================*/

void cascade_bandwidth_allocation(double omega_pos,
                                   double ratio_vel_pos,
                                   double ratio_cur_vel,
                                   double *omega_vel,
                                   double *omega_cur)
{
    /*
     * Cascade control bandwidth separation:
     *
     *   ω_cur ≥ ratio_cur_vel · ω_vel
     *   ω_vel ≥ ratio_vel_pos · ω_pos
     *
     * Typical ratios: 5–10 between each level.
     *
     * Why: Each inner loop must be fast enough that the outer loop
     * sees it as approximately unity gain within the outer bandwidth.
     * Otherwise, the slow inner loop would limit outer loop performance.
     *
     * Ref: Leonhard "Control of Electrical Drives" (2001) §8.5
     */
    if (omega_vel)
        *omega_vel = omega_pos * ratio_vel_pos;
    if (omega_cur)
        *omega_cur = (omega_vel ? *omega_vel : omega_pos * ratio_vel_pos)
                     * ratio_cur_vel;
}

double recommended_sample_rate(double bandwidth_hz)
{
    /*
     * Nyquist-based sampling guideline for digital control:
     *
     *   f_sample ≥ 20 · f_bandwidth
     *
     * Why 20x?:
     *   - Nyquist requires >2x (theoretical minimum)
     *   - Practical digital control needs 10-30x for:
     *     * Phase lag from ZOH ≤ 10° within bandwidth
     *     * Anti-aliasing filter roll-off
     *     * Quantization effects
     *
     * Note: 20x gives ~18° of ZOH phase lag at bandwidth.
     *
     * Ref: Franklin et al. "Digital Control of Dynamic Systems" (1998)
     */
    if (bandwidth_hz <= 0.0) return 0.0;
    return 20.0 * bandwidth_hz;
}

/*===========================================================================*/
/* L5: PI(D) Controllers with Anti-Windup                                    */
/*===========================================================================*/

void pi_antiwindup(double error, double dt,
                    double Kp, double Ki,
                    double u_min, double u_max,
                    double *integral, double *u)
{
    /*
     * PI controller with conditional integration anti-windup:
     *
     *   u_unsaturated = Kp·e + I_accum
     *   if u_unsaturated > u_max:
     *       u = u_max
     *       // Freeze integrator (back-calculation alternative)
     *   elif u_unsaturated < u_min:
     *       u = u_min
     *   else:
     *       u = u_unsaturated
     *       I_accum += Ki·e·dt
     *
     * Integrator windup occurs when the actuator saturates, causing
     * the integrator to accumulate large values that delay recovery
     * when the error changes sign (causing overshoot).
     *
     * Ref: Åström & Hägglund "Advanced PID Control" (2006) §3.5
     */
    if (!integral || !u || dt <= 0.0) return;

    double P_term = Kp * error;
    double u_unsat = P_term + *integral;

    if (u_unsat > u_max) {
        *u = u_max;
        /* Conditional integration: freeze integrator */
    } else if (u_unsat < u_min) {
        *u = u_min;
    } else {
        *u = u_unsat;
        *integral += Ki * error * dt;
    }
}

void pid_filtered_antiwindup(double error, double dt,
                              double Kp, double Ki, double Kd, double tau_f,
                              double u_min, double u_max,
                              double *integral, double *prev_error,
                              double *prev_d, double *u)
{
    /*
     * PID controller with filtered derivative and anti-windup:
     *
     *   P = Kp · e(t)
     *   I = I_prev + Ki·e·dt  (with anti-windup)
     *
     *   Filtered derivative (backward Euler discretization):
     *     D(s) = Kd·s / (τ_f·s + 1) · E(s)
     *     → D[k] = (τ_f/(τ_f+dt))·D[k-1] + (Kd/(τ_f+dt))·(e[k] - e[k-1])
     *
     * The derivative filter τ_f limits high-frequency gain to Kd/τ_f,
     * preventing noise amplification.
     *
     * Ref: Åström & Hägglund (2006) §3.4
     */
    if (!integral || !prev_error || !prev_d || !u || dt <= 0.0) return;

    double P_term = Kp * error;

    /* Filtered derivative */
    double alpha_d = tau_f / (tau_f + dt);
    double beta_d = Kd / (tau_f + dt);
    double d_term = alpha_d * (*prev_d) + beta_d * (error - *prev_error);
    *prev_d = d_term;
    *prev_error = error;

    double u_unsat = P_term + *integral + d_term;

    if (u_unsat > u_max) {
        *u = u_max;
    } else if (u_unsat < u_min) {
        *u = u_min;
    } else {
        *u = u_unsat;
        *integral += Ki * error * dt;
    }
}

/*===========================================================================*/
/* L5: Friction Model                                                        */
/*===========================================================================*/

double friction_torque(double omega, double T_c, double T_s,
                        double omega_s, double B)
{
    /*
     * Combined static + Coulomb + viscous + Stribeck friction model:
     *
     *   T_f(ω) = T_c·sgn(ω) + B·ω + (T_s - T_c)·exp(-|ω/ω_s|)·sgn(ω)
     *
     * Components:
     *   1. Coulomb: constant opposing torque (surface asperity contacts)
     *   2. Viscous: linear with velocity (lubrication)
     *   3. Stribeck: decreasing friction at low velocities (partial lubrication)
     *      — causes stick-slip and limit cycles
     *
     * The Stribeck effect is the transition from static to Coulomb friction
     * and is responsible for the negative damping region at low velocities.
     *
     * Ref: Armstrong-Hélouvry, B. "Control of Machines with Friction" (1994)
     *      §2.3 — Stribeck curve model
     */
    if (fabs(omega) < 1e-30) {
        /* At zero velocity, output zero (stiction handled separately
           by the controller or zero-velocity crossing logic) */
        return 0.0;
    }

    double sign = (omega > 0.0) ? 1.0 : -1.0;
    double stribeck_term = (T_s - T_c) * exp(-fabs(omega / omega_s));

    return T_c * sign + B * omega + stribeck_term * sign;
}

double friction_feedforward(double omega_ref, double T_c, double B)
{
    /*
     * Feedforward friction compensation:
     *
     *   T_ff = T_c · sgn(ω_ref) + B · ω_ref
     *
     * Only Coulomb and viscous terms are fed forward because:
     *   1. Stiction is hard to predict (depends on dwell time)
     *   2. Stribeck is velocity-dependent (feedforward at low velocity
     *      can cause instability if overestimated)
     *
     * The feedback controller handles the residual Stribeck/stiction.
     *
     * Ref: Armstrong-Hélouvry (1994) §4.2
     */
    double sign = (omega_ref > 0.0) ? 1.0 : ((omega_ref < 0.0) ? -1.0 : 0.0);
    return T_c * sign + B * omega_ref;
}

/*===========================================================================*/
/* L5: Backlash Model                                                        */
/*===========================================================================*/

void backlash_model(double theta_in, double theta_in_prev,
                    double theta_out_prev, double backlash,
                    double *theta_out)
{
    (void)theta_in_prev;    /* direction tracking reserved for future use */
    (void)theta_out_prev;   /* uses *theta_out as the current position ref */
    /*
     * Backlash (mechanical play) model:
     *
     * The output shaft only moves when the input moves far enough
     * to take up the backlash gap.
     *
     *   if θ_in - θ_out_prev > b/2:
     *       θ_out = θ_in - b/2   (driving forward)
     *   elif θ_out_prev - θ_in > b/2:
     *       θ_out = θ_in + b/2   (driving reverse)
     *   else:
     *       θ_out = θ_out_prev   (in dead zone, output stuck)
     *
     * where b = total backlash angle.
     *
     * Backlash in gear trains causes:
     *   - Positioning error (lost motion)
     *   - Limit cycles in feedback control
     *   - Impact loads during direction reversals
     *
     * Ref: Nordin & Gutman "Controlling Mechanical Systems with
     *      Backlash — A Survey" (2002) Automatica
     */
    if (!theta_out) return;

    double half_b = backlash * 0.5;

    if (theta_in - *theta_out > half_b) {
        /* Driving forward */
        *theta_out = theta_in - half_b;
    } else if (*theta_out - theta_in > half_b) {
        /* Driving reverse */
        *theta_out = theta_in + half_b;
    }
    /* else: in dead zone, output unchanged */
    (void)theta_in_prev;
}

void backlash_inverse(double theta_des, double theta_out_prev,
                       double backlash, double *theta_ff)
{
    /*
     * Backlash inverse model (for feedforward compensation):
     *
     *   θ_ff = θ_des + sgn(θ_des - θ_out_prev) · b/2
     *
     * This pre-distorts the command to account for the backlash gap,
     * so the actual output matches the desired position.
     *
     * Note: The inverse is only valid when the direction is known
     * (i.e., during continuous motion). Direction reversals still
     * require the backlash to be taken up.
     *
     * Ref: Nordin & Gutman (2002)
     */
    if (!theta_ff) return;

    double diff = theta_des - theta_out_prev;
    double sign = (diff > 0.0) ? 1.0 : ((diff < 0.0) ? -1.0 : 0.0);

    *theta_ff = theta_des + sign * backlash * 0.5;
}

/*===========================================================================*/
/* L6: Complete Servo Axis Design                                            */
/*===========================================================================*/

void design_servo_axis(const DCMotorParams *motor,
                        const GearTrain *gear,
                        const MechatronicLoad *load,
                        const RotarySensor *sensor,
                        double omega_pos_hz,
                        double ratio_vel_pos,
                        double ratio_cur_vel,
                        ServoAxis *axis)
{
    /*
     * Complete servo axis design procedure:
     *
     * 1. Reflect loads to motor shaft
     * 2. Build mechatronic subsystem
     * 3. Allocate cascade bandwidths
     * 4. Design current PI (fastest)
     * 5. Design velocity PI (middle)
     * 6. Design position P (slowest)
     * 7. Set saturation limits
     * 8. Copy sensor parameters
     *
     * This follows the standard industrial servo tuning procedure
     * (cascaded structure, innermost loop first).
     *
     * Ref: Shetty & Kolk "Mechatronics System Design" (2011) §8.1
     */
    if (!motor || !axis) return;

    /* Step 1-2: Build mechatronic subsystem */
    double J_total, B_total;
    compute_reflected_dynamics(motor, gear,
                               load ? load->J_load : 0.0,
                               load ? load->B_load : 0.0,
                               &J_total, &B_total);

    MechatronicSubsystem mech;
    memset(&mech, 0, sizeof(mech));
    mech.motor = *motor;
    if (gear) mech.gear = *gear;
    if (load) mech.load = *load;
    mech.J_total_reflected = J_total;
    mech.B_total_reflected = B_total;
    /* Natural frequency and damping */
    if (gear && gear->K_stiffness > 0.0) {
        double K_eq = gear->K_stiffness / (gear->ratio * gear->ratio);
        mech.omega_n = sqrt(K_eq / J_total);
        mech.zeta = B_total / (2.0 * sqrt(K_eq * J_total));
    }
    mech.tau_sys = (B_total > 0.0) ? J_total / B_total : motor->tau_m;

    *axis = servo_axis_create(mech, sensor ? *sensor
                               : rotary_sensor_create(SENSOR_ENCODER_INCREMENTAL,
                                                      1000, 100e3, 1e-6));

    /* Step 3: Bandwidth allocation */
    double omega_pos_rad = omega_pos_hz * 2.0 * M_PI;
    double omega_vel_rad, omega_cur_rad;
    cascade_bandwidth_allocation(omega_pos_rad, ratio_vel_pos,
                                  ratio_cur_vel, &omega_vel_rad, &omega_cur_rad);

    /* Step 4: Current loop PI */
    double Kp_cur, Ki_cur;
    design_current_pi(motor, omega_cur_rad, &Kp_cur, &Ki_cur);
    axis->Kp_current = Kp_cur;
    axis->Ki_current = Ki_cur;

    /* Step 5: Velocity loop PI */
    double tau_cur = 1.0 / omega_cur_rad;  /* current loop time constant */
    double Kp_vel, Ki_vel;
    design_velocity_pi(motor, J_total, B_total, omega_vel_rad,
                        tau_cur, &Kp_vel, &Ki_vel);
    axis->Kp_velocity = Kp_vel;
    axis->Ki_velocity = Ki_vel;

    /* Step 6: Position loop P */
    double Kp_pos;
    design_position_p(omega_pos_rad, &Kp_pos);
    axis->Kp_position = Kp_pos;

    /* Step 7: Saturation limits */
    axis->velocity_limit = motor->max_speed;
    axis->current_limit = motor->max_current;
    axis->max_accel = motor->max_torque / J_total;

    /* Step 8: Bandwidth reporting */
    axis->bandwidth_cur = omega_cur_rad / (2.0 * M_PI);
    axis->bandwidth_vel = omega_vel_rad / (2.0 * M_PI);
    axis->bandwidth_pos = omega_pos_hz;
}

/*===========================================================================*/
/* L6: Servo Axis Performance Analysis                                       */
/*===========================================================================*/

void servo_axis_step_analysis(const ServoAxis *axis,
                               double step_amplitude,
                               double dt, double T_sim,
                               double *rise_time,
                               double *settling_time,
                               double *overshoot,
                               double *steady_error)
{
    /*
     * Simulate servo axis step response and compute performance metrics.
     *
     * Uses simplified 2nd-order model: position loop + velocity loop.
     * More accurate simulation would use the full 3-loop model with
     * RK4 integration of state-space.
     *
     * Metrics:
     *   Rise time: time from 10% to 90% of final value
     *   Settling time: time to enter and stay within ±2% of final
     *   Overshoot: (peak - final)/final × 100%
     *   Steady-state error: |final - reference|
     *
     * Ref: Ogata (2010) §5.3
     */
    if (!axis) return;

    /* Simplified second-order approximation */
    double omega_n = 2.0 * M_PI * axis->bandwidth_pos;
    double zeta = 0.7;  /* typical servo tuning target */

    if (rise_time) {
        /* Rise time approximation for 2nd-order system:
           t_r ≈ (1.76·ζ³ - 0.417·ζ² + 1.039·ζ + 1) / ω_n
           (empirical fit for 10-90% rise time) */
        double z2 = zeta * zeta;
        double z3 = z2 * zeta;
        *rise_time = (1.76 * z3 - 0.417 * z2 + 1.039 * zeta + 1.0) / omega_n;
    }

    if (settling_time) {
        /* Settling time (2% criterion): t_s ≈ 4 / (ζ·ω_n) */
        *settling_time = 4.0 / (zeta * omega_n);
    }

    if (overshoot) {
        /* Percent overshoot: M_p = 100 · exp(-π·ζ/√(1-ζ²)) */
        if (zeta < 1.0 && zeta >= 0.0) {
            *overshoot = 100.0 * exp(-M_PI * zeta / sqrt(1.0 - zeta * zeta));
        } else {
            *overshoot = 0.0;  /* overdamped or critically damped */
        }
    }

    if (steady_error) {
        /* Steady-state error for P-only position loop with step input */
        /* e_ss = 1 / (1 + Kp) → for nonzero Kp with pure inertia, e_ss=0
           In practice, friction causes small steady error */
        *steady_error = 0.0;  /* P-controlled type-1 system → zero SSE */
    }

    /* Suppress unused parameter warnings */
    (void)step_amplitude;
    (void)dt;
    (void)T_sim;
}

/*===========================================================================*/
/* L6: Ball Screw Axis Sizing                                                */
/*===========================================================================*/

double ball_screw_required_torque(double M_mass, double pitch,
                                   double eta, double mu,
                                   double theta_tilt, double accel,
                                   double J_screw, double F_cutting)
{
    /*
     * Required motor torque for ball screw axis:
     *
     *   T_req = T_accel + T_friction + T_gravity + T_cutting
     *
     * where:
     *   T_accel   = [J_screw + M·(pitch/2π)²] · (accel · 2π/pitch)
     *              + F_cutting · pitch/(2π·η)
     *   T_friction = M·g·cos(θ)·μ · pitch/(2π·η)
     *   T_gravity  = M·g·sin(θ) · pitch/(2π·η)
     *   T_cutting  = F_cutting · pitch/(2π·η)
     *
     * Note: "cutting" force generalizes to any process force.
     *
     * Derivation: Power equality at screw:
     *   T_motor · ω = F_axial · v_axial / η
     *   v_axial = ω · pitch/(2π)
     *
     * Ref: Altintas "Manufacturing Automation" (2012) §4.3
     */
    double g = 9.80665;
    double r = pitch / (2.0 * M_PI);  /* transmission ratio */

    /* Translational inertia reflected to rotary */
    double J_equiv = J_screw + M_mass * r * r;
    double alpha_rotary = accel / r;  /* angular acceleration */

    /* Acceleration torque */
    double T_accel = J_equiv * alpha_rotary;

    /* Friction torque (Coulomb friction between nut and screw,
       plus guide friction from normal force) */
    double F_normal = M_mass * g * cos(theta_tilt);
    double F_friction = F_normal * mu;
    double T_friction = F_friction * r / eta;

    /* Gravity torque (inclined axis) */
    double F_gravity = M_mass * g * sin(theta_tilt);
    double T_gravity = F_gravity * r / eta;

    /* Cutting/process force torque */
    double T_cutting = F_cutting * r / eta;

    return T_accel + T_friction + T_gravity + T_cutting;
}

void linear_to_rotary_conversion(double v_linear, double F_linear,
                                  double pitch, double eta,
                                  double *omega_motor, double *T_motor)
{
    /*
     * Convert linear motion specs to rotary motor specs:
     *
     *   ω_motor = v_linear · (2π / pitch)
     *   T_motor = F_linear · (pitch / (2π·η))
     *
     * These are the fundamental kinematic and kinetic relationships
     * for leadscrew/ballscrew mechanisms.
     *
     * For ball screw: η ≈ 0.90–0.95
     * For lead screw: η ≈ 0.25–0.50
     *
     * Back-driving condition: T_backdrive = F·p·η_back/(2π)
     *   η_back ≈ η² for ball screws (back-driveable)
     *   η_back << 0 for lead screws (self-locking, generally)
     */
    if (pitch <= 0.0 || eta <= 0.0) return;

    double r = pitch / (2.0 * M_PI);

    if (omega_motor) *omega_motor = v_linear / r;
    if (T_motor) *T_motor = F_linear * r / eta;
}

/*===========================================================================*/
/* L6: Robot Joint Design                                                    */
/*===========================================================================*/

double robot_joint_max_torque(double M_link, double L_cm,
                               double J_link, double gear_ratio,
                               double accel_max, double J_motor,
                               double T_friction)
{
    /*
     * Maximum joint torque for single-DOF robot arm:
     *
     *   T_max = T_inertial + T_gravity + T_friction
     *
     *   T_inertial = (J_motor + J_link/N²) · N · α_max
     *              = J_motor·N·α_max + J_link·α_max/N
     *
     *   T_gravity = M_link · g · L_cm · cos(θ) / N
     *   (worst case: horizontal arm, cos(θ) = 1)
     *
     *   T_gravity_motor = M_link·g·L_cm / N
     *
     * Note: motor torque is AFTER the gearbox, so the reflected
     * load torque is divided by N.
     *
     * Ref: Spong et al. "Robot Modeling and Control" (2006) §7.3
     */
    double g = 9.80665;

    /* Inertial torque at motor shaft */
    double T_inertial = (J_motor + J_link / (gear_ratio * gear_ratio))
                        * gear_ratio * accel_max;

    /* Gravity torque at motor shaft (worst case: horizontal) */
    double T_gravity = M_link * g * L_cm / gear_ratio;

    /* Coulomb friction at motor shaft */
    double T_fric = T_friction;

    return T_inertial + T_gravity + T_fric;
}

/*===========================================================================*/
/* L6: Power and Energy Analysis                                             */
/*===========================================================================*/

void motion_profile_energy(const DutyCycle *cycle,
                            const DCMotorParams *motor,
                            const GearTrain *gear,
                            double J_total,
                            double *energy, double *peak_power,
                            double *avg_power)
{
    (void)J_total;  /* reserved for transient accel/decel energy terms */
    /*
     * Compute energy consumption for a motion profile.
     *
     * For each segment:
     *   P_elec(t) = V_a(t)·I_a(t)
     *
     * With V_a = R_a·I_a + K_e·ω and T_em = K_t·I_a = J·α + B·ω + T_load:
     *   P_elec = R_a·(T_em/K_t)² + K_e·ω·T_em/K_t
     *          = R_a/K_t² · T_em² + ω·T_em        (since K_e=K_t)
     *
     * First term: copper loss (heating)
     * Second term: mechanical power delivered to load
     *
     * Ref: Cetinkunt (2007) §6.4
     */
    if (!cycle || !motor) return;

    double total_energy = 0.0;
    double pk_power = 0.0;
    double total_time = 0.0;

    for (int i = 0; i < cycle->n_segments; i++) {
        double T = cycle->segments[i].torque;
        double omega = cycle->segments[i].speed;
        double dt = cycle->segments[i].duration;

        /* Apply gear ratio to reflect to motor */
        double N = gear ? gear->ratio : 1.0;
        double T_motor = T / (N * (gear ? gear->efficiency : 1.0));
        double omega_motor = omega * N;

        /* Electrical power */
        double P_loss = motor->R_a * T_motor * T_motor
                        / (motor->K_t * motor->K_t);
        double P_mech = omega_motor * T_motor;
        double P_elec = P_loss + P_mech;

        total_energy += P_elec * dt;
        total_time += dt;
        if (P_elec > pk_power) pk_power = P_elec;
    }

    if (energy) *energy = total_energy;
    if (peak_power) *peak_power = pk_power;
    if (avg_power && total_time > 0.0) *avg_power = total_energy / total_time;
}

double regenerative_energy(double J_total, double omega_high,
                            double omega_low, double eta_regen)
{
    /*
     * Recoverable kinetic energy during regenerative braking:
     *
     *   E_kin = ½ · J_total · (ω_high² - ω_low²)
     *   E_regen = η_regen · E_kin
     *
     * During deceleration, the motor acts as a generator, converting
     * kinetic energy back to electrical energy. The recoverable amount
     * depends on:
     *   - Motor/generator efficiency
     *   - Drive electronics efficiency (rectification)
     *   - Energy storage capability (battery/capacitor)
     *
     * Typical regeneration efficiency: 50-70%
     *
     * Ref: Lhomme et al. "Energy Savings in Electric Vehicles" (2005)
     */
    double delta_KE = 0.5 * J_total * (omega_high * omega_high
                     - omega_low * omega_low);
    if (delta_KE < 0.0) return 0.0;  /* only deceleration recovers */
    return eta_regen * delta_KE;
}

/*===========================================================================*/
/* L7: Application Analyses                                                  */
/*===========================================================================*/

double hdd_seek_time(double theta_travel, double K_t, double J_total,
                      double I_max, double omega_max)
{
    /*
     * Hard disk drive seek time (bang-bang profile):
     *
     * Maximum torque: T_max = K_t · I_max
     * Maximum acceleration: α_max = T_max / J_total
     *
     * If the travel is short (does not reach ω_max):
     *   T_seek = 2 · √(θ / α_max)
     *
     * If the travel hits speed limit (ω_max):
     *   t_accel = ω_max / α_max
     *   θ_accel = ½ α_max · t_accel² (covered in accel + decel)
     *   θ_cruise = θ_total - 2·θ_accel
     *   T_seek = 2·t_accel + θ_cruise/ω_max
     *
     * Modern HDD: α_max ≈ 5000 rad/s², ω_max ≈ 200 rad/s,
     *              seek time ≈ 5-15 ms.
     *
     * Ref: Abramovitch & Franklin "A Brief History of Disk Drive
     *      Control" (2002) IEEE Control Systems
     */
    double T_max = K_t * I_max;
    double alpha_max = T_max / J_total;

    if (alpha_max <= 0.0) return 1e30;

    /* Distance covered during acceleration to ω_max */
    double theta_accel = 0.5 * omega_max * omega_max / alpha_max;

    if (theta_travel <= 2.0 * theta_accel) {
        /* Short seek: triangular velocity profile */
        double theta_half = theta_travel / 2.0;
        double t_half = sqrt(2.0 * theta_half / alpha_max);
        return 2.0 * t_half;
    } else {
        /* Long seek: trapezoidal velocity profile */
        double t_accel = omega_max / alpha_max;
        double theta_cruise = theta_travel - 2.0 * theta_accel;
        double t_cruise = theta_cruise / omega_max;
        return 2.0 * t_accel + t_cruise;
    }
}

double pick_and_place_cycle_time(double distance, double v_cruise,
                                  double accel, double T_settle,
                                  double T_grip)
{
    /*
     * Pick-and-place cycle time for Cartesian robot:
     *
     * Trapezoidal velocity profile:
     *   t_accel = v_cruise / accel
     *   d_accel = ½·accel·t_accel² = v_cruise² / (2·accel)
     *
     * If distance ≤ 2·d_accel: triangular profile
     *   t_move = 2 · √(distance / accel)
     *
     * If distance > 2·d_accel: trapezoidal profile
     *   t_move = 2·t_accel + (distance - 2·d_accel)/v_cruise
     *
     * Total cycle: T = t_move + T_settle + T_grip
     *
     * Ref: Shetty & Kolk (2011) §9.3
     */
    if (accel <= 0.0) return 1e30;

    double d_accel = (v_cruise > 0.0)
                     ? v_cruise * v_cruise / (2.0 * accel) : 0.0;

    double t_move;
    if (v_cruise <= 0.0 || distance <= 2.0 * d_accel) {
        /* Triangular */
        t_move = 2.0 * sqrt(distance / accel);
    } else {
        /* Trapezoidal */
        double t_accel = v_cruise / accel;
        t_move = 2.0 * t_accel
                 + (distance - 2.0 * d_accel) / v_cruise;
    }

    return t_move + T_settle + T_grip;
}

double antenna_scan_time(double azimuth_range, double elevation_step,
                          double elevation_range, double omega_scan,
                          double T_settle)
{
    /*
     * Raster scan time for antenna pointing system:
     *
     * Number of elevation lines: n_lines = elevation_range / elevation_step
     * Time per azimuth sweep: t_sweep = azimuth_range / omega_scan
     * Settling time per line: T_settle (elevation step + lock)
     * Total time: T = n_lines · (t_sweep + T_settle)
     *
     * This pattern is used in:
     *   - Weather radar scanning
     *   - Satellite ground station tracking
     *   - Radio telescope survey patterns
     *
     * Ref: Skolnik "Radar Handbook" (2008) §7.3
     */
    if (elevation_step <= 0.0 || omega_scan <= 0.0) return 0.0;

    double n_lines = ceil(elevation_range / elevation_step);
    double t_sweep = azimuth_range / omega_scan;

    return n_lines * (t_sweep + T_settle);
}

void printer_carriage_analysis(double F_motor, double M_carriage,
                                double M_ink, double *f_position_bw)
{
    /*
     * Printer carriage dynamics analysis:
     *
     * Maximum acceleration: a_max = F_motor / (M_carriage + M_ink)
     *
     * For a raster scan with sinusoidal velocity profile at
     * carriage reversal, the position loop bandwidth is limited by
     * the achievable acceleration.
     *
     * For simple velocity-reversal S-curve:
     *   f_bw ≈ a_max / (2π · stroke)
     *
     * Or more generally for closed-loop bandwidth:
     *   f_bw ≈ √(F_motor / (M_total · encoder_resolution))
     *
     * Ref: Epson "PrecisionCore Printhead Technology" (2015)
     */
    if (!f_position_bw) return;

    double M_total = M_carriage + M_ink;
    if (M_total <= 0.0) { *f_position_bw = 0.0; return; }

    double a_max = F_motor / M_total;
    /* Bandwidth estimate based on achievable acceleration */
    /* f_bw ≈ √(a_max / (2π · resolution))  with 1μm resolution */
    double resolution = 1e-6;  /* 1 μm typical */
    *f_position_bw = sqrt(a_max / (2.0 * M_PI * resolution)) / (2.0 * M_PI);
}
