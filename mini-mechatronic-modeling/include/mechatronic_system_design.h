#ifndef MECHATRONIC_SYSTEM_DESIGN_H
#define MECHATRONIC_SYSTEM_DESIGN_H

/**
 * @file mechatronic_system_design.h
 * @brief L5-L7 Mechatronic System Design and Integration
 *
 * End-to-end design functions: motor-gear-load matching, cascade controller
 * tuning, bandwidth allocation, axis sizing, and system performance analysis.
 *
 * Key References:
 *  - Shetty & Kolk "Mechatronics System Design" (2011)
 *  - Isermann "Mechatronic Systems: Fundamentals" (2005)
 *  - Altintas "Manufacturing Automation" (2012)
 *
 * 9-Layer Coverage:
 *   L5 — Design methodologies, sizing algorithms
 *   L6 — Canonical problems: servo axis, CNC, robot joint
 *   L7 — Applications: pick-and-place, antenna positioning, HDD servo
 */

#include "mechatronic_definitions.h"
#include "electromechanical_coupling.h"
#include "actuator_modeling.h"
#include "sensor_modeling.h"

/*===========================================================================*/
/* L5: Motor-Gear-Load Matching                                              */
/*===========================================================================*/

/** Compute reflected inertia and damping through a transmission.
 *
 *  J_reflected = J_motor + J_gear_in + (J_gear_out + J_load) / N²
 *  B_reflected = B_motor + B_gear_in + (B_gear_out + B_load) / N²
 *
 *  @param motor    Motor parameters
 *  @param gear     Gear train parameters
 *  @param J_load   Load inertia [kg·m²]
 *  @param B_load   Load viscous damping [N·m·s/rad]
 *  @param J_total  Output: total reflected inertia [kg·m²]
 *  @param B_total  Output: total reflected damping [N·m·s/rad]
 *
 *  Ref: Cetinkunt (2007) Ch.7 */
void compute_reflected_dynamics(const DCMotorParams *motor,
                                const GearTrain *gear,
                                double J_load, double B_load,
                                double *J_total, double *B_total);

/** Inertia ratio check: J_load_reflected / J_motor
 *
 *  Industry guidelines:
 *    • 1:1 to 3:1 — Direct drive or low ratio, high performance
 *    • 3:1 to 5:1 — Standard servo applications
 *    • 5:1 to 10:1 — Acceptable, may need tuning
 *    • >10:1 — Avoid or use advanced control
 *
 *  @param J_load_reflected  Load inertia at motor shaft [kg·m²]
 *  @param J_motor           Motor rotor inertia [kg·m²]
 *  @return Inertia ratio (dimensionless) */
double inertia_ratio(double J_load_reflected, double J_motor);

/** Inertia ratio classification string
 *  @param ratio  Inertia ratio
 *  @return Static string: "Excellent", "Good", "Acceptable", or "Poor" */
const char *inertia_ratio_classification(double ratio);

/*===========================================================================*/
/* L5: Bandwidth Allocation in Cascaded Loops                                */
/*===========================================================================*/

/** Cascade control bandwidth allocation rule
 *
 *  Typical rule: ω_cur ≥ 5-10 × ω_vel ≥ 5-10 × ω_pos
 *
 *  Given the desired position loop bandwidth ω_pos, allocate velocity
 *  and current loop bandwidths using the specified ratios.
 *
 *  @param omega_pos       Desired position loop bandwidth [rad/s]
 *  @param ratio_vel_pos   ω_vel / ω_pos (typical: 5-10)
 *  @param ratio_cur_vel   ω_cur / ω_vel (typical: 5-10)
 *  @param omega_vel       Output: velocity loop bandwidth [rad/s]
 *  @param omega_cur       Output: current loop bandwidth [rad/s] */
void cascade_bandwidth_allocation(double omega_pos,
                                   double ratio_vel_pos,
                                   double ratio_cur_vel,
                                   double *omega_vel,
                                   double *omega_cur);

/** Compute controller sample rates from bandwidths
 *
 *  f_s ≥ 20 · f_bw (rule of thumb for digital control)
 *
 *  @param bandwidth_hz  Control loop bandwidth [Hz]
 *  @return Recommended sample rate [Hz] */
double recommended_sample_rate(double bandwidth_hz);

/*===========================================================================*/
/* L5: Anti-Windup for PI Controllers                                        */
/*===========================================================================*/

/** PI controller with clamping anti-windup
 *
 *  Implements: u = Kp·e + Ki·∫e dt
 *  with back-calculation anti-windup:
 *    if u > u_max: u = u_max; integral term frozen
 *    if u < u_min: u = u_min; integral term frozen
 *
 *  @param error            Current error = reference - measurement
 *  @param dt               Timestep [s]
 *  @param Kp               Proportional gain
 *  @param Ki               Integral gain
 *  @param u_min            Lower saturation limit
 *  @param u_max            Upper saturation limit
 *  @param integral         Integral accumulator (in/out)
 *  @param u                Output: control signal
 *
 *  Ref: Åström & Hägglund "Advanced PID Control" (2006) */
void pi_antiwindup(double error, double dt,
                    double Kp, double Ki,
                    double u_min, double u_max,
                    double *integral, double *u);

/** PID controller with derivative filtering and anti-windup
 *
 *  Filtered derivative: D(s) = Kd·s / (τ_f·s + 1)
 *  Discretized via backward Euler.
 *
 *  @param error       Current error
 *  @param dt          Timestep [s]
 *  @param Kp, Ki, Kd  PID gains
 *  @param tau_f       Derivative filter time constant [s]
 *  @param u_min, u_max Saturation limits
 *  @param integral    Integral accumulator (in/out)
 *  @param prev_error  Previous error (in/out, for derivative)
 *  @param prev_d      Previous filtered derivative (in/out)
 *  @param u           Output: control signal */
void pid_filtered_antiwindup(double error, double dt,
                              double Kp, double Ki, double Kd, double tau_f,
                              double u_min, double u_max,
                              double *integral, double *prev_error,
                              double *prev_d, double *u);

/*===========================================================================*/
/* L5: Friction Compensation                                                 */
/*===========================================================================*/

/** Compute friction torque using Coulomb + viscous + Stribeck model
 *
 *  T_f(ω) = T_c · sign(ω) + B·ω + (T_s - T_c)·exp(-|ω/ω_s|) · sign(ω)
 *
 *  where T_c = Coulomb friction, T_s = static friction (stiction),
 *  ω_s = Stribeck velocity, B = viscous damping
 *
 *  @param omega     Angular velocity [rad/s]
 *  @param T_c       Coulomb friction [N·m]
 *  @param T_s       Static friction / stiction torque [N·m]
 *  @param omega_s   Stribeck velocity [rad/s]
 *  @param B         Viscous damping coefficient [N·m·s/rad]
 *  @return Friction torque [N·m]
 *
 *  Ref: Armstrong-Hélouvry "Control of Machines with Friction" (1994) */
double friction_torque(double omega, double T_c, double T_s,
                        double omega_s, double B);

/** Feedforward friction compensation value
 *
 *  Used in feedforward to cancel predictable friction.
 *  Typically only Coulomb + viscous fed forward (stiction is hard to predict).
 *
 *  @param omega_ref    Reference velocity [rad/s]
 *  @param T_c          Coulomb friction estimate [N·m]
 *  @param B            Viscous friction estimate [N·m·s/rad]
 *  @return Feedforward torque [N·m] */
double friction_feedforward(double omega_ref, double T_c, double B);

/*===========================================================================*/
/* L5: Backlash Model and Compensation                                       */
/*===========================================================================*/

/** Backlash dead-zone model
 *
 *  Output position lags input by backlash angle.
 *
 *  @param theta_in        Input shaft angle [rad]
 *  @param theta_in_prev   Previous input shaft angle [rad] (for direction)
 *  @param theta_out_prev  Current output shaft angle [rad]
 *  @param backlash        Backlash angle [rad]
 *  @param theta_out       Output: output shaft angle [rad]
 *
 *  Ref: Nordin & Gutman "Controlling Mechanical Systems with Backlash" (2002) */
void backlash_model(double theta_in, double theta_in_prev,
                    double theta_out_prev, double backlash,
                    double *theta_out);

/** Backlash inverse compensation (for feedforward)
 *
 *  Pre-distorts input to cancel backlash effect.
 *
 *  @param theta_des      Desired output angle [rad]
 *  @param theta_out_prev Previous actual output [rad]
 *  @param backlash       Estimated backlash [rad]
 *  @param theta_ff       Output: feedforward command [rad] */
void backlash_inverse(double theta_des, double theta_out_prev,
                       double backlash, double *theta_ff);

/*===========================================================================*/
/* L6: Complete Servo Axis Design                                            */
/*===========================================================================*/

/** Design a complete servo axis (motor + gear + load + sensors + controllers)
 *
 *  Steps:
 *  1. Reflect load dynamics to motor shaft
 *  2. Compute total inertia and damping
 *  3. Allocate cascade bandwidths
 *  4. Design current PI controller
 *  5. Design velocity PI controller
 *  6. Design position P controller
 *  7. Set saturation limits
 *  8. Estimate achievable performance
 *
 *  @param motor          Motor parameters
 *  @param gear           Transmission
 *  @param load           Load parameters
 *  @param sensor         Position/velocity sensor
 *  @param omega_pos_hz   Desired position loop bandwidth [Hz]
 *  @param ratio_vel_pos  Velocity/position bandwidth ratio (typ. 5)
 *  @param ratio_cur_vel  Current/velocity bandwidth ratio (typ. 5)
 *  @param axis           Output: fully configured servo axis
 *
 *  Ref: Shetty & Kolk (2011) Ch.8 */
void design_servo_axis(const DCMotorParams *motor,
                        const GearTrain *gear,
                        const MechatronicLoad *load,
                        const RotarySensor *sensor,
                        double omega_pos_hz,
                        double ratio_vel_pos,
                        double ratio_cur_vel,
                        ServoAxis *axis);

/** Evaluate servo axis step response performance
 *
 *  Computes rise time, settling time, overshoot from simulated step.
 *
 *  @param axis            Servo axis configuration
 *  @param step_amplitude  Position step size [rad]
 *  @param dt              Simulation timestep [s]
 *  @param T_sim           Simulation duration [s]
 *  @param rise_time       Output: 10%-90% rise time [s]
 *  @param settling_time   Output: 2% settling time [s]
 *  @param overshoot       Output: percent overshoot [%]
 *  @param steady_error    Output: steady-state error [rad] */
void servo_axis_step_analysis(const ServoAxis *axis,
                               double step_amplitude,
                               double dt, double T_sim,
                               double *rise_time,
                               double *settling_time,
                               double *overshoot,
                               double *steady_error);

/*===========================================================================*/
/* L6: Ball Screw Axis Sizing                                                */
/*===========================================================================*/

/** Compute required motor torque for ball screw axis
 *
 *  T_required = T_accel + T_friction + T_gravity + T_cutting
 *
 *  T_accel = (J_total)·α
 *  T_friction = (M·g·cos(θ)·μ·pitch)/(2π·η)
 *  T_gravity = (M·g·sin(θ)·pitch)/(2π·η)
 *
 *  @param M_mass     Moving mass [kg]
 *  @param pitch      Screw lead [m]
 *  @param eta        Screw efficiency
 *  @param mu         Friction coefficient
 *  @param theta_tilt Tilt angle from horizontal [rad]
 *  @param accel      Acceleration [m/s²]
 *  @param J_screw    Screw inertia [kg·m²]
 *  @param F_cutting  Cutting/process force [N]
 *  @return Required motor torque [N·m] */
double ball_screw_required_torque(double M_mass, double pitch,
                                   double eta, double mu,
                                   double theta_tilt, double accel,
                                   double J_screw, double F_cutting);

/** Convert linear motion requirements to rotary motor requirements
 *
 *  ω_motor = v_linear · 2π / pitch
 *  T_motor = F_linear · pitch / (2π · η)
 *
 *  @param v_linear     Linear velocity [m/s]
 *  @param F_linear     Linear force [N]
 *  @param pitch        Screw lead [m]
 *  @param eta          Screw efficiency
 *  @param omega_motor  Output: motor speed [rad/s]
 *  @param T_motor      Output: motor torque [N·m] */
void linear_to_rotary_conversion(double v_linear, double F_linear,
                                  double pitch, double eta,
                                  double *omega_motor, double *T_motor);

/*===========================================================================*/
/* L6: Robot Joint Design                                                    */
/*===========================================================================*/

/** Compute maximum joint torque for a single-DOF robot arm
 *
 *  T_max = (J_motor + J_gear + J_link/N²) · α_max
 *          + T_gravity_max + T_friction_max
 *
 *  where T_gravity_max = M_link·g·L_cm (horizontal position, worst case)
 *
 *  @param M_link       Link mass [kg]
 *  @param L_cm         Distance from joint to link center of mass [m]
 *  @param J_link       Link inertia about joint [kg·m²]
 *  @param gear_ratio   Gear ratio N
 *  @param accel_max    Maximum angular acceleration [rad/s²]
 *  @param J_motor      Motor inertia [kg·m²]
 *  @param T_friction   Estimated friction [N·m]
 *  @return Maximum required torque at motor shaft [N·m] */
double robot_joint_max_torque(double M_link, double L_cm,
                               double J_link, double gear_ratio,
                               double accel_max, double J_motor,
                               double T_friction);

/*===========================================================================*/
/* L6: Power and Energy Analysis                                             */
/*===========================================================================*/

/** Motor electrical power for a motion profile
 *
 *  P(t) = V_a(t)·I_a(t)
 *
 *  Summed over the duty cycle to compute energy consumption.
 *
 *  @param cycle      Duty cycle
 *  @param motor      Motor parameters
 *  @param gear       Gear train
 *  @param J_total    Total reflected inertia [kg·m²]
 *  @param energy     Output: total energy [J]
 *  @param peak_power Output: peak power [W]
 *  @param avg_power  Output: average power [W] */
void motion_profile_energy(const DutyCycle *cycle,
                            const DCMotorParams *motor,
                            const GearTrain *gear,
                            double J_total,
                            double *energy, double *peak_power,
                            double *avg_power);

/** Regenerative braking energy recovery estimate
 *
 *  When decelerating, the motor acts as a generator.
 *  E_regen = η_regen · ½ · J_total · (ω_high² - ω_low²)
 *
 *  @param J_total      Total reflected inertia [kg·m²]
 *  @param omega_high   Initial speed [rad/s]
 *  @param omega_low    Final speed [rad/s]
 *  @param eta_regen    Regeneration efficiency (0-1)
 *  @return Recoverable energy [J] */
double regenerative_energy(double J_total, double omega_high,
                            double omega_low, double eta_regen);

/*===========================================================================*/
/* L7: Application-Specific Analysis                                         */
/*===========================================================================*/

/** Hard disk drive (HDD) seek time estimation
 *
 *  For a given angular travel, compute minimum seek time under current
 *  and voltage constraints using bang-bang acceleration profile.
 *
 *  T_seek = √(2·θ_travel / α_max) + ω_max/α_max (if speed limit hit)
 *
 *  @param theta_travel  Angular travel [rad]
 *  @param K_t           Torque constant [N·m/A]
 *  @param J_total       Total inertia [kg·m²]
 *  @param I_max         Max current [A]
 *  @param omega_max     Max velocity [rad/s]
 *  @return Seek time [s] */
double hdd_seek_time(double theta_travel, double K_t, double J_total,
                      double I_max, double omega_max);

/** Pick-and-place cycle time for a Cartesian robot
 *
 *  Trapezoidal velocity profile: accel → cruise → decel
 *  T_cycle = T_move + T_settle + T_pick_place
 *
 *  @param distance      Move distance [m]
 *  @param v_cruise      Cruise velocity [m/s]
 *  @param accel         Acceleration [m/s²]
 *  @param T_settle      Settling time after move [s]
 *  @param T_grip        Gripper action time [s]
 *  @return Cycle time per move [s] */
double pick_and_place_cycle_time(double distance, double v_cruise,
                                  double accel, double T_settle,
                                  double T_grip);

/** Scan pattern time for antenna pointing system (azimuth scan)
 *
 *  Raster scan with constant angular velocity and step.
 *
 *  @param azimuth_range  Total azimuth range [rad]
 *  @param elevation_step Elevation step per scan line [rad]
 *  @param elevation_range Total elevation range [rad]
 *  @param omega_scan     Scanning angular velocity [rad/s]
 *  @param T_settle       Settling time per line [s]
 *  @return Total scan time [s] */
double antenna_scan_time(double azimuth_range, double elevation_step,
                          double elevation_range, double omega_scan,
                          double T_settle);

/** Print-head carriage inertia and bandwidth analysis
 *
 *  Computes achievable acceleration given motor force and carriage mass.
 *
 *  @param F_motor        Motor continuous force [N]
 *  @param M_carriage     Carriage mass [kg]
 *  @param M_ink          Ink/payload mass [kg]
 *  @param f_position_bw  Output: position loop bandwidth [Hz] */
void printer_carriage_analysis(double F_motor, double M_carriage,
                                double M_ink, double *f_position_bw);

#endif /* MECHATRONIC_SYSTEM_DESIGN_H */
