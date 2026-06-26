#ifndef ACTUATOR_MODELING_H
#define ACTUATOR_MODELING_H

/**
 * @file actuator_modeling.h
 * @brief L1-L7 Actuator Selection, Modeling, and Sizing
 *
 * Covers DC motor selection, thermal analysis, gear ratio optimization,
 * stepper motor torque-speed curves, voice coil actuators, and
 * piezoelectric actuators.
 *
 * Key References:
 *  - Bolton "Mechatronics" (2015)
 *  - Cetinkunt "Mechatronics" (2007)
 *  - Morar "Stepper Motor Control" (2003)
 *  - Uchino "Piezoelectric Actuators" (2010)
 *
 * 9-Layer Coverage:
 *   L1 — Actuator type definitions, parameter structs
 *   L2 — Selection principles, torque-speed envelopes
 *   L5 — Selection algorithms, thermal models
 *   L7 — Application constraints: robot joint, CNC axis, antenna pointing
 */

#include "mechatronic_definitions.h"
#include "electromechanical_coupling.h"

/*===========================================================================*/
/* L2: Motor Duty Cycle Classification                                       */
/*===========================================================================*/

typedef enum {
    DUTY_CONTINUOUS = 0,  /**< S1: Continuous running duty               */
    DUTY_SHORT_TIME  = 1, /**< S2: Short-time duty                       */
    DUTY_INTERMITTENT = 2 /**< S3: Intermittent periodic duty             */
} MotorDutyType;

/** Define a motion profile segment */
typedef struct {
    double torque;        /**< Required torque [N·m]                      */
    double speed;         /**< Speed during segment [rad/s]               */
    double duration;      /**< Segment duration [s]                       */
} MotionSegment;

/** Complete duty cycle for motor sizing */
typedef struct {
    int          n_segments;
    MotionSegment segments[20];  /**< Up to 20 segments                  */
    double       cycle_time;     /**< Total cycle time [s]               */
    MotorDutyType duty_type;
} DutyCycle;

/*===========================================================================*/
/* L5: Motor Sizing and Selection                                            */
/*===========================================================================*/

/** Compute RMS torque for a duty cycle
 *  T_rms = √(Σ T_i² · Δt_i / T_cycle)
 *
 *  The RMS torque determines the continuous torque requirement.
 *  The motor's rated continuous torque must exceed T_rms.
 *
 *  @param cycle  Duty cycle specification
 *  @return RMS torque [N·m]
 *
 *  Ref: Cetinkunt "Mechatronics" (2007) Ch.7 */
double rms_torque(const DutyCycle *cycle);

/** Compute peak torque requirement from duty cycle
 *  @param cycle  Duty cycle specification
 *  @return Maximum absolute torque [N·m] */
double peak_torque_requirement(const DutyCycle *cycle);

/** Compute average speed over duty cycle (weighted by time)
 *  @param cycle  Duty cycle specification
 *  @return Average speed [rad/s] */
double average_speed(const DutyCycle *cycle);

/** Reflect load inertia and torque through a gear train to motor shaft
 *
 *  J_reflected = J_load / N²
 *  T_reflected = T_load / N (ideal) or T_load / (N·η) for driving
 *
 *  @param J_load     Load inertia [kg·m²]
 *  @param T_load     Load torque [N·m]
 *  @param gear       Gear train parameters
 *  @param direction  0 = motor driving load, 1 = load back-driving motor
 *  @param J_out      Output: reflected inertia [kg·m²]
 *  @param T_out      Output: reflected torque [N·m] */
void reflect_load_through_gear(double J_load, double T_load,
                               const GearTrain *gear, int direction,
                               double *J_out, double *T_out);

/** Optimal gear ratio for minimum motor torque (inertia matching)
 *
 *  Given motor inertia J_m and load inertia J_L,
 *  optimal gear ratio N_opt = √(J_L / J_m)
 *  minimizes the peak acceleration torque.
 *
 *  @param J_m    Motor rotor inertia [kg·m²]
 *  @param J_L    Load inertia [kg·m²]
 *  @return Optimal gear ratio
 *
 *  Ref: Cetinkunt (2007) §7.3 */
double optimal_gear_ratio_inertia_match(double J_m, double J_L);

/** Select DC motor for a given motion profile
 *
 *  Checks: (1) RMS torque ≤ rated continuous torque
 *          (2) Peak torque ≤ stall torque
 *          (3) Max speed ≤ rated speed
 *          (4) Thermal capacity sufficient
 *
 *  @param motor  Candidate motor
 *  @param gear   Transmission
 *  @param cycle  Duty cycle definition
 *  @return MotorSelection with feasibility flag
 *
 *  Ref: Bolton "Mechatronics" (2015) Ch.3 */
MotorSelection select_dc_motor(const DCMotorParams *motor,
                               const GearTrain *gear,
                               const DutyCycle *cycle);

/*===========================================================================*/
/* L5: Motor Thermal Analysis                                                */
/*===========================================================================*/

/** First-order thermal model of DC motor
 *
 *  dT/dt = (P_loss - (T - T_amb) / R_th) / C_th
 *  where P_loss = I²·R_a (copper loss) + core losses
 *
 *  Steady state: T_ss = T_amb + P_loss · R_th
 *
 *  @param motor        Motor parameters
 *  @param I_rms        RMS armature current [A]
 *  @param T_ambient    Ambient temperature [°C]
 *  @param R_th         Thermal resistance [K/W]
 *  @param C_th         Thermal capacitance [J/K]
 *  @param t            Time [s]
 *  @param T_initial    Initial temperature [°C]
 *  @return Temperature at time t [°C] */
double dc_motor_temperature(double I_rms, double R_a,
                            double T_ambient, double R_th, double C_th,
                            double t, double T_initial);

/** Steady-state winding temperature
 *  @param I_rms     RMS current [A]
 *  @param R_a       Resistance at reference temp [Ω]
 *  @param T_ambient Ambient temperature [°C]
 *  @param R_th      Thermal resistance [K/W]
 *  @return Steady-state temperature rise [°C] */
double dc_motor_steady_state_temp(double I_rms, double R_a,
                                   double T_ambient, double R_th);

/** Resistance temperature correction
 *  R(T) = R₀ · [1 + α·(T - T₀)]
 *  Copper: α ≈ 0.00393 K⁻¹ at 20°C
 *
 *  @param R0      Resistance at reference temperature [Ω]
 *  @param T0      Reference temperature [°C]
 *  @param T       Actual temperature [°C]
 *  @param alpha   Temperature coefficient of resistance [1/K]
 *  @return Corrected resistance [Ω] */
double resistance_temp_correction(double R0, double T0, double T,
                                   double alpha);

/*===========================================================================*/
/* L5: Stepper Motor Analysis                                                */
/*===========================================================================*/

/** Stepper motor torque vs. step rate (pull-out torque curve)
 *
 *  Approximate model: T_pullout(speed) = T_hold · ω_c / (ω_c + speed)
 *  where ω_c is the corner frequency of the torque-speed curve.
 *
 *  @param motor       Stepper motor parameters
 *  @param step_rate   Step rate [steps/s]
 *  @param V_supply    Supply voltage [V]
 *  @return Available torque at given step rate [N·m]
 *
 *  Ref: Morar "Stepper Motor Control" (2003) */
double stepper_torque_at_rate(const StepperMotorParams *motor,
                               double step_rate, double V_supply);

/** Maximum start-stop step rate (no ramping)
 *  Rate at which the motor can start without losing steps under load.
 *
 *  @param motor        Stepper parameters
 *  @param T_load       Load torque [N·m]
 *  @param J_load       Load inertia [kg·m²]
 *  @return Maximum start-stop rate [steps/s] */
double stepper_max_start_rate(const StepperMotorParams *motor,
                               double T_load, double J_load);

/** Maximum slew rate with acceleration ramp
 *  @param motor        Stepper parameters
 *  @param T_load       Load torque [N·m]
 *  @return Maximum slew rate [steps/s] */
double stepper_max_slew_rate(const StepperMotorParams *motor,
                               double T_load);

/** Required acceleration ramp for stepper motor
 *  Linear ramp: α = (ω_target - ω_start) / t_ramp
 *
 *  @param motor        Stepper parameters
 *  @param target_rate  Target step rate [steps/s]
 *  @param T_load       Load torque [N·m]
 *  @param J_total      Total inertia [kg·m²]
 *  @param t_ramp       Output: required ramp time [s]
 *  @return 0 if feasible, -1 if target exceeds motor capability */
int stepper_ramp_time(const StepperMotorParams *motor,
                       double target_rate, double T_load, double J_total,
                       double *t_ramp);

/*===========================================================================*/
/* L6: Voice Coil Actuator (VCA)                                             */
/*===========================================================================*/

/** Voice coil actuator parameters */
typedef struct {
    double K_f;           /**< Force constant [N/A]                        */
    double R_coil;        /**< Coil resistance [Ω]                        */
    double L_coil;        /**< Coil inductance [H]                        */
    double stroke;        /**< Total stroke [m]                            */
    double M_moving;      /**< Moving mass [kg]                            */
    double K_spring;      /**< Return spring stiffness [N/m] (0 if none)   */
    double max_current;   /**< Max continuous current [A]                  */
    double max_force;     /**< Peak force [N]                              */
} VoiceCoilActuator;

/** Voice coil force: F = K_f · I
 *  @param vca   VCA parameters
 *  @param I     Coil current [A]
 *  @return Output force [N] */
double vca_force(const VoiceCoilActuator *vca, double I);

/** Voice coil electrical dynamics: V = R·I + L·dI/dt + K_f·v
 *  (back-EMF term from motion: K_f·v)
 *
 *  @param vca   VCA parameters
 *  @param I     Current [A]
 *  @param v     Velocity [m/s]
 *  @param V     Applied voltage [V]
 *  @return dI/dt [A/s] */
double vca_current_derivative(const VoiceCoilActuator *vca,
                               double I, double v, double V);

/** VCA bandwidth estimation
 *  ω_bw ≈ min( R/L,  K_f / √(M·K_spring) )
 *
 *  @param vca   VCA parameters
 *  @return Bandwidth [rad/s] */
double vca_bandwidth(const VoiceCoilActuator *vca);

/*===========================================================================*/
/* L8: Piezoelectric Actuator                                                */
/*===========================================================================*/

typedef struct {
    double d33;           /**< Piezo strain coefficient [m/V] or [C/N]    */
    double K_stiffness;   /**< Stiffness [N/m]                            */
    double C_cap;         /**< Capacitance [F]                            */
    double max_voltage;   /**< Maximum operating voltage [V]               */
    double max_stroke;    /**< Free stroke at max voltage [m]              */
    double blocking_force;/**< Blocked force at max voltage [N]            */
    double resonance_freq;/**< First mechanical resonance [Hz]             */
    double M_eff;         /**< Effective mass [kg]                         */
} PiezoActuator;

/** Piezo free stroke: ΔL = d33 · n_layers · V
 *  @param piezo  Piezo parameters
 *  @param V      Applied voltage [V]
 *  @return Free displacement [m] */
double piezo_free_stroke(const PiezoActuator *piezo, double V);

/** Piezo blocked force: F_block = K_stiffness · ΔL_free
 *  @param piezo  Piezo parameters
 *  @param V      Applied voltage [V]
 *  @return Blocked force [N] */
double piezo_blocked_force(const PiezoActuator *piezo, double V);

/** Piezo actuator with elastic load
 *  Effective stroke: ΔL_eff = ΔL_free · K_piezo / (K_piezo + K_load)
 *  @param piezo   Piezo parameters
 *  @param V       Applied voltage [V]
 *  @param K_load  Load stiffness [N/m]
 *  @param F_ext   External force [N]
 *  @return Actual displacement [m] */
double piezo_stroke_under_load(const PiezoActuator *piezo,
                                double V, double K_load, double F_ext);

/** Piezo electrical bandwidth (RC limited)
 *  f_bw = 1 / (2π · R_drive · C_piezo)
 *  @param piezo    Piezo parameters
 *  @param R_drive  Driver output resistance [Ω]
 *  @return Electrical bandwidth [Hz] */
double piezo_electrical_bandwidth(const PiezoActuator *piezo,
                                   double R_drive);

/** Piezo power consumption at frequency
 *  P = π · f · C · V_pp² · tan(δ)
 *  @param piezo   Piezo parameters
 *  @param f       Operating frequency [Hz]
 *  @param V_pp    Peak-to-peak voltage [V]
 *  @param tan_delta Dielectric loss tangent
 *  @return Power dissipation [W] */
double piezo_power_dissipation(const PiezoActuator *piezo,
                                double f, double V_pp, double tan_delta);

#endif /* ACTUATOR_MODELING_H */
