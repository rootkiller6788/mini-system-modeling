#ifndef MECHATRONIC_DEFINITIONS_H
#define MECHATRONIC_DEFINITIONS_H

/**
 * @file mechatronic_definitions.h
 * @brief L1 Core Definitions — Mechatronic Systems
 *
 * Defines all fundamental data structures for mechatronic modeling:
 * motors, gears, sensors, actuators, loads, and integrated subsystems.
 *
 * Reference:
 *  - Bolton, W. "Mechatronics: Electronic Control Systems" (2015)
 *  - Cetinkunt, S. "Mechatronics" (2007)
 *  - Shetty & Kolk "Mechatronics System Design" (2011)
 *
 * 9-Layer Coverage: L1 (Definitions) — Complete
 */

#include <stddef.h>
#include <stdint.h>

/*---------------------------------------------------------------------------*/
/* L1.1: DC Motor Parameters (brushed, permanent magnet)                     */
/*---------------------------------------------------------------------------*/

/** Physical parameters for a brushed DC motor */
typedef struct {
    double R_a;           /**< Armature resistance [Ω]                    */
    double L_a;           /**< Armature inductance [H]                    */
    double K_e;           /**< Back-EMF constant [V·s/rad]               */
    double K_t;           /**< Torque constant [N·m/A], = K_e in SI      */
    double J_m;           /**< Rotor inertia [kg·m²]                     */
    double B_m;           /**< Viscous damping [N·m·s/rad]               */
    double tau_e;         /**< Electrical time constant L_a/R_a [s]       */
    double tau_m;         /**< Mechanical time constant J_m/B_m [s]      */
    double max_voltage;   /**< Rated voltage [V]                          */
    double max_current;   /**< Rated current [A]                          */
    double max_speed;     /**< Max no-load speed [rad/s]                  */
    double max_torque;    /**< Stall torque at rated current [N·m]        */
    double P_rated;       /**< Rated power [W]                            */
    double pole_pairs;    /**< Number of pole pairs                       */
} DCMotorParams;

/*---------------------------------------------------------------------------*/
/* L1.2: Brushless DC Motor (BLDC) Parameters                                */
/*---------------------------------------------------------------------------*/

typedef struct {
    double R_ph;          /**< Phase resistance [Ω]                       */
    double L_ph;          /**< Phase inductance [H]                       */
    double K_e_ll;        /**< Line-to-line back-EMF constant [V·s/rad]   */
    double K_t;           /**< Torque constant [N·m/A]                    */
    double J_r;           /**< Rotor inertia [kg·m²]                     */
    double B_r;           /**< Viscous damping [N·m·s/rad]               */
    double max_voltage_dc;/**< DC bus voltage [V]                         */
    double max_current;   /**< Peak phase current [A]                     */
    double max_speed;     /**< Max speed [rad/s]                          */
    double pole_pairs;    /**< Number of pole pairs                       */
    int    commutation;   /**< 0=sine, 1=six-step, 2=FOC                 */
} BLDCMotorParams;

/*---------------------------------------------------------------------------*/
/* L1.3: Stepper Motor Parameters                                            */
/*---------------------------------------------------------------------------*/

typedef struct {
    double R_coil;        /**< Coil resistance [Ω]                        */
    double L_coil;        /**< Coil inductance [H]                        */
    double holding_torque;/**< Holding torque [N·m]                       */
    double detent_torque; /**< Detent torque [N·m]                        */
    double step_angle;    /**< Full-step angle [rad]                      */
    int    steps_per_rev; /**< Steps per revolution                       */
    double max_slew_rate; /**< Max step rate [steps/s]                    */
    double J_r;           /**< Rotor inertia [kg·m²]                     */
} StepperMotorParams;

/*---------------------------------------------------------------------------*/
/* L1.4: Gear Train                                                          */
/*---------------------------------------------------------------------------*/

typedef struct {
    double ratio;         /**< Gear ratio N (>1 = speed reduction)        */
    double efficiency;    /**< Power transmission efficiency (0-1)        */
    double J_input;       /**< Input gear inertia [kg·m²]                */
    double J_output;      /**< Output gear inertia [kg·m²]               */
    double backlash;      /**< Backlash angle [rad]                       */
    double K_stiffness;   /**< Torsional stiffness [N·m/rad]             */
    double B_damping;     /**< Torsional damping [N·m·s/rad]             */
} GearTrain;

/*---------------------------------------------------------------------------*/
/* L1.5: Lead Screw / Ball Screw                                             */
/*---------------------------------------------------------------------------*/

typedef struct {
    double pitch;         /**< Screw lead (axial travel per rev) [m]      */
    double efficiency;    /**< Forward efficiency (0-1)                   */
    double back_efficiency;/**< Back-driving efficiency (0-1)             */
    double diameter;      /**< Nominal diameter [m]                       */
    double J_screw;       /**< Screw inertia [kg·m²]                     */
    double K_axial;       /**< Axial stiffness [N/m]                     */
    double mu_v;          /**< Coefficient of friction                    */
    double max_load;      /**< Maximum axial load [N]                     */
} LeadScrew;

/*---------------------------------------------------------------------------*/
/* L1.6: Belt Drive                                                          */
/*---------------------------------------------------------------------------*/

typedef struct {
    double ratio;         /**< Pulley ratio (driven/driver)               */
    double efficiency;    /**< Transmission efficiency (0-1)              */
    double belt_stiffness;/**< Belt stiffness [N/m]                       */
    double belt_damping;  /**< Belt viscoelastic damping [N·s/m]          */
    double pulley1_radius;/**< Driver pulley radius [m]                   */
    double pulley2_radius;/**< Driven pulley radius [m]                   */
    double J_pulley1;     /**< Driver pulley inertia [kg·m²]            */
    double J_pulley2;     /**< Driven pulley inertia [kg·m²]            */
    double preload;       /**< Belt preload force [N]                     */
} BeltDrive;

/*---------------------------------------------------------------------------*/
/* L1.7: Rotary Sensor (Encoder, Resolver, Tachometer)                       */
/*---------------------------------------------------------------------------*/

typedef enum {
    SENSOR_ENCODER_INCREMENTAL = 0,
    SENSOR_ENCODER_ABSOLUTE   = 1,
    SENSOR_RESOLVER           = 2,
    SENSOR_TACHOMETER         = 3,
    SENSOR_HALL_EFFECT        = 4,
    SENSOR_POTENTIOMETER      = 5
} RotarySensorType;

typedef struct {
    RotarySensorType type;
    double resolution;    /**< Encoder: rad/count; Resolver: electrical cycles/rev */
    int    ppr;           /**< Pulses per revolution (encoder)              */
    double max_speed;     /**< Maximum tracking speed [rad/s]               */
    double bandwidth;     /**< Measurement bandwidth [Hz]                   */
    double noise_density; /**< Noise spectral density [(rad)²/Hz]          */
    double latency;       /**< Measurement latency [s]                      */
} RotarySensor;

/*---------------------------------------------------------------------------*/
/* L1.8: Linear Sensor (LVDT, Strain Gauge, Laser)                            */
/*---------------------------------------------------------------------------*/

typedef enum {
    LINEAR_SENSOR_LVDT         = 0,
    LINEAR_SENSOR_STRAIN_GAUGE = 1,
    LINEAR_SENSOR_LASER        = 2,
    LINEAR_SENSOR_CAPACITIVE   = 3,
    LINEAR_SENSOR_EDDY_CURRENT = 4
} LinearSensorType;

typedef struct {
    LinearSensorType type;
    double range;         /**< Measurement range [m]                        */
    double resolution;    /**< Resolution [m]                               */
    double bandwidth;     /**< Bandwidth [Hz]                               */
    double linearity;     /**< Linearity error (% full scale / 100)        */
    double noise_density; /**< Noise [m/√Hz]                               */
    double sensitivity;   /**< Output [V/m] (for analog sensors)            */
} LinearSensor;

/*---------------------------------------------------------------------------*/
/* L1.9: Inertial Measurement Unit (IMU)                                     */
/*---------------------------------------------------------------------------*/

typedef struct {
    double accel_range;       /**< Accelerometer range [m/s²]             */
    double accel_noise;       /**< Accelerometer noise [m/s²/√Hz]        */
    double accel_bias_stab;   /**< Accel bias stability [m/s²]            */
    double gyro_range;        /**< Gyroscope range [rad/s]                 */
    double gyro_noise;        /**< Gyro noise [(rad/s)/√Hz]              */
    double gyro_bias_stab;    /**< Gyro bias stability [rad/s]             */
    double bandwidth;         /**< Sensor bandwidth [Hz]                   */
    double update_rate;       /**< Output data rate [Hz]                   */
} IMUParams;

/*---------------------------------------------------------------------------*/
/* L1.10: Mechatronic Load Model                                             */
/*---------------------------------------------------------------------------*/

typedef struct {
    double J_load;        /**< Load inertia [kg·m²]                       */
    double M_load;        /**< Load mass (linear) [kg]                    */
    double B_load;        /**< Viscous friction [N·m·s/rad] or [N·s/m]   */
    double T_coulomb;     /**< Coulomb friction [N·m] or [N]              */
    double T_stribeck;    /**< Stribeck friction amplitude [N·m]          */
    double v_stribeck;    /**< Stribeck velocity [rad/s] or [m/s]        */
    double T_static;      /**< Static friction (stiction) [N·m]          */
    double K_elastic;     /**< Elastic stiffness seen by motor [N·m/rad] */
    double D_elastic;     /**< Structural damping [N·m·s/rad]            */
} MechatronicLoad;

/*---------------------------------------------------------------------------*/
/* L1.11: Motor Selection Result                                             */
/*---------------------------------------------------------------------------*/

typedef struct {
    DCMotorParams motor;
    double reflected_J;   /**< Load inertia reflected to motor [kg·m²]   */
    double root_mean_sq_torque;  /**< RMS torque over duty cycle [N·m]   */
    double peak_torque;   /**< Required peak torque [N·m]                 */
    double continuous_speed;/**< Required continuous speed [rad/s]         */
    double torque_margin; /**< Torque margin (avail/req'd - 1) in %       */
    double thermal_margin;/**< Thermal margin                             */
    int    feasible;       /**< 1 if motor satisfies constraints          */
} MotorSelection;

/*---------------------------------------------------------------------------*/
/* L1.12: Mechatronic Subsystem — Integrated Motor + Transmission + Load     */
/*---------------------------------------------------------------------------*/

typedef struct {
    DCMotorParams    motor;
    GearTrain        gear;
    MechatronicLoad  load;
    /* Reflected parameters to motor shaft */
    double J_total_reflected;  /**< Total inertia at motor shaft [kg·m²] */
    double B_total_reflected;  /**< Total damping at motor shaft [N·m·s/rad] */
    double omega_n;            /**< Natural frequency [rad/s]            */
    double zeta;               /**< Damping ratio                        */
    double tau_sys;            /**< Overall time constant [s]            */
} MechatronicSubsystem;

/*---------------------------------------------------------------------------*/
/* L1.13: Servo Axis (complete single-axis drive)                            */
/*---------------------------------------------------------------------------*/

typedef struct {
    MechatronicSubsystem mech;
    RotarySensor         sensor;
    double Kp_position;  /**< Position loop proportional gain [1/s]       */
    double Kp_velocity;  /**< Velocity loop gain [(N·m)/(rad/s)]         */
    double Ki_velocity;  /**< Velocity loop integral gain                 */
    double Kp_current;   /**< Current loop gain [V/A]                     */
    double Ki_current;   /**< Current loop integral gain [V/(A·s)]       */
    double velocity_limit;/**< Velocity saturation [rad/s]                */
    double current_limit; /**< Current / torque saturation [A]            */
    double max_accel;     /**< Acceleration limit [rad/s²]               */
    double bandwidth_pos; /**< Position loop bandwidth [Hz]               */
    double bandwidth_vel; /**< Velocity loop bandwidth [Hz]               */
    double bandwidth_cur; /**< Current loop bandwidth [Hz]                */
} ServoAxis;

/*---------------------------------------------------------------------------*/
/* L1.14: Multi-Axis Manipulator                                             */
/*---------------------------------------------------------------------------*/

#define MAX_MANIPULATOR_JOINTS 7

typedef struct {
    int    n_joints;
    ServoAxis joints[MAX_MANIPULATOR_JOINTS];
    double DH_a[MAX_MANIPULATOR_JOINTS];      /**< Denavit-Hartenberg a [m]  */
    double DH_alpha[MAX_MANIPULATOR_JOINTS];  /**< DH alpha [rad]            */
    double DH_d[MAX_MANIPULATOR_JOINTS];      /**< DH d [m]                  */
    double DH_theta_offset[MAX_MANIPULATOR_JOINTS]; /**< DH theta offset     */
    double payload_mass;     /**< End-effector payload [kg]                 */
    double max_reach;        /**< Maximum reach [m]                         */
    double repeatability;    /**< Positioning repeatability [m]              */
    double max_payload;      /**< Rated payload capacity [kg]               */
} Manipulator;

/*---------------------------------------------------------------------------*/
/* L1.15: PWM Drive Parameters                                               */
/*---------------------------------------------------------------------------*/

typedef struct {
    double V_dc;          /**< DC bus voltage [V]                          */
    double f_pwm;         /**< PWM switching frequency [Hz]                */
    double f_sample;      /**< Current loop sampling frequency [Hz]        */
    double dead_time;     /**< Dead-time [s]                               */
    double R_ds_on;       /**< MOSFET on-resistance [Ω]                   */
    double V_diode;       /**< Freewheeling diode forward voltage [V]     */
    double max_duty;      /**< Maximum duty cycle (0-1)                    */
    double min_duty;      /**< Minimum duty cycle (0-1)                    */
    int    is_h_bridge;   /**< 1 = H-bridge, 0 = half-bridge               */
} PWMDrive;

/*===========================================================================*/
/* L1 API: Constructor / Initializer Functions                               */
/*===========================================================================*/

DCMotorParams      dc_motor_create(double R_a, double L_a, double K_e,
                                   double J_m, double B_m,
                                   double max_voltage, double max_current);

BLDCMotorParams    bldc_motor_create(double R_ph, double L_ph, double K_e_ll,
                                     double J_r, double B_r,
                                     double max_vdc, double max_I, double pp);

StepperMotorParams stepper_motor_create(double R_coil, double L_coil,
                                        double hold_torque, double step_angle,
                                        double J_r);

GearTrain          gear_train_create(double ratio, double efficiency,
                                     double J_in, double J_out,
                                     double backlash, double K_stiff);

LeadScrew          lead_screw_create(double pitch, double efficiency,
                                     double diameter, double J_screw,
                                     double mu_v);

BeltDrive          belt_drive_create(double ratio, double efficiency,
                                     double r1, double r2,
                                     double belt_K, double belt_D);

RotarySensor       rotary_sensor_create(RotarySensorType type, int ppr,
                                        double bandwidth, double latency);

LinearSensor       linear_sensor_create(LinearSensorType type, double range,
                                        double resolution, double sensitivity);

MechatronicLoad    mechatronic_load_create(double J, double M,
                                           double B, double T_coulomb,
                                           double T_static);

MechatronicSubsystem mechatronic_subsystem_create(DCMotorParams motor,
                                                   GearTrain gear,
                                                   MechatronicLoad load);

ServoAxis          servo_axis_create(MechatronicSubsystem mech,
                                     RotarySensor sensor);

PWMDrive           pwm_drive_create(double V_dc, double f_pwm,
                                    double dead_time, int is_h_bridge);

#endif /* MECHATRONIC_DEFINITIONS_H */
