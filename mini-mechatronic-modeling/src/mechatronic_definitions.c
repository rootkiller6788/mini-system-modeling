/**
 * @file mechatronic_definitions.c
 * @brief L1 Constructor implementations for mechatronic data types.
 *
 * Every struct in mechatronic_definitions.h has a corresponding constructor
 * here. Each constructor validates input ranges and precomputes derived
 * parameters (time constants, natural frequencies).
 *
 * 9-Layer Coverage: L1 (Definitions) — Complete
 *
 * References:
 *   Bolton (2015), Cetinkunt (2007), Shetty & Kolk (2011)
 */

#include "mechatronic_definitions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================*/
/* L1.1: DC Motor Constructor                                                */
/*===========================================================================*/

DCMotorParams dc_motor_create(double R_a, double L_a, double K_e,
                               double J_m, double B_m,
                               double max_voltage, double max_current)
{
    /*
     * Source: Krause et al. "Analysis of Electric Machinery" (2013) Ch.2
     *
     * The DC motor is characterized by:
     *   Electrical: R_a (armature resistance), L_a (armature inductance)
     *   Electromechanical: K_e = K_t (back-EMF/torque constant in SI)
     *   Mechanical: J_m (rotor inertia), B_m (viscous friction)
     *
     * Derived:
     *   τ_e = L_a / R_a  — electrical time constant (current response)
     *   τ_m = J_m / B_m  — mechanical time constant (speed response)
     *   P_rated = max_voltage · max_current
     *   max_speed = max_voltage / K_e (ideal no-load)
     *   max_torque = K_t · max_current (stall at rated current)
     */
    DCMotorParams m;
    memset(&m, 0, sizeof(m));

    m.R_a = R_a > 0.0 ? R_a : 1.0;
    m.L_a = L_a > 0.0 ? L_a : 0.001;
    m.K_e = K_e > 0.0 ? K_e : 0.01;
    m.K_t = m.K_e;  /* SI equivalence: K_t = K_e */
    m.J_m = J_m > 0.0 ? J_m : 1e-6;
    m.B_m = B_m > 0.0 ? B_m : 1e-6;
    m.tau_e = m.L_a / m.R_a;
    m.tau_m = m.J_m / m.B_m;
    m.max_voltage = max_voltage > 0.0 ? max_voltage : 12.0;
    m.max_current = max_current > 0.0 ? max_current : 1.0;
    m.max_speed = m.max_voltage / m.K_e;
    m.max_torque = m.K_t * m.max_current;
    m.P_rated = m.max_voltage * m.max_current;
    m.pole_pairs = 1;

    return m;
}

/*===========================================================================*/
/* L1.2: BLDC Motor Constructor                                              */
/*===========================================================================*/

BLDCMotorParams bldc_motor_create(double R_ph, double L_ph, double K_e_ll,
                                   double J_r, double B_r,
                                   double max_vdc, double max_I, double pp)
{
    /*
     * Source: Krishnan "Permanent Magnet Synchronous and Brushless
     *         DC Motor Drives" (2010)
     *
     * For sinusoidal commutation:
     *   K_t = (3/2) · (K_e_ll / √3)  — torque constant from line-line BEMF
     *
     * K_e_ll is the line-to-line back-EMF constant in V·s/rad.
     * For six-step commutation, K_t = K_e_ll (different convention).
     */
    BLDCMotorParams m;
    memset(&m, 0, sizeof(m));

    m.R_ph = R_ph > 0.0 ? R_ph : 0.5;
    m.L_ph = L_ph > 0.0 ? L_ph : 0.002;
    m.K_e_ll = K_e_ll > 0.0 ? K_e_ll : 0.01;
    /* For FOC: K_t per-phase = K_e_ll / √3 */
    m.K_t = m.K_e_ll / 1.7320508075688772;
    m.J_r = J_r > 0.0 ? J_r : 1e-6;
    m.B_r = B_r > 0.0 ? B_r : 1e-6;
    m.max_voltage_dc = max_vdc > 0.0 ? max_vdc : 24.0;
    m.max_current = max_I > 0.0 ? max_I : 2.0;
    m.max_speed = m.max_voltage_dc / m.K_e_ll;
    m.pole_pairs = pp > 0.0 ? pp : 2.0;
    m.commutation = 2;  /* Default: FOC */

    return m;
}

/*===========================================================================*/
/* L1.3: Stepper Motor Constructor                                           */
/*===========================================================================*/

StepperMotorParams stepper_motor_create(double R_coil, double L_coil,
                                         double hold_torque, double step_angle,
                                         double J_r)
{
    /*
     * Source: Kenjo & Sugawara "Stepping Motors and Their
     *         Microprocessor Controls" (1994)
     *
     * Hybrid stepper motor typical parameters:
     *   Step angle = 1.8° (200 steps/rev) — most common
     *   Holding torque: 0.01–10 N·m range
     *
     * Detent torque ≈ 5-10% of holding torque (permanent magnet residual)
     */
    StepperMotorParams m;
    memset(&m, 0, sizeof(m));

    m.R_coil = R_coil > 0.0 ? R_coil : 2.0;
    m.L_coil = L_coil > 0.0 ? L_coil : 0.005;
    m.holding_torque = hold_torque > 0.0 ? hold_torque : 0.5;
    m.step_angle = step_angle > 0.0 ? step_angle : 0.03141592653589793; /* 1.8° */
    m.steps_per_rev = (int)(2.0 * M_PI / m.step_angle + 0.5);
    m.detent_torque = 0.05 * m.holding_torque;
    m.max_slew_rate = 2000.0;  /* typical: 2000 steps/s */
    m.J_r = J_r > 0.0 ? J_r : 1e-6;

    return m;
}

/*===========================================================================*/
/* L1.4: Gear Train Constructor                                              */
/*===========================================================================*/

GearTrain gear_train_create(double ratio, double efficiency,
                             double J_in, double J_out,
                             double backlash, double K_stiff)
{
    /*
     * Source: Dudley "Handbook of Practical Gear Design" (1994)
     *
     * Gear ratio N = ω_in / ω_out (N > 1 = speed reduction)
     * Efficiency η = P_out / P_in
     * Reflected inertia: J_ref = J_in + J_out / N²
     * The gear ratio amplifies torque by N and reduces speed by N.
     */
    GearTrain g;
    memset(&g, 0, sizeof(g));

    g.ratio = ratio > 0.0 ? ratio : 1.0;
    g.efficiency = (efficiency > 0.0 && efficiency <= 1.0) ? efficiency : 0.9;
    g.J_input = J_in >= 0.0 ? J_in : 1e-6;
    g.J_output = J_out >= 0.0 ? J_out : 1e-6;
    g.backlash = backlash >= 0.0 ? backlash : 0.0;
    g.K_stiffness = K_stiff > 0.0 ? K_stiff : 1e6;
    g.B_damping = 0.01 * g.K_stiffness; /* approximate from stiffness */

    return g;
}

/*===========================================================================*/
/* L1.5: Lead Screw Constructor                                              */
/*===========================================================================*/

LeadScrew lead_screw_create(double pitch, double efficiency,
                             double diameter, double J_screw,
                             double mu_v)
{
    /*
     * Source: Budynas & Nisbett "Shigley's Mechanical Engineering Design"
     *         (2015) Ch.8
     *
     * Lead screw converts rotary motion to linear:
     *   x_linear = θ_rotary · pitch / (2π)
     *   F_axial = T_rotary · 2π · η / pitch
     *
     * Ball screws: η ≈ 0.90–0.95 (rolling contact)
     * Lead screws: η ≈ 0.30–0.50 (sliding contact)
     */
    LeadScrew s;
    memset(&s, 0, sizeof(s));

    s.pitch = pitch > 0.0 ? pitch : 0.005;
    s.efficiency = (efficiency > 0.0 && efficiency <= 1.0) ? efficiency : 0.9;
    s.back_efficiency = 0.7 * s.efficiency;  /* ball screws are back-driveable */
    s.diameter = diameter > 0.0 ? diameter : 0.016;
    s.J_screw = J_screw > 0.0 ? J_screw : 1e-5;
    s.K_axial = 1e8;  /* typical ball screw axial stiffness */
    s.mu_v = mu_v > 0.0 ? mu_v : 0.003;
    s.max_load = 5000.0;  /* typical for 16mm ball screw */

    return s;
}

/*===========================================================================*/
/* L1.6: Belt Drive Constructor                                              */
/*===========================================================================*/

BeltDrive belt_drive_create(double ratio, double efficiency,
                             double r1, double r2,
                             double belt_K, double belt_D)
{
    /*
     * Source: Gates Rubber Company "Belt Drive Design Manual"
     *
     * Timing belts provide positive (no-slip) power transmission
     * with typical efficiency 0.95–0.98.
     * Stiffness depends on belt width and cord material.
     */
    BeltDrive b;
    memset(&b, 0, sizeof(b));

    b.ratio = ratio > 0.0 ? ratio : 1.0;
    b.efficiency = (efficiency > 0.0 && efficiency <= 1.0) ? efficiency : 0.95;
    b.pulley1_radius = r1 > 0.0 ? r1 : 0.02;
    b.pulley2_radius = r2 > 0.0 ? r2 : 0.02;
    b.belt_stiffness = belt_K > 0.0 ? belt_K : 1e5;
    b.belt_damping = belt_D > 0.0 ? belt_D : 100.0;
    b.J_pulley1 = 1e-5;
    b.J_pulley2 = 1e-5;
    b.preload = 200.0;  /* typical preload [N] */

    return b;
}

/*===========================================================================*/
/* L1.7-8: Sensor Constructors                                               */
/*===========================================================================*/

RotarySensor rotary_sensor_create(RotarySensorType type, int ppr,
                                   double bandwidth, double latency)
{
    /*
     * Source: Fraden "Handbook of Modern Sensors" (2016) Ch.5
     *
     * Incremental encoder: ppr = pulses/rev, resolution = 2π/(ppr·4) with quad
     * Resolver: typically 2-4 pole, accuracy 5-20 arcmin
     * Tachometer: analog DC output, K_tach in V·s/rad
     */
    RotarySensor s;
    memset(&s, 0, sizeof(s));

    s.type = type;
    s.ppr = ppr > 0 ? ppr : 1000;

    switch (type) {
    case SENSOR_ENCODER_INCREMENTAL:
        s.resolution = (2.0 * M_PI) / (s.ppr * 4);  /* with 4x quadrature */
        s.bandwidth = bandwidth > 0.0 ? bandwidth : 100e3;
        s.noise_density = 1e-8;
        break;
    case SENSOR_ENCODER_ABSOLUTE:
        s.resolution = (2.0 * M_PI) / (1 << 17);  /* 17-bit typical */
        s.bandwidth = bandwidth > 0.0 ? bandwidth : 10e3;
        s.noise_density = 1e-9;
        break;
    case SENSOR_RESOLVER:
        s.resolution = 2.0 * M_PI / 65536.0;  /* 16-bit RDC typical */
        s.bandwidth = bandwidth > 0.0 ? bandwidth : 1000.0;
        s.noise_density = 1e-7;
        break;
    case SENSOR_TACHOMETER:
        s.resolution = 0.01;  /* rad/s per mV typical */
        s.bandwidth = bandwidth > 0.0 ? bandwidth : 500.0;
        s.noise_density = 1e-5;
        break;
    case SENSOR_HALL_EFFECT:
        s.resolution = (2.0 * M_PI) / (s.ppr * 6);  /* 6-step commutation */
        s.bandwidth = bandwidth > 0.0 ? bandwidth : 10000.0;
        s.noise_density = 1e-6;
        break;
    case SENSOR_POTENTIOMETER:
        s.resolution = 1e-4;
        s.bandwidth = bandwidth > 0.0 ? bandwidth : 100.0;
        s.noise_density = 1e-6;
        break;
    default:
        s.resolution = 0.001;
        s.bandwidth = 1000.0;
        s.noise_density = 1e-6;
    }

    s.max_speed = (2.0 * M_PI * s.bandwidth) / s.ppr;
    s.bandwidth = bandwidth > 0.0 ? bandwidth : 1000.0;
    s.latency = latency >= 0.0 ? latency : 1e-6;

    return s;
}

LinearSensor linear_sensor_create(LinearSensorType type, double range,
                                   double resolution, double sensitivity)
{
    /*
     * Source: Nyce "Linear Position Sensors: Theory and Application" (2004)
     *
     * LVDT: infinite resolution in theory, sub-μm in practice
     * Strain gauge: resolution ~0.1 με (microstrain)
     * Laser interferometer: resolution ~1 nm
     */
    LinearSensor s;
    memset(&s, 0, sizeof(s));

    s.type = type;
    s.range = range > 0.0 ? range : 0.1;
    s.resolution = resolution > 0.0 ? resolution : 1e-6;
    s.sensitivity = sensitivity > 0.0 ? sensitivity : 1.0;

    switch (type) {
    case LINEAR_SENSOR_LVDT:
        s.bandwidth = 1000.0;
        s.linearity = 0.0025;  /* 0.25% typical */
        s.noise_density = 1e-9;
        break;
    case LINEAR_SENSOR_STRAIN_GAUGE:
        s.bandwidth = 100.0;
        s.linearity = 0.001;
        s.noise_density = 1e-9;
        break;
    case LINEAR_SENSOR_LASER:
        s.bandwidth = 10000.0;
        s.linearity = 0.0001;
        s.noise_density = 1e-12;
        break;
    case LINEAR_SENSOR_CAPACITIVE:
        s.bandwidth = 5000.0;
        s.linearity = 0.001;
        s.noise_density = 1e-10;
        break;
    case LINEAR_SENSOR_EDDY_CURRENT:
        s.bandwidth = 50000.0;
        s.linearity = 0.002;
        s.noise_density = 1e-9;
        break;
    default:
        s.bandwidth = 1000.0;
        s.linearity = 0.01;
        s.noise_density = 1e-8;
    }

    return s;
}

/*===========================================================================*/
/* L1.9-10: Load and Subsystem Constructors                                  */
/*===========================================================================*/

MechatronicLoad mechatronic_load_create(double J, double M,
                                         double B, double T_coulomb,
                                         double T_static)
{
    /*
     * Source: Armstrong-Hélouvry "Control of Machines with Friction" (1994)
     *
     * Friction model includes:
     *   Coulomb friction: constant opposing torque
     *   Stiction: higher static friction at zero velocity
     *   Stribeck effect: exponential decay from stiction to Coulomb
     *   Viscous damping: linear with velocity
     */
    MechatronicLoad L;
    memset(&L, 0, sizeof(L));

    L.J_load = J >= 0.0 ? J : 0.001;
    L.M_load = M >= 0.0 ? M : 1.0;
    L.B_load = B >= 0.0 ? B : 0.001;
    L.T_coulomb = T_coulomb >= 0.0 ? T_coulomb : 0.01;
    L.T_stribeck = 0.3 * L.T_coulomb;
    L.v_stribeck = 0.1;  /* rad/s typical */
    L.T_static = T_static >= 0.0 ? T_static : 1.5 * L.T_coulomb;
    L.K_elastic = 1e6;   /* typical coupling stiffness */
    L.D_elastic = 0.01 * L.K_elastic;

    return L;
}

MechatronicSubsystem mechatronic_subsystem_create(DCMotorParams motor,
                                                   GearTrain gear,
                                                   MechatronicLoad load)
{
    /*
     * Source: Cetinkunt "Mechatronics" (2007) §7.4
     *
     * Computes total reflected dynamics from motor through transmission
     * to load. Uses the equivalence principle for mechanical systems.
     *
     * Reflected inertia: J_total = J_m + J_g1 + (J_g2 + J_L)/N²
     * Reflected damping: B_total = B_m + B_g1 + (B_g2 + B_L)/N²
     */
    MechatronicSubsystem s;
    memset(&s, 0, sizeof(s));

    s.motor = motor;
    s.gear = gear;
    s.load = load;

    double N = gear.ratio;
    double N2 = N * N;

    /* Reflected dynamics to motor shaft */
    s.J_total_reflected = motor.J_m + gear.J_input
                          + (gear.J_output + load.J_load) / N2;
    s.B_total_reflected = motor.B_m + gear.B_damping
                          + (gear.B_damping + load.B_load) / N2;

    /* Natural frequency and damping ratio */
    double K_eff = gear.K_stiffness / N2;  /* effective stiffness at motor */
    if (s.J_total_reflected > 0.0 && K_eff > 0.0) {
        s.omega_n = sqrt(K_eff / s.J_total_reflected);
    } else {
        s.omega_n = 0.0;
    }

    if (s.omega_n > 0.0 && s.J_total_reflected > 0.0) {
        s.zeta = s.B_total_reflected / (2.0 * s.J_total_reflected * s.omega_n);
    } else {
        s.zeta = 1.0;  /* default: critically damped assumption */
    }

    /* Overall system time constant */
    if (s.B_total_reflected > 0.0) {
        s.tau_sys = s.J_total_reflected / s.B_total_reflected;
    } else {
        s.tau_sys = motor.tau_m;  /* fall back to motor mechanical */
    }

    return s;
}

/*===========================================================================*/
/* L1.13: Servo Axis Constructor                                             */
/*===========================================================================*/

ServoAxis servo_axis_create(MechatronicSubsystem mech,
                             RotarySensor sensor)
{
    /*
     * Source: Shetty & Kolk "Mechatronics System Design" (2011) Ch.8
     *
     * Initializes a servo axis with default cascaded PID structure:
     *   Inner loop: current (torque) control
     *   Middle loop: velocity control
     *   Outer loop: position control
     */
    ServoAxis axis;
    memset(&axis, 0, sizeof(axis));

    axis.mech = mech;
    axis.sensor = sensor;

    /* Default conservative gains */
    axis.Kp_position = 10.0;     /* [1/s] */
    axis.Kp_velocity = 0.1;      /* [(N·m)/(rad/s)] */
    axis.Ki_velocity = 1.0;
    axis.Kp_current = 10.0;      /* [V/A] */
    axis.Ki_current = 100.0;
    axis.velocity_limit = mech.motor.max_speed;
    axis.current_limit = mech.motor.max_current;
    axis.max_accel = mech.motor.max_torque / mech.J_total_reflected;

    /* Bandwidth estimates */
    axis.bandwidth_cur = 1000.0;  /* Hz — current loop */
    axis.bandwidth_vel = 100.0;   /* Hz — velocity loop */
    axis.bandwidth_pos = 10.0;    /* Hz — position loop */

    return axis;
}

/*===========================================================================*/
/* L1.15: PWM Drive Constructor                                              */
/*===========================================================================*/

PWMDrive pwm_drive_create(double V_dc, double f_pwm,
                           double dead_time, int is_h_bridge)
{
    /*
     * Source: Mohan et al. "Power Electronics" (2003) Ch.8
     *
     * PWM drive for DC motor:
     *   Average voltage = V_dc · duty_cycle
     *   Current ripple ΔI = V_dc·duty·(1-duty)/(f_pwm·L_a)
     *
     * Dead time prevents shoot-through in half-bridge.
     * H-bridge enables bidirectional current (4-quadrant operation).
     */
    PWMDrive d;
    memset(&d, 0, sizeof(d));

    d.V_dc = V_dc > 0.0 ? V_dc : 24.0;
    d.f_pwm = f_pwm > 0.0 ? f_pwm : 20000.0;
    d.f_sample = d.f_pwm;  /* same as PWM for most drives */
    d.dead_time = dead_time > 0.0 ? dead_time : 1e-6;
    d.R_ds_on = 0.01;   /* typical power MOSFET */
    d.V_diode = 0.7;     /* typical silicon diode */
    d.max_duty = 0.95;
    d.min_duty = 0.05;
    d.is_h_bridge = is_h_bridge;

    return d;
}
