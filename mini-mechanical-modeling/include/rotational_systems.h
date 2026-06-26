/**
 * @file rotational_systems.h
 * @brief Rotational mechanical systems - gears, torsion, bearings, rotating unbalance
 *
 * Knowledge Coverage:
 *   L1 - Definitions: moment of inertia, torsional stiffness, angular velocity/acceleration,
 *        gear ratio, reflected inertia, torsional natural frequency
 *   L2 - Core Concepts: rotational Newton"s 2nd law, torsional Hooke"s law,
 *        angular momentum, power transmission, gear efficiency
 *   L3 - Mathematical Structures: rotational ODEs, torsional impedance
 *   L4 - Fundamental Laws: conservation of angular momentum, torque balance
 *   L5 - Engineering Methods: gear train analysis, reflected parameters,
 *        rotational-translational conversion (rack-pinion, lead screw)
 *
 * Reference:
 *   MIT 2.003 Dynamics and Control I
 *   Stanford ENGR105
 *   Thomson & Dahleh, Theory of Vibration (1998)
 *   Norton, Machine Design (2019)
 */

#ifndef ROTATIONAL_SYSTEMS_H
#define ROTATIONAL_SYSTEMS_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* ===== L1: Core Data Types ===== */

typedef struct {
    double inertia;
    double angle;
    double angular_vel;
    double angular_accel;
} Rotor;

typedef struct {
    double stiffness;
    double free_angle;
    double max_torque;
} TorsionSpring;

typedef struct {
    double damping_coeff;
    double max_torque;
} RotDamper;

typedef struct {
    int teeth_driver;
    int teeth_driven;
    double efficiency;
    double backlash_rad;
} GearPair;

typedef struct {
    int n_stages;
    GearPair *stages;
    double total_ratio;
    double total_efficiency;
} GearTrain;

/* ===== L2: Rotational Newton"s 2nd Law ===== */

double rotational_torque(double inertia, double angular_accel);
double rotational_accel(double torque, double inertia);
double angular_momentum(double inertia, double angular_vel);

/* ===== L2: Torsional Spring ===== */

double torsional_spring_torque(double stiffness, double deflection);
double torsional_potential_energy(double stiffness, double deflection);
double torsional_spring_series(const double *k_vals, size_t n);
double torsional_spring_parallel(const double *k_vals, size_t n);

/* ===== L2: Rotational Damping ===== */

double rotational_damping_torque(double coeff, double angular_vel);

/* ===== L2: Rotational Kinetic Energy ===== */

double rotational_kinetic_energy(double inertia, double angular_vel);
double total_kinetic_energy(double m, double v, double J, double w);

/* ===== L3: Torsional Vibration Parameters ===== */

double torsional_natural_freq(double k, double J);
double torsional_natural_freq_hz(double k, double J);
double torsional_damping_ratio(double b, double k, double J);
double torsional_damped_freq(double k, double J, double b);

/* ===== L5: Gear Ratio and Reflected Parameters ===== */

double gear_ratio_calc(int teeth_driver, int teeth_driven);
double gear_output_speed(double omega_in, double ratio);
double gear_output_torque(double T_in, double ratio);
double gear_output_torque_eff(double T_in, double ratio, double efficiency);
double reflected_inertia(double J_load, double ratio);
double reflected_damping(double b_load, double ratio);
double reflected_stiffness(double k_load, double ratio);
double total_reflected_inertia(double J_motor, const double *J_loads,
                                const double *ratios, int n);
double optimal_gear_ratio(double J_motor, double J_load);
double gear_train_total_ratio(const GearTrain *train);
int gear_train_build(const int *t_driver, const int *t_driven,
                     const double *eff, int n, GearTrain *train);
void gear_train_free(GearTrain *train);

/* ===== L5: Bearing Friction ===== */

double bearing_friction_torque(double T_static, double b_visc, double w);
double bearing_power_loss(double T_friction, double w);

/* ===== L5: Rotational-to-Translational ===== */

double rack_pinion_position(double r, double theta);
double rack_pinion_velocity(double r, double w);
double rack_pinion_force(double r, double T);
double rack_pinion_equiv_mass(double J_pinion, double r);
double leadscrew_position(double lead, double theta);
double leadscrew_velocity(double lead, double w);
double leadscrew_force(double lead, double T, double eta);
double leadscrew_reflected_inertia(double m, double lead);

/* ===== L6: Rotating Unbalance ===== */

double centrifugal_force(double m_u, double e, double w);
double unbalance_displacement_amp(double m_u, double e, double M_total,
                                   double r, double zeta);
double normalized_unbalance_response(double r, double zeta);
double critical_speed_rad_s(double k, double m_rotor);
double critical_speed_rpm(double k, double m_rotor);

/* ===== L6: Rotor Dynamics ===== */

double jeffcott_rotor_nat_freq(double k, double m);
double rotor_amplification_factor(double zeta);
double rotor_orbit_radius(double m_u, double e, double M, double r, double zeta);

/* ===== L7: Applications ===== */

double flywheel_energy(double J, double w);
double flywheel_energy_density(double J, double m, double w);
double clutch_engagement_time(double J_load, double T_clutch, double w_target);
double belt_drive_natural_freq(double belt_stiffness, double pulley_inertia1,
                                double pulley_inertia2);
double shaft_whirling_speed(double E_mod, double area_I, double m_per_length, double L,
                             int end_condition);

#endif /* ROTATIONAL_SYSTEMS_H */
