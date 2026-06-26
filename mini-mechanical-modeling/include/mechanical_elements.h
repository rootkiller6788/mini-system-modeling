/**
 * @file mechanical_elements.h
 * @brief Fundamental mechanical element definitions - mass, spring, damper, friction
 *
 * Knowledge Coverage:
 *   L1 - Definitions: mass, stiffness, damping coefficient, inertia, natural frequency,
 *        damping ratio, logarithmic decrement, transmissibility, mechanical impedance
 *   L2 - Core Concepts: Newton"s laws, Hooke"s law, viscous damping, Coulomb friction,
 *        Stribeck friction, kinetic/potential energy, resonance, critical damping
 *   L3 - Mathematical Structures: 2nd-order ODE parameters, complex impedance,
 *        eigenvalue-based damping, static deflection, half-power bandwidth
 *   L4 - Fundamental Laws: Newton"s 2nd (F=ma), Hooke"s (F=kx), viscous (F=bv),
 *        conservation of mechanical energy
 *
 * Reference:
 *   MIT 2.003 Dynamics and Control I
 *   Stanford ENGR105 Feedback Control
 *   Thomson & Dahleh, Theory of Vibration with Applications (1998)
 *   Rao, Mechanical Vibrations (2017)
 *
 * All functions documented with physical theorem source and O(.) complexity.
 */

#ifndef MECHANICAL_ELEMENTS_H
#define MECHANICAL_ELEMENTS_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* ===== L1: Core Data Types ===== */

/** Mass element [kg] - stores kinetic energy E = 0.5*m*v^2 */
typedef struct {
    double mass;
    double position;
    double velocity;
    double acceleration;
} MassEl;

/** Linear spring - stores potential energy E = 0.5*k*x^2. Obeys Hooke"s law. */
typedef struct {
    double stiffness;
    double free_length;
    double max_deflection;
    double internal_damping;
} SpringEl;

/** Linear damper (dashpot) - dissipates energy via F = b*v */
typedef struct {
    double damping_coeff;
    double max_force;
    double stroke_length;
} DamperEl;

/** Friction model types */
typedef enum {
    FRICTION_NONE = 0,
    FRICTION_VISCOUS = 1,
    FRICTION_COULOMB = 2,
    FRICTION_STRIBECK = 3,
    FRICTION_STATIC = 4
} FrictionType;

/** Friction model parameters */
typedef struct {
    FrictionType type;
    double static_coeff;
    double coulomb_coeff;
    double viscous_coeff;
    double stribeck_vel;
    double normal_force;
} FrictionModel;

/* ===== L2: Newton"s Second Law (F = ma) ===== */

/** Compute force from mass and acceleration. O(1).
 *  Theorem: Newton (1687) Principia. */
double newton_second_law_force(double mass, double acceleration);

/** Compute acceleration from force and mass. O(1).
 *  Guard: mass <= 0 returns DBL_MAX or -DBL_MAX. */
double newton_second_law_acceleration(double force, double mass);

/* ===== L2: Hooke"s Law (F = k*dx) ===== */

double hookes_law_force(double stiffness, double displacement);
double hookes_law_displacement(double stiffness, double force);
double spring_potential_energy(double stiffness, double displacement);
double spring_series_equivalent(const double *stiffness, size_t n);
double spring_parallel_equivalent(const double *stiffness, size_t n);
double spring_combined_two(double k1, double k2, int is_series);

/* ===== L2: Viscous Damping (F = b*v) ===== */

double viscous_damping_force(double damping_coeff, double velocity);
double damper_series_equivalent(const double *damping_coeffs, size_t n);
double damper_parallel_equivalent(const double *damping_coeffs, size_t n);

/* ===== L2: Friction Models ===== */

double coulomb_friction_magnitude(double kinetic_coeff, double normal_force);
double coulomb_friction_force(double kinetic_coeff, double normal_force, double velocity);
double stribeck_friction_force(const FrictionModel *fm, double velocity);
int static_friction_holds(const FrictionModel *fm, double applied_force);
double equivalent_viscous_from_coulomb(double coulomb_force, double velocity_amplitude, double freq);

/* ===== L2: Kinetic and Potential Energy ===== */

double kinetic_energy_translational(double mass, double velocity);
double gravitational_potential_energy(double mass, double height, double gravity);
double total_mechanical_energy(double mass, double velocity, double stiffness,
                                double displacement, double height, double gravity);
double energy_dissipated_per_cycle(double damping_coeff, double angular_freq, double amplitude);

/* ===== L3: Natural Frequency and Damping Parameters ===== */

double natural_frequency_undamped(double stiffness, double mass);
double natural_frequency_hz(double stiffness, double mass);
double damping_ratio(double damping_coeff, double stiffness, double mass);
double critical_damping_coefficient(double stiffness, double mass);
double damped_natural_frequency(double stiffness, double mass, double damping_coeff);
double logarithmic_decrement(double amplitude_0, double amplitude_n, int num_cycles);
double damping_ratio_from_logdec(double log_dec);
double quality_factor_from_damping(double damping_ratio);
double damping_ratio_from_q_factor(double q_factor);
double loss_factor_from_damping(double damping_ratio);

/* ===== L3: Static Deflection ===== */

double static_deflection(double mass, double stiffness, double gravity);
double stiffness_from_static_deflection(double mass, double deflection, double gravity);
double natural_freq_from_static_deflection(double deflection, double gravity);

/* ===== L3: Mechanical Impedance ===== */

double complex mechanical_impedance_sdof(double mass, double damping, double stiffness,
                                          double angular_freq);
double mechanical_impedance_magnitude(double mass, double damping, double stiffness,
                                       double angular_freq);
double mechanical_impedance_phase(double mass, double damping, double stiffness,
                                   double angular_freq);
double antiresonant_frequency(double stiffness, double mass);
double complex mechanical_mobility(double mass, double damping, double stiffness,
                                    double angular_freq);

/* ===== L5: Vibration Transmissibility ===== */

double force_transmissibility(double freq_ratio, double damping_ratio);
double displacement_transmissibility(double freq_ratio, double damping_ratio);
double peak_transmissibility_freq_ratio(double damping_ratio);
double maximum_transmissibility(double damping_ratio);
double isolation_efficiency(double freq_ratio, double damping_ratio);

/* ===== L3: Mass Moments of Inertia (rigid bodies) ===== */

double inertia_point_mass(double mass, double distance);
double inertia_rod_end(double mass, double length);
double inertia_rod_center(double mass, double length);
double inertia_solid_cylinder(double mass, double radius);
double inertia_hollow_cylinder(double mass, double r_outer, double r_inner);
double inertia_solid_sphere(double mass, double radius);
double inertia_rectangular_plate(double mass, double side_a, double side_b);
double parallel_axis_theorem(double i_cm, double mass, double distance);
double radius_of_gyration(double moment_of_inertia, double mass);

/* ===== L5: Time-Domain Response Parameters ===== */

double peak_time_underdamped(double damping_ratio, double natural_freq);
double percent_overshoot(double damping_ratio);
double settling_time_2percent(double damping_ratio, double natural_freq);
double settling_time_5percent(double damping_ratio, double natural_freq);
double rise_time_underdamped(double damping_ratio, double natural_freq);
double delay_time(double damping_ratio, double natural_freq);
double sdof_step_response_displacement(double static_deflection_val, double natural_freq,
                                        double damping_ratio_val, double time);
double resonance_amplitude_force(double force_amplitude, double stiffness, double damping_ratio);
double phase_angle_sdof(double freq_ratio, double damping_ratio);
double half_power_bandwidth(double damping_ratio, double natural_freq);
double damping_from_half_power(double omega1, double omega2, double omega_n);

/* ===== L5: Equivalent Systems ===== */

double effective_mass_sdof_lumped(double mass, double mode_shape_at_point);
double effective_stiffness_from_omega(double natural_freq, double effective_mass);
double mass_normalized_stiffness(double mass, double stiffness);

/* ===== L5: Structural Stiffness Formulas ===== */

double shaft_torsional_stiffness(double shear_modulus, double polar_moment, double length);
double cantilever_beam_stiffness(double E_modulus, double inertia_area, double length);
double simply_supported_beam_stiffness(double E_modulus, double inertia_area, double length);
double area_moment_rectangular(double width, double height);
double area_moment_circular(double diameter);
double polar_moment_circular(double diameter);

#endif /* MECHANICAL_ELEMENTS_H */
