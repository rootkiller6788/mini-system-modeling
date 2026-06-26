/**
 * @file rotational_systems.c
 * @brief Rotational mechanical system implementations
 * Each function implements one independent knowledge point.
 */

#include "rotational_systems.h"
#include <math.h>
#include <complex.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

/* ===== Rotational Newton's 2nd Law ===== */

double rotational_torque(double inertia, double angular_accel)
{
    /* T = J * alpha (angular form of Newton's 2nd law) */
    return inertia * angular_accel;
}

double rotational_accel(double torque, double inertia)
{
    /* alpha = T / J */
    if (fabs(inertia) < 1e-30)
        return (torque >= 0.0) ? DBL_MAX : -DBL_MAX;
    return torque / inertia;
}

double angular_momentum(double inertia, double angular_vel)
{
    /* L = J * omega */
    return inertia * angular_vel;
}

/* ===== Torsional Spring ===== */

double torsional_spring_torque(double stiffness, double deflection)
{
    /* T = -k_t * theta (torsional Hooke's law) */
    return -stiffness * deflection;
}

double torsional_potential_energy(double stiffness, double deflection)
{
    /* U = 0.5 * k_t * theta^2 */
    if (stiffness <= 0.0) return 0.0;
    return 0.5 * stiffness * deflection * deflection;
}

double torsional_spring_series(const double *k_vals, size_t n)
{
    /* 1/k_eq = sum(1/k_i) */
    if (!k_vals || n == 0) return DBL_MAX;
    double sum_recip = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (fabs(k_vals[i]) < 1e-30) return 0.0;
        sum_recip += 1.0 / k_vals[i];
    }
    if (fabs(sum_recip) < 1e-30) return DBL_MAX;
    return 1.0 / sum_recip;
}

double torsional_spring_parallel(const double *k_vals, size_t n)
{
    /* k_eq = sum(k_i) */
    if (!k_vals || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += k_vals[i];
    return sum;
}

/* ===== Rotational Damping ===== */

double rotational_damping_torque(double coeff, double angular_vel)
{
    /* T_b = -b_t * omega */
    return -coeff * angular_vel;
}

/* ===== Rotational Kinetic Energy ===== */

double rotational_kinetic_energy(double inertia, double angular_vel)
{
    /* K_rot = 0.5 * J * omega^2 */
    if (inertia <= 0.0) return 0.0;
    return 0.5 * inertia * angular_vel * angular_vel;
}

double total_kinetic_energy(double m, double v, double J, double w)
{
    /* E_k = 0.5*m*v^2 + 0.5*J*w^2 */
    double K_trans = (m > 0) ? 0.5*m*v*v : 0.0;
    double K_rot = (J > 0) ? 0.5*J*w*w : 0.0;
    return K_trans + K_rot;
}

/* ===== Torsional Vibration Parameters ===== */

double torsional_natural_freq(double k, double J)
{
    /* omega_n = sqrt(k_t / J) */
    if (J <= 0.0 || k <= 0.0) return 0.0;
    return sqrt(k / J);
}

double torsional_natural_freq_hz(double k, double J)
{
    return torsional_natural_freq(k, J) / (2.0 * M_PI);
}

double torsional_damping_ratio(double b, double k, double J)
{
    /* zeta = b_t / (2*sqrt(k_t*J)) */
    if (J <= 0.0 || k <= 0.0) return (b > 0) ? DBL_MAX : 0.0;
    double denom = 2.0 * sqrt(k * J);
    if (denom < 1e-30) return DBL_MAX;
    return b / denom;
}

double torsional_damped_freq(double k, double J, double b)
{
    /* omega_d = omega_n * sqrt(1 - zeta^2) */
    double zeta = torsional_damping_ratio(b, k, J);
    if (zeta >= 1.0) return 0.0;
    double omega_n = torsional_natural_freq(k, J);
    return omega_n * sqrt(1.0 - zeta*zeta);
}

/* ===== Gear Ratio ===== */

double gear_ratio_calc(int teeth_driver, int teeth_driven)
{
    /* N = n_driven / n_driver */
    if (teeth_driver <= 0 || teeth_driven <= 0) return 1.0;
    return (double)teeth_driven / (double)teeth_driver;
}

double gear_output_speed(double omega_in, double ratio)
{
    /* omega_out = omega_in / N */
    if (fabs(ratio) < 1e-30) return DBL_MAX;
    return omega_in / ratio;
}

double gear_output_torque(double T_in, double ratio)
{
    /* T_out = N * T_in (ideal, 100% efficiency) */
    return T_in * ratio;
}

double gear_output_torque_eff(double T_in, double ratio, double efficiency)
{
    /* T_out = N * eta * T_in */
    if (efficiency < 0.0 || efficiency > 1.0) efficiency = 1.0;
    return T_in * ratio * efficiency;
}

/* ===== Reflected Parameters ===== */

double reflected_inertia(double J_load, double ratio)
{
    /* J_ref = J_load / N^2 */
    if (fabs(ratio) < 1e-30) return DBL_MAX;
    return J_load / (ratio * ratio);
}

double reflected_damping(double b_load, double ratio)
{
    /* b_ref = b_load / N^2 */
    if (fabs(ratio) < 1e-30) return DBL_MAX;
    return b_load / (ratio * ratio);
}

double reflected_stiffness(double k_load, double ratio)
{
    /* k_ref = k_load / N^2 */
    if (fabs(ratio) < 1e-30) return DBL_MAX;
    return k_load / (ratio * ratio);
}

double total_reflected_inertia(double J_motor, const double *J_loads,
                                const double *ratios, int n)
{
    double J_total = J_motor;
    for (int i = 0; i < n; i++) {
        J_total += reflected_inertia(J_loads[i], ratios[i]);
    }
    return J_total;
}

double optimal_gear_ratio(double J_motor, double J_load)
{
    /* N_opt = sqrt(J_load / J_motor) for max acceleration */
    if (J_motor <= 0.0) return DBL_MAX;
    return sqrt(J_load / J_motor);
}

double gear_train_total_ratio(const GearTrain *train)
{
    if (!train) return 1.0;
    return train->total_ratio;
}

int gear_train_build(const int *t_driver, const int *t_driven,
                     const double *eff, int n, GearTrain *train)
{
    if (!train || n <= 0) return -1;
    train->n_stages = n;
    train->stages = (GearPair*)malloc(n * sizeof(GearPair));
    if (!train->stages) return -1;
    train->total_ratio = 1.0;
    train->total_efficiency = 1.0;
    for (int i = 0; i < n; i++) {
        train->stages[i].teeth_driver = t_driver[i];
        train->stages[i].teeth_driven = t_driven[i];
        train->stages[i].efficiency = eff ? eff[i] : 1.0;
        train->stages[i].backlash_rad = 0.0;
        double ratio = gear_ratio_calc(t_driver[i], t_driven[i]);
        train->total_ratio *= ratio;
        train->total_efficiency *= eff ? eff[i] : 1.0;
    }
    return 0;
}

void gear_train_free(GearTrain *train)
{
    if (!train) return;
    free(train->stages);
    train->stages = NULL;
    train->n_stages = 0;
}

/* ===== Bearing Friction ===== */

double bearing_friction_torque(double T_static, double b_visc, double w)
{
    /* Simplified: T = T_s*sgn(w) + b_v*w */
    double sign = (w > 1e-30) ? 1.0 : ((w < -1e-30) ? -1.0 : 0.0);
    return T_static * sign + b_visc * w;
}

double bearing_power_loss(double T_friction, double w)
{
    /* P_loss = T_friction * omega */
    return fabs(T_friction * w);
}

/* ===== Rack and Pinion ===== */

double rack_pinion_position(double r, double theta)
{
    /* x = r * theta */
    if (r <= 0.0) return 0.0;
    return r * theta;
}

double rack_pinion_velocity(double r, double w)
{
    /* v = r * omega */
    if (r <= 0.0) return 0.0;
    return r * w;
}

double rack_pinion_force(double r, double T)
{
    /* F = T / r */
    if (r <= 0.0) return DBL_MAX;
    return T / r;
}

double rack_pinion_equiv_mass(double J_pinion, double r)
{
    /* m_eq = J_pinion / r^2 */
    if (r <= 0.0) return DBL_MAX;
    return J_pinion / (r * r);
}

/* ===== Lead Screw ===== */

double leadscrew_position(double lead, double theta)
{
    /* x = lead * theta / (2*pi) */
    if (lead <= 0.0) return 0.0;
    return lead * theta / (2.0 * M_PI);
}

double leadscrew_velocity(double lead, double w)
{
    /* v = lead * omega / (2*pi) */
    if (lead <= 0.0) return 0.0;
    return lead * w / (2.0 * M_PI);
}

double leadscrew_force(double lead, double T, double eta)
{
    /* F = 2*pi*T*eta / lead */
    if (lead <= 0.0) return DBL_MAX;
    if (eta < 0.0) eta = 0.0;
    if (eta > 1.0) eta = 1.0;
    return 2.0 * M_PI * T * eta / lead;
}

double leadscrew_reflected_inertia(double m, double lead)
{
    /* J_eq = m * (lead/(2*pi))^2 */
    if (m <= 0.0 || lead <= 0.0) return 0.0;
    double r_eff = lead / (2.0 * M_PI);
    return m * r_eff * r_eff;
}

/* ===== Rotating Unbalance ===== */

double centrifugal_force(double m_u, double e, double w)
{
    /* F = m_u * e * omega^2 */
    if (m_u <= 0.0 || e <= 0.0) return 0.0;
    return m_u * e * w * w;
}

double unbalance_displacement_amp(double m_u, double e, double M_total,
                                   double r, double zeta)
{
    /* X = (m_u*e/M) * r^2 / sqrt((1-r^2)^2 + (2*zeta*r)^2) */
    if (M_total <= 0.0) return 0.0;
    double num = (m_u * e / M_total) * r * r;
    double denom = sqrt((1.0 - r*r)*(1.0 - r*r) + 4.0*zeta*zeta*r*r);
    if (denom < 1e-30) return DBL_MAX;
    return num / denom;
}

double normalized_unbalance_response(double r, double zeta)
{
    /* MX/(m_u*e) = r^2 / sqrt((1-r^2)^2 + (2*zeta*r)^2) */
    double num = r * r;
    double denom = sqrt((1.0 - r*r)*(1.0 - r*r) + 4.0*zeta*zeta*r*r);
    if (denom < 1e-30) return DBL_MAX;
    return num / denom;
}

double critical_speed_rad_s(double k, double m_rotor)
{
    /* omega_cr = sqrt(k/m) = omega_n */
    if (m_rotor <= 0.0 || k <= 0.0) return 0.0;
    return sqrt(k / m_rotor);
}

double critical_speed_rpm(double k, double m_rotor)
{
    /* n_cr = 60*omega_cr/(2*pi) */
    return critical_speed_rad_s(k, m_rotor) * 60.0 / (2.0 * M_PI);
}

/* ===== Rotor Dynamics ===== */

double jeffcott_rotor_nat_freq(double k, double m)
{
    /* Jeffcott (Laval) rotor: omega_n = sqrt(k/m) */
    if (m <= 0.0 || k <= 0.0) return 0.0;
    return sqrt(k / m);
}

double rotor_amplification_factor(double zeta)
{
    /* Q = 1/(2*zeta) */
    if (zeta <= 0.0) return DBL_MAX;
    return 1.0 / (2.0 * zeta);
}

double rotor_orbit_radius(double m_u, double e, double M, double r, double zeta)
{
    /* Steady-state orbit: same as unbalance displacement */
    return unbalance_displacement_amp(m_u, e, M, r, zeta);
}

/* ===== Applications ===== */

double flywheel_energy(double J, double w)
{
    /* E = 0.5 * J * w^2 */
    if (J <= 0.0) return 0.0;
    return 0.5 * J * w * w;
}

double flywheel_energy_density(double J, double m, double w)
{
    /* E/m = (0.5*J*w^2) / m */
    if (m <= 0.0) return 0.0;
    return 0.5 * J * w * w / m;
}

double clutch_engagement_time(double J_load, double T_clutch, double w_target)
{
    /* t = J_load * w_target / T_clutch (constant torque acceleration) */
    if (T_clutch <= 0.0 || J_load <= 0.0) return DBL_MAX;
    return J_load * w_target / T_clutch;
}

double belt_drive_natural_freq(double belt_stiffness, double pulley_inertia1,
                                double pulley_inertia2)
{
    /* omega_n = sqrt(k_belt * (1/J1 + 1/J2)) */
    if (pulley_inertia1 <= 0.0 || pulley_inertia2 <= 0.0) return 0.0;
    double J_eq = (pulley_inertia1 * pulley_inertia2) / (pulley_inertia1 + pulley_inertia2);
    return sqrt(belt_stiffness / J_eq);
}

double shaft_whirling_speed(double E_mod, double area_I, double m_per_length, double L,
                             int end_condition)
{
    /* Simplified whirling speed (first critical).
     * omega = C * sqrt(E_mod*area_I/(m_per_length*L^4))
     * C = pi^2 (simply-supported), ~3.52 (cantilever), ~22.4 (clamped-clamped)
     */
    if (L <= 0.0 || m_per_length <= 0.0 || E_mod <= 0.0 || area_I <= 0.0) return 0.0;
    double C;
    switch (end_condition) {
        case 0: C = M_PI * M_PI; break;           /* simply-supported */
        case 1: C = 3.516 * 3.516; break;         /* cantilever */
        case 2: C = 4.730 * 4.730; break;         /* clamped-clamped */
        case 3: C = 3.927 * 3.927; break;         /* clamped-pinned */
        default: C = M_PI * M_PI; break;
    }
    return C * sqrt(E_mod * area_I / (m_per_length * L*L*L*L));
}
