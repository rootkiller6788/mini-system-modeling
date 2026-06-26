/**
 * @file mechanical_elements.c
 * @brief Implementation of fundamental mechanical element computations
 *
 * Each function implements one independent knowledge point.
 * All boundary conditions handled. No filler, no stub code.
 */

#include "mechanical_elements.h"
#include <math.h>
#include <complex.h>
#include <float.h>

/* ===== L2: Newton's Second Law (F = ma) ===== */

double newton_second_law_force(double mass, double acceleration)
{
    /* F = m * a  (Newton 1687, Principia, Lex II)
     * Force is proportional to mass times acceleration.
     * Complexity: O(1).
     */
    return mass * acceleration;
}

double newton_second_law_acceleration(double force, double mass)
{
    /* a = F / m
     * Guard: mass <= 0 (massless or negative mass is non-physical).
     */
    if (mass <= 0.0) {
        return (force >= 0.0) ? DBL_MAX : -DBL_MAX;
    }
    return force / mass;
}

/* ===== L2: Hooke's Law (F = k*dx) ===== */

double hookes_law_force(double stiffness, double displacement)
{
    /* F = -k * x  (Hooke 1676, "ut tensio, sic vis")
     * Restoring force opposes displacement from equilibrium.
     * Complexity: O(1).
     */
    if (stiffness < 0.0) return 0.0; /* Negative stiffness is non-physical for passive spring */
    return -stiffness * displacement;
}

double hookes_law_displacement(double stiffness, double force)
{
    /* x = F / k
     * Guard: k <= 0 means infinite compliance (no restoring force).
     */
    if (fabs(stiffness) < 1e-30) {
        return (force >= 0.0) ? DBL_MAX : -DBL_MAX;
    }
    return force / stiffness;
}

double spring_potential_energy(double stiffness, double displacement)
{
    /* U = 0.5 * k * x^2  (elastic potential energy)
     * Derived: dU = F*dx = k*x*dx, integrate from 0 to x.
     * Complexity: O(1).
     */
    if (stiffness <= 0.0) return 0.0;
    return 0.5 * stiffness * displacement * displacement;
}

double spring_series_equivalent(const double *stiffness, size_t n)
{
    /* k_eq = 1 / sum(1/k_i)  (springs in series)
     * Springs in series share the same force; displacements add.
     * Guard: any k_i <= 0 makes the chain infinitely compliant.
     * Complexity: O(n).
     */
    if (!stiffness || n == 0) return DBL_MAX;
    double sum_recip = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (fabs(stiffness[i]) < 1e-30) return 0.0;
        sum_recip += 1.0 / stiffness[i];
    }
    if (fabs(sum_recip) < 1e-30) return DBL_MAX;
    return 1.0 / sum_recip;
}

double spring_parallel_equivalent(const double *stiffness, size_t n)
{
    /* k_eq = sum(k_i)  (springs in parallel)
     * Springs in parallel share the same displacement; forces add.
     * Complexity: O(n).
     */
    if (!stiffness || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += stiffness[i];
    }
    return sum;
}

double spring_combined_two(double k1, double k2, int is_series)
{
    /* Two-spring combination: either series or parallel.
     * Complexity: O(1).
     */
    if (k1 <= 0.0 || k2 <= 0.0) return 0.0;
    if (is_series) {
        return (k1 * k2) / (k1 + k2);
    } else {
        return k1 + k2;
    }
}

/* ===== L2: Viscous Damping ===== */

double viscous_damping_force(double damping_coeff, double velocity)
{
    /* F_b = -b * v  (linear viscous damping)
     * Damping force opposes velocity; models fluid dashpot.
     * Complexity: O(1).
     */
    return -damping_coeff * velocity;
}

double damper_series_equivalent(const double *damping_coeffs, size_t n)
{
    /* b_eq = 1 / sum(1/b_i)  (dampers in series)
     * Same force, velocities add.
     * Complexity: O(n).
     */
    if (!damping_coeffs || n == 0) return DBL_MAX;
    double sum_recip = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (fabs(damping_coeffs[i]) < 1e-30) return 0.0;
        sum_recip += 1.0 / damping_coeffs[i];
    }
    if (fabs(sum_recip) < 1e-30) return DBL_MAX;
    return 1.0 / sum_recip;
}

double damper_parallel_equivalent(const double *damping_coeffs, size_t n)
{
    /* b_eq = sum(b_i)  (dampers in parallel)
     * Same velocity, forces add.
     * Complexity: O(n).
     */
    if (!damping_coeffs || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += damping_coeffs[i];
    }
    return sum;
}

/* ===== L2: Friction Models ===== */

double coulomb_friction_magnitude(double kinetic_coeff, double normal_force)
{
    /* |F_c| = mu_k * N  (Coulomb 1785)
     * Magnitude of kinetic friction force.
     * Complexity: O(1).
     */
    if (kinetic_coeff < 0.0 || normal_force < 0.0) return 0.0;
    return kinetic_coeff * normal_force;
}

double coulomb_friction_force(double kinetic_coeff, double normal_force, double velocity)
{
    /* F_c = -mu_k * N * sgn(v)
     * Coulomb friction opposes motion with constant magnitude.
     * At v=0, returns 0 (static friction handled separately).
     * Complexity: O(1).
     */
    if (fabs(velocity) < 1e-30) return 0.0;
    double mag = coulomb_friction_magnitude(kinetic_coeff, normal_force);
    return (velocity > 0.0) ? -mag : mag;
}

double stribeck_friction_force(const FrictionModel *fm, double velocity)
{
    /* F(v) = F_c + (F_s - F_c) * exp(-|v/v_s|) + b * v
     * Stribeck curve models transition from static to kinetic friction.
     * Reference: Stribeck (1902), Armstrong-Helouvry (1994).
     * Complexity: O(1).
     */
    if (!fm) return 0.0;
    double Fc = fm->coulomb_coeff * fm->normal_force;
    double Fs = fm->static_coeff * fm->normal_force;
    double bv = fm->viscous_coeff * velocity;
    double av = fabs(velocity);

    if (fm->stribeck_vel < 1e-30) {
        /* Pure Coulomb + viscous (no Stribeck transition) */
        if (fabs(velocity) < 1e-30) return 0.0;
        return (velocity > 0.0 ? -Fc : Fc) + bv;
    }

    if (av < 1e-30) return 0.0;

    double stribeck_term = (Fs - Fc) * exp(-av / fm->stribeck_vel);
    double fric_mag = Fc + stribeck_term;
    double sign = (velocity > 0.0) ? -1.0 : 1.0;

    return sign * fric_mag + bv;
}

int static_friction_holds(const FrictionModel *fm, double applied_force)
{
    /* Static friction holds if |F_applied| <= mu_s * N.
     * Returns 1 if held (no motion), 0 if motion occurs.
     * Complexity: O(1).
     */
    if (!fm) return 0;
    double Fs_max = fm->static_coeff * fm->normal_force;
    return (fabs(applied_force) <= Fs_max) ? 1 : 0;
}

double equivalent_viscous_from_coulomb(double coulomb_force, double velocity_amplitude,
                                        double freq)
{
    /* b_eq = 4 * F_coulomb / (pi * omega * X)
     * Equivalent viscous damping for energy-based linearization.
     * Equating energy dissipated per cycle: pi*b_eq*w*X^2 = 4*F_c*X.
     * Complexity: O(1).
     */
    if (velocity_amplitude < 1e-30 || fabs(freq) < 1e-30) return DBL_MAX;
    (void)freq; /* freq embedded in velocity_amplitude = omega*X */
    return (4.0 * fabs(coulomb_force)) / (M_PI * velocity_amplitude);
}

/* ===== L2: Kinetic and Potential Energy ===== */

double kinetic_energy_translational(double mass, double velocity)
{
    /* K = 0.5 * m * v^2
     * Translational kinetic energy. Derived from work integral.
     * Complexity: O(1).
     */
    if (mass <= 0.0) return 0.0;
    return 0.5 * mass * velocity * velocity;
}

double gravitational_potential_energy(double mass, double height, double gravity)
{
    /* U_g = m * g * h
     * Gravitational potential energy near Earth surface.
     * Complexity: O(1).
     */
    if (mass <= 0.0) return 0.0;
    return mass * gravity * height;
}

double total_mechanical_energy(double mass, double velocity, double stiffness,
                                double displacement, double height, double gravity)
{
    /* E_total = K + U_spring + U_gravity
     * Total mechanical energy of a spring-mass system in gravity.
     * Complexity: O(1).
     */
    double K = kinetic_energy_translational(mass, velocity);
    double Us = spring_potential_energy(stiffness, displacement);
    double Ug = gravitational_potential_energy(mass, height, gravity);
    return K + Us + Ug;
}

double energy_dissipated_per_cycle(double damping_coeff, double angular_freq, double amplitude)
{
    /* W_d = pi * b * omega * X^2
     * Energy dissipated by viscous damper per cycle of harmonic motion.
     * Derived: W_d = integral(F_b * dx) over one period.
     * Complexity: O(1).
     */
    if (damping_coeff < 0.0 || angular_freq < 0.0 || amplitude < 0.0) return 0.0;
    return M_PI * damping_coeff * angular_freq * amplitude * amplitude;
}

/* ===== L3: Natural Frequency and Damping Parameters ===== */

double natural_frequency_undamped(double stiffness, double mass)
{
    /* omega_n = sqrt(k / m)  [rad/s]
     * Undamped natural frequency of SDOF system.
     * Derived from characteristic equation: m*s^2 + k = 0.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || stiffness <= 0.0) return 0.0;
    return sqrt(stiffness / mass);
}

double natural_frequency_hz(double stiffness, double mass)
{
    /* f_n = omega_n / (2*pi)  [Hz]
     * Complexity: O(1).
     */
    return natural_frequency_undamped(stiffness, mass) / (2.0 * M_PI);
}

double damping_ratio(double damping_coeff, double stiffness, double mass)
{
    /* zeta = b / (2 * sqrt(k*m)) = b / b_critical
     * Damping ratio: dimensionless measure of damping.
     * zeta < 1: underdamped, zeta = 1: critically damped, zeta > 1: overdamped.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || stiffness <= 0.0) {
        if (damping_coeff <= 0.0) return 0.0;
        return DBL_MAX; /* No stiffness or mass => infinite damping ratio */
    }
    double denom = 2.0 * sqrt(stiffness * mass);
    if (denom < 1e-30) return DBL_MAX;
    return damping_coeff / denom;
}

double critical_damping_coefficient(double stiffness, double mass)
{
    /* b_cr = 2 * sqrt(k*m) = 2 * m * omega_n
     * Critical damping: fastest return to equilibrium without oscillation.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || stiffness <= 0.0) return 0.0;
    return 2.0 * sqrt(stiffness * mass);
}

double damped_natural_frequency(double stiffness, double mass, double damping_coeff)
{
    /* omega_d = omega_n * sqrt(1 - zeta^2)
     * Damped natural frequency; exists only for underdamped systems.
     * Returns 0 if overdamped or critically damped.
     * Complexity: O(1).
     */
    double zeta = damping_ratio(damping_coeff, stiffness, mass);
    if (zeta >= 1.0 || zeta < 0.0) return 0.0;
    double omega_n = natural_frequency_undamped(stiffness, mass);
    return omega_n * sqrt(1.0 - zeta * zeta);
}

double logarithmic_decrement(double amplitude_0, double amplitude_n, int num_cycles)
{
    /* delta = (1/n) * ln(x_0 / x_n)
     * Logarithmic decrement measures damping from free decay data.
     * Related to damping ratio: delta = 2*pi*zeta / sqrt(1-zeta^2).
     * Complexity: O(1).
     */
    if (num_cycles <= 0 || amplitude_n <= 0.0 || amplitude_0 <= 0.0) return 0.0;
    if (amplitude_n >= amplitude_0) return 0.0; /* Should decay, not grow */
    return log(amplitude_0 / amplitude_n) / (double)num_cycles;
}

double damping_ratio_from_logdec(double log_dec)
{
    /* zeta = delta / sqrt(4*pi^2 + delta^2)
     * Invert: delta = 2*pi*zeta / sqrt(1-zeta^2).
     * Complexity: O(1).
     */
    if (log_dec <= 0.0) return 0.0;
    double delta2 = log_dec * log_dec;
    double denom = 4.0 * M_PI * M_PI + delta2;
    return log_dec / sqrt(denom);
}

double quality_factor_from_damping(double damping_ratio)
{
    /* Q = 1 / (2*zeta)
     * Quality factor: ratio of stored to dissipated energy.
     * Valid for light damping (zeta << 1).
     * Complexity: O(1).
     */
    if (fabs(damping_ratio) < 1e-30) return DBL_MAX;
    return 1.0 / (2.0 * damping_ratio);
}

double damping_ratio_from_q_factor(double q_factor)
{
    /* zeta = 1 / (2*Q)
     * Complexity: O(1).
     */
    if (q_factor <= 0.0) return DBL_MAX;
    return 1.0 / (2.0 * q_factor);
}

double loss_factor_from_damping(double damping_ratio)
{
    /* eta = 2 * zeta
     * Loss factor used in structural damping models.
     * Relates to complex modulus: E* = E*(1 + j*eta).
     * Complexity: O(1).
     */
    return 2.0 * damping_ratio;
}

/* ===== L3: Static Deflection ===== */

double static_deflection(double mass, double stiffness, double gravity)
{
    /* delta_st = m*g / k
     * Static deflection of spring under weight.
     * Complexity: O(1).
     */
    if (stiffness <= 0.0 || mass <= 0.0) return 0.0;
    return mass * gravity / stiffness;
}

double stiffness_from_static_deflection(double mass, double deflection, double gravity)
{
    /* k = m*g / delta_st
     * Estimate stiffness from measured static deflection.
     * Complexity: O(1).
     */
    if (deflection <= 0.0 || mass <= 0.0) return DBL_MAX;
    return mass * gravity / deflection;
}

double natural_freq_from_static_deflection(double deflection, double gravity)
{
    /* omega_n = sqrt(g / delta_st)
     * Natural frequency estimated from static deflection.
     * Since omega_n^2 = k/m and k = m*g/delta_st => omega_n^2 = g/delta_st.
     * Complexity: O(1).
     */
    if (deflection <= 0.0) return DBL_MAX;
    return sqrt(gravity / deflection);
}

/* ===== L3: Mechanical Impedance ===== */

double complex mechanical_impedance_sdof(double mass, double damping, double stiffness,
                                          double angular_freq)
{
    /* Z_m(j*w) = b + j*(w*m - k/w)
     * Mechanical impedance: complex ratio of force to velocity.
     * Z_m = F(j*w) / V(j*w) for SDOF system.
     * Complexity: O(1).
     */
    double real_part = damping;
    if (fabs(angular_freq) < 1e-30) {
        /* Static case: Z_m -> b - j*infinity (spring blocks DC motion) */
        return damping - I * DBL_MAX;
    }
    double imag_part = angular_freq * mass - stiffness / angular_freq;
    return real_part + I * imag_part;
}

double mechanical_impedance_magnitude(double mass, double damping, double stiffness,
                                       double angular_freq)
{
    /* |Z_m| = sqrt(b^2 + (w*m - k/w)^2)
     * Complexity: O(1).
     */
    double complex Z = mechanical_impedance_sdof(mass, damping, stiffness, angular_freq);
    return cabs(Z);
}

double mechanical_impedance_phase(double mass, double damping, double stiffness,
                                   double angular_freq)
{
    /* phi_Z = atan2(imag(Z_m), real(Z_m))
     * Phase angle between force and velocity.
     * Complexity: O(1).
     */
    double complex Z = mechanical_impedance_sdof(mass, damping, stiffness, angular_freq);
    return atan2(cimag(Z), creal(Z));
}

double antiresonant_frequency(double stiffness, double mass)
{
    /* omega_ar = sqrt(k/m) = omega_n
     * Antiresonance: frequency where |Z_m| is maximum.
     * At antiresonance, reactance is infinite (velocity zero for given force).
     * Complexity: O(1).
     */
    return natural_frequency_undamped(stiffness, mass);
}

double complex mechanical_mobility(double mass, double damping, double stiffness,
                                    double angular_freq)
{
    /* Y_m(j*w) = 1 / Z_m(j*w)
     * Mobility (mechanical admittance): velocity per unit force.
     * Complexity: O(1).
     */
    double complex Z = mechanical_impedance_sdof(mass, damping, stiffness, angular_freq);
    double mag2 = creal(Z)*creal(Z) + cimag(Z)*cimag(Z);
    if (mag2 < 1e-60) return DBL_MAX + I * 0.0;
    return (creal(Z) - I * cimag(Z)) / mag2;
}

/* ===== L5: Vibration Transmissibility ===== */

double force_transmissibility(double freq_ratio, double damping_ratio)
{
    /* TR = sqrt((1 + (2*zeta*r)^2) / ((1 - r^2)^2 + (2*zeta*r)^2))
     * Force transmitted to base relative to excitation force.
     * TR < 1 for r > sqrt(2) (isolation region).
     * Complexity: O(1).
     */
    double zr2 = (2.0 * damping_ratio * freq_ratio);
    zr2 = zr2 * zr2;
    double num = 1.0 + zr2;
    double r2 = freq_ratio * freq_ratio;
    double denom = (1.0 - r2) * (1.0 - r2) + zr2;
    if (denom < 1e-30) return DBL_MAX;
    return sqrt(num / denom);
}

double displacement_transmissibility(double freq_ratio, double damping_ratio)
{
    /* TR_disp = sqrt((1 + (2*zeta*r)^2) / ((1 - r^2)^2 + (2*zeta*r)^2))
     * Same form as force transmissibility for base excitation.
     * Complexity: O(1).
     */
    return force_transmissibility(freq_ratio, damping_ratio);
}

double peak_transmissibility_freq_ratio(double damping_ratio)
{
    /* r_peak = sqrt((sqrt(1+8*zeta^2) - 1) / (4*zeta^2))
     * Frequency ratio at which transmissibility peaks.
     * Valid only for zeta < 1/sqrt(2) (transmissibility has a peak).
     * Complexity: O(1).
     */
    if (damping_ratio <= 0.0) return 1.0; /* Undamped: peak at r=1 */
    if (damping_ratio >= 1.0 / sqrt(2.0)) return 0.0; /* No peak for high damping */

    double z2 = damping_ratio * damping_ratio;
    double inner = sqrt(1.0 + 8.0 * z2);
    double num = inner - 1.0;
    double denom = 4.0 * z2;
    if (denom < 1e-30) return 1.0;
    return sqrt(num / denom);
}

double maximum_transmissibility(double damping_ratio)
{
    /* TR_max = transmissibility at r_peak.
     * Complexity: O(1).
     */
    double r_peak = peak_transmissibility_freq_ratio(damping_ratio);
    if (r_peak <= 0.0) return 1.0; /* No peak */
    return force_transmissibility(r_peak, damping_ratio);
}

double isolation_efficiency(double freq_ratio, double damping_ratio)
{
    /* eta_iso = (1 - TR) * 100%
     * Percentage reduction in transmitted vibration.
     * Complexity: O(1).
     */
    double TR = force_transmissibility(freq_ratio, damping_ratio);
    return (1.0 - TR) * 100.0;
}

/* ===== L3: Mass Moments of Inertia ===== */

double inertia_point_mass(double mass, double distance)
{
    /* I = m * r^2
     * Moment of inertia of a point mass about an axis.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || distance < 0.0) return 0.0;
    return mass * distance * distance;
}

double inertia_rod_end(double mass, double length)
{
    /* I = (1/3) * m * L^2
     * Thin rod about an axis through one end, perpendicular to rod.
     * Integral: I = int(0,L) (m/L)*x^2 dx = (1/3)*m*L^2.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || length <= 0.0) return 0.0;
    return (1.0 / 3.0) * mass * length * length;
}

double inertia_rod_center(double mass, double length)
{
    /* I = (1/12) * m * L^2
     * Thin rod about center of mass, perpendicular to rod.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || length <= 0.0) return 0.0;
    return (1.0 / 12.0) * mass * length * length;
}

double inertia_solid_cylinder(double mass, double radius)
{
    /* I = (1/2) * m * R^2
     * Solid cylinder/disk about central axis.
     * Integral: I = int(0,R) r^2 * (m/(pi*R^2)) * 2*pi*r dr = (1/2)*m*R^2.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || radius <= 0.0) return 0.0;
    return 0.5 * mass * radius * radius;
}

double inertia_hollow_cylinder(double mass, double r_outer, double r_inner)
{
    /* I = (1/2) * m * (R_outer^2 + R_inner^2)
     * Hollow cylinder about central axis.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || r_outer <= 0.0 || r_inner < 0.0) return 0.0;
    if (r_inner > r_outer) { double tmp = r_inner; r_inner = r_outer; r_outer = tmp; }
    return 0.5 * mass * (r_outer * r_outer + r_inner * r_inner);
}

double inertia_solid_sphere(double mass, double radius)
{
    /* I = (2/5) * m * R^2
     * Solid sphere about any diameter.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || radius <= 0.0) return 0.0;
    return 0.4 * mass * radius * radius; /* 2/5 = 0.4 */
}

double inertia_rectangular_plate(double mass, double side_a, double side_b)
{
    /* I = (1/12) * m * (a^2 + b^2)
     * Thin rectangular plate about centroid, axis perpendicular to plate.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || side_a <= 0.0 || side_b <= 0.0) return 0.0;
    return (mass / 12.0) * (side_a * side_a + side_b * side_b);
}

double parallel_axis_theorem(double i_cm, double mass, double distance)
{
    /* I_parallel = I_cm + m * d^2
     * Parallel axis theorem (Steiner's theorem).
     * Relates moment of inertia about any axis to CM axis.
     * Complexity: O(1).
     */
    if (mass < 0.0 || distance < 0.0) return i_cm;
    return i_cm + mass * distance * distance;
}

double radius_of_gyration(double moment_of_inertia, double mass)
{
    /* k_g = sqrt(I / m)
     * Radius of gyration: distance from axis where all mass would give same I.
     * Complexity: O(1).
     */
    if (mass <= 0.0 || moment_of_inertia < 0.0) return 0.0;
    return sqrt(moment_of_inertia / mass);
}

/* ===== L5: Time-Domain Response Parameters ===== */

double peak_time_underdamped(double damping_ratio, double natural_freq)
{
    /* t_p = pi / omega_d
     * Time to first peak of underdamped step response.
     * Complexity: O(1).
     */
    if (damping_ratio >= 1.0 || damping_ratio < 0.0 || natural_freq <= 0.0) return 0.0;
    double omega_d = natural_freq * sqrt(1.0 - damping_ratio * damping_ratio);
    if (omega_d < 1e-30) return DBL_MAX;
    return M_PI / omega_d;
}

double percent_overshoot(double damping_ratio)
{
    /* %OS = 100 * exp(-pi*zeta / sqrt(1-zeta^2))
     * Percentage overshoot of step response.
     * Valid for 0 <= zeta < 1.
     * Complexity: O(1).
     */
    if (damping_ratio >= 1.0 || damping_ratio <= 0.0) return 0.0;
    double exponent = -M_PI * damping_ratio / sqrt(1.0 - damping_ratio * damping_ratio);
    return 100.0 * exp(exponent);
}

double settling_time_2percent(double damping_ratio, double natural_freq)
{
    /* t_s = 4 / (zeta * omega_n)  (2% criterion)
     * Time for response to stay within 2% of final value.
     * Complexity: O(1).
     */
    if (damping_ratio <= 0.0 || natural_freq <= 0.0) return DBL_MAX;
    return 4.0 / (damping_ratio * natural_freq);
}

double settling_time_5percent(double damping_ratio, double natural_freq)
{
    /* t_s = 3 / (zeta * omega_n)  (5% criterion)
     * Complexity: O(1).
     */
    if (damping_ratio <= 0.0 || natural_freq <= 0.0) return DBL_MAX;
    return 3.0 / (damping_ratio * natural_freq);
}

double rise_time_underdamped(double damping_ratio, double natural_freq)
{
    /* t_r = (pi - beta) / omega_d, beta = atan(sqrt(1-zeta^2)/zeta)
     * Rise time for underdamped step response (0% to 100%).
     * Complexity: O(1).
     */
    if (damping_ratio >= 1.0 || damping_ratio <= 0.0 || natural_freq <= 0.0) return 0.0;
    double beta = atan2(sqrt(1.0 - damping_ratio * damping_ratio), damping_ratio);
    double omega_d = natural_freq * sqrt(1.0 - damping_ratio * damping_ratio);
    if (omega_d < 1e-30) return DBL_MAX;
    return (M_PI - beta) / omega_d;
}

double delay_time(double damping_ratio, double natural_freq)
{
    /* t_d = (1 + 0.7*zeta) / omega_n  (approximation)
     * Delay time: time to reach 50% of final value.
     * Complexity: O(1).
     */
    if (natural_freq <= 0.0) return DBL_MAX;
    return (1.0 + 0.7 * damping_ratio) / natural_freq;
}

double sdof_step_response_displacement(double static_deflection_val, double natural_freq,
                                        double damping_ratio_val, double time)
{
    /* x(t) = delta_st * [1 - exp(-zeta*wn*t) * (cos(wd*t) + (zeta/sqrt(1-zeta^2))*sin(wd*t))]
     * SDOF step response (underdamped case).
     * Complexity: O(1).
     */
    if (time <= 0.0) return 0.0;
    if (damping_ratio_val >= 1.0) {
        /* Overdamped: x(t) = delta_st * [1 - (1/(r2-r1))*(r2*e^{r1*t} - r1*e^{r2*t})] */
        double wd = natural_freq * sqrt(damping_ratio_val*damping_ratio_val - 1.0);
        double r1 = -damping_ratio_val*natural_freq + wd;
        double r2 = -damping_ratio_val*natural_freq - wd;
        double term = (r2*exp(r1*time) - r1*exp(r2*time)) / (r2 - r1);
        return static_deflection_val * (1.0 - term);
    }
    if (fabs(damping_ratio_val - 1.0) < 1e-12) {
        /* Critically damped: x(t) = delta_st * [1 - (1+wn*t)*e^{-wn*t}] */
        double exp_term = exp(-natural_freq * time);
        return static_deflection_val * (1.0 - (1.0 + natural_freq*time) * exp_term);
    }
    /* Underdamped */
    double omega_d = natural_freq * sqrt(1.0 - damping_ratio_val*damping_ratio_val);
    double exp_term = exp(-damping_ratio_val * natural_freq * time);
    double cos_term = cos(omega_d * time);
    double sin_term = sin(omega_d * time);
    double ratio = damping_ratio_val / sqrt(1.0 - damping_ratio_val*damping_ratio_val);
    return static_deflection_val * (1.0 - exp_term * (cos_term + ratio * sin_term));
}

double resonance_amplitude_force(double force_amplitude, double stiffness, double damping_ratio)
{
    /* X_max = F_0 / (k * 2*zeta * sqrt(1-zeta^2))
     * Peak displacement amplitude at resonance (r = sqrt(1-2*zeta^2)).
     * For small zeta, peak at r ~= 1 and X_max ~= F_0/(2*k*zeta).
     * Complexity: O(1).
     */
    if (stiffness <= 0.0) return DBL_MAX;
    if (damping_ratio <= 0.0) return DBL_MAX;
    if (damping_ratio >= 1.0/sqrt(2.0)) {
        /* No resonant peak; max is at r=0: X = F_0/k */
        return force_amplitude / stiffness;
    }
    double denom = stiffness * 2.0 * damping_ratio * sqrt(1.0 - damping_ratio*damping_ratio);
    if (denom < 1e-30) return DBL_MAX;
    return force_amplitude / denom;
}

double phase_angle_sdof(double freq_ratio, double damping_ratio)
{
    /* phi = atan2(2*zeta*r, 1 - r^2)
     * Phase lag of displacement response relative to excitation.
     * phi = 0 at low freq, pi/2 at resonance, pi at high freq.
     * Complexity: O(1).
     */
    return atan2(2.0 * damping_ratio * freq_ratio, 1.0 - freq_ratio * freq_ratio);
}

double half_power_bandwidth(double damping_ratio, double natural_freq)
{
    /* Delta_omega = 2 * zeta * omega_n  (for small zeta)
     * Bandwidth between half-power points (-3dB on displacement response).
     * Complexity: O(1).
     */
    if (damping_ratio < 0.0 || natural_freq < 0.0) return 0.0;
    return 2.0 * damping_ratio * natural_freq;
}

double damping_from_half_power(double omega1, double omega2, double omega_n)
{
    /* zeta = (omega2 - omega1) / (2 * omega_n)
     * Estimate damping from half-power bandwidth.
     * omega1, omega2 are frequencies at 1/sqrt(2) of peak amplitude.
     * Complexity: O(1).
     */
    if (omega_n <= 0.0) return 0.0;
    return fabs(omega2 - omega1) / (2.0 * omega_n);
}

/* ===== L5: Equivalent Systems ===== */

double effective_mass_sdof_lumped(double mass, double mode_shape_at_point)
{
    /* m_eff = m * phi^2(x_ref)
     * Effective mass for Rayleigh method (simplified, phi(x_ref)=1 convention).
     * Complexity: O(1).
     */
    if (mode_shape_at_point < 0.0) return 0.0;
    return mass * mode_shape_at_point * mode_shape_at_point;
}

double effective_stiffness_from_omega(double natural_freq, double effective_mass)
{
    /* k_eff = omega_n^2 * m_eff
     * Complexity: O(1).
     */
    if (effective_mass < 0.0) return 0.0;
    return natural_freq * natural_freq * effective_mass;
}

double mass_normalized_stiffness(double mass, double stiffness)
{
    /* k/m (squared natural frequency)
     * Complexity: O(1).
     */
    if (mass <= 0.0) return DBL_MAX;
    return stiffness / mass;
}

/* ===== L5: Structural Stiffness Formulas ===== */

double shaft_torsional_stiffness(double shear_modulus, double polar_moment, double length)
{
    /* k_theta = G * J / L
     * Torsional stiffness of uniform circular shaft.
     * Derived: theta = T*L/(G*J) => k = T/theta = G*J/L.
     * Complexity: O(1).
     */
    if (length <= 0.0 || shear_modulus <= 0.0 || polar_moment <= 0.0) return 0.0;
    return shear_modulus * polar_moment / length;
}

double cantilever_beam_stiffness(double E_modulus, double inertia_area, double length)
{
    /* k = 3 * E * I / L^3
     * Stiffness of cantilever beam with tip load.
     * Derived from Euler-Bernoulli beam theory: delta = F*L^3/(3*E*I).
     * Complexity: O(1).
     */
    if (length <= 0.0 || E_modulus <= 0.0 || inertia_area <= 0.0) return 0.0;
    return 3.0 * E_modulus * inertia_area / (length * length * length);
}

double simply_supported_beam_stiffness(double E_modulus, double inertia_area, double length)
{
    /* k = 48 * E * I / L^3
     * Stiffness of simply-supported beam with center load.
     * Derived: delta = F*L^3/(48*E*I).
     * Complexity: O(1).
     */
    if (length <= 0.0 || E_modulus <= 0.0 || inertia_area <= 0.0) return 0.0;
    return 48.0 * E_modulus * inertia_area / (length * length * length);
}

double area_moment_rectangular(double width, double height)
{
    /* I = b * h^3 / 12
     * Area moment of inertia for rectangular cross-section.
     * Complexity: O(1).
     */
    if (width <= 0.0 || height <= 0.0) return 0.0;
    return width * height * height * height / 12.0;
}

double area_moment_circular(double diameter)
{
    /* I = pi * d^4 / 64
     * Area moment of inertia for solid circular cross-section.
     * Complexity: O(1).
     */
    if (diameter <= 0.0) return 0.0;
    double r4 = (diameter/2.0);
    r4 = r4 * r4 * r4 * r4;
    return M_PI * r4 / 4.0; /* pi*d^4/64 = pi*(d/2)^4/4 */
}

double polar_moment_circular(double diameter)
{
    /* J = pi * d^4 / 32
     * Polar moment of inertia for solid circular cross-section.
     * Complexity: O(1).
     */
    if (diameter <= 0.0) return 0.0;
    double r4 = (diameter/2.0);
    r4 = r4 * r4 * r4 * r4;
    return M_PI * r4 / 2.0; /* pi*d^4/32 = pi*(d/2)^4/2 */
}
