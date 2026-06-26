/**
 * mini-fluid-thermal-modeling — Fluid Mechanics Implementation
 *
 * Conservation laws (Bernoulli, continuity), Navier-Stokes
 * analytical solutions, pipe flow (Darcy-Weisbach, Moody),
 * boundary layer theory, dimensional analysis.
 *
 * Knowledge Levels: L3-L5
 */

#include "fluid_mechanics.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =====================================================================
 * L3: Engineering Quantities — Typical Value Ranges
 * ===================================================================== */

ReynoldsRange typical_reynolds_range(const char *application)
{
    ReynoldsRange r = {0, 0, ""};
    if (!application) return r;

    if (strcmp(application, "pipe_water") == 0) {
        r.min_re = 1000.0; r.max_re = 1.0e6;
        r.description = "Water pipe (domestic to industrial)";
    } else if (strcmp(application, "pipe_air") == 0) {
        r.min_re = 5000.0; r.max_re = 5.0e5;
        r.description = "Air duct (HVAC)";
    } else if (strcmp(application, "open_channel") == 0) {
        r.min_re = 5000.0; r.max_re = 1.0e8;
        r.description = "Open channel flow (rivers, canals)";
    } else if (strcmp(application, "aircraft_wing") == 0) {
        r.min_re = 5.0e5; r.max_re = 2.0e8;
        r.description = "Aircraft wing (takeoff to cruise)";
    } else if (strcmp(application, "automobile") == 0) {
        r.min_re = 1.0e5; r.max_re = 2.0e7;
        r.description = "Automobile aerodynamics";
    } else if (strcmp(application, "microchannel") == 0) {
        r.min_re = 0.1; r.max_re = 2000.0;
        r.description = "Microchannel flow (lab-on-chip)";
    } else if (strcmp(application, "blood_vessel") == 0) {
        r.min_re = 0.001; r.max_re = 4000.0;
        r.description = "Blood flow (capillaries to aorta)";
    } else if (strcmp(application, "heat_exchanger_tube") == 0) {
        r.min_re = 3000.0; r.max_re = 50000.0;
        r.description = "Heat exchanger tube-side";
    } else if (strcmp(application, "wind_turbine") == 0) {
        r.min_re = 1.0e5; r.max_re = 2.0e7;
        r.description = "Wind turbine blade";
    }
    return r;
}

double typical_prandtl(const char *fluid_name)
{
    if (!fluid_name) return NAN;
    if (strcmp(fluid_name, "air") == 0) return 0.71;
    if (strcmp(fluid_name, "water") == 0) return 5.83;
    if (strcmp(fluid_name, "oil_engine") == 0) return 6400.0;
    if (strcmp(fluid_name, "glycerin") == 0) return 12600.0;
    if (strcmp(fluid_name, "mercury") == 0) return 0.025;
    if (strcmp(fluid_name, "ethanol") == 0) return 15.7;
    if (strcmp(fluid_name, "r134a") == 0) return 0.82;
    return NAN;
}

void typical_h_range(const char *scenario, double *h_min, double *h_max)
{
    if (!scenario || !h_min || !h_max) return;
    *h_min = 0.0; *h_max = 0.0;
    if (strcmp(scenario, "natural_air") == 0) {
        *h_min = 2.0; *h_max = 25.0;
    } else if (strcmp(scenario, "forced_air") == 0) {
        *h_min = 25.0; *h_max = 250.0;
    } else if (strcmp(scenario, "forced_water") == 0) {
        *h_min = 500.0; *h_max = 20000.0;
    } else if (strcmp(scenario, "boiling_water") == 0) {
        *h_min = 2500.0; *h_max = 100000.0;
    } else if (strcmp(scenario, "condensing_steam") == 0) {
        *h_min = 5000.0; *h_max = 100000.0;
    } else if (strcmp(scenario, "oil_forced") == 0) {
        *h_min = 50.0; *h_max = 2000.0;
    }
}

/* =====================================================================
 * L4: Conservation Laws
 * ===================================================================== */

int continuity_equation(double density1, double area1, double velocity1,
                        double density2, double area2, double *velocity2)
{
    if (!velocity2) return -1;
    if (area2 == 0.0) { *velocity2 = INFINITY; return -1; }
    if (area1 == 0.0) { *velocity2 = 0.0; return -1; }
    if (density2 == 0.0) { *velocity2 = INFINITY; return -1; }

    *velocity2 = (density1 * area1 * velocity1) / (density2 * area2);
    return 0;
}

double bernoulli_solve(double p1, double p2, double u1, double u2,
                       double z1, double z2, double rho, double g)
{
    int nan_count = 0;
    if (isnan(p1)) nan_count++;
    if (isnan(p2)) nan_count++;
    if (isnan(u1)) nan_count++;
    if (isnan(u2)) nan_count++;
    if (isnan(z1)) nan_count++;
    if (isnan(z2)) nan_count++;

    if (nan_count != 1) return NAN; /* must have exactly one unknown */
    if (rho <= 0.0) return NAN;

    /* Known sum: K = p_known + 0.5·ρ·u²_known + ρ·g·z_known */
    double K_known = 0.0;
    int side_known = 0; /* which side has all known values? */

    /* Try left side */
    if (!isnan(p1) && !isnan(u1) && !isnan(z1)) {
        K_known = p1 + 0.5 * rho * u1 * u1 + rho * g * z1;
        side_known = 1;
    }
    /* Try right side */
    if (!isnan(p2) && !isnan(u2) && !isnan(z2)) {
        if (side_known == 0) {
            K_known = p2 + 0.5 * rho * u2 * u2 + rho * g * z2;
            side_known = 2;
        } else {
            return NAN; /* both sides fully known → not one unknown */
        }
    }
    if (side_known == 0) return NAN;

    if (side_known == 1) {
        /* Solve for unknown on right side */
        if (isnan(p2)) {
            return K_known - 0.5 * rho * u2 * u2 - rho * g * z2;
        } else if (isnan(u2)) {
            double arg = 2.0 * (K_known - p2 - rho * g * z2) / rho;
            if (arg < 0.0) return NAN;
            return sqrt(arg);
        } else if (isnan(z2)) {
            return (K_known - p2 - 0.5 * rho * u2 * u2) / (rho * g);
        }
    } else {
        /* Solve for unknown on left side */
        if (isnan(p1)) {
            return K_known - 0.5 * rho * u1 * u1 - rho * g * z1;
        } else if (isnan(u1)) {
            double arg = 2.0 * (K_known - p1 - rho * g * z1) / rho;
            if (arg < 0.0) return NAN;
            return sqrt(arg);
        } else if (isnan(z1)) {
            return (K_known - p1 - 0.5 * rho * u1 * u1) / (rho * g);
        }
    }
    return NAN;
}

double bernoulli_head_loss(double p1, double p2, double u1, double u2,
                           double z1, double z2, double gamma, double g)
{
    if (gamma == 0.0) return INFINITY;
    double hL = (p1 - p2) / gamma
              + (u1 * u1 - u2 * u2) / (2.0 * g)
              + (z1 - z2);
    return hL;
}

/* =====================================================================
 * L4: Navier-Stokes — Simplified Analytical Solutions
 * ===================================================================== */

double plane_poiseuille_velocity(double dp_dx, double mu,
                                 double h, double y)
{
    if (mu == 0.0) return INFINITY;
    double half_h = h / 2.0;
    if (fabs(y) > half_h) return 0.0; /* outside channel */
    double dp_effective = -dp_dx; /* dp_dx is negative for +x flow */
    return dp_effective / (2.0 * mu) * (half_h * half_h - y * y);
}

double plane_poiseuille_flow_rate(double dp_dx, double mu, double h)
{
    if (mu == 0.0) return INFINITY;
    double dp_effective = -dp_dx;
    return dp_effective * h * h * h / (12.0 * mu);
}

double hagen_poiseuille_velocity(double dp_dx, double mu,
                                 double R, double r)
{
    if (mu == 0.0) return INFINITY;
    if (fabs(r) > R) return 0.0;
    double dp_effective = -dp_dx;
    return dp_effective / (4.0 * mu) * (R * R - r * r);
}

double hagen_poiseuille_flow_rate(double dp_dx, double mu, double R)
{
    if (mu == 0.0) return INFINITY;
    double dp_effective = -dp_dx;
    return M_PI * R * R * R * R * dp_effective / (8.0 * mu);
}

double couette_velocity(double U, double h, double y)
{
    if (h == 0.0) return INFINITY;
    if (y < 0.0) return 0.0;
    if (y > h) return U;
    return U * y / h;
}

double couette_poiseuille_velocity(double U, double dp_dx,
                                   double mu, double h, double y)
{
    if (mu == 0.0 || h == 0.0) return INFINITY;
    if (y < 0.0) return 0.0;
    if (y > h) return U;
    double dp_effective = -dp_dx;
    double shear_term = U * y / h;
    double pressure_term = (h * h / (2.0 * mu)) * dp_effective
                           * (y / h) * (1.0 - y / h);
    return shear_term + pressure_term;
}

double stokes_first_problem(double U, double nu, double y, double t)
{
    if (t <= 0.0) return 0.0;
    if (nu <= 0.0) return U; /* inviscid limit */
    /* erfc(z) = 1 - erf(z) */
    /* Since we don't have standard erfc, implement via erf */
    double z = y / (2.0 * sqrt(nu * t));

    /* Rational approximation for complementary error function
     * Abramowitz & Stegun 7.1.26, max error 1.5e-7 */
    double erf_z;
    double abs_z = fabs(z);
    double t_denom = 1.0 / (1.0 + 0.3275911 * abs_z);
    double poly = t_denom * (0.254829592
                   + t_denom * (-0.284496736
                   + t_denom * (1.421413741
                   + t_denom * (-1.453152027
                   + t_denom * 1.061405429))));
    erf_z = 1.0 - poly * exp(-abs_z * abs_z);
    if (z < 0) erf_z = -erf_z;

    double erfc_z = 1.0 - erf_z;
    return U * erfc_z;
}

/* =====================================================================
 * L5: Engineering Methods — Pipe Flow
 * ===================================================================== */

double darcy_friction_laminar(double reynolds_number)
{
    if (reynolds_number <= 0.0) return INFINITY;
    return 64.0 / reynolds_number;
}

double colebrook_friction_factor(double reynolds_number,
                                 double roughness_ratio,
                                 int max_iter, double tolerance)
{
    if (reynolds_number <= 0.0) return INFINITY;

    /* Initial guess: Swamee-Jain */
    double f = swamee_jain_friction_factor(reynolds_number, roughness_ratio);
    if (f <= 0.0) f = 0.02;

    /* Newton-Raphson iteration on g(f) = 1/√f + 2·log₁₀(ε/D/3.7 + 2.51/(Re·√f)) = 0
     * Let x = 1/√f, then f = 1/x²
     * g(x) = x + 2·log₁₀(ε/D/3.7 + 2.51·x/Re) = 0
     * g'(x) = 1 + (2/ln10)·(2.51/Re)/(ε/D/3.7 + 2.51·x/Re)
     */
    double x = 1.0 / sqrt(f);
    double eps_D = roughness_ratio;
    double Re = reynolds_number;
    double ln10 = log(10.0);

    for (int iter = 0; iter < max_iter; iter++) {
        double denom = eps_D / 3.7 + 2.51 * x / Re;
        if (denom <= 0.0) break;
        double g = x + 2.0 * log10(denom);
        double dg_dx = 1.0 + (2.0 / ln10) * (2.51 / Re) / denom;

        if (fabs(dg_dx) < 1e-15) break;
        double dx = -g / dg_dx;
        x += dx;

        if (fabs(dx) < tolerance) {
            double f_new = 1.0 / (x * x);
            if (f_new > 0.0) return f_new;
        }
    }
    f = 1.0 / (x * x);
    return f > 0.0 ? f : 0.02;
}

double swamee_jain_friction_factor(double reynolds_number,
                                   double roughness_ratio)
{
    if (reynolds_number <= 0.0) return INFINITY;
    double Re = reynolds_number;
    double eps_D = roughness_ratio;
    double arg = eps_D / 3.7 + 5.74 / pow(Re, 0.9);
    if (arg <= 0.0) return INFINITY;
    double denom = log10(arg);
    if (denom == 0.0) return INFINITY;
    return 0.25 / (denom * denom);
}

double darcy_weisbach_head_loss(double f, double L, double D,
                                double u, double g)
{
    if (D == 0.0 || g == 0.0) return INFINITY;
    return f * (L / D) * (u * u / (2.0 * g));
}

double darcy_weisbach_pressure_drop(double f, double L, double D,
                                    double rho, double u)
{
    if (D == 0.0) return INFINITY;
    return f * (L / D) * (0.5 * rho * u * u);
}

double minor_loss_coefficient(const char *fitting_type)
{
    if (!fitting_type) return 0.0;

    if (strcmp(fitting_type, "elbow_90_standard") == 0) return 0.75;
    if (strcmp(fitting_type, "elbow_45") == 0) return 0.35;
    if (strcmp(fitting_type, "tee_straight") == 0) return 0.4;
    if (strcmp(fitting_type, "tee_branch") == 0) return 1.0;
    if (strcmp(fitting_type, "gate_valve_open") == 0) return 0.15;
    if (strcmp(fitting_type, "globe_valve_open") == 0) return 10.0;
    if (strcmp(fitting_type, "ball_valve_open") == 0) return 0.05;
    if (strcmp(fitting_type, "sudden_expansion") == 0) return 1.0;
    if (strcmp(fitting_type, "sudden_contraction") == 0) return 0.5;
    if (strcmp(fitting_type, "pipe_entrance") == 0) return 0.5;
    if (strcmp(fitting_type, "pipe_exit") == 0) return 1.0;
    if (strcmp(fitting_type, "bend_90_r/d=1") == 0) return 0.35;
    if (strcmp(fitting_type, "bend_90_r/d=2") == 0) return 0.19;

    return 0.0; /* unknown fitting */
}

double hydraulic_diameter(double area, double perimeter)
{
    if (perimeter == 0.0) return INFINITY;
    return 4.0 * area / perimeter;
}

/* =====================================================================
 * L5: Boundary Layer Theory
 * ===================================================================== */

double blasius_boundary_layer_thickness(double x, double rex)
{
    if (rex <= 0.0) return 0.0;
    return 5.0 * x / sqrt(rex);
}

double blasius_skin_friction_local(double rex)
{
    if (rex <= 0.0) return INFINITY;
    return 0.664 / sqrt(rex);
}

double blasius_skin_friction_average(double re_l)
{
    if (re_l <= 0.0) return INFINITY;
    return 1.328 / sqrt(re_l);
}

double turbulent_boundary_layer_thickness(double x, double rex)
{
    if (rex <= 0.0) return 0.0;
    return 0.37 * x / pow(rex, 0.2);
}

double turbulent_skin_friction_local(double rex)
{
    if (rex <= 0.0) return INFINITY;
    return 0.0592 / pow(rex, 0.2);
}

double turbulent_skin_friction_average(double re_l)
{
    if (re_l <= 0.0) return INFINITY;
    return 0.074 / pow(re_l, 0.2);
}

double estimate_transition_reynolds(double turbulence_intensity)
{
    /* Simplified correlation for transition onset Reynolds number
     * based on freestream turbulence intensity Tu [%].
     *
     * Re_tr ≈ 5×10⁵ · exp(-0.5·Tu)
     * For Tu = 0.1% → Re_tr ≈ 4.75×10⁵
     * For Tu = 1.0% → Re_tr ≈ 3.03×10⁵
     * For Tu = 3.0% → Re_tr ≈ 1.12×10⁵
     */
    if (turbulence_intensity < 0.0) turbulence_intensity = 0.0;
    return 5.0e5 * exp(-0.5 * turbulence_intensity);
}

double displacement_thickness_laminar(double x, double rex)
{
    if (rex <= 0.0) return 0.0;
    return 1.7208 * x / sqrt(rex);
}

double momentum_thickness_laminar(double x, double rex)
{
    if (rex <= 0.0) return 0.0;
    return 0.664 * x / sqrt(rex);
}

double shape_factor_laminar(void)
{
    /* H = delta_star/theta = 1.7208/0.664 approx 2.59 for Blasius */
    return 2.5915662650602407;
}

/* =====================================================================
 * L5: Dimensional Analysis — Buckingham Pi Theorem
 * ===================================================================== */

int buckingham_pi_count(int n_variables, int n_dimensions)
{
    if (n_variables < n_dimensions) return 0;
    return n_variables - n_dimensions;
}

double dimensionless_group(const char *group_name, const double params[])
{
    if (!group_name || !params) return NAN;

    if (strcmp(group_name, "Re") == 0) {
        return compute_reynolds(params[0], params[1], params[2], params[3]);
    }
    if (strcmp(group_name, "Pr") == 0) {
        return compute_prandtl(params[0], params[1], params[2]);
    }
    if (strcmp(group_name, "Nu") == 0) {
        return compute_nusselt(params[0], params[1], params[2]);
    }
    if (strcmp(group_name, "Gr") == 0) {
        return compute_grashof(params[0], params[1], params[2],
                               params[3], params[4]);
    }
    if (strcmp(group_name, "Ra") == 0) {
        double gr = compute_grashof(params[0], params[1], params[2],
                                     params[3], params[4]);
        return compute_rayleigh(gr, params[5]);
    }
    if (strcmp(group_name, "Pe") == 0) {
        double re = compute_reynolds(params[0], params[1], params[2], params[3]);
        return compute_peclet(re, params[4]);
    }
    if (strcmp(group_name, "Ma") == 0) {
        return compute_mach(params[0], params[1]);
    }
    if (strcmp(group_name, "Fr") == 0) {
        return compute_froude(params[0], params[1], params[2]);
    }
    if (strcmp(group_name, "We") == 0) {
        return compute_weber(params[0], params[1], params[2], params[3]);
    }
    if (strcmp(group_name, "Eu") == 0) {
        return compute_euler(params[0], params[1], params[2]);
    }
    if (strcmp(group_name, "St") == 0) {
        return compute_strouhal(params[0], params[1], params[2]);
    }
    return NAN;
}

double compute_froude(double velocity, double g, double char_length)
{
    if (g == 0.0 || char_length == 0.0) return INFINITY;
    return velocity / sqrt(g * char_length);
}

double compute_weber(double density, double velocity,
                     double char_length, double surface_tension)
{
    if (surface_tension == 0.0) return INFINITY;
    return density * velocity * velocity * char_length / surface_tension;
}

double compute_euler(double pressure_drop, double density, double velocity)
{
    if (density == 0.0 || velocity == 0.0) return INFINITY;
    return pressure_drop / (density * velocity * velocity);
}

double compute_strouhal(double frequency, double char_length, double velocity)
{
    if (velocity == 0.0) return INFINITY;
    return frequency * char_length / velocity;
}
