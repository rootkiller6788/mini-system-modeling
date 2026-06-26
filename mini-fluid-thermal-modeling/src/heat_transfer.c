/**
 * mini-fluid-thermal-modeling — Heat Transfer Implementation
 *
 * Conduction (Fourier's law, resistance networks), convection
 * (Newton's law, Nusselt correlations), radiation (Stefan-Boltzmann,
 * view factors), fin analysis, transient conduction,
 * and heat exchanger design methods (LMTD, ε-NTU).
 *
 * Knowledge Levels: L3-L6
 */

#include "heat_transfer.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double SIGMA_SB = 5.670374419e-8; /* W/(m²·K⁴) */

/* =====================================================================
 * L3: Typical Material Properties
 * ===================================================================== */

static double get_fluid_k(const char *material);
static double get_solid_k(const char *material);

double typical_thermal_conductivity(const char *material)
{
    if (!material) return NAN;
    double k = get_solid_k(material);
    if (!isnan(k)) return k;
    return get_fluid_k(material);
}

static double get_solid_k(const char *material)
{
    if (strcmp(material, "copper") == 0) return 401.0;
    if (strcmp(material, "aluminum") == 0) return 237.0;
    if (strcmp(material, "steel_carbon") == 0) return 60.5;
    if (strcmp(material, "stainless_steel") == 0) return 14.9;
    if (strcmp(material, "glass") == 0) return 1.4;
    if (strcmp(material, "concrete") == 0) return 1.4;
    if (strcmp(material, "brick") == 0) return 0.72;
    if (strcmp(material, "wood_pine") == 0) return 0.12;
    if (strcmp(material, "fiberglass_insulation") == 0) return 0.038;
    if (strcmp(material, "silicon") == 0) return 148.0;
    if (strcmp(material, "diamond") == 0) return 2300.0;
    if (strcmp(material, "titanium") == 0) return 21.9;
    if (strcmp(material, "tungsten") == 0) return 174.0;
    if (strcmp(material, "platinum") == 0) return 71.6;
    if (strcmp(material, "gold") == 0) return 317.0;
    if (strcmp(material, "silver") == 0) return 429.0;
    if (strcmp(material, "brass") == 0) return 110.0;
    if (strcmp(material, "iron_cast") == 0) return 52.0;
    return NAN;
}

static double get_fluid_k(const char *material)
{
    if (strcmp(material, "water") == 0) return 0.613;
    if (strcmp(material, "air") == 0) return 0.0263;
    if (strcmp(material, "oil_engine") == 0) return 0.145;
    if (strcmp(material, "glycerin") == 0) return 0.286;
    if (strcmp(material, "mercury") == 0) return 8.54;
    if (strcmp(material, "ethanol") == 0) return 0.171;
    return NAN;
}

double typical_specific_heat(const char *material)
{
    if (!material) return NAN;
    if (strcmp(material, "copper") == 0) return 385.0;
    if (strcmp(material, "aluminum") == 0) return 903.0;
    if (strcmp(material, "steel_carbon") == 0) return 434.0;
    if (strcmp(material, "stainless_steel") == 0) return 500.0;
    if (strcmp(material, "glass") == 0) return 750.0;
    if (strcmp(material, "concrete") == 0) return 880.0;
    if (strcmp(material, "brick") == 0) return 835.0;
    if (strcmp(material, "wood_pine") == 0) return 2800.0;
    if (strcmp(material, "silicon") == 0) return 712.0;
    if (strcmp(material, "water") == 0) return 4179.0;
    if (strcmp(material, "air") == 0) return 1007.0;
    if (strcmp(material, "oil_engine") == 0) return 1900.0;
    return NAN;
}

double typical_density(const char *material)
{
    if (!material) return NAN;
    if (strcmp(material, "copper") == 0) return 8933.0;
    if (strcmp(material, "aluminum") == 0) return 2702.0;
    if (strcmp(material, "steel_carbon") == 0) return 7854.0;
    if (strcmp(material, "stainless_steel") == 0) return 7900.0;
    if (strcmp(material, "glass") == 0) return 2500.0;
    if (strcmp(material, "concrete") == 0) return 2300.0;
    if (strcmp(material, "brick") == 0) return 1920.0;
    if (strcmp(material, "wood_pine") == 0) return 510.0;
    if (strcmp(material, "fiberglass_insulation") == 0) return 24.0;
    if (strcmp(material, "silicon") == 0) return 2330.0;
    if (strcmp(material, "water") == 0) return 997.0;
    if (strcmp(material, "air") == 0) return 1.1614;
    if (strcmp(material, "oil_engine") == 0) return 884.0;
    return NAN;
}

/* =====================================================================
 * L4: Conduction — Fourier's Law
 * ===================================================================== */

double fourier_conduction(double k, double A, double T1,
                          double T2, double L)
{
    if (L == 0.0) return INFINITY;
    return k * A * (T1 - T2) / L;
}

double thermal_resistance_plane_wall(double L, double k, double A)
{
    if (k == 0.0 || A == 0.0) return INFINITY;
    return L / (k * A);
}

double thermal_resistance_cylinder(double r1, double r2,
                                   double k, double L)
{
    if (k == 0.0 || L == 0.0 || r1 <= 0.0 || r2 <= r1) return INFINITY;
    return log(r2 / r1) / (2.0 * M_PI * k * L);
}

double thermal_resistance_sphere(double r1, double r2, double k)
{
    if (k == 0.0 || r1 <= 0.0 || r2 <= r1) return INFINITY;
    return (1.0 / r1 - 1.0 / r2) / (4.0 * M_PI * k);
}

double thermal_resistance_convection(double h, double A)
{
    if (h == 0.0 || A == 0.0) return INFINITY;
    return 1.0 / (h * A);
}

double contact_resistance(double specific_resistance, double area)
{
    if (area == 0.0) return INFINITY;
    return specific_resistance / area;
}

/* =====================================================================
 * L4: Convection — Newton's Law of Cooling
 * ===================================================================== */

double newton_cooling(double h, double A, double Ts, double Tinf)
{
    return h * A * (Ts - Tinf);
}

/* =====================================================================
 * L5: Nusselt Number Correlations
 * ===================================================================== */

double dittus_boelter_nusselt(double re, double pr, int heating)
{
    if (re < 10000.0 || pr < 0.6 || pr > 160.0) return NAN;
    double n = heating ? 0.4 : 0.3;
    return 0.023 * pow(re, 0.8) * pow(pr, n);
}

double sieder_tate_nusselt(double re, double pr,
                           double mu_bulk, double mu_wall)
{
    if (re < 10000.0 || mu_wall == 0.0) return NAN;
    return 0.027 * pow(re, 0.8) * pow(pr, 1.0/3.0)
           * pow(mu_bulk / mu_wall, 0.14);
}

double churchill_chu_vertical_plate(double ra, double pr)
{
    if (ra < 0.0) return NAN;
    double num = 0.387 * pow(ra, 1.0/6.0);
    double denom = pow(1.0 + pow(0.492 / pr, 9.0/16.0), 8.0/27.0);
    if (denom == 0.0) return INFINITY;
    double nu_sqrt = 0.825 + num / denom;
    return nu_sqrt * nu_sqrt;
}

double mcadams_horizontal_plate(double ra, int heated_upward)
{
    if (ra < 0.0) return NAN;

    if (heated_upward) {
        if (ra <= 1.0e7) {
            return 0.54 * pow(ra, 0.25);
        } else {
            return 0.15 * pow(ra, 1.0/3.0);
        }
    } else {
        /* Heated downward: only laminar correlation */
        return 0.27 * pow(ra, 0.25);
    }
}

double morgan_horizontal_cylinder(double ra)
{
    /* Morgan (1975) correlation for natural convection
     * from a horizontal cylinder.
     * Nu = C · Ra^n with piecewise C, n.
     */
    if (ra <= 0.0) return NAN;

    if (ra < 1.0e-5) {
        return 0.675 * pow(ra, 0.058);
    } else if (ra < 1.0e-3) {
        return 0.889 * pow(ra, 0.088);
    } else if (ra < 1.0) {
        return 0.85 * pow(ra, 0.148);
    } else if (ra < 1.0e2) {
        return 0.48 * pow(ra, 0.188);
    } else if (ra < 1.0e4) {
        return 0.85 * pow(ra, 0.188);
    } else if (ra < 1.0e7) {
        return 0.48 * pow(ra, 0.25);
    } else {
        return 0.125 * pow(ra, 1.0/3.0);
    }
}

double flat_plate_laminar_nusselt_local(double re_x, double pr)
{
    if (re_x <= 0.0 || pr < 0.6) return NAN;
    return 0.332 * sqrt(re_x) * pow(pr, 1.0/3.0);
}

double flat_plate_laminar_nusselt_average(double re_l, double pr)
{
    if (re_l <= 0.0 || pr < 0.6) return NAN;
    return 0.664 * sqrt(re_l) * pow(pr, 1.0/3.0);
}

double flat_plate_turbulent_nusselt_average(double re_l, double pr)
{
    if (re_l < 5.0e5 || re_l > 1.0e7
        || pr < 0.6 || pr > 60.0) return NAN;
    return (0.037 * pow(re_l, 0.8) - 871.0) * pow(pr, 1.0/3.0);
}

double zhukauskas_tube_bank_nusselt(double re_max, double pr,
                                    double pr_s, int n_rows)
{
    if (re_max <= 0.0 || pr <= 0.0) return NAN;

    /* Zhukauskas correlation constants C, m depend on Re_max and
     * tube arrangement. Using inline arrangement defaults. */
    double C, m;
    if (re_max < 100.0) {
        C = 0.80; m = 0.40;
    } else if (re_max < 1000.0) {
        C = 0.51; m = 0.50;
    } else if (re_max < 2.0e5) {
        C = 0.27; m = 0.63;
    } else {
        C = 0.021; m = 0.84;
    }

    /* Row count correction for n_rows < 20 */
    double row_factor;
    if (n_rows >= 20) {
        row_factor = 1.0;
    } else if (n_rows >= 10) {
        row_factor = 0.98;
    } else if (n_rows >= 5) {
        row_factor = 0.94;
    } else if (n_rows >= 3) {
        row_factor = 0.89;
    } else {
        row_factor = 0.81;
    }

    double pr_ratio = pr_s > 0.0 ? pr / pr_s : 1.0;
    return row_factor * C * pow(re_max, m)
           * pow(pr, 0.36) * pow(pr_ratio, 0.25);
}

/* =====================================================================
 * L4: Radiation — Stefan-Boltzmann Law
 * ===================================================================== */

double stefan_boltzmann_emission(double T)
{
    return SIGMA_SB * T * T * T * T;
}

double blackbody_net_radiation(double sigma, double A1, double F12,
                               double T1, double T2)
{
    return sigma * A1 * F12 * (T1*T1*T1*T1 - T2*T2*T2*T2);
}

double gray_net_radiation(double T1, double T2, double A1, double A2,
                          double eps1, double eps2, double F12)
{
    /* Surface resistance method:
     * q = σ·(T₁⁴ - T₂⁴) / ((1-ε₁)/(ε₁·A₁) + 1/(A₁·F₁₂) + (1-ε₂)/(ε₂·A₂))
     */
    if (eps1 <= 0.0 || eps2 <= 0.0 || (A1 * F12) == 0.0) return NAN;
    double R1 = (1.0 - eps1) / (eps1 * A1);
    double R12 = 1.0 / (A1 * F12);
    double R2 = (1.0 - eps2) / (eps2 * A2);
    double R_total = R1 + R12 + R2;
    if (R_total == 0.0) return INFINITY;
    return SIGMA_SB * (T1*T1*T1*T1 - T2*T2*T2*T2) / R_total;
}

double view_factor_parallel_disks(double r1, double r2, double L)
{
    if (L <= 0.0 || r1 <= 0.0 || r2 <= 0.0) return NAN;
    double R1 = r1 / L;
    double R2 = r2 / L;
    double X = 1.0 + (1.0 + R2 * R2) / (R1 * R1);
    double disc = X * X - 4.0 * (r2 * r2) / (r1 * r1);
    if (disc < 0.0) return 0.0;
    return 0.5 * (X - sqrt(disc));
}

double view_factor_perpendicular_rectangles(double w, double h, double d)
{
    /* View factor between two perpendicular rectangles
     * with a common edge of length h.
     * Rectangle 1: width w, height h
     * Rectangle 2: width d, height h
     * F₁₂ from Hottel's string method.
     *
     * Using Hamilton-Morgan formula (NACA TN 2836) simplified:
     * F₁₂ = 1/(π·W) [ W·tan⁻¹(1/W) + H·tan⁻¹(1/H) - √(H²+W²)·tan⁻¹(1/√(H²+W²))
     *   + 0.25·ln( ((1+W²)(1+H²))/(1+W²+H²) · (W²(1+W²+H²)/((1+W²)(W²+H²)))^{W²}
     *               · (H²(1+W²+H²)/((1+H²)(W²+H²)))^{H²} ) ]
     * where W = w/d, H = h/d
     */
    if (d <= 0.0) return NAN;
    double W = w / d;
    double H_rect = h / d;
    double W2 = W * W;
    double H2 = H_rect * H_rect;
    double WH2 = W2 + H2;

    if (W == 0.0 || H_rect == 0.0) return 0.0;

    double term1 = W * atan(1.0 / W);
    double term2 = H_rect * atan(1.0 / H_rect);
    double term3 = sqrt(1.0 + WH2) * atan(1.0 / sqrt(1.0 + WH2));

    double arg1 = ((1.0 + W2) * (1.0 + H2)) / (1.0 + WH2);
    double arg2_num = W2 * (1.0 + WH2);
    double arg2_den = (1.0 + W2) * (W2 + H2);
    double arg3_num = H2 * (1.0 + WH2);
    double arg3_den = (1.0 + H2) * (W2 + H2);

    if (arg1 <= 0.0 || arg2_num <= 0.0 || arg2_den <= 0.0
        || arg3_num <= 0.0 || arg3_den <= 0.0) return NAN;

    double log_term = 0.25 * log(arg1)
                    + 0.25 * W2 * log(arg2_num / arg2_den)
                    + 0.25 * H2 * log(arg3_num / arg3_den);

    return (term1 + term2 - term3 + log_term) / (M_PI * W);
}

double view_factor_parallel_rectangles(double a, double b, double c)
{
    /* View factor between two aligned parallel rectangles
     * of dimensions a×b separated by distance c.
     *
     * Using Howell's catalog formula.
     *
     * X = a/c, Y = b/c
     * F₁₂ = 2/(π·X·Y) [ ln(√((1+X²)(1+Y²)/(1+X²+Y²)))
     *   + X·√(1+Y²)·tan⁻¹(X/√(1+Y²))
     *   + Y·√(1+X²)·tan⁻¹(Y/√(1+X²))
     *   - X·tan⁻¹(X) - Y·tan⁻¹(Y) ]
     */
    if (c <= 0.0 || a <= 0.0 || b <= 0.0) return NAN;
    double X = a / c;
    double Y = b / c;
    double X2 = X * X;
    double Y2 = Y * Y;
    double one_plus_X2 = 1.0 + X2;
    double one_plus_Y2 = 1.0 + Y2;
    double one_plus_X2_plus_Y2 = 1.0 + X2 + Y2;

    double arg = sqrt(one_plus_X2 * one_plus_Y2 / one_plus_X2_plus_Y2);
    if (arg <= 0.0) return NAN;

    double term1 = log(arg);
    double term2 = X * sqrt(one_plus_Y2) * atan(X / sqrt(one_plus_Y2));
    double term3 = Y * sqrt(one_plus_X2) * atan(Y / sqrt(one_plus_X2));
    double term4 = -X * atan(X) - Y * atan(Y);

    return 2.0 / (M_PI * X * Y) * (term1 + term2 + term3 + term4);
}

/* =====================================================================
 * L5: Fin Analysis
 * ===================================================================== */

double fin_temperature_distribution(double m, double L, double x)
{
    if (m == 0.0) return 1.0; /* isothermal fin */
    if (x < 0.0) x = 0.0;
    if (x > L) x = L;
    return cosh(m * (L - x)) / cosh(m * L);
}

double fin_heat_transfer_rate(double h, double P, double k,
                              double Ac, double theta_b, double L)
{
    if (k == 0.0 || Ac == 0.0) return 0.0;
    double m = sqrt(h * P / (k * Ac));
    if (m == 0.0) return h * P * L * theta_b; /* uniform temperature */
    return sqrt(h * P * k * Ac) * theta_b * tanh(m * L);
}

double fin_efficiency(double m, double L)
{
    if (m == 0.0 || L == 0.0) return 1.0;
    return tanh(m * L) / (m * L);
}

double fin_effectiveness(double q_f, double h, double Ac, double theta_b)
{
    double q_no_fin = h * Ac * theta_b;
    if (q_no_fin == 0.0) return INFINITY;
    return q_f / q_no_fin;
}

/* =====================================================================
 * L5: Transient Conduction
 * ===================================================================== */

double lumped_capacitance_temp(double T_inf, double T_i, double h,
                               double A, double rho, double V, double cp,
                               double t)
{
    if (V == 0.0 || cp == 0.0 || rho == 0.0) return T_i;
    double tau = rho * V * cp / (h * A);
    if (tau <= 0.0) return T_inf; /* zero resistance */
    return T_inf + (T_i - T_inf) * exp(-t / tau);
}

double lumped_capacitance_time_constant(double rho, double V,
                                        double cp, double h, double A)
{
    if (h == 0.0 || A == 0.0) return INFINITY;
    return rho * V * cp / (h * A);
}

double one_term_approximation_slab(double Bi, double Fo, double x_L)
{
    /* One-term approximation for infinite slab (Fo > 0.2).
     *
     * ζ₁·tan(ζ₁) = Bi  (eigenvalue equation)
     * C₁ = 4·sin(ζ₁) / (2·ζ₁ + sin(2·ζ₁))
     *
     * Tabulated values for common Bi:
     */
    double zeta1, C1;

    if (Bi <= 0.0) {
        /* Perfectly insulated: uniform temperature = initial */
        zeta1 = 0.0;
        C1 = 1.0;
    } else if (Bi < 0.01) {
        /* Lumped approximation: ζ₁ ≈ √(Bi), C₁ ≈ 1 */
        zeta1 = sqrt(Bi);
        C1 = 1.0;
    } else if (Bi < 0.1) {
        zeta1 = 0.3111;
        C1 = 1.0161;
    } else if (Bi < 0.5) {
        zeta1 = 0.6533;
        C1 = 1.0701;
    } else if (Bi < 1.0) {
        zeta1 = 0.8603;
        C1 = 1.1191;
    } else if (Bi < 5.0) {
        zeta1 = 1.3138;
        C1 = 1.2402;
    } else if (Bi < 10.0) {
        zeta1 = 1.4289;
        C1 = 1.2620;
    } else if (Bi < 50.0) {
        zeta1 = 1.5400;
        C1 = 1.2698;
    } else {
        /* Bi → ∞, constant surface temperature */
        zeta1 = M_PI / 2.0;
        C1 = 4.0 / M_PI;
    }

    if (x_L < 0.0) x_L = 0.0;
    if (x_L > 1.0) x_L = 1.0;
    double arg1 = -zeta1 * zeta1 * Fo;
    double arg2 = zeta1 * x_L;
    if (arg1 < -50.0) return 0.0; /* avoid underflow */
    return C1 * exp(arg1) * cos(arg2);
}

/* =====================================================================
 * L5: Combined Heat Transfer
 * ===================================================================== */

double overall_heat_transfer_coefficient_plane(double h1, double L,
                                                double k, double h2)
{
    double R_total = 0.0;
    if (h1 > 0.0) R_total += 1.0 / h1;
    if (k > 0.0 && L > 0.0) R_total += L / k;
    if (h2 > 0.0) R_total += 1.0 / h2;
    if (R_total == 0.0) return INFINITY;
    return 1.0 / R_total;
}

double overall_heat_transfer_coefficient_cylinder(
    double hi, double ho, double ri, double ro, double k)
{
    if (ri <= 0.0 || ro <= ri || k <= 0.0) return NAN;
    double R_total = 0.0;
    if (hi > 0.0) R_total += 1.0 / hi;
    if (k > 0.0) R_total += ri * log(ro / ri) / k;
    if (ho > 0.0) R_total += ri / (ro * ho);
    if (R_total == 0.0) return INFINITY;
    return 1.0 / R_total;
}

double critical_insulation_radius(double k, double h_outer)
{
    if (h_outer == 0.0) return INFINITY;
    return k / h_outer;
}

double critical_insulation_radius_sphere(double k, double h_outer)
{
    if (h_outer == 0.0) return INFINITY;
    return 2.0 * k / h_outer;
}

/* =====================================================================
 * L6: Heat Exchanger Design Methods
 * ===================================================================== */

double lmtd(double delta_T1, double delta_T2)
{
    if (delta_T1 <= 0.0 || delta_T2 <= 0.0) return NAN;
    if (delta_T1 == delta_T2) return delta_T1; /* limit case */
    return (delta_T1 - delta_T2) / log(delta_T1 / delta_T2);
}

double heat_exchanger_duty_lmtd(double U, double A, double F,
                                double dT_lm)
{
    return U * A * F * dT_lm;
}

double heat_exchanger_effectiveness(double ntu, double Cr,
                                    const char *arrangement)
{
    if (!arrangement || ntu < 0.0 || Cr < 0.0 || Cr > 1.0) return NAN;

    if (strcmp(arrangement, "parallel") == 0) {
        double denom = 1.0 + Cr;
        if (denom == 0.0) return 1.0 - exp(-ntu);
        return (1.0 - exp(-ntu * (1.0 + Cr))) / denom;
    }
    if (strcmp(arrangement, "counter") == 0) {
        if (Cr == 0.0) return 1.0 - exp(-ntu);
        if (Cr == 1.0) return ntu / (1.0 + ntu);
        double exp_neg = exp(-ntu * (1.0 - Cr));
        return (1.0 - exp_neg) / (1.0 - Cr * exp_neg);
    }
    if (strcmp(arrangement, "shell_tube_1") == 0) {
        /* One shell pass, 2·n tube passes */
        if (Cr == 0.0) return 1.0 - exp(-ntu);
        double e1 = exp(-ntu * sqrt(1.0 + Cr * Cr));
        double num = 2.0;
        double denom = 1.0 + Cr + sqrt(1.0 + Cr * Cr) * (1.0 + e1) / (1.0 - e1);
        return num / denom;
    }
    if (strcmp(arrangement, "cross_both_unmixed") == 0) {
        /* Cross-flow, both fluids unmixed
         * ε = 1 - exp((1/Cr)·NTU^0.22·(exp(-Cr·NTU^0.78)-1))
         */
        if (Cr == 0.0) return 1.0 - exp(-ntu);
        double inner = exp(-Cr * pow(ntu, 0.78)) - 1.0;
        return 1.0 - exp((1.0 / Cr) * pow(ntu, 0.22) * inner);
    }
    return NAN;
}

double ntu_from_effectiveness_parallel(double epsilon, double Cr)
{
    if (epsilon >= 1.0 || epsilon < 0.0) return NAN;
    double denom = 1.0 + Cr;
    if (denom == 0.0) return -log(1.0 - epsilon);
    return -log(1.0 - epsilon * denom) / denom;
}

double ntu_from_effectiveness_counter(double epsilon, double Cr)
{
    if (epsilon >= 1.0 || epsilon < 0.0) return NAN;
    if (Cr == 0.0) return -log(1.0 - epsilon);
    if (Cr == 1.0) return epsilon / (1.0 - epsilon);
    return log((epsilon - 1.0) / (epsilon * Cr - 1.0)) / (Cr - 1.0);
}

double heat_capacity_rate(double mass_flow_rate, double specific_heat)
{
    return mass_flow_rate * specific_heat;
}

void heat_exchanger_outlet_temps(double Th_in, double Tc_in,
                                 double Ch, double Cc,
                                 double epsilon,
                                 double *Th_out, double *Tc_out)
{
    if (!Th_out || !Tc_out) return;
    double C_min = (Ch < Cc) ? Ch : Cc;
    if (C_min == 0.0) {
        *Th_out = Th_in; *Tc_out = Tc_in; return;
    }
    double q = epsilon * C_min * (Th_in - Tc_in);
    if (Ch > 0.0) *Th_out = Th_in - q / Ch;
    else *Th_out = Th_in;
    if (Cc > 0.0) *Tc_out = Tc_in + q / Cc;
    else *Tc_out = Tc_in;
}

double hx_correction_factor_shell_tube(double P, double R)
{
    /* Correction factor F for 1 shell, 2 tube pass configuration.
     *
     * P = (Tc_out - Tc_in) / (Th_in - Tc_in) — temperature effectiveness
     * R = (Th_in - Th_out) / (Tc_out - Tc_in) — capacity ratio
     *
     * Source: Bowman, Mueller, Nagle (1940)
     */
    if (P <= 0.0 || P >= 1.0) return NAN;

    if (R == 1.0) {
        /* Limit case R=1 */
        double S = P / (1.0 - P);
        return sqrt(2.0) * S / log((1.0 + S * sqrt(2.0)) /
                                   (1.0 - S * sqrt(2.0)));
    }

    /* Standard form: arg = (1-P)/(1-PR) — see Incropera Eq 11.32 */
    double arg1_num = 1.0 - P;
    double arg1_den = 1.0 - P * R;
    double arg1 = arg1_num / arg1_den;
    if (arg1 <= 0.0) return NAN;

    double inner = (2.0 - P * (R + 1.0 - sqrt(R*R + 1.0)))
                 / (2.0 - P * (R + 1.0 + sqrt(R*R + 1.0)));
    if (inner <= 0.0) return NAN;

    double numerator = sqrt(R*R + 1.0) * log(arg1);
    double denominator = (R - 1.0) * log(inner);

    if (denominator == 0.0) return 1.0;
    return numerator / denominator;
}
