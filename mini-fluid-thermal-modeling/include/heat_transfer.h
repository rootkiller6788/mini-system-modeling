/**
 * mini-fluid-thermal-modeling — Heat Transfer
 *
 * Conduction, convection, radiation, and combined-mode analysis.
 * Steady and transient problems, fin analysis, and thermal resistance networks.
 *
 * Reference:
 *   Incropera & DeWitt, "Fundamentals of Heat and Mass Transfer" (2007)
 *   Holman, "Heat Transfer" (2010)
 *   Cengel, "Heat and Mass Transfer: Fundamentals and Applications" (2014)
 *
 * Knowledge Levels:
 *   L3: Engineering quantities (h ranges, Nu correlations, k values)
 *   L4: Conservation laws (Fourier, Newton cooling, Stefan-Boltzmann)
 *   L5: Engineering methods (LMTD, ε-NTU, fin efficiency, resistance networks)
 */

#ifndef HEAT_TRANSFER_H
#define HEAT_TRANSFER_H

#include "fluid_thermal_core.h"

/* ---------------------------------------------------------------------------
 * L3: Typical Material Properties
 * ------------------------------------------------------------------------- */

/** Get thermal conductivity [W/(m·K)] for common materials at ~300K.
 *
 *  Supported: "copper", "aluminum", "steel_carbon", "stainless_steel",
 *             "glass", "water", "air", "oil_engine", "concrete",
 *             "brick", "wood_pine", "fiberglass_insulation", "silicon"
 */
double typical_thermal_conductivity(const char *material);

/** Get specific heat cp [J/(kg·K)] for common materials at ~300K. */
double typical_specific_heat(const char *material);

/** Get density [kg/m³] for common materials at ~300K. */
double typical_density(const char *material);

/* ---------------------------------------------------------------------------
 * L4: Conduction — Fourier's Law
 * ------------------------------------------------------------------------- */

/**
 * Fourier's law of heat conduction (1D steady).
 * q = -k · dT/dx
 *
 * Discrete form: q = k · A · (T₁ - T₂) / L
 *
 * @param k   Thermal conductivity [W/(m·K)]
 * @param A   Cross-sectional area [m²]
 * @param T1  Temperature at x=0 [K]
 * @param T2  Temperature at x=L [K]
 * @param L   Distance [m]
 * @return Heat transfer rate q [W]
 *
 * Source: Fourier (1822), "Théorie analytique de la chaleur", Article 67
 */
double fourier_conduction(double k, double A, double T1, double T2, double L);

/**
 * Thermal resistance for plane wall conduction.
 * R_cond = L / (k·A)
 *
 * @return Thermal resistance [K/W]
 */
double thermal_resistance_plane_wall(double L, double k, double A);

/**
 * Thermal resistance for cylindrical conduction (radial).
 * R_cond = ln(r₂/r₁) / (2π·k·L)
 *
 * @param r1  Inner radius [m]
 * @param r2  Outer radius [m]
 * @param k   Thermal conductivity [W/(m·K)]
 * @param L   Cylinder length [m]
 * @return Thermal resistance [K/W]
 */
double thermal_resistance_cylinder(double r1, double r2, double k, double L);

/**
 * Thermal resistance for spherical conduction (radial).
 * R_cond = (1/r₁ - 1/r₂) / (4π·k)
 */
double thermal_resistance_sphere(double r1, double r2, double k);

/**
 * Thermal resistance for convection.
 * R_conv = 1 / (h·A)
 */
double thermal_resistance_convection(double h, double A);

/**
 * Interface contact resistance model.
 * R_t,c = R"_t,c / A
 *
 * Typical R"_t,c values [m²·K/W]:
 *   Silicon-to-silicon (grease): 0.000025
 *   Aluminum-to-aluminum (air): 0.000275
 *   Stainless-to-stainless (air): 0.0009
 */
double contact_resistance(double specific_resistance, double area);

/* ---------------------------------------------------------------------------
 * L4: Convection — Newton's Law of Cooling
 * ------------------------------------------------------------------------- */

/**
 * Newton's law of cooling.
 * q = h · A · (Ts - T∞)
 *
 * @param h    Heat transfer coefficient [W/(m²·K)]
 * @param A    Surface area [m²]
 * @param Ts   Surface temperature [K]
 * @param Tinf Fluid temperature [K]
 * @return Heat transfer rate [W] (positive = from surface to fluid)
 *
 * Source: Newton (1701), "Scala graduum Caloris"
 */
double newton_cooling(double h, double A, double Ts, double Tinf);

/* ---------------------------------------------------------------------------
 * L5: Nusselt Number Correlations
 * ------------------------------------------------------------------------- */

/**
 * Dittus-Boelter correlation for turbulent pipe flow (heating).
 * Nu = 0.023 · Re^0.8 · Pr^n
 * n = 0.4 for heating (Ts > Tb), n = 0.3 for cooling (Ts < Tb)
 *
 * Valid: 0.6 ≤ Pr ≤ 160, Re ≥ 10000, L/D ≥ 10
 *
 * Source: Dittus & Boelter (1930), UC Berkeley Pub. Eng. 2:443
 */
double dittus_boelter_nusselt(double re, double pr, int heating);

/**
 * Sieder-Tate correlation.
 * Nu = 0.027 · Re^0.8 · Pr^{1/3} · (μ/μs)^0.14
 *
 * Accounts for property variation; valid 0.7 ≤ Pr ≤ 16700, Re ≥ 10000
 *
 * Source: Sieder & Tate (1936), Ind. Eng. Chem. 28:1429
 */
double sieder_tate_nusselt(double re, double pr,
                           double mu_bulk, double mu_wall);

/**
 * Churchill-Chu correlation for natural convection (vertical plate).
 * Nu = {0.825 + 0.387·Ra^{1/6} / [1 + (0.492/Pr)^{9/16}]^{8/27}}²
 *
 * Valid for all RaL (laminar + turbulent)
 *
 * Source: Churchill & Chu (1975), Int. J. Heat Mass Transfer 18:1323
 */
double churchill_chu_vertical_plate(double ra, double pr);

/**
 * McAdams correlation for natural convection (horizontal plate, heated up).
 * Nu = 0.54 · Ra^{1/4}  (10⁴ ≤ Ra ≤ 10⁷, laminar)
 * Nu = 0.15 · Ra^{1/3}  (10⁷ ≤ Ra ≤ 10¹⁰, turbulent)
 */
double mcadams_horizontal_plate(double ra, int heated_upward);

/**
 * Morgan correlation for natural convection (horizontal cylinder).
 * Nu = C · Ra^n
 * C and n from tabulated values.
 */
double morgan_horizontal_cylinder(double ra);

/**
 * Forced convection over a flat plate (laminar, constant Ts).
 * Nu_x = 0.332 · Re_x^{1/2} · Pr^{1/3}
 * Valid: Pr ≥ 0.6
 *
 * Source: Pohlhausen (1921)
 */
double flat_plate_laminar_nusselt_local(double re_x, double pr);

/**
 * Forced convection over a flat plate (laminar, average).
 * Nu_L = 0.664 · Re_L^{1/2} · Pr^{1/3}
 */
double flat_plate_laminar_nusselt_average(double re_l, double pr);

/**
 * Forced convection over a flat plate (turbulent, average).
 * Nu_L = (0.037 · Re_L^{4/5} - 871) · Pr^{1/3}
 *
 * Valid: 0.6 ≤ Pr ≤ 60, 5×10⁵ ≤ Re_L ≤ 10⁷
 */
double flat_plate_turbulent_nusselt_average(double re_l, double pr);

/**
 * Zhukauskas correlation for flow over tube banks.
 * Nu = C · Re_max^m · Pr^0.36 · (Pr/Pr_s)^{1/4}
 *
 * Source: Zhukauskas (1972), Adv. Heat Transfer 8:93
 */
double zhukauskas_tube_bank_nusselt(double re_max, double pr,
                                    double pr_s, int n_rows);

/* ---------------------------------------------------------------------------
 * L4: Radiation — Stefan-Boltzmann Law
 * ------------------------------------------------------------------------- */

/**
 * Stefan-Boltzmann law for blackbody radiation.
 * E_b = σ · T⁴
 * σ = 5.670374419 × 10⁻⁸ W/(m²·K⁴)
 *
 * @param T  Absolute temperature [K]
 * @return Emissive power [W/m²]
 */
double stefan_boltzmann_emission(double T);

/**
 * Net radiation between two black surfaces A₁ and A₂.
 * q₁₂ = σ · A₁ · F₁₂ · (T₁⁴ - T₂⁴)
 *
 * @param sigma  Stefan-Boltzmann constant
 * @param A1     Area of surface 1 [m²]
 * @param F12    View factor from 1 to 2
 * @param T1,T2  Temperatures [K]
 * @return Net radiative heat transfer [W]
 */
double blackbody_net_radiation(double sigma, double A1, double F12,
                               double T1, double T2);

/**
 * Net radiation between two gray diffuse surfaces.
 * q₁₂ = σ · (T₁⁴ - T₂⁴) / R_rad
 * R_rad = (1-ε₁)/(ε₁·A₁) + 1/(A₁·F₁₂) + (1-ε₂)/(ε₂·A₂)
 *
 * @return Net radiative heat transfer [W]
 */
double gray_net_radiation(double T1, double T2, double A1, double A2,
                          double eps1, double eps2, double F12);

/**
 * View factor between two parallel coaxial disks of radii r₁, r₂.
 *
 * F₁₂ = ½ {X - [X² - 4(r₂/r₁)²]^{1/2}}
 * where X = 1 + (1+R₂²)/R₁², R₁=r₁/L, R₂=r₂/L
 *
 * Source: Howell (2010), "A Catalog of Radiation Heat Transfer Configurations"
 */
double view_factor_parallel_disks(double r1, double r2, double L);

/**
 * View factor between two perpendicular rectangles with common edge.
 * Uses contour integration solution.
 *
 * Source: Hamilton & Morgan (1952), NACA TN 2836
 */
double view_factor_perpendicular_rectangles(double w, double h, double d);

/**
 * View factor between two parallel identical rectangles.
 *
 * Source: Howell (2010)
 */
double view_factor_parallel_rectangles(double a, double b, double c);

/* ---------------------------------------------------------------------------
 * L5: Fin Analysis
 * ------------------------------------------------------------------------- */

/**
 * Temperature distribution in a fin of uniform cross-section.
 *
 * θ(x)/θ_b = cosh[m·(L-x)] / cosh(m·L)  (adiabatic tip)
 *
 * where m² = h·P/(k·Ac), θ = T - T∞, θ_b = Tb - T∞
 *
 * @param m    Fin parameter √(hP/kAc) [1/m]
 * @param L    Fin length [m]
 * @param x    Position from base [m]
 * @return θ(x)/θ_b [dimensionless]
 */
double fin_temperature_distribution(double m, double L, double x);

/**
 * Fin heat transfer rate (adiabatic tip).
 * q_f = √(h·P·k·Ac) · θ_b · tanh(m·L)
 */
double fin_heat_transfer_rate(double h, double P, double k, double Ac,
                              double theta_b, double L);

/**
 * Fin efficiency (adiabatic tip).
 * η_f = tanh(m·L) / (m·L)
 */
double fin_efficiency(double m, double L);

/**
 * Fin effectiveness.
 * ε_f = q_f / (h·Ac·θ_b)
 */
double fin_effectiveness(double q_f, double h, double Ac, double theta_b);

/* ---------------------------------------------------------------------------
 * L5: Transient Conduction
 * ------------------------------------------------------------------------- */

/**
 * Lumped capacitance method: temperature vs time.
 *
 * T(t) = T∞ + (Ti - T∞) · exp(-Bi·Fo)
 *      = T∞ + (Ti - T∞) · exp(-h·A·t/(ρ·V·cp))
 *
 * Valid when Biot number Bi < 0.1.
 *
 * @param T_inf   Ambient temperature [K]
 * @param T_i     Initial temperature [K]
 * @param h       Convection coefficient [W/(m²·K)]
 * @param A       Surface area [m²]
 * @param rho     Density [kg/m³]
 * @param V       Volume [m³]
 * @param cp      Specific heat [J/(kg·K)]
 * @param t       Time [s]
 * @return Temperature T(t) [K]
 */
double lumped_capacitance_temp(double T_inf, double T_i, double h, double A,
                               double rho, double V, double cp, double t);

/**
 * Lumped capacitance time constant.
 * τ = ρ·V·cp / (h·A)
 */
double lumped_capacitance_time_constant(double rho, double V,
                                        double cp, double h, double A);

/**
 * One-term approximation for infinite slab (Fo > 0.2).
 *
 * θ*(x,t) = C₁ · exp(-ζ₁²·Fo) · cos(ζ₁·x/L)
 *
 * C₁, ζ₁ from tabulated coefficients for given Biot number.
 *
 * @param Bi  Biot number hL/k
 * @param Fo  Fourier number αt/L²
 * @param x_L Normalized position x/L ∈ [0,1]
 * @return Dimensionless temperature θ*
 */
double one_term_approximation_slab(double Bi, double Fo, double x_L);

/* ---------------------------------------------------------------------------
 * L5: Combined Heat Transfer
 * ------------------------------------------------------------------------- */

/**
 * Overall heat transfer coefficient (plane wall).
 * U = 1 / (1/h₁ + L/k + 1/h₂)
 *
 * @return U [W/(m²·K)]
 */
double overall_heat_transfer_coefficient_plane(double h1, double L,
                                                double k, double h2);

/**
 * Overall heat transfer coefficient (cylindrical wall).
 * U_i = 1 / (1/h_i + r_i·ln(r_o/r_i)/k + r_i/(r_o·h_o))
 * Based on inner area.
 *
 * @return U_i [W/(m²·K)]
 */
double overall_heat_transfer_coefficient_cylinder(
    double hi, double ho, double ri, double ro, double k);

/**
 * Critical radius of insulation (cylinder).
 * r_cr = k / h_o
 *
 * Adding insulation increases heat loss if r_o < r_cr.
 */
double critical_insulation_radius(double k, double h_outer);

/**
 * Critical radius of insulation (sphere).
 * r_cr = 2·k / h_o
 */
double critical_insulation_radius_sphere(double k, double h_outer);

/* ---------------------------------------------------------------------------
 * L6: Heat Exchanger Design Methods
 * ------------------------------------------------------------------------- */

/**
 * Log Mean Temperature Difference (LMTD).
 *
 * ΔT_lm = (ΔT₁ - ΔT₂) / ln(ΔT₁/ΔT₂)
 *
 * For parallel flow:
 *   ΔT₁ = Th,i - Tc,i
 *   ΔT₂ = Th,o - Tc,o
 * For counter flow:
 *   ΔT₁ = Th,i - Tc,o
 *   ΔT₂ = Th,o - Tc,i
 *
 * @return LMTD [K]. Returns NaN if ΔT₁×ΔT₂ ≤ 0 (invalid).
 */
double lmtd(double delta_T1, double delta_T2);

/**
 * Heat exchanger duty via LMTD method.
 * q = U·A·F·ΔT_lm
 *
 * @param U     Overall heat transfer coefficient [W/(m²·K)]
 * @param A     Heat transfer area [m²]
 * @param F     Correction factor (F=1 for pure counter/parallel)
 * @param dT_lm LMTD [K]
 * @return Heat duty q [W]
 */
double heat_exchanger_duty_lmtd(double U, double A, double F, double dT_lm);

/**
 * Effectiveness-NTU method.
 *
 * ε = q / q_max = f(NTU, C_r, flow arrangement)
 *
 * NTU = U·A / C_min
 * C_r = C_min / C_max
 *
 * @param ntu   Number of Transfer Units
 * @param Cr    Capacity rate ratio
 * @param arrangement "parallel", "counter", "shell_tube_1", "cross_both_unmixed"
 * @return Effectiveness ε [0, 1]
 */
double heat_exchanger_effectiveness(double ntu, double Cr,
                                    const char *arrangement);

/**
 * NTU from effectiveness (inverse relation, parallel flow).
 * NTU = -ln(1 - ε·(1+Cr)) / (1+Cr)
 */
double ntu_from_effectiveness_parallel(double epsilon, double Cr);

/**
 * NTU from effectiveness (counter flow, Cr < 1).
 * NTU = (1/(Cr-1)) · ln((ε-1)/(ε·Cr-1))
 */
double ntu_from_effectiveness_counter(double epsilon, double Cr);

/**
 * Heat capacity rate C [W/K].
 * C = ṁ · cp
 */
double heat_capacity_rate(double mass_flow_rate, double specific_heat);

/**
 * Heat exchanger outlet temperature computation.
 * Given hot and cold inlet temperatures, C_h, C_c, ε.
 */
void heat_exchanger_outlet_temps(double Th_in, double Tc_in,
                                 double Ch, double Cc,
                                 double epsilon,
                                 double *Th_out, double *Tc_out);

/**
 * Correction factor F for shell-and-tube (1 shell, 2 tube passes).
 *
 * F = (√(R²+1)·ln[(1-P)/(1-P·R)]) / ((R-1)·ln[(2-P(R+1-√(R²+1)))/
 *                                         (2-P(R+1+√(R²+1)))])
 *
 * Source: Bowman, Mueller, Nagle (1940), Trans. ASME 62:283
 */
double hx_correction_factor_shell_tube(double P, double R);

#endif /* HEAT_TRANSFER_H */
