/**
 * mini-fluid-thermal-modeling — Coupled Thermal-Fluid Systems
 *
 * Combined fluid-thermal analysis: pipe flow with heat transfer,
 * pump/fan system curves, system-level thermal management,
 * HVAC fundamentals, electronics cooling, and natural circulation loops.
 *
 * Reference:
 *   Kays & London, "Compact Heat Exchangers" (1984)
 *   Cengel & Cimbala, "Fluid Mechanics: Fundamentals and Applications" (2013)
 *   Nellis & Klein, "Heat Transfer" (2009)
 *
 * Knowledge Levels:
 *   L6: Engineering problems (pipe sizing + heating, cooling system design)
 *   L7: Applications (electronics cooling, HVAC, automotive thermal)
 *   L8: Advanced methods (system coupling, conjugate analysis)
 */

#ifndef THERMAL_FLUID_SYSTEMS_H
#define THERMAL_FLUID_SYSTEMS_H

#include "fluid_thermal_core.h"

/* ---------------------------------------------------------------------------
 * L6: Pipe Flow with Heat Transfer
 * ------------------------------------------------------------------------- */

/**
 * Temperature rise of fluid flowing through a heated pipe.
 *
 * T_out = T_in + q_w · π · D · L / (ṁ · cp)
 *
 * where q_w is the uniform wall heat flux [W/m²].
 *
 * @param T_in   Inlet temperature [K]
 * @param q_w    Wall heat flux [W/m²]
 * @param D      Pipe diameter [m]
 * @param L      Pipe length [m]
 * @param m_dot  Mass flow rate [kg/s]
 * @param cp     Specific heat [J/(kg·K)]
 * @return Outlet temperature [K]
 */
double pipe_heating_outlet_temp(double T_in, double q_w,
                                double D, double L,
                                double m_dot, double cp);

/**
 * Wall temperature in constant-heat-flux pipe (fully developed).
 *
 * Ts(x) = T_in + (q_w·π·D·L/(ṁ·cp))·(x/L) + q_w/h
 *
 * @param x   Distance from inlet [m]
 * @param h   Local heat transfer coefficient [W/(m²·K)]
 * @return Wall temperature at x [K]
 */
double pipe_wall_temperature_chf(double T_in, double q_w, double D,
                                 double L, double m_dot, double cp,
                                 double x, double h);

/**
 * Outlet temperature for constant-wall-temperature pipe.
 *
 * T_out = Ts - (Ts - T_in)·exp(-π·D·L·h/(ṁ·cp))
 *
 * @param Ts   Wall temperature [K] (constant)
 * @param T_in Inlet temperature [K]
 * @param h    Average heat transfer coefficient [W/(m²·K)]
 * @return Outlet temperature [K]
 */
double pipe_outlet_temp_cwt(double Ts, double T_in, double D, double L,
                            double h, double m_dot, double cp);

/**
 * Required pipe length for given temperature change (constant wall temp).
 *
 * L = -(ṁ·cp/(π·D·h)) · ln((T_out - Ts)/(T_in - Ts))
 */
double pipe_length_for_temp_change(double Ts, double T_in, double T_out,
                                   double D, double h,
                                   double m_dot, double cp);

/**
 * Nusselt number for fully-developed laminar flow in circular pipe.
 *
 * Constant heat flux:  Nu = 4.36
 * Constant wall temp:   Nu = 3.66
 *
 * Valid: Re < 2300, hydrodynamically and thermally fully developed.
 */
double nusselt_laminar_pipe(int constant_heat_flux);

/* ---------------------------------------------------------------------------
 * L6: Pump and Fan System Curves
 * ------------------------------------------------------------------------- */

/**
 * Pump hydraulic power.
 * P_hyd = γ · Q · H = ρ · g · Q · H
 *
 * @param rho  Density [kg/m³]
 * @param g    Gravity [m/s²]
 * @param Q    Volumetric flow rate [m³/s]
 * @param H    Head [m]
 * @return Hydraulic power [W]
 */
double pump_hydraulic_power(double rho, double g, double Q, double H);

/**
 * Pump shaft power (incorporating efficiency).
 * P_shaft = P_hyd / η
 */
double pump_shaft_power(double hydraulic_power, double efficiency);

/**
 * System curve: head required at given flow rate.
 * H_sys = H_static + K_sys · Q²
 *
 * @param H_static  Static head (elevation + pressure difference) [m]
 * @param K_sys     System loss coefficient [(m)/(m³/s)²]
 * @param Q         Volumetric flow rate [m³/s]
 * @return Total system head [m]
 */
double system_curve_head(double H_static, double K_sys, double Q);

/**
 * Pump curve (quadratic model).
 * H_pump = H_shutoff - K_pump · Q²
 *
 * @param H_shutoff  Head at zero flow [m]
 * @param K_pump     Pump loss coefficient
 * @param Q          Volumetric flow rate [m³/s]
 * @return Pump head [m]
 */
double pump_curve_head(double H_shutoff, double K_pump, double Q);

/**
 * Find operating point where pump curve intersects system curve.
 * H_shutoff - K_pump·Q² = H_static + K_sys·Q²
 * → Q = √((H_shutoff - H_static) / (K_pump + K_sys))
 *
 * @param H_shutoff  Pump shutoff head [m]
 * @param H_static   System static head [m]
 * @param K_pump     Pump curve coefficient
 * @param K_sys      System curve coefficient
 * @param Q_op       Output: operating flow rate [m³/s]
 * @param H_op       Output: operating head [m]
 * @return 0 on success, -1 if no intersection
 */
int pump_system_operating_point(double H_shutoff, double H_static,
                                double K_pump, double K_sys,
                                double *Q_op, double *H_op);

/**
 * Net Positive Suction Head Available (NPSHa).
 * NPSHa = (p_atm - p_vapor)/γ + z_suction - hL_suction
 *
 * @param p_atm    Atmospheric pressure [Pa]
 * @param p_vapor  Vapor pressure at fluid temperature [Pa]
 * @param gamma    Specific weight [N/m³]
 * @param z_suction Suction height (negative if pump below source) [m]
 * @param hL_suction Suction line head loss [m]
 * @return NPSHa [m]
 */
double npsh_available(double p_atm, double p_vapor, double gamma,
                      double z_suction, double hL_suction);

/* ---------------------------------------------------------------------------
 * L6: Natural Circulation Loops
 * ------------------------------------------------------------------------- */

/**
 * Driving pressure head for natural circulation.
 * Δp_drive = ρ₀ · g · β · ∮ ΔT · dz
 *
 * Simplified: Δp_drive = ρ₀ · g · β · ΔT · H_loop
 *
 * @param rho0   Reference density [kg/m³]
 * @param g      Gravity [m/s²]
 * @param beta   Thermal expansion coefficient [1/K]
 * @param deltaT Temperature difference (hot - cold) [K]
 * @param H_loop Height of loop [m]
 * @return Driving pressure [Pa]
 */
double natural_circulation_driving_pressure(double rho0, double g,
                                            double beta, double deltaT,
                                            double H_loop);

/**
 * Flow rate in natural circulation loop (turbulent friction).
 * ṁ = [ (2·D·A²·ρ²·g·β·ΔT·H) / (f·L_eff) ]^{1/2}
 *
 * @param D       Pipe diameter [m]
 * @param A       Flow area [m²]
 * @param rho     Density [kg/m³]
 * @param g       Gravity [m/s²]
 * @param beta    Expansion coefficient [1/K]
 * @param deltaT  Hot-cold difference [K]
 * @param H       Height [m]
 * @param f       Darcy friction factor
 * @param L_eff   Effective loop length [m]
 * @return Mass flow rate [kg/s]
 */
double natural_circulation_flow_rate(double D, double A, double rho,
                                     double g, double beta, double deltaT,
                                     double H, double f, double L_eff);

/* ---------------------------------------------------------------------------
 * L7: Electronics Cooling Applications
 * ------------------------------------------------------------------------- */

/**
 * Heat sink thermal resistance model.
 *
 * R_total = R_spreading + R_fin + R_convection
 *
 * R_spreading ≈ 1/(2·k·√A) (semi-infinite approximation)
 * R_convection = 1/(η_f·h·A_fin_total)
 *
 * @param k_base      Base thermal conductivity [W/(m·K)]
 * @param A_base      Base area [m²]
 * @param eta_fin     Overall fin efficiency
 * @param h           Convection coefficient [W/(m²·K)]
 * @param A_fin_total Total fin surface area [m²]
 * @return Total thermal resistance [K/W]
 */
double heat_sink_thermal_resistance(double k_base, double A_base,
                                    double eta_fin, double h,
                                    double A_fin_total);

/**
 * Junction temperature in electronics package.
 * T_junction = T_ambient + q · (R_jc + R_cs + R_sa)
 *
 * @param T_ambient Ambient temperature [K]
 * @param q         Heat dissipation [W]
 * @param R_jc      Junction-to-case resistance [K/W]
 * @param R_cs      Case-to-sink resistance [K/W]
 * @param R_sa      Sink-to-ambient resistance [K/W]
 * @return Junction temperature [K]
 */
double junction_temperature(double T_ambient, double q,
                            double R_jc, double R_cs, double R_sa);

/**
 * Required airflow for electronics cooling.
 *
 * Q_air = q_total / (ρ_air · cp_air · ΔT_allowable)
 *
 * @param q_total    Total heat load [W]
 * @param rho_air    Air density [kg/m³]
 * @param cp_air     Air specific heat [J/(kg·K)]
 * @param deltaT     Allowable temperature rise [K]
 * @return Required volumetric airflow [m³/s]
 */
double required_airflow_cooling(double q_total, double rho_air,
                                double cp_air, double deltaT);

/**
 * Chip temperature in forced convection (3D IC thermal analysis).
 *
 * Uses thermal wake function for stacked dies.
 *
 * T_die_n = T_in + Σ_{i=1}^{n} q_i/(ṁ·cp) + q_n/(h_n·A_n)
 *
 * @param T_in      Coolant inlet temperature [K]
 * @param n_dies    Number of stacked dies
 * @param q_power   Array of die powers [W] (length n_dies)
 * @param m_dot     Coolant mass flow rate [kg/s]
 * @param cp        Coolant specific heat [J/(kg·K)]
 * @param h         Array of heat transfer coefficients [W/(m²·K)]
 * @param A         Array of die areas [m²]
 * @param T_die     Output: die temperatures [K] (length n_dies)
 */
void stacked_die_temperatures(double T_in, int n_dies,
                              const double q_power[],
                              double m_dot, double cp,
                              const double h[], const double A[],
                              double T_die[]);

/* ---------------------------------------------------------------------------
 * L7: HVAC and Building Thermal Analysis
 * ------------------------------------------------------------------------- */

/**
 * Sensible heating/cooling load.
 * q_sensible = ṁ · cp · (T_supply - T_room)
 */
double sensible_heat_load(double m_dot, double cp,
                          double T_supply, double T_room);

/**
 * Latent heat load (dehumidification).
 * q_latent = ṁ · h_fg · (ω_room - ω_supply)
 *
 * @param h_fg  Latent heat of vaporization [J/kg]
 * @param w     Humidity ratio [kg water / kg dry air]
 */
double latent_heat_load(double m_dot, double h_fg,
                        double w_room, double w_supply);

/**
 * Room air change rate.
 * ACH = Q_air · 3600 / V_room
 *
 * @param Q_air   Supply airflow [m³/s]
 * @param V_room  Room volume [m³]
 * @return Air changes per hour [1/hr]
 */
double air_changes_per_hour(double Q_air, double V_room);

/**
 * Mean radiant temperature (MRT) from surrounding surfaces.
 * MRT ≈ Σ (T_i · A_i) / Σ A_i  (area-weighted average)
 *
 * @param T_surfaces  Array of surface temperatures [K]
 * @param A_surfaces  Array of surface areas [m²]
 * @param n_surfaces  Number of surfaces
 * @return Mean radiant temperature [K]
 */
double mean_radiant_temperature(const double T_surfaces[],
                                const double A_surfaces[],
                                int n_surfaces);

/* ---------------------------------------------------------------------------
 * L7: Automotive / Industrial Cooling
 * ------------------------------------------------------------------------- */

/**
 * Radiator heat rejection model (ε-NTU method for compact HX).
 *
 * q_rad = ε · C_min · (Th_in - Tc_in)
 *
 * Typical automotive radiator: cross-flow, both unmixed.
 */
double radiator_heat_rejection(double T_coolant_in, double T_air_in,
                               double m_dot_coolant, double cp_coolant,
                               double m_dot_air, double cp_air,
                               double UA);

/**
 * Coolant pump power scaling.
 * P_pump ∝ Q³ (affinity laws)
 *
 * P₂ = P₁ · (N₂/N₁)³
 */
double pump_power_scaling(double P_ref, double N_new, double N_ref);

/**
 * Fan affinity laws.
 *
 * Q₂/Q₁ = N₂/N₁
 * H₂/H₁ = (N₂/N₁)²
 * P₂/P₁ = (N₂/N₁)³
 */
void fan_affinity_laws(double N_ratio,
                       double Q_ref, double H_ref, double P_ref,
                       double *Q_new, double *H_new, double *P_new);

/* ---------------------------------------------------------------------------
 * L8: System-Level Coupling
 * ------------------------------------------------------------------------- */

/**
 * Coupled fluid-thermal solver convergence check.
 * Residual R = |T_new - T_old|/T_old for each node.
 *
 * @param T_new    Current iteration temperatures [K]
 * @param T_old    Previous iteration temperatures [K]
 * @param n        Number of nodes
 * @return Maximum relative residual
 */
double thermal_convergence_residual(const double T_new[],
                                    const double T_old[], int n);

/**
 * Conjugate heat transfer Nusselt number correction.
 *
 * For gas flows with finite wall conductivity, the true Nusselt
 * number is reduced due to axial conduction in the wall.
 *
 * Nu_true / Nu = 1 / (1 + k_wall·A_wall/(ṁ·cp·L) · Nu)
 *
 * @param Nu      Nusselt number with isothermal wall assumption
 * @param k_wall  Wall thermal conductivity [W/(m·K)]
 * @param A_wall  Wall cross-section area [m²]
 * @param m_dot   Mass flow rate [kg/s]
 * @param cp      Fluid specific heat [J/(kg·K)]
 * @param L       Characteristic length [m]
 * @return Corrected Nusselt number
 */
double conjugate_nusselt_correction(double Nu, double k_wall,
                                    double A_wall, double m_dot,
                                    double cp, double L);

/**
 * Temperature effectiveness of a counterflow recuperator.
 *
 * ε = NTU / (1 + NTU)  when Cr = 1 (balanced flow)
 */
double recuperator_effectiveness_balanced(double ntu);

/**
 * Staggered tube bank pressure drop (Zukauskas correlation).
 *
 * Δp = N_L · χ · (ρ·u_max²/2) · f
 *
 * @param N_L    Number of tube rows
 * @param chi    Correction factor
 * @param rho    Density [kg/m³]
 * @param u_max  Maximum velocity [m/s]
 * @param f      Friction factor
 * @return Pressure drop [Pa]
 */
double tube_bank_pressure_drop(int N_L, double chi, double rho,
                               double u_max, double f);

/* ---------------------------------------------------------------------------
 * L8: Multi-Phase Flow Basics
 * ------------------------------------------------------------------------- */

/**
 * Void fraction (homogeneous model).
 *
 * α = x / [x + (1-x)·(ρ_v/ρ_l)·S]
 *
 * where x = quality, S = slip ratio.
 * For homogeneous model, S = 1.
 *
 * @param quality   Vapor mass fraction x
 * @param rho_v     Vapor density [kg/m³]
 * @param rho_l     Liquid density [kg/m³]
 * @param slip_ratio Slip ratio S (use 1.0 for homogeneous)
 * @return Void fraction α [dimensionless]
 */
double void_fraction_homogeneous(double quality, double rho_v,
                                 double rho_l, double slip_ratio);

/**
 * Two-phase frictional multiplier (Lockhart-Martinelli).
 *
 * φ_l² = 1 + C/X + 1/X²
 *
 * where X = (dp/dx)_l/(dp/dx)_v is the Martinelli parameter,
 * C depends on flow regime (C ≈ 20 for turbulent-turbulent).
 *
 * @param X     Martinelli parameter
 * @param C     Chisholm coefficient (regime-dependent)
 * @return Liquid-phase multiplier φ_l²
 */
double two_phase_multiplier(double X, double C);

/**
 * Lockhart-Martinelli parameter.
 *
 * X = ( (1-x)/x )^0.9 · (ρ_v/ρ_l)^0.5 · (μ_l/μ_v)^0.1
 *
 * @param quality  Vapor mass fraction
 * @param rho_v    Vapor density [kg/m³]
 * @param rho_l    Liquid density [kg/m³]
 * @param mu_v     Vapor viscosity [Pa·s]
 * @param mu_l     Liquid viscosity [Pa·s]
 * @return Martinelli parameter X
 */
double martinelli_parameter(double quality, double rho_v, double rho_l,
                            double mu_v, double mu_l);

#endif /* THERMAL_FLUID_SYSTEMS_H */
