/**
 * mini-fluid-thermal-modeling — Fluid Mechanics
 *
 * Incompressible and compressible flow analysis:
 * Bernoulli, Navier-Stokes simplifications, pipe flow,
 * boundary layer theory, drag/lift, and dimensional analysis.
 *
 * Reference:
 *   White, "Fluid Mechanics" (2016), Chapters 3-7
 *   Munson, Young, Okiishi, "Fundamentals of Fluid Mechanics" (2013)
 *   Schlichting & Gersten, "Boundary Layer Theory" (2016)
 *
 * Knowledge Levels:
 *   L3: Engineering quantities (Re/Pr/Nu ranges, friction factors)
 *   L4: Conservation laws (mass, momentum, energy — Bernoulli, NS)
 *   L5: Engineering methods (Darcy-Weisbach, Moody, dimensional analysis)
 */

#ifndef FLUID_MECHANICS_H
#define FLUID_MECHANICS_H

#include "fluid_thermal_core.h"

/* ---------------------------------------------------------------------------
 * L3: Engineering Quantities — Typical Value Ranges
 * ------------------------------------------------------------------------- */

/** Reynolds number ranges for common engineering flows */
typedef struct {
    double min_re;
    double max_re;
    const char *description;
} ReynoldsRange;

/** Get typical Reynolds number range for a given application.
 *
 *  Supported: "pipe_water", "pipe_air", "open_channel",
 *             "aircraft_wing", "automobile", "microchannel",
 *             "blood_vessel", "heat_exchanger_tube", "wind_turbine"
 */
ReynoldsRange typical_reynolds_range(const char *application);

/** Get typical Prandtl number for common fluids at ~300K */
double typical_prandtl(const char *fluid_name);

/** Get typical convective heat transfer coefficient range [W/(m²·K)].
 *  "natural_air", "forced_air", "forced_water", "boiling_water",
 *  "condensing_steam", "oil_forced"
 */
void typical_h_range(const char *scenario, double *h_min, double *h_max);

/* ---------------------------------------------------------------------------
 * L4: Conservation Laws
 * ------------------------------------------------------------------------- */

/**
 * Conservation of mass (continuity equation, 1D steady).
 * ρ₁·A₁·u₁ = ρ₂·A₂·u₂
 *
 * For incompressible flow (ρ₁=ρ₂): A₁·u₁ = A₂·u₂
 *
 * @param density1  ρ₁ [kg/m³]
 * @param area1     A₁ [m²]
 * @param velocity1 u₁ [m/s]
 * @param density2  ρ₂ [kg/m³]
 * @param area2     A₂ [m²]
 * @param velocity2 Output: u₂ [m/s]
 * @return 0 on success, -1 if area1=0
 *
 * Source: Euler (1757), "Principes généraux du mouvement des fluides"
 */
int continuity_equation(double density1, double area1, double velocity1,
                        double density2, double area2, double *velocity2);

/**
 * Bernoulli equation along a streamline (steady, inviscid, incompressible).
 * p₁ + ½ρ·u₁² + ρ·g·z₁ = p₂ + ½ρ·u₂² + ρ·g·z₂
 *
 * Given 6 of 7 parameters, solve for the unknown.
 * Set the unknown parameter to NaN to solve for it.
 *
 * @param p1,p2    Static pressures [Pa]
 * @param u1,u2    Velocities [m/s]
 * @param z1,z2    Elevations [m]
 * @param rho      Density [kg/m³]
 * @param g        Gravitational acceleration [m/s²]
 * @return Solved value on success, NaN on failure
 *
 * Source: Daniel Bernoulli (1738), "Hydrodynamica"
 */
double bernoulli_solve(double p1, double p2, double u1, double u2,
                       double z1, double z2, double rho, double g);

/**
 * Bernoulli equation with head loss.
 * p₁/γ + u₁²/(2g) + z₁ = p₂/γ + u₂²/(2g) + z₂ + hL
 *
 * @return Head loss hL [m] needed to satisfy the equation
 */
double bernoulli_head_loss(double p1, double p2, double u1, double u2,
                           double z1, double z2, double gamma, double g);

/* ---------------------------------------------------------------------------
 * L4: Navier-Stokes — Simplified Analytical Solutions
 * ------------------------------------------------------------------------- */

/**
 * Plane Poiseuille flow: steady, 2D, fully-developed laminar flow
 * between parallel plates (gap height h).
 *
 * u(y) = (1/(2μ)) · (-dp/dx) · (h²/4 - y²)
 * where y ∈ [-h/2, h/2]
 *
 * Maximum velocity at centerline y=0: u_max = (h²/(8μ)) · (-dp/dx)
 * Volumetric flow rate per unit width: q = (h³/(12μ)) · (-dp/dx)
 *
 * @param dp_dx     Pressure gradient [Pa/m] (negative for flow in +x)
 * @param mu        Dynamic viscosity [Pa·s]
 * @param h         Gap height [m]
 * @param y         Position from centerline [m] ∈ [-h/2, h/2]
 * @return Velocity u(y) [m/s]
 *
 * Source: Navier (1823), Poiseuille (1840)
 */
double plane_poiseuille_velocity(double dp_dx, double mu, double h, double y);

/**
 * Plane Poiseuille volumetric flow rate per unit width.
 */
double plane_poiseuille_flow_rate(double dp_dx, double mu, double h);

/**
 * Hagen-Poiseuille flow: steady, fully-developed laminar flow in a circular pipe.
 *
 * u(r) = (1/(4μ)) · (-dp/dx) · (R² - r²)
 * Q = (π·R⁴/(8μ)) · (-dp/dx)
 * u_avg = Q/(πR²) = (R²/(8μ)) · (-dp/dx)
 *
 * @param dp_dx  Pressure gradient [Pa/m]
 * @param mu     Dynamic viscosity [Pa·s]
 * @param R      Pipe radius [m]
 * @param r      Radial position [m] ∈ [0, R]
 * @return Velocity u(r) [m/s]
 *
 * Source: Hagen (1839), Poiseuille (1840)
 */
double hagen_poiseuille_velocity(double dp_dx, double mu, double R, double r);

/**
 * Hagen-Poiseuille volumetric flow rate.
 */
double hagen_poiseuille_flow_rate(double dp_dx, double mu, double R);

/**
 * Couette flow: steady flow between two parallel plates,
 * lower plate stationary, upper plate moving at velocity U.
 * Zero pressure gradient.
 *
 * u(y) = U · y / h
 *
 * @param U  Upper plate velocity [m/s]
 * @param h  Gap height [m]
 * @param y  Position [m] ∈ [0, h]
 * @return Velocity [m/s]
 */
double couette_velocity(double U, double h, double y);

/**
 * Couette-Poiseuille flow: combined pressure-driven and shear-driven.
 *
 * u(y) = U·y/h + (h²/(2μ))·(-dp/dx)·(y/h)·(1 - y/h)
 *
 * @param U     Upper plate velocity [m/s]
 * @param dp_dx Pressure gradient [Pa/m]
 * @param mu    Dynamic viscosity [Pa·s]
 * @param h     Gap height [m]
 * @param y     Position [m] ∈ [0, h]
 * @return Velocity [m/s]
 */
double couette_poiseuille_velocity(double U, double dp_dx,
                                   double mu, double h, double y);

/**
 * Stokes first problem: impulsively started flat plate.
 * Fluid initially at rest, plate suddenly moves at constant velocity U.
 *
 * u(y, t) = U · erfc(y / (2·√(ν·t)))
 *
 * @param U   Plate velocity [m/s]
 * @param nu  Kinematic viscosity [m²/s]
 * @param y   Distance from plate [m]
 * @param t   Time since impulse [s]
 * @return Velocity u(y,t) [m/s]
 *
 * Source: Stokes (1851)
 */
double stokes_first_problem(double U, double nu, double y, double t);

/* ---------------------------------------------------------------------------
 * L5: Engineering Methods — Pipe Flow
 * ------------------------------------------------------------------------- */

/**
 * Darcy friction factor for laminar flow (Re < 2300).
 * f = 64 / Re
 *
 * Source: Darcy (1857), Weisbach (1845)
 */
double darcy_friction_laminar(double reynolds_number);

/**
 * Colebrook-White equation for turbulent friction factor.
 * 1/√f = -2.0·log₁₀( (ε/D)/3.7 + 2.51/(Re·√f) )
 *
 * Solved iteratively via Newton-Raphson.
 *
 * @param reynolds_number  Re
 * @param roughness_ratio  ε/D (relative roughness)
 * @param max_iter         Maximum iterations (suggest 20)
 * @param tolerance        Convergence tolerance (suggest 1e-10)
 * @return Darcy friction factor f
 *
 * Source: Colebrook & White (1937), "Turbulent flow in pipes"
 */
double colebrook_friction_factor(double reynolds_number,
                                 double roughness_ratio,
                                 int max_iter, double tolerance);

/**
 * Explicit approximation of Colebrook-White (Swamee-Jain, 1976).
 * Valid for: 5000 ≤ Re ≤ 10⁸, 10⁻⁶ ≤ ε/D ≤ 0.05
 *
 * f = 0.25 / [log₁₀(ε/D/3.7 + 5.74/Re^0.9)]²
 *
 * Source: Swamee & Jain (1976), "Explicit equations for pipe-flow problems"
 */
double swamee_jain_friction_factor(double reynolds_number,
                                   double roughness_ratio);

/**
 * Darcy-Weisbach head loss.
 * hL = f · (L/D) · (u²/(2g))
 *
 * @param f     Darcy friction factor
 * @param L     Pipe length [m]
 * @param D     Pipe diameter [m]
 * @param u     Average velocity [m/s]
 * @param g     Gravity [m/s²]
 * @return Head loss [m]
 *
 * Source: Weisbach (1845), "Lehrbuch der Ingenieur- und Maschinen-Mechanik"
 */
double darcy_weisbach_head_loss(double f, double L, double D,
                                double u, double g);

/**
 * Darcy-Weisbach pressure drop.
 * Δp = f · (L/D) · (½ρ·u²)
 */
double darcy_weisbach_pressure_drop(double f, double L, double D,
                                    double rho, double u);

/**
 * Minor loss coefficient for common fittings.
 *
 * Supported fittings:
 *   "elbow_90_standard", "elbow_45", "tee_straight", "tee_branch",
 *   "gate_valve_open", "globe_valve_open", "ball_valve_open",
 *   "sudden_expansion", "sudden_contraction", "pipe_entrance",
 *   "pipe_exit", "bend_90_r/d=1", "bend_90_r/d=2"
 *
 * Values from Crane Technical Paper No. 410.
 */
double minor_loss_coefficient(const char *fitting_type);

/**
 * Hydraulic diameter for non-circular ducts.
 * Dh = 4 · A / P
 *
 * @param area    Cross-sectional area [m²]
 * @param perimeter Wetted perimeter [m]
 * @return Hydraulic diameter [m]
 */
double hydraulic_diameter(double area, double perimeter);

/* ---------------------------------------------------------------------------
 * L5: Boundary Layer Theory
 * ------------------------------------------------------------------------- */

/**
 * Blasius solution: laminar flat-plate boundary layer thickness.
 * δ = 5.0 · x / √Re_x
 *
 * @param x   Distance from leading edge [m]
 * @param rex Reynolds number at x: Re_x = u∞·x/ν
 * @return Boundary layer thickness δ [m]
 *
 * Source: Blasius (1908), "Grenzschichten in Flüssigkeiten"
 */
double blasius_boundary_layer_thickness(double x, double rex);

/**
 * Blasius laminar skin friction coefficient (local).
 * cf = 0.664 / √Re_x
 */
double blasius_skin_friction_local(double rex);

/**
 * Blasius laminar skin friction coefficient (average over plate length L).
 * Cf = 1.328 / √Re_L
 */
double blasius_skin_friction_average(double re_l);

/**
 * Turbulent flat-plate boundary layer thickness (1/7 power law).
 * δ = 0.37 · x / Re_x^{1/5}
 *
 * Valid for: 5×10⁵ ≤ Re_x ≤ 10⁷
 *
 * Source: Prandtl (1927), von Kármán (1930)
 */
double turbulent_boundary_layer_thickness(double x, double rex);

/**
 * Turbulent skin friction coefficient (local).
 * cf = 0.0592 / Re_x^{1/5}
 *
 * Source: Prandtl-Schlichting formula
 */
double turbulent_skin_friction_local(double rex);

/**
 * Turbulent skin friction coefficient (average).
 * Cf = 0.074 / Re_L^{1/5}
 *
 * Valid for: 5×10⁵ ≤ Re_L ≤ 10⁷
 */
double turbulent_skin_friction_average(double re_l);

/**
 * Transition Reynolds number estimation.
 * Uses Michel's criterion for transition onset.
 *
 * Re_θ,tr ≈ 1.174 · (1 + 22400/Re_x,tr) · Re_x,tr^0.46
 *
 * Simplified engineering estimate: Re_tr ≈ 5×10⁵ for smooth plates.
 */
double estimate_transition_reynolds(double turbulence_intensity);

/**
 * Displacement thickness δ* for laminar flat plate.
 * δ* = 1.7208 · x / √Re_x
 */
double displacement_thickness_laminar(double x, double rex);

/**
 * Momentum thickness θ for laminar flat plate.
 * θ = 0.664 · x / √Re_x
 */
double momentum_thickness_laminar(double x, double rex);

/**
 * Shape factor H = delta_star / theta for laminar flat plate.
 * H approx 2.59 (Blasius)
 */
double shape_factor_laminar(void);

/* ---------------------------------------------------------------------------
 * L5: Dimensional Analysis — Buckingham Pi Theorem
 * ------------------------------------------------------------------------- */

/**
 * Buckingham Pi theorem: determine number of dimensionless groups.
 * For n variables with k fundamental dimensions, expect (n - k) Pi groups.
 *
 * @param n_variables  Number of physical variables
 * @param n_dimensions Number of fundamental dimensions involved
 * @return Expected number of dimensionless Pi groups
 */
int buckingham_pi_count(int n_variables, int n_dimensions);

/**
 * Compute common dimensionless group given its name.
 *
 * Supported groups:
 *   "Re" (Reynolds), "Pr" (Prandtl), "Nu" (Nusselt),
 *   "Gr" (Grashof), "Ra" (Rayleigh), "Pe" (Peclet),
 *   "Ma" (Mach), "Fr" (Froude), "We" (Weber),
 *   "Ca" (Capillary), "Eu" (Euler), "St" (Strouhal)
 *
 * @param group_name Group identifier
 * @param params     Array of physical parameters (order depends on group)
 * @return Computed value, NaN if unknown group
 */
double dimensionless_group(const char *group_name, const double params[]);

/**
 * Froude number: Fr = u / √(g·L)
 * Ratio of inertial to gravitational forces.
 */
double compute_froude(double velocity, double g, double char_length);

/**
 * Weber number: We = ρ·u²·L / σ
 * Ratio of inertial to surface tension forces.
 */
double compute_weber(double density, double velocity,
                     double char_length, double surface_tension);

/**
 * Euler number: Eu = Δp / (ρ·u²)
 * Ratio of pressure to inertial forces.
 */
double compute_euler(double pressure_drop, double density, double velocity);

/**
 * Strouhal number: St = f·L / u
 * Ratio of unsteady to steady inertial forces.
 */
double compute_strouhal(double frequency, double char_length, double velocity);

#endif /* FLUID_MECHANICS_H */
