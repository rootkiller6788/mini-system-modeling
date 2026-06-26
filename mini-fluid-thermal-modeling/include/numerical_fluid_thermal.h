/**
 * mini-fluid-thermal-modeling — Numerical Methods
 *
 * Finite difference solvers for conduction, iterative methods
 * for pipe networks, and numerical integration for boundary layers.
 *
 * Reference:
 *   Patankar, "Numerical Heat Transfer and Fluid Flow" (1980)
 *   Anderson, "Computational Fluid Dynamics: The Basics" (1995)
 *   Versteeg & Malalasekera, "An Introduction to CFD: FVM" (2007)
 *
 * Knowledge Levels:
 *   L5: Finite difference methods, relaxation solvers
 *   L6: Steady/transient conduction problems
 *   L8: CFD basics, SIMPLE algorithm outline
 */

#ifndef NUMERICAL_FLUID_THERMAL_H
#define NUMERICAL_FLUID_THERMAL_H

#include "fluid_thermal_core.h"

/* ---------------------------------------------------------------------------
 * L5: 1D Steady Conduction — Finite Difference
 * ------------------------------------------------------------------------- */

/**
 * 1D steady conduction with internal heat generation.
 * Solves: d²T/dx² + q̇/k = 0  on [0, L]
 *
 * Dirichlet BC: T(0) = T_left, T(L) = T_right
 *
 * Finite difference: T_i = (T_{i-1} + T_{i+1} + q̇·Δx²/k) / 2
 *
 * Uses Gauss-Seidel iteration.
 *
 * @param n          Number of grid points (≥ 3)
 * @param L          Domain length [m]
 * @param k          Thermal conductivity [W/(m·K)]
 * @param q_dot      Volumetric heat generation [W/m³]
 * @param T_left     Left boundary temperature [K]
 * @param T_right    Right boundary temperature [K]
 * @param T          Input: initial guess / Output: solution [K], size n
 * @param max_iter   Maximum iterations
 * @param tol        Convergence tolerance
 * @return Number of iterations performed, -1 if did not converge
 *
 * Source: Patankar (1980), Chapter 3
 */
int solve_1d_steady_conduction(int n, double L, double k, double q_dot,
                               double T_left, double T_right,
                               double T[], int max_iter, double tol);

/**
 * 1D transient conduction (explicit scheme).
 * Solves: ∂T/∂t = α·∂²T/∂x²
 *
 * T_i^{n+1} = Fo·(T_{i-1}^n + T_{i+1}^n) + (1 - 2·Fo)·T_i^n
 *
 * Stability: Fo = α·Δt/Δx² ≤ 0.5
 *
 * @param nx     Number of spatial grid points
 * @param nt     Number of time steps
 * @param L      Domain length [m]
 * @param alpha  Thermal diffusivity [m²/s]
 * @param dt     Time step [s]
 * @param T_left Left boundary temperature [K] (constant)
 * @param T_right Right boundary temperature [K] (constant)
 * @param T0     Initial temperature field [K], size nx
 * @param T      Output: temperature after nt steps [K], size nx
 * @return 0 on success, -1 if stability violated (Fo > 0.5)
 */
int solve_1d_transient_conduction_explicit(int nx, int nt, double L,
                                           double alpha, double dt,
                                           double T_left, double T_right,
                                           const double T0[], double T[]);

/**
 * 1D transient conduction (implicit/Crank-Nicolson scheme).
 *
 * Unconditionally stable. Uses Thomas algorithm for tridiagonal solve.
 *
 * @param nx     Number of spatial grid points
 * @param nt     Number of time steps
 * @param L      Domain length [m]
 * @param alpha  Thermal diffusivity [m²/s]
 * @param dt     Time step [s]
 * @param T_left Left boundary temperature [K]
 * @param T_right Right boundary temperature [K]
 * @param T0     Initial temperature [K], size nx
 * @param T      Output: temperature after nt steps [K], size nx
 * @return 0 on success, -1 on invalid input
 */
int solve_1d_transient_conduction_implicit(int nx, int nt, double L,
                                           double alpha, double dt,
                                           double T_left, double T_right,
                                           const double T0[], double T[]);

/* ---------------------------------------------------------------------------
 * L5: Thomas Algorithm (Tridiagonal Solver)
 * ------------------------------------------------------------------------- */

/**
 * Thomas algorithm for solving tridiagonal linear systems.
 *
 * a_i·x_{i-1} + b_i·x_i + c_i·x_{i+1} = d_i  for i=0..n-1
 *
 * a[0] and c[n-1] are unused.
 *
 * @param n  System size
 * @param a  Sub-diagonal (size n, a[0] unused)
 * @param b  Main diagonal (size n)
 * @param c  Super-diagonal (size n, c[n-1] unused)
 * @param d  RHS vector (size n)
 * @param x  Solution vector (size n)
 * @return 0 on success, -1 if singular
 */
int thomas_algorithm(int n, const double a[], const double b[],
                     const double c[], const double d[], double x[]);

/* ---------------------------------------------------------------------------
 * L5: Gauss-Seidel Relaxation
 * ------------------------------------------------------------------------- */

/**
 * Gauss-Seidel solver for 2D Laplace equation (steady conduction).
 * Solves: ∂²T/∂x² + ∂²T/∂y² = 0
 *
 * @param nx,ny       Grid dimensions
 * @param T           Temperature field [nx × ny], row-major
 *                     Boundary values must be set before calling.
 * @param max_iter    Maximum iterations
 * @param tol         Convergence tolerance
 * @return Number of iterations, -1 if did not converge
 */
int gauss_seidel_2d_laplace(int nx, int ny, double T[],
                            int max_iter, double tol);

/**
 * Successive Over-Relaxation (SOR) for 2D Laplace.
 *
 * T_new = (1 - ω)·T_old + ω·(T_{i+1,j}+T_{i-1,j}+T_{i,j+1}+T_{i,j-1})/4
 *
 * @param omega  Relaxation factor (1 < ω < 2 for over-relaxation)
 */
int sor_2d_laplace(int nx, int ny, double T[], double omega,
                   int max_iter, double tol);

/* ---------------------------------------------------------------------------
 * L6: Pipe Network Solver (Hardy Cross Method)
 * ------------------------------------------------------------------------- */

/**
 * Hardy Cross method for pipe network analysis.
 *
 * Iteratively adjusts flow in each loop to satisfy:
 *   Σ_{loop} hL(Q) = 0
 *
 * Correction: ΔQ = -Σ hL / (n·Σ |hL/Q|)
 * where n = 2 for Darcy-Weisbach (turbulent).
 *
 * @param n_pipes       Number of pipes
 * @param n_loops       Number of independent loops
 * @param lengths       Pipe lengths [m] (size n_pipes)
 * @param diameters     Pipe diameters [m] (size n_pipes)
 * @param roughness     Pipe roughness [m] (size n_pipes)
 * @param loop_members  Row-major loop-pipe incidence matrix (+1/-1/0)
 *                       (size n_loops × n_pipes)
 * @param Q             Input: initial guess / Output: flow rates [m³/s] (size n_pipes)
 * @param max_iter      Maximum iterations
 * @param tol           Convergence tolerance [m³/s]
 * @return Number of iterations, -1 if did not converge
 */
int hardy_cross_solver(int n_pipes, int n_loops,
                       const double lengths[], const double diameters[],
                       const double roughness[],
                       const int loop_members[],
                       double Q[], int max_iter, double tol);

/* ---------------------------------------------------------------------------
 * L5: Numerical Integration for Boundary Layer Quantities
 * ------------------------------------------------------------------------- */

/**
 * Compute displacement thickness δ* from velocity profile data.
 *
 * δ* = ∫₀^δ (1 - u/U∞) dy
 *
 * Uses composite Simpson's rule.
 *
 * @param n        Number of data points (must be odd ≥ 3)
 * @param y        Distance from wall [m], size n (uniform spacing)
 * @param u        Velocity profile [m/s], size n
 * @param U_inf    Freestream velocity [m/s]
 * @return Displacement thickness [m]
 */
double integrate_displacement_thickness(int n, const double y[],
                                        const double u[], double U_inf);

/**
 * Compute momentum thickness θ from velocity profile data.
 *
 * θ = ∫₀^δ (u/U∞)·(1 - u/U∞) dy
 *
 * Uses composite Simpson's rule.
 */
double integrate_momentum_thickness(int n, const double y[],
                                    const double u[], double U_inf);

/**
 * Integrate using composite Simpson's rule.
 *
 * ∫ₐᵇ f(x) dx ≈ (h/3)·[f₀+f_n + 4(f₁+f₃+...) + 2(f₂+f₄+...)]
 *
 * @param n  Number of points (must be odd)
 * @param h  Uniform spacing
 * @param f  Function values, size n
 * @return Integral approximation
 */
double simpson_integral(int n, double h, const double f[]);

/**
 * 4th-order Runge-Kutta for ODE initial value problems.
 * dy/dx = f(x, y), y(x0) = y0
 *
 * Used for solving boundary layer ODEs (Falkner-Skan) numerically.
 *
 * @param x0     Initial x
 * @param y0     Initial y
 * @param h      Step size
 * @param n_steps Number of steps
 * @param f      Derivative function: f(x, y, params)
 * @param params User parameters passed to f
 * @param x_out  Output x values, size n_steps+1
 * @param y_out  Output y values, size n_steps+1
 */
typedef double (*ode_function)(double x, double y, void *params);

void rk4_solve(double x0, double y0, double h, int n_steps,
               ode_function f, void *params,
               double x_out[], double y_out[]);

/* ---------------------------------------------------------------------------
 * L8: SIMPLE Algorithm Pressure Correction (CFD Basics)
 * ------------------------------------------------------------------------- */

/**
 * Simplified pressure correction step (SIMPLE-like).
 *
 * Given current velocity field u*, compute pressure correction p'.
 *
 * p'_P = p'_P - ω·b_P/a_P  (under-relaxation applied)
 *
 * This is a 1D demonstration of the SIMPLE algorithm logic.
 * Full 2D/3D implementation requires staggered grids.
 *
 * @param n_cells       Number of 1D cells
 * @param u_star        Current velocity guess [m/s], size n_cells+1 (faces)
 * @param p             Current pressure [Pa], size n_cells (centers) — updated in place
 * @param rho           Density [kg/m³]
 * @param dx            Cell size [m]
 * @param omega_p       Pressure under-relaxation factor (0 < ω < 1)
 * @param max_iter      Max iterations
 * @param tol           Continuity residual tolerance
 * @return Number of iterations, -1 if diverged
 */
int simple_pressure_correction_1d(int n_cells, const double u_star[],
                                  double p[], double rho, double dx,
                                  double omega_p, int max_iter, double tol);

#endif /* NUMERICAL_FLUID_THERMAL_H */
