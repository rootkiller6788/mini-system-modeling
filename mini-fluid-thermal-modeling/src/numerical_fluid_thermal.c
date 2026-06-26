/**
 * mini-fluid-thermal-modeling — Numerical Methods Implementation
 *
 * Finite difference solvers (1D steady/transient conduction),
 * Thomas algorithm for tridiagonal systems, Gauss-Seidel/SOR
 * for 2D Laplace, Hardy Cross pipe network solver,
 * numerical integration for boundary layers, RK4 ODE solver,
 * and SIMPLE pressure correction (1D demonstration).
 *
 * Knowledge Levels: L5-L6, L8
 */

#include "numerical_fluid_thermal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =====================================================================
 * L5: 1D Steady Conduction — Finite Difference
 * ===================================================================== */

int solve_1d_steady_conduction(int n, double L, double k, double q_dot,
                               double T_left, double T_right,
                               double T[], int max_iter, double tol)
{
    if (!T || n < 3 || L <= 0.0 || k <= 0.0) return -1;

    double dx = L / (n - 1);
    double dx2 = dx * dx;
    double source = q_dot * dx2 / k;

    /* Set boundary conditions */
    T[0] = T_left;
    T[n - 1] = T_right;

    /* Gauss-Seidel iteration */
    for (int iter = 0; iter < max_iter; iter++) {
        double max_change = 0.0;
        for (int i = 1; i < n - 1; i++) {
            double T_old = T[i];
            T[i] = 0.5 * (T[i-1] + T[i+1] + source);
            double change = fabs(T[i] - T_old);
            if (change > max_change) max_change = change;
        }
        if (max_change < tol) return iter + 1;
    }
    return -1; /* did not converge */
}

/* =====================================================================
 * L5: 1D Transient Conduction
 * ===================================================================== */

int solve_1d_transient_conduction_explicit(int nx, int nt, double L,
                                           double alpha, double dt,
                                           double T_left, double T_right,
                                           const double T0[], double T[])
{
    if (!T0 || !T || nx < 3 || nt < 1 || L <= 0.0
        || alpha <= 0.0 || dt <= 0.0) return -1;

    double dx = L / (nx - 1);
    double Fo = alpha * dt / (dx * dx);

    /* Stability check */
    if (Fo > 0.5) return -1;

    /* Internal work arrays */
    double *T_curr = (double *)malloc(nx * sizeof(double));
    double *T_next = (double *)malloc(nx * sizeof(double));
    if (!T_curr || !T_next) { free(T_curr); free(T_next); return -1; }

    memcpy(T_curr, T0, nx * sizeof(double));
    T_curr[0] = T_left;
    T_curr[nx - 1] = T_right;

    for (int step = 0; step < nt; step++) {
        T_next[0] = T_left;
        T_next[nx - 1] = T_right;
        for (int i = 1; i < nx - 1; i++) {
            T_next[i] = Fo * (T_curr[i-1] + T_curr[i+1])
                        + (1.0 - 2.0 * Fo) * T_curr[i];
        }
        double *tmp = T_curr; T_curr = T_next; T_next = tmp;
    }

    memcpy(T, T_curr, nx * sizeof(double));
    free(T_curr);
    free(T_next);
    return 0;
}

int solve_1d_transient_conduction_implicit(int nx, int nt, double L,
                                           double alpha, double dt,
                                           double T_left, double T_right,
                                           const double T0[], double T[])
{
    if (!T0 || !T || nx < 3 || nt < 1 || L <= 0.0
        || alpha <= 0.0 || dt <= 0.0) return -1;

    double dx = L / (nx - 1);
    double Fo = alpha * dt / (dx * dx);

    /* Build tridiagonal system for Crank-Nicolson:
     * -Fo/2·T^{n+1}_{i-1} + (1+Fo)·T^{n+1}_i - Fo/2·T^{n+1}_{i+1}
     *   = Fo/2·T^n_{i-1} + (1-Fo)·T^n_i + Fo/2·T^n_{i+1}
     */
    int n_int = nx - 2; /* interior points */
    double *a = (double *)malloc(n_int * sizeof(double));
    double *b = (double *)malloc(n_int * sizeof(double));
    double *c = (double *)malloc(n_int * sizeof(double));
    double *d = (double *)malloc(n_int * sizeof(double));
    double *x = (double *)malloc(n_int * sizeof(double));
    double *T_curr = (double *)malloc(nx * sizeof(double));
    double *T_out = (double *)malloc(nx * sizeof(double));

    if (!a || !b || !c || !d || !x || !T_curr || !T_out) {
        free(a); free(b); free(c); free(d); free(x);
        free(T_curr); free(T_out); return -1;
    }

    memcpy(T_curr, T0, nx * sizeof(double));
    T_curr[0] = T_left;
    T_curr[nx - 1] = T_right;

    /* Coefficients for tridiagonal (same for all steps) */
    for (int i = 0; i < n_int; i++) {
        a[i] = -Fo / 2.0;
        b[i] = 1.0 + Fo;
        c[i] = -Fo / 2.0;
    }
    /* Adjust boundaries: a[0] and c[n_int-1] already correct */

    for (int step = 0; step < nt; step++) {
        /* Build RHS */
        for (int i = 1; i < nx - 1; i++) {
            int j = i - 1;
            d[j] = (Fo / 2.0) * T_curr[i-1]
                 + (1.0 - Fo) * T_curr[i]
                 + (Fo / 2.0) * T_curr[i+1];
        }
        /* Boundary contributions */
        d[0] += (Fo / 2.0) * T_left;
        d[n_int - 1] += (Fo / 2.0) * T_right;

        if (thomas_algorithm(n_int, a, b, c, d, x) != 0) {
            free(a); free(b); free(c); free(d); free(x);
            free(T_curr); free(T_out); return -1;
        }

        /* Update T_curr */
        for (int i = 1; i < nx - 1; i++) {
            T_curr[i] = x[i - 1];
        }
    }

    memcpy(T, T_curr, nx * sizeof(double));
    free(a); free(b); free(c); free(d); free(x);
    free(T_curr); free(T_out);
    return 0;
}

/* =====================================================================
 * L5: Thomas Algorithm (Tridiagonal Solver)
 * ===================================================================== */

int thomas_algorithm(int n, const double a[], const double b[],
                     const double c[], const double d[], double x[])
{
    if (!a || !b || !c || !d || !x || n <= 0) return -1;

    /* Allocate scratch arrays */
    double *c_prime = (double *)malloc(n * sizeof(double));
    double *d_prime = (double *)malloc(n * sizeof(double));
    if (!c_prime || !d_prime) {
        free(c_prime); free(d_prime); return -1;
    }

    if (fabs(b[0]) < 1e-15) {
        free(c_prime); free(d_prime); return -1; /* singular */
    }

    /* Forward sweep */
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double denom = b[i] - a[i] * c_prime[i - 1];
        if (fabs(denom) < 1e-15) {
            free(c_prime); free(d_prime); return -1;
        }
        if (i < n - 1) {
            c_prime[i] = c[i] / denom;
        }
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom;
    }

    /* Back substitution */
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    free(c_prime);
    free(d_prime);
    return 0;
}

/* =====================================================================
 * L5: Gauss-Seidel Relaxation
 * ===================================================================== */

int gauss_seidel_2d_laplace(int nx, int ny, double T[],
                            int max_iter, double tol)
{
    if (!T || nx < 3 || ny < 3) return -1;

    for (int iter = 0; iter < max_iter; iter++) {
        double max_change = 0.0;
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                int idx = j * nx + i;
                double T_old = T[idx];
                T[idx] = 0.25 * (T[j * nx + (i + 1)]
                               + T[j * nx + (i - 1)]
                               + T[(j + 1) * nx + i]
                               + T[(j - 1) * nx + i]);
                double change = fabs(T[idx] - T_old);
                if (change > max_change) max_change = change;
            }
        }
        if (max_change < tol) return iter + 1;
    }
    return -1;
}

int sor_2d_laplace(int nx, int ny, double T[], double omega,
                   int max_iter, double tol)
{
    if (!T || nx < 3 || ny < 3 || omega <= 0.0 || omega >= 2.0)
        return -1;

    for (int iter = 0; iter < max_iter; iter++) {
        double max_change = 0.0;
        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                int idx = j * nx + i;
                double T_old = T[idx];
                double T_gs = 0.25 * (T[j * nx + (i + 1)]
                                    + T[j * nx + (i - 1)]
                                    + T[(j + 1) * nx + i]
                                    + T[(j - 1) * nx + i]);
                T[idx] = (1.0 - omega) * T_old + omega * T_gs;
                double change = fabs(T[idx] - T_old);
                if (change > max_change) max_change = change;
            }
        }
        if (max_change < tol) return iter + 1;
    }
    return -1;
}

/* =====================================================================
 * L6: Pipe Network Solver (Hardy Cross Method)
 * ===================================================================== */

int hardy_cross_solver(int n_pipes, int n_loops,
                       const double lengths[], const double diameters[],
                       const double roughness[],
                       const int loop_members[],
                       double Q[], int max_iter, double tol)
{
    if (!lengths || !diameters || !roughness || !loop_members || !Q)
        return -1;
    if (n_pipes <= 0 || n_loops <= 0) return -1;
    if (max_iter <= 0) max_iter = 100;

    double g = 9.81; /* gravity [m/s²] */
    double nu = 1.0e-6; /* water kinematic viscosity [m²/s] (for Re calc) */

    for (int iter = 0; iter < max_iter; iter++) {
        double max_correction = 0.0;
        int converged = 1;

        for (int loop = 0; loop < n_loops; loop++) {
            double sum_hL = 0.0;
            double sum_hL_div_Q = 0.0;

            for (int p = 0; p < n_pipes; p++) {
                int sign = loop_members[loop * n_pipes + p];
                if (sign == 0) continue;

                double D = diameters[p];
                if (D <= 0.0) continue;
                double A = M_PI * D * D / 4.0;
                double u = fabs(Q[p]) / A;
                double Re = u * D / nu;
                double eps_D = roughness[p] / D;

                double f;
                if (Re < 2300.0) {
                    f = 64.0 / Re;
                } else {
                    /* Swamee-Jain approximation */
                    f = 0.25 / pow(log10(eps_D / 3.7 + 5.74 / pow(Re, 0.9)), 2.0);
                }

                double hL = f * (lengths[p] / D) * (u * u / (2.0 * g));
                /* Head loss follows sign convention */
                sum_hL += sign * (Q[p] > 0 ? hL : -hL);
                sum_hL_div_Q += fabs(hL / (Q[p] + 1e-15));
            }

            if (fabs(sum_hL_div_Q) < 1e-15) continue;

            double dQ = -sum_hL / (2.0 * sum_hL_div_Q);

            if (fabs(dQ) > max_correction) {
                max_correction = fabs(dQ);
            }
            if (fabs(dQ) > tol) converged = 0;

            /* Apply correction to all pipes in loop */
            for (int p = 0; p < n_pipes; p++) {
                int sign = loop_members[loop * n_pipes + p];
                if (sign != 0) {
                    Q[p] += sign * dQ;
                }
            }
        }

        if (converged) return iter + 1;
        if (max_correction < tol) return iter + 1;
    }
    return -1;
}

/* =====================================================================
 * L5: Numerical Integration for Boundary Layer Quantities
 * ===================================================================== */

double simpson_integral(int n, double h, const double f[])
{
    if (!f || n < 3 || (n % 2) == 0) return NAN; /* n must be odd */
    double sum = f[0] + f[n - 1];
    for (int i = 1; i < n - 1; i++) {
        if (i % 2 == 1) sum += 4.0 * f[i];  /* odd indices */
        else sum += 2.0 * f[i];              /* even indices */
    }
    return h * sum / 3.0;
}

double integrate_displacement_thickness(int n, const double y[],
                                        const double u[], double U_inf)
{
    if (!y || !u || n < 3 || (n % 2) == 0 || U_inf == 0.0) return NAN;
    double h = y[1] - y[0];
    /* Check uniform spacing */
    for (int i = 1; i < n; i++) {
        if (fabs((y[i] - y[i-1]) - h) > 1e-10 * h) return NAN;
    }

    double *integrand = (double *)malloc(n * sizeof(double));
    if (!integrand) return NAN;
    for (int i = 0; i < n; i++) {
        integrand[i] = 1.0 - u[i] / U_inf;
    }
    double result = simpson_integral(n, h, integrand);
    free(integrand);
    return result;
}

double integrate_momentum_thickness(int n, const double y[],
                                    const double u[], double U_inf)
{
    if (!y || !u || n < 3 || (n % 2) == 0 || U_inf == 0.0) return NAN;
    double h = y[1] - y[0];
    for (int i = 1; i < n; i++) {
        if (fabs((y[i] - y[i-1]) - h) > 1e-10 * h) return NAN;
    }

    double *integrand = (double *)malloc(n * sizeof(double));
    if (!integrand) return NAN;
    for (int i = 0; i < n; i++) {
        double ratio = u[i] / U_inf;
        integrand[i] = ratio * (1.0 - ratio);
    }
    double result = simpson_integral(n, h, integrand);
    free(integrand);
    return result;
}

/* =====================================================================
 * L5: 4th-Order Runge-Kutta ODE Solver
 * ===================================================================== */

void rk4_solve(double x0, double y0, double h, int n_steps,
               ode_function f, void *params,
               double x_out[], double y_out[])
{
    if (!f || !x_out || !y_out || n_steps <= 0) return;

    double x = x0;
    double y = y0;
    x_out[0] = x;
    y_out[0] = y;

    for (int i = 0; i < n_steps; i++) {
        double k1 = f(x, y, params);
        double k2 = f(x + h / 2.0, y + h * k1 / 2.0, params);
        double k3 = f(x + h / 2.0, y + h * k2 / 2.0, params);
        double k4 = f(x + h, y + h * k3, params);

        y = y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        x = x + h;

        x_out[i + 1] = x;
        y_out[i + 1] = y;
    }
}

/* =====================================================================
 * L8: SIMPLE Algorithm Pressure Correction (CFD Basics)
 * ===================================================================== */

int simple_pressure_correction_1d(int n_cells, const double u_star[],
                                  double p[], double rho, double dx,
                                  double omega_p, int max_iter, double tol)
{
    if (!u_star || !p || n_cells < 2 || rho <= 0.0 || dx <= 0.0)
        return -1;
    if (omega_p <= 0.0 || omega_p > 1.0) omega_p = 0.3;

    for (int iter = 0; iter < max_iter; iter++) {
        double max_residual = 0.0;
        double *p_new = (double *)malloc(n_cells * sizeof(double));
        if (!p_new) return -1;

        /* Copy current pressure */
        for (int i = 0; i < n_cells; i++) p_new[i] = p[i];

        /* Pressure correction based on continuity:
         * For 1D staggered grid:
         * u_e = u_star_e + (p_P - p_E) * dx / (rho * dx)  simplified
         *
         * The pressure correction equation:
         * a_P·p'_P = a_W·p'_W + a_E·p'_E + b
         *
         * Simplified: p'_P = (p'_W + p'_E)/2 + (rho·dx/(2dt))*(u_w - u_e)
         *
         * For steady flow: b = ρ·(u_w* - u_e*)
         *
         * Using Jacobi-like iteration with under-relaxation.
         */
        for (int i = 1; i < n_cells - 1; i++) {
            /* Mass imbalance in cell i */
            double u_w = u_star[i];     /* velocity at west face */
            double u_e = u_star[i + 1]; /* velocity at east face */
            double mass_imbalance = rho * (u_w - u_e);

            /* Pressure correction */
            double p_correction = 0.5 * (p[i-1] + p[i+1])
                                 + mass_imbalance * dx / (2.0 * rho)
                                 - p[i];

            p_new[i] = p[i] + omega_p * p_correction;

            double residual = fabs(mass_imbalance);
            if (residual > max_residual) max_residual = residual;
        }

        /* Boundary cells unchanged */
        p_new[0] = p[0];
        p_new[n_cells - 1] = p[n_cells - 1];

        /* Update pressure */
        for (int i = 0; i < n_cells; i++) p[i] = p_new[i];
        free(p_new);

        if (max_residual < tol) return iter + 1;
    }
    return -1;
}
