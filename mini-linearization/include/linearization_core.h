/**
 * @file linearization_core.h
 * @brief Core linearization theory
 */
#ifndef LINEARIZATION_CORE_H
#define LINEARIZATION_CORE_H
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>

typedef struct {
    double *state;
    double *input;
    int n_states;
    int n_inputs;
    bool is_hyperbolic;
} EquilibriumPoint;

typedef struct {
    double *data;
    int n_rows;
    int n_cols;
    double condition;
} JacobianMatrix;

typedef struct {
    double *A, *B, *C, *D;
    int n_states, n_inputs, n_outputs;
    EquilibriumPoint eq;
} LinearizedSystem;

typedef struct {
    double *real_parts, *imag_parts, *magnitudes;
    int n_eigenvalues;
    bool is_stable;
    double max_real_part;
} LinearizationSpectrum;

typedef enum {
    LINEARIZE_ANALYTICAL,
    LINEARIZE_NUMERICAL_FORWARD,
    LINEARIZE_NUMERICAL_CENTRAL,
    LINEARIZE_COMPLEX_STEP
} LinearizationMethod;

typedef enum {
    ODE_EULER_FORWARD,
    ODE_RK4,
    ODE_MATRIX_EXPONENTIAL
} LinearizedODEMethod;

typedef struct {
    bool is_conclusive, is_stable, is_unstable;
    int n_unstable_modes, n_center_modes;
} LyapunovIndirectResult;

/* Core API */
EquilibriumPoint *create_equilibrium(int n_states, int n_inputs);
void free_equilibrium(EquilibriumPoint *eq);
void set_equilibrium_values(EquilibriumPoint *eq, const double *x_eq, const double *u_eq);
bool check_equilibrium_hyperbolicity(const EquilibriumPoint *eq, const JacobianMatrix *A);

JacobianMatrix *create_jacobian(int n_rows, int n_cols);
void free_jacobian(JacobianMatrix *J);
void jacobian_set(JacobianMatrix *J, int i, int j, double value);
double jacobian_get(const JacobianMatrix *J, int i, int j);

JacobianMatrix *numerical_jacobian_central(
    void (*f)(const double *x, const double *u, double *result, int n_rows, void *user_data),
    const double *x, const double *u, int n_rows, int n_cols, double epsilon, void *user_data);

JacobianMatrix *numerical_jacobian_complex_step(
    void (*f_cmplx)(const double _Complex *x, const double _Complex *u,
                    double _Complex *result, int n_rows, void *user_data),
    const double *x, const double *u, int n_rows, int n_cols, double epsilon, void *user_data);

LinearizedSystem *linearize_system(
    void (*f_dyn)(const double *x, const double *u, double *dx, int n, void *data),
    void (*h_out)(const double *x, const double *u, double *y, int p, void *data),
    const double *x_eq, const double *u_eq,
    int n_states, int n_inputs, int n_outputs,
    LinearizationMethod method, double epsilon, void *user_data);

void free_linearized_system(LinearizedSystem *sys);
LyapunovIndirectResult *apply_lyapunov_indirect(const LinearizedSystem *sys);
void free_lyapunov_indirect_result(LyapunovIndirectResult *result);

void simulate_linearized(const LinearizedSystem *sys, const double *delta_x0,
    void (*u_func)(double t, double *delta_u, int m, void *data),
    double t_start, double t_end, int n_steps, LinearizedODEMethod method,
    int n_output, double *t_out, double *y_out, void *user_data);

bool validate_linearization(
    void (*f_dyn)(const double *x, const double *u, double *dx, int n, void *data),
    const LinearizedSystem *sys, const double *delta_x, double tolerance, void *user_data);

double *compute_controllability_matrix(const LinearizedSystem *sys, int *rank);
double *compute_observability_matrix(const LinearizedSystem *sys, int *rank);
double estimate_linearization_domain(
    void (*f_dyn)(const double *x, const double *u, double *dx, int n, void *data),
    const LinearizedSystem *sys, double tolerance, void *user_data);
int compute_eigenvalues_qr(const double *A, int n, double *real_parts, double *imag_parts);

 /* LINEARIZATION_CORE_H */

/**
 * L5: Find equilibrium point via Newton-Raphson iteration.
 *
 * Solve f(x, u) = 0 for x given fixed u.
 * Iteration: x_{k+1} = x_k - J^{-1}(x_k) * f(x_k, u)
 *
 * @param f_dyn      Dynamics function
 * @param u_eq       Fixed input
 * @param x0         Initial guess (length n), overwritten with solution
 * @param n          State dimension
 * @param max_iter   Maximum iterations
 * @param tol        Convergence tolerance
 * @param user_data  Arbitrary data
 * @return           Number of iterations used (0 = failure)
 *
 * Complexity: O(max_iter * n^3)
 * Reference: Nocedal & Wright (2006), section 2.4
 */
int find_equilibrium_newton(
    void (*f_dyn)(const double *x, const double *u, double *dx, int n, void *data),
    const double *u_eq, double *x0, int n, int max_iter, double tol, void *user_data);

/**
 * L8: Continuation method for finding equilibrium along parameter path.
 *
 * When an equilibrium at one parameter value is known, traces the
 * equilibrium manifold as parameters vary using predictor-corrector steps.
 *
 * @param f_dyn       Dynamics f(x, p) where p is the continuation parameter
 * @param x           Known equilibrium at p0 (length n), overwritten at p1
 * @param p0          Starting parameter value
 * @param p1          Target parameter value
 * @param n           State dimension
 * @param n_steps     Number of continuation steps
 * @param user_data   Arbitrary data
 *
 * Complexity: O(n_steps * n^3)
 * Reference: Allgower & Georg, Numerical Continuation Methods (1990)
 */
void continuation_equilibrium(
    void (*f_dyn)(const double *x, double p, double *dx, int n, void *data),
    double *x, double p0, double p1, int n, int n_steps, void *user_data);

/**
 * L5: Linearize system using analytical Jacobian (user-provided).
 *
 * When the Jacobian is known analytically, this function assembles
 * the linearized system directly without numerical differentiation.
 *
 * @param f_dyn       Dynamics f(x,u)
 * @param J_f         Analytical Jacobian df/d[x;u] (n x (n+m), row-major)
 * @param J_h         Analytical Jacobian dh/d[x;u] (p x (n+m), row-major)
 * @param x_eq, u_eq  Equilibrium point
 * @param n_states, n_inputs, n_outputs
 * @param sys         Output linearized system (caller allocates)
 *
 * Complexity: O(n * (n+m))
 * Reference: Ogata (2010), section 12.2
 */
void linearize_analytical(
    void (*f_dyn)(const double *x, const double *u, double *dx, int n, void *data),
    const double *J_f, const double *J_h,
    const double *x_eq, const double *u_eq,
    int n_states, int n_inputs, int n_outputs,
    LinearizedSystem *sys, void *user_data);

/**
 * L5: Compute linearized system's transfer function matrix.
 *
 * G(s) = C * (sI - A)^{-1} * B + D
 *
 * For SISO: single transfer function.
 * For MIMO: transfer function matrix G_{ij}(s).
 *
 * @param sys          Linearized system
 * @param s            Complex frequency variable
 * @param G            Output transfer matrix (n_outputs x n_inputs, row-major)
 *
 * Complexity: O(n^3) for matrix inversion
 * Reference: Ogata (2010), section 9.3
 */
void linearized_transfer_function(
    const LinearizedSystem *sys, double complex s, double complex *G);

/**
 * L3: Compute Gramian matrices for the linearized system.
 *
 * Controllability Gramian: W_c = integral_0^inf e^{At} B B^T e^{A^T t} dt
 * Observability Gramian:  W_o = integral_0^inf e^{A^T t} C^T C e^{At} dt
 *
 * Solve Lyapunov equations: A*W_c + W_c*A^T + B*B^T = 0
 *                           A^T*W_o + W_o*A + C^T*C = 0
 *
 * @param sys          Stabilizable/detectable linearized system
 * @param W_c          Output controllability Gramian (n x n, caller allocates)
 * @param W_o          Output observability Gramian (n x n, caller allocates)
 *
 * Complexity: O(n^3)
 * Reference: Dullerud & Paganini, A Course in Robust Control (2000)
 */
void compute_gramians(
    const LinearizedSystem *sys, double *W_c, double *W_o);

/**
 * L8: Balanced truncation model reduction of linearized system.
 *
 * Uses Hankel singular values from Gramians to produce a reduced-order
 * model that preserves input-output behavior.
 *
 * @param sys               Full-order linearized system
 * @param reduced_order     Desired order k (< n)
 * @param reduced_sys       Output reduced system (caller allocates)
 *
 * Complexity: O(n^3)
 * Reference: Moore (1981), IEEE TAC; Antoulas (2005)
 */
void balanced_truncation(
    const LinearizedSystem *sys, int reduced_order,
    LinearizedSystem *reduced_sys);

/**
 * L3: Compute matrix exponential e^{A*t} using Pade approximation.
 *
 * e^{A*t} approx N_{pq}(A) * D_{pq}^{-1}(A)
 *
 * Uses [8/8] Pade approximant with scaling-and-squaring method.
 *
 * @param A           Square matrix (n x n, row-major)
 * @param n           Dimension
 * @param t           Time
 * @param expAt       Output e^{A*t} (n x n, row-major, caller allocates)
 *
 * Complexity: O(n^3)
 * Reference: Higham (2005), SIAM J. Matrix Anal.
 */
void matrix_exponential_pade(const double *A, int n, double t, double *expAt);




#endif /* LINEARIZATION_CORE_H */
