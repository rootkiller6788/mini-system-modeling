/**
 * @file jacobian_methods.h
 * @brief Jacobian matrix computation methods.
 */
#ifndef JACOBIAN_METHODS_H
#define JACOBIAN_METHODS_H
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>
#include "linearization_core.h"

typedef enum {
    FD_FORWARD_O1,
    FD_CENTRAL_O2,
    FD_FORWARD_O2,
    FD_CENTRAL_O4
} FiniteDifferenceOrder;

JacobianMatrix *compute_jacobian_fd(
    void (*f)(const double *x, const double *u, double *result, int n_funcs, void *data),
    const double *x, const double *u, int n_vars, int n_funcs,
    FiniteDifferenceOrder order, double eps, void *user_data);

double optimal_fd_step(int p, double f_scale);

JacobianMatrix *compute_jacobian_complex_step(
    void (*f_cmplx)(const double complex *x, const double complex *u,
                    double complex *result, int n_funcs, void *data),
    const double *x, const double *u, int n_vars, int n_funcs,
    double eps, void *user_data);

void jacobian_vector_product(
    void (*f)(const double *x, double *result, int m, void *data),
    const double *x, const double *v, int n, double eps,
    double *Jv, int m, void *user_data);

void compute_hessian_fd(
    double (*f_scalar)(const double *x, int n, void *data),
    const double *x, int n, double *H, double eps, void *user_data);

double matrix_condition_number(const double *A, int m, int n);

void compute_sparse_jacobian(
    void (*f)(const double *x, double *result, int m, void *data),
    const double *x, int n, int m,
    const int *sparsity_pattern, const int *col_ptr, int nnz,
    double eps, double *J, void *user_data);

double verify_jacobian_convergence(
    void (*f)(const double *x, double *result, int n_funcs, void *data),
    const JacobianMatrix *J, const double *x,
    int n_vars, int n_funcs, void *user_data);

 /* JACOBIAN_METHODS_H */

/* Additional Jacobian computation utilities and advanced methods */

/**
 * L5: Compute Jacobian using forward-mode automatic differentiation (AD).
 *
 * Uses dual numbers: J_{ij} = dual_part[f_i(x + eps*e_j)] where
 * each variable carries a tangent component.
 *
 * Complexity: O(n * cost(f_dual))
 * Reference: Griewank & Walther, Evaluating Derivatives (2008)
 */
JacobianMatrix *compute_jacobian_autodiff_forward(
    void (*f_dual)(const double *x_val, const double *x_dot,
                   double *y_val, double *y_dot, int n, void *data),
    const double *x, int n, void *user_data);

/**
 * L8: Use Krylov-subspace method to compute J*v products efficiently
 * for large-scale systems where forming the full Jacobian is prohibitive.
 *
 * Complexity: O(krylov_dim * cost(f))
 * Reference: Knoll & Keyes (2004), J. Comp. Physics
 */
void jacobian_krylov_action(
    void (*f)(const double *x, double *result, int m, void *data),
    const double *x, int n, int m,
    const double *V, int krylov_dim, double *JV,
    double eps, void *user_data);

/**
 * L8: Estimate Jacobian structure via sparse differencing with graph coloring.
 *
 * Partition columns of J into structurally orthogonal groups.
 * Each group can be evaluated with a single function evaluation.
 *
 * Complexity: O(n_colors * cost(f))
 * Reference: Coleman & More (1984), SIAM J. Numer. Anal.
 */
int jacobian_coloring_groups(
    const int *sparsity_pattern, const int *col_ptr, int n,
    int *color_assignment);

/**
 * L5: Compute Gauss-Newton Jacobian for least-squares problems.
 *
 * For r(x) = [r_1(x), ..., r_m(x)]^T with objective f = ||r(x)||^2,
 * the Gauss-Newton approximation to the Hessian is J_r^T * J_r.
 *
 * Complexity: O(m * n^2)
 * Reference: Nocedal & Wright (2006), section 10.2
 */
void gauss_newton_hessian_approx(const JacobianMatrix *J,
    double *H_approx, int n);

/**
 * L8: Compute Jacobian via ensemble method for Monte Carlo sensitivity.
 *
 * Uses simultaneous perturbations for all parameters.
 *
 * Complexity: O(n_ensemble * cost(f))
 * Reference: Spall (2003), Introduction to Stochastic Search
 */
JacobianMatrix *ensemble_jacobian(
    void (*f)(const double *x, double *result, int m, void *data),
    const double *x, int n, int m, int n_ensemble,
    double eps, void *user_data);

/**
 * L3: Compute singular values of Jacobian for sensitivity analysis.
 *
 * sigma_i = sqrt(eigenvalue_i(J^T * J))
 * Condition number = sigma_max / sigma_min
 *
 * Complexity: O(min(m,n)^3)
 * Reference: Golub & Van Loan (2013), section 2.5
 */
void jacobian_singular_values(const JacobianMatrix *J,
    double *singular_values, int *n_sv);




#endif /* JACOBIAN_METHODS_H */
