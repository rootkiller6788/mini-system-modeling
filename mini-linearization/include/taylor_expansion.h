/**
 * @file taylor_expansion.h
 * @brief Taylor series expansion for multi-variable nonlinear systems.
 *
 * Knowledge Coverage:
 *   L1 - Taylor series, Maclaurin series, remainder forms (Lagrange, Cauchy)
 *   L2 - Truncation error, convergence radius, order of approximation
 *   L3 - Multi-index notation, factorial, partial derivative tensors
 *   L4 - Taylor's Theorem with remainder
 *
 * Reference:
 *   Rudin, Principles of Mathematical Analysis (1976), Theorem 5.15
 *   Nocedal & Wright, Numerical Optimization (2006), section 2.2
 */

#ifndef TAYLOR_EXPANSION_H
#define TAYLOR_EXPANSION_H

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/**
 * L1: Taylor expansion order specification.
 * f(a+h) = sum_{|alpha|<=k} (D^alpha f(a) / alpha!) * h^alpha + R_k(h)
 */
typedef struct {
    int order;
    int n_variables;
    double expansion_point[32];
    double *coefficients;
    int n_coefficients;
} TaylorExpansion;

/**
 * L3: Multi-index for partial derivatives.
 * alpha = (alpha_1, ..., alpha_n), |alpha| = sum alpha_i
 */
typedef struct {
    int *indices;
    int n_variables;
    int total_order;
} MultiIndex;

/**
 * L4: Taylor remainder forms.
 * Lagrange: R_k(h) = D^{k+1}f(c) * h^{k+1} / (k+1)!
 * Cauchy:   R_k(h) = D^{k+1}f(c) * (h-c)^k * h / k!
 * Integral: R_k(h) = (1/k!) * integral_0^1 (1-t)^k * D^{k+1}f(a+th) * h^{k+1} dt
 */
typedef enum {
    TAYLOR_REMAINDER_LAGRANGE,
    TAYLOR_REMAINDER_CAUCHY,
    TAYLOR_REMAINDER_INTEGRAL
} TaylorRemainderForm;

/* =====================================================================
 * L1-L2: Scalar Taylor Expansion (1D)
 * ===================================================================== */

/**
 * L1: Compute k-th order Taylor coefficient for 1D function.
 * f(a+h) = sum_{j=0}^k f^{(j)}(a) * h^j / j! + O(h^{k+1})
 *
 * Complexity: O(2^order)
 * Reference: Fornberg (1988), Mathematics of Computation
 */
double taylor_coefficient_1d(
    double (*f)(double x, void *user_data),
    double a, int order, double h, void *user_data);

/**
 * L1: Evaluate 1D Taylor polynomial using Horner's method.
 * P_k(x) = sum_{j=0}^k f^{(j)}(a) * (x-a)^j / j!
 *
 * Complexity: O(k)
 * Reference: Higham, Accuracy and Stability (2002), section 5.2
 */
double evaluate_taylor_polynomial_1d(
    const double *coeffs, double a, double x, int k);

/**
 * L4: Estimate Lagrange remainder for 1D Taylor expansion.
 * |R_k(x)| <= M_{k+1} * |x-a|^{k+1} / (k+1)!
 *
 * Complexity: O(k)
 * Reference: Apostol, Calculus Vol.I (1967), Theorem 7.7
 */
double taylor_remainder_lagrange_1d(
    double max_deriv, double x, double a, int k);

/* =====================================================================
 * L2-L3: Multi-Variable Taylor Expansion
 * ===================================================================== */

/**
 * L2: Count terms in k-th order Taylor expansion of n variables.
 * N(n,k) = C(n+k, k) = (n+k)! / (n! * k!)
 *
 * Complexity: O(min(n,k))
 */
int taylor_term_count(int n_variables, int order);

/**
 * L3: Generate all multi-indices with |alpha| <= max_order.
 * For n=2, k=2: (0,0), (1,0), (0,1), (2,0), (1,1), (0,2)
 *
 * Complexity: O(C(n+k,k) * n)
 */
MultiIndex *generate_multi_indices(int n_variables, int max_order, int *n_out);

/**
 * L3: Evaluate multi-variable monomial h^alpha.
 * h^alpha = h_1^{alpha_1} * h_2^{alpha_2} * ... * h_n^{alpha_n}
 *
 * Complexity: O(|alpha|)
 */
double evaluate_monomial(const double *h, const MultiIndex *mi);

/**
 * L1: Compute multivariate factorial alpha! = alpha_1! * ... * alpha_n!
 *
 * Complexity: O(|alpha| + n)
 */
double multi_index_factorial(const MultiIndex *mi);

/**
 * L2: Compute quadratic form h^T * H * h for Hessian-based expansion.
 *
 * Complexity: O(n^2)
 */
double quadratic_form(const double *H, const double *h, int n);

/**
 * L2: Evaluate full multi-variable Taylor polynomial.
 * P_k(a+h) = sum_{|alpha|<=k} (D^alpha f(a) / alpha!) * h^alpha
 *
 * Complexity: O(n_coeffs * max_order)
 */
double evaluate_taylor_polynomial_mv(
    const double *coeffs, const MultiIndex *mis, int n_coeffs,
    const double *a, const double *h, int n);

/**
 * L4: Estimate remainder bound for multivariate Taylor expansion.
 * ||R_k(h)|| <= (||h||^{k+1} / (k+1)!) * max_{|alpha|=k+1} ||D^alpha f||
 *
 * Complexity: O(k)
 * Reference: Nocedal & Wright (2006), Theorem 2.1
 */
double taylor_remainder_bound_mv(double max_deriv_norm, double h_norm, int k);

/**
 * L5: Approximate Hessian matrix using finite differences.
 * H_{ij} = [f(x+hi*ei+hj*ej) - f(x+hi*ei-hj*ej)
 *           - f(x-hi*ei+hj*ej) + f(x-hi*ei-hj*ej)] / (4*hi*hj)
 *
 * Complexity: O(n^2 * cost(f))
 * Reference: Dennis & Schnabel (1996), section 5.6
 */
void numerical_hessian(
    double (*f)(const double *x, int n, void *user_data),
    const double *x, int n, double *H, double h, void *user_data);

/**
 * L5: Compute Taylor coefficients for nonlinear ODE system up to order k.
 * For dx/dt = f(x):
 *   x^{(0)} = x(0),  x^{(1)} = f(x),  x^{(2)} = J(x)*f(x), ...
 *
 * Complexity: O(k * n^2 * cost(f))
 * Reference: Hairer, Norsett, Wanner, Solving ODE I (1993), section II.14
 */
void taylor_coefficients_ode(
    void (*f)(const double *x, double *dx, int n, void *user_data),
    void (*f_jac)(const double *x, double *J, int n, void *user_data),
    const double *x0, int n, int k,
    double *coeffs, double h, void *user_data);

/** Free TaylorExpansion. O(1). */
void free_taylor_expansion(TaylorExpansion *te);

/** Free MultiIndex array. O(n). */
void free_multi_indices(MultiIndex *mis, int n);


/* L1: Taylor series for elementary functions */
double taylor_series_sin(double x, int k);
double taylor_series_cos(double x, int k);
double taylor_series_exp(double x, int k);
double taylor_series_log1p(double x, int k);
double taylor_series_arctan(double x, int k);

#endif /* TAYLOR_EXPANSION_H */
