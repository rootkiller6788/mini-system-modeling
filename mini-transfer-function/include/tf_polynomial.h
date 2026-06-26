/**
 * @file tf_polynomial.h
 * @brief Polynomial Algebra for Transfer Function Analysis
 *
 * Knowledge Coverage:
 *   L1: Polynomial struct, coefficient representation, degree
 *   L2: Polynomial arithmetic (add, sub, mul, scale, power, derivative, integral)
 *   L3: Horner evaluation, complex evaluation, division algorithm, GCD, LCM
 *   L4: Root finding (Newton-Raphson with deflation, companion matrix method)
 *   L5: Polynomial composition, Chebyshev representation, numerical conditioning
 *
 * Polynomials form the algebraic foundation of transfer functions.
 * Every LTI transfer function G(s) = N(s)/D(s) is a ratio of two
 * polynomials in the complex variable s (or z for discrete-time).
 *
 * Coefficient ordering: coeff[0] + coeff[1]*s + ... + coeff[order]*s^order
 * i.e., ascending powers, constant term at index 0.
 *
 * Key theorems used:
 *   - Fundamental Theorem of Algebra: polynomial of degree n has n complex roots
 *   - Division Algorithm: unique q, r exist with deg(r) < deg(b)
 *   - Euclidean Algorithm: computes GCD of two polynomials
 *   - Companion matrix: eigenvalues = polynomial roots
 *   - Horner's Method: optimal polynomial evaluation in O(n)
 *
 * References:
 *   W.H. Press et al., "Numerical Recipes", 3rd ed., Ch. 5, 9
 *   N.J. Higham, "Accuracy and Stability of Numerical Algorithms", 2nd ed.
 *   G.H. Golub, C.F. Van Loan, "Matrix Computations", 4th ed., Ch. 7
 *   MIT 6.241 (Dynamic Systems) / Stanford ENGR207B / Berkeley ME132
 */
#ifndef TF_POLYNOMIAL_H
#define TF_POLYNOMIAL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Polynomial: coeff[0] + coeff[1]*s + coeff[2]*s^2 + ... + coeff[order]*s^order.
 *  Coefficients stored in ascending power order. */
typedef struct {
    double *coeff;   /**< Coefficient array [order+1] */
    int     order;   /**< Polynomial degree */
} Polynomial;

/* ================================================================
 * L1/L2: Lifecycle and Construction
 * ================================================================ */

/**
 * Create polynomial from coefficient array.
 * @param coeff  Array of [order+1] coefficients in ascending power order
 * @param order  Polynomial degree (>=0)
 * @return Allocated Polynomial, or NULL on error
 * Complexity: O(order) time and space
 */
Polynomial* poly_create(const double *coeff, int order);

/** Create zero polynomial P(s) = 0. Complexity: O(1) */
Polynomial* poly_zero(void);

/** Create unit polynomial P(s) = 1. Complexity: O(1) */
Polynomial* poly_one(void);

/** Create identity polynomial P(s) = s. Complexity: O(1) */
Polynomial* poly_s(void);

/** Create constant polynomial P(s) = c. Complexity: O(1) */
Polynomial* poly_constant(double c);

/** Create polynomial from roots: P(s) = gain * prod(s - root_i).
 *  @param roots   Array of n roots
 *  @param n       Number of roots
 *  @param gain    Leading coefficient multiplier
 *  @return Expanded polynomial
 *  Complexity: O(n^2)
 */
Polynomial* poly_from_roots(const double *roots, int n, double gain);

/**
 * Deep-copy polynomial.
 * Complexity: O(order)
 */
Polynomial* poly_clone(const Polynomial *p);

/**
 * Free polynomial and all memory. NULL-safe.
 * Complexity: O(1)
 */
void poly_destroy(Polynomial *p);

/* ================================================================
 * L2: Polynomial Arithmetic
 * ================================================================ */

/**
 * Addition: r(s) = a(s) + b(s).
 * Theorem: deg(r) <= max(deg(a), deg(b)).
 * Complexity: O(max(deg_a, deg_b))
 */
Polynomial* poly_add(const Polynomial *a, const Polynomial *b);

/**
 * Subtraction: r(s) = a(s) - b(s).
 * Complexity: O(max(deg_a, deg_b))
 */
Polynomial* poly_subtract(const Polynomial *a, const Polynomial *b);

/**
 * Multiplication (convolution): r(s) = a(s) * b(s).
 * Theorem: deg(r) = deg(a) + deg(b) (unless leading terms cancel).
 * Implements discrete convolution: r_k = sum_{i+j=k} a_i * b_j.
 * Complexity: O(deg_a * deg_b)
 */
Polynomial* poly_multiply(const Polynomial *a, const Polynomial *b);

/**
 * Scalar multiplication: r(s) = k * p(s).
 * Complexity: O(order)
 */
Polynomial* poly_scale(const Polynomial *p, double k);

/**
 * Power: r(s) = p(s)^n via repeated squaring.
 * @param n  Non-negative integer exponent
 * @return p^n, or NULL on error
 * Complexity: O(log(n) * n * order^2)
 *
 * Algorithm: Binary exponentiation reduces multiply operations from n-1
 * to O(log n) using: p^n = (p^2)^{n/2} if n even, p*(p^2)^{(n-1)/2} if n odd.
 */
Polynomial* poly_power(const Polynomial *p, int n);

/**
 * Derivative: r(s) = d/ds P(s).
 * Theorem: d/ds (sum a_k*s^k) = sum k*a_k*s^{k-1}.
 * Complexity: O(order)
 *
 * Time-domain correspondence: s*G(s) <-> d/dt g(t).
 */
Polynomial* poly_derivative(const Polynomial *p);

/**
 * Indefinite integral: r(s) = integral P(s) ds (constant term = 0).
 * Complexity: O(order)
 * Time-domain: G(s)/s <-> integral_0^t g(tau) dtau.
 */
Polynomial* poly_integral(const Polynomial *p);

/**
 * Polynomial composition: r(s) = a(b(s)).
 * Substitutes the polynomial b(s) for s in a(s).
 * Used for frequency shift (s -> s+a), bilinear transform, etc.
 * Complexity: O(deg_a^2 * deg_b) using Horner's method for composition.
 */
Polynomial* poly_compose(const Polynomial *a, const Polynomial *b);

/* ================================================================
 * L3: Evaluation Methods
 * ================================================================ */

/**
 * Evaluate polynomial at real x using Horner's method.
 * @return P(x)
 * Complexity: O(order), optimal numerical stability
 *
 * Horner's method (also called nested multiplication):
 *   p = a_n
 *   for i = n-1 downto 0: p = p*x + a_i
 * Requires n multiplications and n additions. Backward error is
 * bounded by 2n*eps/(1-2n*eps) * |P|(|x|) (Higham §5.1).
 */
double poly_eval(const Polynomial *p, double x);

/**
 * Evaluate polynomial at complex argument s = re + j*im.
 * @param ore  Output: real part of P(s)
 * @param oim  Output: imaginary part of P(s)
 * Complexity: O(order)
 *
 * Uses complex Horner's method: 4n multiplications, 3n additions.
 */
void poly_eval_complex(const Polynomial *p, double re, double im,
                        double *ore, double *oim);

/**
 * Evaluate P(x) and P'(x) simultaneously via Horner.
 * @param dval  Output: P'(x)
 * @return P(x)
 * Complexity: O(order), more efficient than two separate evaluations.
 *
 * Combined form: p = a_n, dp = 0; for i=n-1 downto 0:
 *   dp = dp*x + p; p = p*x + a_i
 * (Higham §5.2)
 */
double poly_eval_derivative(const Polynomial *p, double x, double *dval);

/**
 * Evaluate polynomial at matrix argument (for polynomial matrix functions).
 * @param X     Input square matrix [n*n], row-major
 * @param n     Matrix dimension
 * @param PX    Output P(X) matrix [n*n], row-major
 * @return 0 on success
 * Complexity: O(order * n^3) using Horner with matrix multiply
 *
 * Application: Compute exp(A), (sI-A)^{-1} via Cayley-Hamilton.
 */
int poly_eval_matrix(const Polynomial *p, const double *X, int n, double *PX);

/* ================================================================
 * L3: Polynomial Queries
 * ================================================================ */

/** True degree (highest index with |coeff| >= 1e-15).
 *  Returns -1 for zero polynomial. */
int poly_degree(const Polynomial *p);

/** Check if identically zero (within eps). */
int poly_is_zero(const Polynomial *p);

/** Check if constant (degree 0 or zero). */
int poly_is_constant(const Polynomial *p);

/** Check if monic (leading coefficient = 1 within eps). */
int poly_is_monic(const Polynomial *p);

/** Get leading coefficient. */
double poly_leading_coeff(const Polynomial *p);

/** Count sign changes in coefficient sequence (for Descartes' rule of signs). */
int poly_sign_changes(const Polynomial *p);

/* ================================================================
 * L3: Polynomial Division and GCD
 * ================================================================ */

/**
 * Polynomial long division: a(s) = q(s)*b(s) + r(s).
 * @param quotient   Output quotient (pass NULL to skip)
 * @param remainder  Output remainder (pass NULL to skip)
 * @return 0 on success, -1 if b is zero
 * Complexity: O((deg_a - deg_b + 1) * deg_b)
 *
 * Theorem (Division Algorithm for polynomials over a field):
 * For a, b in F[s] with b != 0, there exist unique q, r in F[s]
 * such that a = q*b + r and deg(r) < deg(b).
 */
int poly_divide(const Polynomial *a, const Polynomial *b,
                Polynomial **quotient, Polynomial **remainder);

/**
 * Polynomial GCD via Euclidean algorithm.
 * Complexity: O(deg_a * deg_b)
 *
 * Algorithm: While b != 0: (a, b) = (b, a mod b). Final a = GCD.
 * Uses pseudo-remainder to handle floating-point.
 *
 * Application: Detect pole-zero cancellation, compute coprime factorization.
 */
Polynomial* poly_gcd(const Polynomial *a, const Polynomial *b);

/**
 * Polynomial LCM: lcm(a,b) = (a*b) / gcd(a,b).
 * Complexity: O(deg_a * deg_b)
 */
Polynomial* poly_lcm(const Polynomial *a, const Polynomial *b);

/**
 * Extended Euclidean algorithm: compute u, v, g such that
 * a*u + b*v = g where g = gcd(a,b).
 * Bezout coefficients u, v enable Diophantine equation solutions.
 * Complexity: O(deg_a * deg_b)
 */
int poly_extended_gcd(const Polynomial *a, const Polynomial *b,
                      Polynomial **u, Polynomial **v, Polynomial **g);

/**
 * Polynomial remainder: r = a mod b (remainder of division).
 * Complexity: O((deg_a - deg_b + 1) * deg_b)
 */
Polynomial* poly_remainder(const Polynomial *a, const Polynomial *b);

/* ================================================================
 * L4: Root Finding
 * ================================================================ */

/**
 * Find real roots via multi-start Newton-Raphson with deflation.
 *
 * Algorithm:
 *   1. Try 20+ starting points evenly spaced in [-10, 10]
 *   2. For each, run Newton iteration up to 50 steps
 *   3. Keep best root found (smallest |P(root)|)
 *   4. Deflate: divide by (s - root)
 *   5. Repeat for remaining degree
 *
 * @param coeff  Coefficient array [order+1]
 * @param order  Polynomial degree
 * @param roots  Output real roots [maxr]
 * @param maxr   Maximum roots to find
 * @return Number of real roots found
 * Complexity: O(maxr * 20 * 50 * order)
 *
 * Numerical stability: Uses central-difference derivative for robustness.
 * Deflation uses forward-synthetic division.
 * Ref: Press et al., Numerical Recipes §9.5
 */
int poly_find_roots_real(const double *coeff, int order,
                          double *roots, int maxr);

/**
 * Find all roots (real and complex) via companion matrix eigenvalues.
 *
 * The companion matrix for monic polynomial s^n + a_{n-1}s^{n-1} + ... + a_0:
 *   C = [[0, 0, ..., 0, -a_0],
 *        [1, 0, ..., 0, -a_1],
 *        [0, 1, ..., 0, -a_2],
 *        ...
 *        [0, 0, ..., 1, -a_{n-1}]]
 *
 * det(sI - C) = s^n + a_{n-1}s^{n-1} + ... + a_0 = P(s)
 *
 * So eigenvalues of C = roots of P(s). QR algorithm gives all roots.
 *
 * @param roots_real  Output: real parts [order]
 * @param roots_imag  Output: imag parts [order]
 * @return Number of roots found (=order for well-conditioned polynomials)
 * Complexity: O(order^3) via Hessenberg QR
 * Ref: Golub & Van Loan §7.4.6
 */
int poly_find_roots_complex(const double *coeff, int order,
                             double *roots_real, double *roots_imag);

/* ================================================================
 * L5: Polynomial Utilities
 * ================================================================ */

/** Print polynomial in human-readable form with given variable name. */
void poly_print(const Polynomial *p, const char *var);

/**
 * Convert to monic form: divide all coefficients by leading coefficient.
 * @return Original leading coefficient (0 if zero polynomial).
 * Complexity: O(order)
 */
double poly_make_monic(Polynomial *p);

/**
 * Trim negligible coefficients (|coeff| < 1e-15) from both ends.
 * Complexity: O(order)
 */
void poly_trim(Polynomial *p);

/**
 * Check if two polynomials are equal within tolerance.
 * Complexity: O(max(deg_a, deg_b))
 */
int poly_equals(const Polynomial *a, const Polynomial *b, double tol);

/**
 * Compute polynomial norm (Euclidean / 2-norm of coefficient vector).
 * Complexity: O(order)
 */
double poly_norm(const Polynomial *p);

#ifdef __cplusplus
}
#endif
#endif /* TF_POLYNOMIAL_H */
