/**
 * sfg_complex.h — Complex Number Utilities for Signal Flow Graphs
 *
 * Signal flow graph gains are fundamentally complex-valued (transmittances
 * can represent magnitude and phase, or real and imaginary parts of
 * frequency-dependent gains). This module provides complex arithmetic
 * operations needed throughout the SFG library.
 *
 * L3 Math Structure: The complex field C with operations +, -, *, /,
 * magnitude, phase, conjugation.
 *
 * While C99 <complex.h> provides basic complex support, this module
 * adds convenience functions, printing, comparison, and polynomial
 * evaluation optimized for SFG computation patterns.
 *
 * References:
 *  - C99 Standard, ISO/IEC 9899:1999, Section 7.3 (Complex arithmetic)
 *  - Churchill, R.V., Brown, J.W. "Complex Variables and Applications"
 *    McGraw-Hill, 9th Ed., 2013.
 */

#ifndef SFG_COMPLEX_H
#define SFG_COMPLEX_H

#include <complex.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L3: Complex Number Creation / Access
 * ================================================================ */

/**
 * sfg_cmplx — Create a complex number from real and imaginary parts
 */
#define sfg_cmplx(re, im) ((re) + (im) * I)

/**
 * sfg_real_part — Extract real part
 */
#define sfg_real_part(z) (creal(z))

/**
 * sfg_imag_part — Extract imaginary part
 */
#define sfg_imag_part(z) (cimag(z))

/**
 * sfg_cmplx_from_polar — Create complex from magnitude and phase
 *
 * L3: Polar representation: z = r · e^(jθ)
 */
static inline double complex sfg_cmplx_from_polar(double r, double theta) {
    return r * (cos(theta) + I * sin(theta));
}

/* ================================================================
 * L3: Complex Arithmetic Helpers
 * ================================================================ */

/**
 * sfg_cmplx_magnitude — |z| = sqrt(re² + im²)
 */
static inline double sfg_cmplx_magnitude(double complex z) {
    return cabs(z);
}

/**
 * sfg_cmplx_phase — arg(z) = atan2(im, re) in radians
 */
static inline double sfg_cmplx_phase(double complex z) {
    return carg(z);
}

/**
 * sfg_cmplx_conjugate — z̄ = re - j·im
 */
static inline double complex sfg_cmplx_conjugate(double complex z) {
    return conj(z);
}

/**
 * sfg_cmplx_sqmag — |z|² = re² + im² (avoids sqrt, faster)
 */
static inline double sfg_cmplx_sqmag(double complex z) {
    double re = creal(z), im = cimag(z);
    return re * re + im * im;
}

/**
 * sfg_cmplx_reciprocal — 1/z = z̄ / |z|²
 *
 * Handles near-zero case.
 */
static inline double complex sfg_cmplx_reciprocal(double complex z) {
    double sq = sfg_cmplx_sqmag(z);
    if (sq < 1e-30) return 0.0;
    return conj(z) / sq;
}

/**
 * sfg_cmplx_add — z = a + b
 */
static inline double complex sfg_cmplx_add(double complex a, double complex b) {
    return a + b;
}

/**
 * sfg_cmplx_sub — z = a - b
 */
static inline double complex sfg_cmplx_sub(double complex a, double complex b) {
    return a - b;
}

/**
 * sfg_cmplx_mul — z = a * b
 */
static inline double complex sfg_cmplx_mul(double complex a, double complex b) {
    return a * b;
}

/**
 * sfg_cmplx_div — z = a / b (with protection against division by zero)
 */
static inline double complex sfg_cmplx_div(double complex a, double complex b) {
    double sq = sfg_cmplx_sqmag(b);
    if (sq < 1e-30) return 0.0;
    return (a * conj(b)) / sq;
}

/**
 * sfg_cmplx_scale — z = k * z (scale by real constant)
 */
static inline double complex sfg_cmplx_scale(double complex z, double k) {
    return k * z;
}

/**
 * sfg_cmplx_negate — z = -z
 */
static inline double complex sfg_cmplx_negate(double complex z) {
    return -z;
}

/* ================================================================
 * L3: Complex Comparison and Properties
 * ================================================================ */

/**
 * sfg_cmplx_eq — Fuzzy equality check (within epsilon)
 */
static inline int sfg_cmplx_eq(double complex a, double complex b, double eps) {
    return cabs(a - b) < eps;
}

/**
 * sfg_cmplx_is_zero — Check if complex number is approximately zero
 */
static inline int sfg_cmplx_is_zero(double complex z, double eps) {
    return cabs(z) < eps;
}

/**
 * sfg_cmplx_is_real — Check if imaginary part is negligible
 */
static inline int sfg_cmplx_is_real(double complex z, double eps) {
    return fabs(cimag(z)) < eps;
}

/**
 * sfg_cmplx_is_pure_imag — Check if real part is negligible
 */
static inline int sfg_cmplx_is_pure_imag(double complex z, double eps) {
    return fabs(creal(z)) < eps;
}

/* ================================================================
 * L3: Complex Polynomial Evaluation
 * ================================================================ */

/**
 * sfg_poly_eval — Evaluate polynomial with complex coefficients at z
 *
 * P(z) = c[0] + c[1]·z + c[2]·z² + ... + c[n]·z^n
 *
 * Uses Horner's method for numerical stability.
 *
 * Complexity: O(n)
 */
static inline double complex sfg_poly_eval(const double complex* coef,
                                            int n, double complex z) {
    if (n < 0) return 0.0;
    double complex result = coef[n];
    for (int i = n - 1; i >= 0; i--) {
        result = result * z + coef[i];
    }
    return result;
}

/**
 * sfg_poly_eval_real — Evaluate polynomial with real coefficients at real x
 *
 * Horner's method for real-valued polynomials.
 *
 * Complexity: O(n)
 */
static inline double sfg_poly_eval_real(const double* coef, int n, double x) {
    if (n < 0) return 0.0;
    double result = coef[n];
    for (int i = n - 1; i >= 0; i--) {
        result = result * x + coef[i];
    }
    return result;
}

/**
 * sfg_poly_derivative — Compute coefficients of polynomial derivative
 *
 * If P(x) = Σ c_i x^i, then P'(x) = Σ i·c_i x^{i-1}
 *
 * Complexity: O(n)
 */
void sfg_poly_derivative(const double* coef, int n, double* deriv);

/**
 * sfg_poly_multiply — Multiply two polynomials
 *
 * (Σ a_i x^i) · (Σ b_j x^j) = Σ c_k x^k
 *
 * Complexity: O(n·m)
 * Returns: degree of result (n + m)
 */
int sfg_poly_multiply(const double* a, int n,
                       const double* b, int m,
                       double* result);

/* ================================================================
 * L3: Complex Display
 * ================================================================ */

/**
 * sfg_cmplx_print — Print complex number as "re + j·im"
 *
 * @param z      Complex number to print
 * @param label  Optional label (NULL for none)
 */
void sfg_cmplx_print(double complex z, const char* label);

/**
 * sfg_cmplx_print_polar — Print as "r ∠ θ°" (magnitude / angle in degrees)
 */
void sfg_cmplx_print_polar(double complex z, const char* label);

/**
 * sfg_cmplx_format — Format complex number into string
 *
 * @param buf    Output buffer (must be ≥ 64 chars)
 * @param size   Buffer size
 * @param z      Complex number
 */
void sfg_cmplx_format(char* buf, int size, double complex z);

/* ================================================================
 * L3: Complex Matrix Utilities
 * ================================================================ */

/**
 * sfg_cmplx_matrix_print — Print an M×N complex matrix
 */
void sfg_cmplx_matrix_print(const double complex* mat, int M, int N,
                              const char* label);

/**
 * sfg_cmplx_matrix_mul_vec — Matrix-vector multiply: y = A x
 *
 * y[i] = Σ_j A[i*N + j] * x[j]
 *
 * Complexity: O(M·N)
 */
void sfg_cmplx_matrix_mul_vec(const double complex* A,
                               const double complex* x,
                               double complex* y, int M, int N);

/**
 * sfg_cmplx_vector_dot — Dot product of two complex vectors
 *
 * Σ x_i · y_i  (no conjugation — algebraic dot product)
 *
 * Complexity: O(N)
 */
double complex sfg_cmplx_vector_dot(const double complex* x,
                                      const double complex* y, int N);

#ifdef __cplusplus
}
#endif

#endif /* SFG_COMPLEX_H */
