/**
 * sfg_complex.c — Complex Number Utilities Implementation
 *
 * Provides non-inline complex number functions: polynomial operations,
 * printing, formatting, and matrix-vector operations.
 *
 * Knowledge coverage:
 *   L3: Complex polynomial algebra, matrix operations
 */

#include "sfg_complex.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * L3: Polynomial Operations
 * ================================================================ */

void sfg_poly_derivative(const double* coef, int n, double* deriv) {
    /* deriv has degree n-1:
     * P(x) = c_0 + c_1 x + c_2 x^2 + ... + c_n x^n
     * P'(x) = c_1 + 2 c_2 x + 3 c_3 x^2 + ... + n c_n x^{n-1}
     */
    if (n <= 0) {
        deriv[0] = 0.0;
        return;
    }
    for (int i = 0; i < n; i++) {
        deriv[i] = (i + 1) * coef[i + 1];
    }
}

int sfg_poly_multiply(const double* a, int n,
                       const double* b, int m,
                       double* result) {
    /*
     * Convolution: (a_0 + ... + a_n x^n) * (b_0 + ... + b_m x^m)
     *            = c_0 + ... + c_{n+m} x^{n+m}
     * where c_k = Σ_{i+j=k} a_i * b_j
     */
    if (n < 0 || m < 0) return -1;

    int deg = n + m;
    for (int k = 0; k <= deg; k++) {
        result[k] = 0.0;
        for (int i = 0; i <= n; i++) {
            int j = k - i;
            if (j >= 0 && j <= m) {
                result[k] += a[i] * b[j];
            }
        }
    }
    return deg;
}

/* ================================================================
 * L3: Complex Display Functions
 * ================================================================ */

void sfg_cmplx_print(double complex z, const char* label) {
    double re = creal(z);
    double im = cimag(z);
    if (label && label[0]) {
        printf("%s = ", label);
    }
    if (fabs(im) < 1e-12) {
        printf("%.6g", re);
    } else if (fabs(re) < 1e-12) {
        printf("j%.6g", im);
    } else if (im > 0) {
        printf("%.6g + j%.6g", re, im);
    } else {
        printf("%.6g - j%.6g", re, -im);
    }
    printf("\n");
}

void sfg_cmplx_print_polar(double complex z, const char* label) {
    double r   = cabs(z);
    double ang = carg(z) * 180.0 / 3.14159265358979323846;
    if (label && label[0]) {
        printf("%s = ", label);
    }
    printf("%.6g ∠ %.2f°\n", r, ang);
}

void sfg_cmplx_format(char* buf, int size, double complex z) {
    if (!buf || size < 1) return;
    double re = creal(z);
    double im = cimag(z);
    if (fabs(im) < 1e-12) {
        snprintf(buf, size, "%.6g", re);
    } else if (fabs(re) < 1e-12) {
        snprintf(buf, size, "j%.6g", im);
    } else if (im > 0) {
        snprintf(buf, size, "%.6g + j%.6g", re, im);
    } else {
        snprintf(buf, size, "%.6g - j%.6g", re, -im);
    }
}

/* ================================================================
 * L3: Complex Matrix Utilities
 * ================================================================ */

void sfg_cmplx_matrix_print(const double complex* mat, int M, int N,
                              const char* label) {
    if (!mat) return;
    if (label && label[0]) {
        printf("%s [%d x %d]:\n", label, M, N);
    }
    for (int i = 0; i < M; i++) {
        printf("  ");
        for (int j = 0; j < N; j++) {
            double re = creal(mat[i * N + j]);
            double im = cimag(mat[i * N + j]);
            if (im >= 0) {
                printf("%7.3f+%7.3fj ", re, im);
            } else {
                printf("%7.3f-%7.3fj ", re, -im);
            }
        }
        printf("\n");
    }
}

void sfg_cmplx_matrix_mul_vec(const double complex* A,
                               const double complex* x,
                               double complex* y, int M, int N) {
    if (!A || !x || !y) return;
    for (int i = 0; i < M; i++) {
        y[i] = 0.0;
        for (int j = 0; j < N; j++) {
            y[i] += A[i * N + j] * x[j];
        }
    }
}

double complex sfg_cmplx_vector_dot(const double complex* x,
                                      const double complex* y, int N) {
    if (!x || !y) return 0.0;
    double complex sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += x[i] * y[i];
    }
    return sum;
}
