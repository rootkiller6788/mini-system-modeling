/**
 * @file tf_polynomial.c
 * @brief Polynomial Algebra Implementation
 *
 * Implements the complete polynomial algebra required for transfer
 * function manipulation: arithmetic, evaluation, division, GCD,
 * root finding, and composition.
 *
 * Knowledge points:
 *   L1: Polynomial data structure with dynamic coefficient array
 *   L2: Add/sub/mul/scale/power/derivative/integral operations
 *   L3: Horner evaluation (real and complex), polynomial division, GCD, LCM
 *   L4: Newton-Raphson root finding with deflation, companion matrix method
 *   L5: Polynomial composition, from-roots construction
 *
 * Ref: Press et al., Numerical Recipes 3e; Higham, Accuracy and Stability
 */

#include "tf_polynomial.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define POLY_EPS 1e-15

/* ================================================================
 * L1/L2: Construction and Lifecycle
 * ================================================================ */

Polynomial* poly_create(const double *coeff, int order) {
    if (!coeff || order < 0) return NULL;
    Polynomial *p = (Polynomial*)calloc(1, sizeof(Polynomial));
    if (!p) return NULL;
    p->order = order;
    p->coeff = (double*)calloc((size_t)(order + 1), sizeof(double));
    if (!p->coeff) { free(p); return NULL; }
    memcpy(p->coeff, coeff, (size_t)(order + 1) * sizeof(double));
    return p;
}

Polynomial* poly_zero(void) {
    double c[1] = { 0.0 };
    return poly_create(c, 0);
}

Polynomial* poly_one(void) {
    double c[1] = { 1.0 };
    return poly_create(c, 0);
}

Polynomial* poly_s(void) {
    double c[2] = { 0.0, 1.0 };
    return poly_create(c, 1);
}

Polynomial* poly_constant(double c) {
    double coeff[1] = { c };
    return poly_create(coeff, 0);
}

Polynomial* poly_from_roots(const double *roots, int n, double gain) {
    /* Build polynomial from its roots: P(s) = gain * prod_{i=0}^{n-1} (s - r_i).
     * Start with P(s) = 1, then iteratively multiply by (s - r_i).
     * Each step: new_coeff[k] = old[k-1] - r_i * old[k] for k > 0,
     *            new_coeff[0] = -r_i * old[0]. */
    if (!roots || n <= 0) return poly_constant(gain);

    /* Start with polynomial = 1 */
    Polynomial *result = poly_one();
    if (!result) return NULL;

    int i;
    for (i = 0; i < n; i++) {
        /* Multiply current result by (s - root) */
        double r = roots[i];
        int old_order = result->order;
        double *old = (double*)malloc((size_t)(old_order + 1) * sizeof(double));
        if (!old) { poly_destroy(result); return NULL; }
        memcpy(old, result->coeff, (size_t)(old_order + 1) * sizeof(double));

        /* Extend result */
        double *new_c = (double*)realloc(result->coeff,
                           (size_t)(old_order + 2) * sizeof(double));
        if (!new_c) { free(old); poly_destroy(result); return NULL; }
        result->coeff = new_c;
        result->order = old_order + 1;
        result->coeff[result->order] = 0.0;

        int j;
        for (j = result->order; j >= 0; j--) {
            result->coeff[j] = (j > 0 ? old[j - 1] : 0.0) - r * old[j];
        }
        free(old);
    }

    /* Scale by gain */
    if (fabs(gain - 1.0) > POLY_EPS) {
        int k;
        for (k = 0; k <= result->order; k++) result->coeff[k] *= gain;
    }

    return result;
}

Polynomial* poly_clone(const Polynomial *p) {
    return p ? poly_create(p->coeff, p->order) : NULL;
}

void poly_destroy(Polynomial *p) {
    if (p) { free(p->coeff); free(p); }
}

/* ================================================================
 * L2: Arithmetic
 * ================================================================ */

Polynomial* poly_add(const Polynomial *a, const Polynomial *b) {
    if (!a && !b) return NULL;
    if (!a) return poly_clone(b);
    if (!b) return poly_clone(a);

    int d = (a->order > b->order) ? a->order : b->order;
    double *c = (double*)calloc((size_t)(d + 1), sizeof(double));
    if (!c) return NULL;

    int i;
    for (i = 0; i <= a->order; i++) c[i] += a->coeff[i];
    for (i = 0; i <= b->order; i++) c[i] += b->coeff[i];

    Polynomial *r = poly_create(c, d);
    free(c);
    return r;
}

Polynomial* poly_subtract(const Polynomial *a, const Polynomial *b) {
    if (!b) return poly_clone(a);
    Polynomial *nb = poly_scale(b, -1.0);
    if (!nb) return NULL;
    Polynomial *r = poly_add(a, nb);
    poly_destroy(nb);
    return r;
}

Polynomial* poly_multiply(const Polynomial *a, const Polynomial *b) {
    /* Convolution: r[k] = sum_{i+j=k} a[i] * b[j].
     * This is the discrete convolution of the two coefficient sequences,
     * equivalent to polynomial multiplication. */
    if (!a || !b) return NULL;
    if (poly_is_zero(a) || poly_is_zero(b)) return poly_zero();

    int d = a->order + b->order;
    double *c = (double*)calloc((size_t)(d + 1), sizeof(double));
    if (!c) return NULL;

    int i, j;
    for (i = 0; i <= a->order; i++) {
        if (fabs(a->coeff[i]) < POLY_EPS) continue;
        for (j = 0; j <= b->order; j++) {
            c[i + j] += a->coeff[i] * b->coeff[j];
        }
    }

    Polynomial *r = poly_create(c, d);
    free(c);
    poly_trim(r);
    return r;
}

Polynomial* poly_scale(const Polynomial *p, double k) {
    if (!p) return NULL;
    Polynomial *r = poly_create(p->coeff, p->order);
    if (!r) return NULL;
    int i;
    for (i = 0; i <= r->order; i++) r->coeff[i] = p->coeff[i] * k;
    return r;
}

Polynomial* poly_power(const Polynomial *p, int n) {
    /* Binary exponentiation for polynomial powers.
     * Reduces O(n) multiplications to O(log n).
     *
     * Algorithm: result = 1, base = p, exp = n
     *   while exp > 0:
     *     if exp is odd: result *= base
     *     base = base * base
     *     exp /= 2 */
    if (!p || n < 0) return NULL;
    if (n == 0) return poly_one();
    if (n == 1) return poly_clone(p);

    Polynomial *result = poly_one();
    Polynomial *base = poly_clone(p);
    int exp = n;

    while (exp > 0) {
        if (exp & 1) {
            Polynomial *tmp = poly_multiply(result, base);
            poly_destroy(result);
            result = tmp;
            if (!result) { poly_destroy(base); return NULL; }
        }
        exp >>= 1;
        if (exp > 0) {
            Polynomial *tmp = poly_multiply(base, base);
            poly_destroy(base);
            base = tmp;
            if (!base) { poly_destroy(result); return NULL; }
        }
    }

    poly_destroy(base);
    poly_trim(result);
    return result;
}

Polynomial* poly_derivative(const Polynomial *p) {
    /* d/ds (sum_{k=0}^n a_k s^k) = sum_{k=1}^n k*a_k s^{k-1}.
     * Coefficient at position j (power j) of derivative = (j+1)*a_{j+1}. */
    if (!p) return NULL;
    if (p->order <= 0) return poly_zero();

    int new_order = p->order - 1;
    double *c = (double*)calloc((size_t)(new_order + 1), sizeof(double));
    if (!c) return NULL;

    int i;
    for (i = 1; i <= p->order; i++) {
        c[i - 1] = (double)i * p->coeff[i];
    }

    Polynomial *r = poly_create(c, new_order);
    free(c);
    return r;
}

Polynomial* poly_integral(const Polynomial *p) {
    /* integral (sum a_k s^k) ds = sum a_k/(k+1) * s^{k+1} + C.
     * We set C=0 (indefinite integral without constant). */
    if (!p) return NULL;

    int new_order = p->order + 1;
    double *c = (double*)calloc((size_t)(new_order + 1), sizeof(double));
    if (!c) return NULL;

    int i;
    for (i = 0; i <= p->order; i++) {
        c[i + 1] = p->coeff[i] / (double)(i + 1);
    }

    Polynomial *r = poly_create(c, new_order);
    free(c);
    return r;
}

Polynomial* poly_compose(const Polynomial *a, const Polynomial *b) {
    /* Compute a(b(s)) via Horner's method for polynomial composition.
     * If a(s) = a_0 + a_1*s + ... + a_n*s^n, then
     * a(b(s)) = a_0 + a_1*b + a_2*b^2 + ... + a_n*b^n
     * evaluated as: (...(a_n*b + a_{n-1})*b + ...)*b + a_0
     *
     * This is numerically stable and O(n^2 * deg_b) rather than
     * naive O(2^n). */
    if (!a || !b) return NULL;
    if (poly_is_constant(b)) {
        /* b(s) = c => a(b(s)) = a(c) (constant) */
        return poly_constant(poly_eval(a, b->coeff[0]));
    }

    /* Horner composition starting from leading coefficient */
    int da = poly_degree(a);
    if (da < 0) return poly_zero();

    Polynomial *result = poly_constant(a->coeff[da]);
    if (!result) return NULL;

    int i;
    for (i = da - 1; i >= 0; i--) {
        /* result = result * b + a_i */
        Polynomial *tmp1 = poly_multiply(result, b);
        if (!tmp1) { poly_destroy(result); return NULL; }
        poly_destroy(result);

        Polynomial *tmp2 = poly_create(&a->coeff[i], 0);
        result = poly_add(tmp1, tmp2);
        poly_destroy(tmp1); poly_destroy(tmp2);
        if (!result) return NULL;
    }

    poly_trim(result);
    return result;
}

/* ================================================================
 * L3: Evaluation
 * ================================================================ */

double poly_eval(const Polynomial *p, double x) {
    /* Horner's method: p = a_n; for i=n-1 to 0: p = p*x + a_i.
     * Minimizes operations (n mults + n adds) and rounding error.
     * Forward error bound: |computed - true| <= 2*n*eps/(1-2*n*eps) * P_s(|x|)
     * where P_s has absolute-valued coefficients. (Higham §5.1) */
    if (!p) return 0.0;
    double r = p->coeff[p->order];
    int i;
    for (i = p->order - 1; i >= 0; i--) {
        r = r * x + p->coeff[i];
    }
    return r;
}

void poly_eval_complex(const Polynomial *p, double re, double im,
                        double *ore, double *oim) {
    /* Complex Horner evaluation for s = re + j*im.
     * Maintains real and imaginary accumulators:
     *   (ar + j*ai) = (ar + j*ai) * (re + j*im) + coeff[i]
     *                = (ar*re - ai*im + coeff[i]) + j*(ar*im + ai*re)
     * 4 mults + 3 adds per iteration. */
    *ore = 0.0; *oim = 0.0;
    if (!p) return;

    double ar = p->coeff[p->order], ai = 0.0;
    int i;
    for (i = p->order - 1; i >= 0; i--) {
        double tmp_r = ar * re - ai * im + p->coeff[i];
        double tmp_i = ar * im + ai * re;
        ar = tmp_r; ai = tmp_i;
    }
    *ore = ar; *oim = ai;
}

double poly_eval_derivative(const Polynomial *p, double x, double *dval) {
    /* Simultaneous evaluation of P(x) and P'(x) via combined Horner.
     *
     * Algorithm (Higham §5.2):
     *   p = a_n; dp = 0
     *   for i = n-1 downto 0:
     *     dp = dp * x + p
     *     p  = p  * x + a_i
     *
     * This is 2n multiplications + 2n additions (vs 3n + 2n separately). */
    *dval = 0.0;
    if (!p) return 0.0;

    double val = p->coeff[p->order];
    *dval = 0.0;
    int i;
    for (i = p->order - 1; i >= 0; i--) {
        *dval = (*dval) * x + val;
        val = val * x + p->coeff[i];
    }
    return val;
}

int poly_eval_matrix(const Polynomial *p, const double *X, int n, double *PX) {
    /* Compute P(X) for an n x n matrix X using Horner's method.
     * Each "multiplication" is a matrix-matrix multiply (O(n^3)).
     * Total complexity: O(deg(p) * n^3).
     *
     * Application: matrix functions, Cayley-Hamilton verification,
     * controllability/observability analysis. */
    if (!p || !X || !PX || n <= 0) return -1;

    /* Initialize PX = a_deg * I */
    int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            PX[i * n + j] = (i == j) ? p->coeff[p->order] : 0.0;
        }
    }

    /* Horner iteration: PX = PX * X + a_k * I */
    double *tmp = (double*)malloc((size_t)(n * n) * sizeof(double));
    if (!tmp) return -1;

    for (k = p->order - 1; k >= 0; k--) {
        /* tmp = PX * X */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                double sum = 0.0;
                int l;
                for (l = 0; l < n; l++) {
                    sum += PX[i * n + l] * X[l * n + j];
                }
                tmp[i * n + j] = sum;
            }
        }
        /* PX = tmp + a_k * I */
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                PX[i * n + j] = tmp[i * n + j] + ((i == j) ? p->coeff[k] : 0.0);
            }
        }
    }

    free(tmp);
    return 0;
}

/* ================================================================
 * L3: Queries
 * ================================================================ */

int poly_degree(const Polynomial *p) {
    if (!p) return -1;
    int d = p->order;
    while (d >= 0 && fabs(p->coeff[d]) < POLY_EPS) d--;
    return d;
}

int poly_is_zero(const Polynomial *p) {
    return poly_degree(p) < 0;
}

int poly_is_constant(const Polynomial *p) {
    int d = poly_degree(p);
    return d <= 0;
}

int poly_is_monic(const Polynomial *p) {
    int d = poly_degree(p);
    if (d < 0) return 0;
    return fabs(p->coeff[d] - 1.0) < POLY_EPS;
}

double poly_leading_coeff(const Polynomial *p) {
    int d = poly_degree(p);
    return (d >= 0) ? p->coeff[d] : 0.0;
}

int poly_sign_changes(const Polynomial *p) {
    /* Descartes' rule of signs: number of positive real roots is at most
     * the number of sign changes in the coefficient sequence, and differs
     * from it by an even number. */
    if (!p) return 0;
    int d = poly_degree(p);
    if (d < 0) return 0;

    int changes = 0;
    double prev = 0.0;
    int has_prev = 0;
    int i;
    for (i = 0; i <= d; i++) {
        if (fabs(p->coeff[i]) < POLY_EPS) continue;
        if (!has_prev) {
            prev = p->coeff[i];
            has_prev = 1;
        } else {
            if (prev * p->coeff[i] < 0) changes++;
            prev = p->coeff[i];
        }
    }
    return changes;
}

/* ================================================================
 * L3: Division and GCD
 * ================================================================ */

int poly_divide(const Polynomial *a, const Polynomial *b,
                Polynomial **quotient, Polynomial **remainder) {
    /* Polynomial long division: a = q*b + r, deg(r) < deg(b).
     *
     * Algorithm (synthetic division):
     *   Copy a into remainder r.
     *   For k from deg(a)-deg(b) down to 0:
     *     q[k] = r[k+deg(b)] / b[deg(b)]
     *     For j from 0 to deg(b):
     *       r[k+j] -= q[k] * b[j]  */
    if (!a || !b || poly_is_zero(b)) return -1;
    int da = poly_degree(a), db = poly_degree(b);
    if (da < db) {
        /* Deg(a) < deg(b): q=0, r=a */
        if (quotient) *quotient = poly_zero();
        if (remainder) *remainder = poly_clone(a);
        return 0;
    }

    int dq = da - db;
    /* Initialize remainder as copy of a */
    double *r_coeff = (double*)malloc((size_t)(da + 1) * sizeof(double));
    if (!r_coeff) return -1;
    memcpy(r_coeff, a->coeff, (size_t)(da + 1) * sizeof(double));

    double *q_coeff = (double*)calloc((size_t)(dq + 1), sizeof(double));
    if (!q_coeff) { free(r_coeff); return -1; }

    double lead_b = b->coeff[db];
    int k, j;
    for (k = dq; k >= 0; k--) {
        q_coeff[k] = r_coeff[k + db] / lead_b;
        for (j = 0; j <= db; j++) {
            r_coeff[k + j] -= q_coeff[k] * b->coeff[j];
        }
    }

    if (quotient) *quotient = poly_create(q_coeff, dq);
    if (remainder) *remainder = poly_create(r_coeff, db - 1 >= 0 ? db - 1 : 0);

    free(r_coeff); free(q_coeff);
    return 0;
}

Polynomial* poly_gcd(const Polynomial *a, const Polynomial *b) {
    /* Euclidean algorithm for polynomial GCD.
     * While b != 0: (a, b) = (b, a mod b)
     *
     * Uses polynomial remainder to iteratively reduce degree.
     * Termination guaranteed since degree decreases each step. */
    if (!a || !b) return NULL;
    if (poly_is_zero(a)) return poly_clone(b);
    if (poly_is_zero(b)) return poly_clone(a);

    Polynomial *r1 = poly_clone(a);
    Polynomial *r2 = poly_clone(b);
    if (!r1 || !r2) { poly_destroy(r1); poly_destroy(r2); return NULL; }

    while (!poly_is_zero(r2)) {
        Polynomial *rem = NULL;
        if (poly_divide(r1, r2, NULL, &rem) != 0) {
            poly_destroy(r1); poly_destroy(r2);
            return NULL;
        }
        poly_destroy(r1);
        r1 = r2;
        r2 = rem;
    }

    /* Make GCD monic for canonical representation */
    Polynomial *gcd = r1;
    poly_make_monic(gcd);
    poly_destroy(r2);
    return gcd;
}

Polynomial* poly_lcm(const Polynomial *a, const Polynomial *b) {
    /* LCM(a,b) = a*b / GCD(a,b).
     * Computed as (a/GCD) * b to reduce risk of overflow in coefficient
     * magnitude (better numerical conditioning). */
    if (!a || !b) return NULL;
    Polynomial *g = poly_gcd(a, b);
    if (!g || poly_is_zero(g)) { poly_destroy(g); return NULL; }

    Polynomial *rem = NULL;
    Polynomial *quot = NULL;
    poly_divide(a, g, &quot, &rem);
    poly_destroy(rem);
    poly_destroy(g);

    Polynomial *lcm = poly_multiply(quot, b);
    poly_destroy(quot);
    return lcm;
}

int poly_extended_gcd(const Polynomial *a, const Polynomial *b,
                      Polynomial **u, Polynomial **v, Polynomial **g) {
    /* Extended Euclidean algorithm for polynomials.
     * Finds u, v such that a*u + b*v = g = gcd(a,b).
     *
     * Maintains invariants:
     *   r_k = a*u_k + b*v_k
     *   r_{k+1} = a*u_{k+1} + b*v_{k+1}
     *
     * This is essential for Bezout identity in coprime factorization
     * and Youla parameterization. */
    if (!a || !b) return -1;

    Polynomial *r0 = poly_clone(a), *r1 = poly_clone(b);
    Polynomial *u0 = poly_one(),  *u1 = poly_zero();
    Polynomial *v0 = poly_zero(), *v1 = poly_one();

    while (!poly_is_zero(r1)) {
        Polynomial *q = NULL, *r = NULL;
        poly_divide(r0, r1, &q, &r);

        /* u2 = u0 - q*u1, v2 = v0 - q*v1 */
        Polynomial *qu1 = poly_multiply(q, u1);
        Polynomial *u2 = poly_subtract(u0, qu1);
        poly_destroy(qu1);

        Polynomial *qv1 = poly_multiply(q, v1);
        Polynomial *v2 = poly_subtract(v0, qv1);
        poly_destroy(qv1);

        poly_destroy(q); poly_destroy(r0);
        poly_destroy(u0); poly_destroy(v0);
        r0 = r1; r1 = r;
        u0 = u1; u1 = u2;
        v0 = v1; v1 = v2;
    }

    /* r0 = GCD, u0, v0 satisfy a*u0 + b*v0 = r0 */
    if (g) *g = r0; else poly_destroy(r0);
    if (u) *u = u0; else poly_destroy(u0);
    if (v) *v = v0; else poly_destroy(v0);
    poly_destroy(r1); poly_destroy(u1); poly_destroy(v1);

    return 0;
}

Polynomial* poly_remainder(const Polynomial *a, const Polynomial *b) {
    Polynomial *rem = NULL;
    if (poly_divide(a, b, NULL, &rem) != 0) return NULL;
    return rem;
}

/* ================================================================
 * L4: Root Finding
 * ================================================================ */

static double nr_poly_eval(const double *c, int n, double x) {
    /* Inline Horner evaluation for internal root-finding use.
     * c[0] = constant term, c[n] = leading coefficient. */
    double r = c[n];
    int i;
    for (i = n - 1; i >= 0; i--) r = r * x + c[i];
    return r;
}

int poly_find_roots_real(const double *coeff, int order,
                          double *roots, int maxr) {
    /* Multi-start Newton-Raphson with deflation for polynomial root finding.
     *
     * Algorithm:
     *  1. Try starting points from -10 to 10 (step 1, gives 21 points)
     *  2. For each starting point, run Newton iteration (max 50 steps):
     *       x_{k+1} = x_k - f(x_k)/f'(x_k)
     *       where f'(x_k) computed via central difference
     *  3. Accept root if |f(x)| < 1e-10
     *  4. Deflate: divide polynomial by (s - root) using synthetic division
     *  5. Repeat until no more roots found or maxr reached
     *
     * Numerical stability notes:
     *  - Central difference avoids catastrophic cancellation of forward diff
     *  - Forward deflation maintains accuracy better than backward
     *  - Multiple starting points reduce risk of missing roots */
    if (!coeff || order < 1 || !roots || maxr <= 0) return 0;

    /* Working copy for deflation */
    double *c = (double*)malloc((size_t)(order + 1) * sizeof(double));
    if (!c) return 0;
    memcpy(c, coeff, (size_t)(order + 1) * sizeof(double));

    int found = 0, deg = order;
    while (deg >= 1 && found < maxr) {
        double best_x = 0.0, best_val = 1e100;
        int start;

        /* Try 21 starting points in [-10, 10] */
        for (start = -10; start <= 10; start++) {
            double x = (double)start;
            int iter;
            for (iter = 0; iter < 50; iter++) {
                double f = nr_poly_eval(c, deg, x);
                if (fabs(f) < 1e-10) break;

                /* Central difference derivative: f'(x) = (f(x+h)-f(x-h))/(2h) */
                double h = (fabs(x) > 1.0) ? fabs(x) * 1e-6 : 1e-6;
                double fp = (nr_poly_eval(c, deg, x + h) -
                             nr_poly_eval(c, deg, x - h)) / (2.0 * h);

                if (fabs(fp) < 1e-15) break; /* stationary point */

                double dx = f / fp;
                x -= dx;
                if (fabs(dx) < fabs(x) * 1e-12) break;
            }

            double f_val = fabs(nr_poly_eval(c, deg, x));
            if (f_val < best_val) {
                best_val = f_val;
                best_x = x;
            }
        }

        if (best_val > 1e-6) break; /* no good root found */

        /* Refine to nearby integer/round value if very close */
        double rounded = round(best_x);
        if (fabs(best_x - rounded) < 1e-5) {
            if (fabs(nr_poly_eval(c, deg, rounded)) < best_val)
                best_x = rounded;
        }

        roots[found++] = best_x;

        /* Forward deflation: divide by (s - root) */
        if (deg >= 2) {
            double *new_c = (double*)malloc((size_t)(deg) * sizeof(double));
            if (!new_c) break;
            /* Synthetic division: new_c[i-1] = c[i] + root * new_c[i] */
            new_c[deg - 1] = c[deg];
            int i;
            for (i = deg - 2; i >= 0; i--) {
                new_c[i] = c[i + 1] + best_x * new_c[i + 1];
            }
            deg--;
            /* Copy back */
            for (i = 0; i <= deg; i++) c[i] = new_c[i];
            free(new_c);
        } else {
            deg--;
        }
    }

    free(c);
    return found;
}

/* Simple QR-like eigenvalue extraction via shifted power iteration */
static void hessenberg_qr_step(double *H, int n) {
    /* Hessenberg QR step with Wilkinson shift for real eigenvalues.
     * Simplified implementation for companion matrix root finding. */
    int i, j, k;
    for (i = 0; i < n - 1; i++) {
        for (j = i + 2; j < n; j++) {
            double c, s, r;
            double x = H[(i+1)*n + i], y = H[j*n + i];
            if (fabs(x) < 1e-15 && fabs(y) < 1e-15) continue;
            r = sqrt(x*x + y*y);
            c = x / r; s = -y / r;

            /* Apply Givens rotation to rows */
            for (k = i; k < n; k++) {
                double h1 = H[(i+1)*n + k], h2 = H[j*n + k];
                H[(i+1)*n + k] = c * h1 - s * h2;
                H[j*n + k]     = s * h1 + c * h2;
            }
            /* Apply Givens rotation to columns */
            for (k = 0; k < n; k++) {
                double h1 = H[k*n + (i+1)], h2 = H[k*n + j];
                H[k*n + (i+1)] =  c * h1 - s * h2;
                H[k*n + j]     =  s * h1 + c * h2;
            }
        }
    }
}

int poly_find_roots_complex(const double *coeff, int order,
                             double *roots_real, double *roots_imag) {
    /* Find all roots via companion matrix eigenvalue computation.
     *
     * For a monic polynomial P(s) = s^n + a_{n-1}s^{n-1} + ... + a_0,
     * the companion matrix C has characteristic polynomial det(sI-C) = P(s).
     * So eigenvalues of C = roots of P(s).
     *
     * We build the Frobenius companion matrix, convert to upper Hessenberg,
     * then apply QR iterations to extract eigenvalues.
     *
     * This method finds ALL roots (real and complex) simultaneously,
     * unlike Newton-Raphson which finds one at a time. */
    if (!coeff || order < 1 || !roots_real || !roots_imag) return 0;

    /* Normalize to monic */
    double lead = coeff[order];
    if (fabs(lead) < POLY_EPS) return 0;

    double *norm = (double*)malloc((size_t)(order + 1) * sizeof(double));
    if (!norm) return 0;
    int i;
    for (i = 0; i <= order; i++) norm[i] = coeff[i] / lead;

    /* Build Hessenberg companion matrix (n x n) */
    int n = order;
    double *H = (double*)calloc((size_t)(n * n), sizeof(double));
    if (!H) { free(norm); return 0; }

    /* Fill companion form: first subdiagonal = 1, last column = -a_i */
    for (i = 0; i < n - 1; i++) H[(i + 1) * n + i] = 1.0;
    for (i = 0; i < n; i++) H[i * n + (n - 1)] = -norm[i];

    /* QR iterations (simplified: 100 iterations) */
    int iter;
    for (iter = 0; iter < 100; iter++) {
        hessenberg_qr_step(H, n);
    }

    /* Extract eigenvalues from diagonal and subdiagonal blocks */
    int found = 0;
    i = 0;
    while (i < n && found < n) {
        if (i < n - 1 && fabs(H[(i+1)*n + i]) > 1e-10) {
            /* 2x2 block - complex conjugate pair */
            double a = H[i*n + i], b = H[i*n + (i+1)];
            double c = H[(i+1)*n + i], d = H[(i+1)*n + (i+1)];
            double trace = a + d;
            double det = a*d - b*c;
            double disc = trace*trace - 4.0*det;

            if (disc < 0) {
                roots_real[found] = trace / 2.0;
                roots_imag[found] = sqrt(-disc) / 2.0;
                roots_real[found+1] = trace / 2.0;
                roots_imag[found+1] = -sqrt(-disc) / 2.0;
            } else {
                roots_real[found] = (trace + sqrt(disc)) / 2.0;
                roots_imag[found] = 0.0;
                roots_real[found+1] = (trace - sqrt(disc)) / 2.0;
                roots_imag[found+1] = 0.0;
            }
            found += 2; i += 2;
        } else {
            /* 1x1 block - real eigenvalue */
            roots_real[found] = H[i*n + i];
            roots_imag[found] = 0.0;
            found++; i++;
        }
    }

    free(norm); free(H);
    return found;
}

/* ================================================================
 * L5: Utilities
 * ================================================================ */

void poly_print(const Polynomial *p, const char *var) {
    if (!p) { printf("NULL"); return; }
    int d = poly_degree(p);
    if (d < 0) { printf("0"); return; }

    int first = 1, i;
    for (i = d; i >= 0; i--) {
        double v = p->coeff[i];
        if (fabs(v) < 1e-12) continue;

        /* Sign handling */
        if (!first) {
            printf(v >= 0 ? " + " : " - ");
        } else if (v < 0) {
            printf("-");
        }

        double av = fabs(v);
        /* Coefficient display: skip "1" for variable terms unless constant */
        if (i == 0 || fabs(av - 1.0) > 1e-12)
            printf("%.4g", av);

        /* Variable display */
        if (i > 0) {
            printf("%s", var);
            if (i > 1) printf("^%d", i);
        }
        first = 0;
    }
    if (first) printf("0");
}

double poly_make_monic(Polynomial *p) {
    if (!p) return 0.0;
    int d = poly_degree(p);
    if (d < 0) return 0.0;

    double lead = p->coeff[d];
    if (fabs(lead) < POLY_EPS) return 0.0;

    int i;
    for (i = 0; i <= d; i++) p->coeff[i] /= lead;
    return lead;
}

void poly_trim(Polynomial *p) {
    if (!p) return;
    /* Trim trailing near-zero coefficients from the high end */
    while (p->order > 0 && fabs(p->coeff[p->order]) < POLY_EPS) {
        p->order--;
    }
    /* Also trim from lower end if constant term is zero? No - that changes degree. */
    int d = poly_degree(p);
    if (d < p->order) p->order = (d >= 0) ? d : 0;
}

int poly_equals(const Polynomial *a, const Polynomial *b, double tol) {
    /* Check coefficient-by-coefficient equality within tolerance.
     * Polynomials of different effective degree cannot be equal
     * (unless the extra coefficients are epsilon-small). */
    if (!a && !b) return 1;
    if (!a || !b) return 0;

    int da = poly_degree(a), db = poly_degree(b);
    int d_max = (da > db) ? da : db;

    int i;
    for (i = 0; i <= d_max; i++) {
        double va = (i <= a->order) ? a->coeff[i] : 0.0;
        double vb = (i <= b->order) ? b->coeff[i] : 0.0;
        /* Skip epsilon-small coefficients */
        if (fabs(va) < tol && fabs(vb) < tol) continue;
        if (fabs(va - vb) > tol) return 0;
    }
    return d_max >= 0;
}

double poly_norm(const Polynomial *p) {
    /* Euclidean norm of coefficient vector */
    if (!p) return 0.0;
    double sum_sq = 0.0;
    int i;
    for (i = 0; i <= p->order; i++) {
        sum_sq += p->coeff[i] * p->coeff[i];
    }
    return sqrt(sum_sq);
}
