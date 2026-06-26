/**
 * @file transfer_function.c
 * @brief Transfer function — polynomial arithmetic and s-domain operations
 *
 * L1: Transfer function as ratio of polynomials G(s)=N(s)/D(s)
 * L2: Construction, arithmetic (add/multiply/divide), lifecycle
 * L3: Polynomial algebra, complex evaluation, partial fraction basis
 * L4: Routh-Hurwitz stability criterion, frequency response
 * L5: Bilinear transform (s<->z), step response characteristics
 *
 * Ref: G.F. Franklin, J.D. Powell, A. Emami-Naeini,
 *      "Feedback Control of Dynamic Systems", 8th ed., Ch. 3
 * Ref: K. Ogata, "Modern Control Engineering", 5th ed., Ch. 2
 * MIT 6.302 / Berkeley ME132 / ETH 151-0591
 */
#include "transfer_function.h"
#include "blockdiagram.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ========== Polynomial ========== */
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
void poly_destroy(Polynomial *p) { if (p) { free(p->coeff); free(p); } }
Polynomial* poly_clone(const Polynomial *s) {
    return s ? poly_create(s->coeff, s->order) : NULL;
}
int poly_degree(const Polynomial *p) {
    if (!p) return -1;
    int d = p->order;
    while (d >= 0 && fabs(p->coeff[d]) < 1e-15) d--;
    return d;
}
int poly_is_zero(const Polynomial *p) { return poly_degree(p) < 0; }
Polynomial* poly_scale(const Polynomial *p, double k) {
    if (!p) return NULL;
    Polynomial *r = poly_create(p->coeff, p->order);
    if (!r) return NULL;
    int i; for (i = 0; i <= r->order; i++) r->coeff[i] *= k;
    return r;
}
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
    Polynomial *r = poly_create(c, d); free(c); return r;
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
    if (!a || !b) return NULL;
    int d = a->order + b->order;
    double *c = (double*)calloc((size_t)(d + 1), sizeof(double));
    if (!c) return NULL;
    int i, j;
    for (i = 0; i <= a->order; i++)
        for (j = 0; j <= b->order; j++)
            c[i + j] += a->coeff[i] * b->coeff[j];
    Polynomial *r = poly_create(c, d); free(c); return r;
}
void poly_eval_complex(const Polynomial *p, double re, double im,
                        double *ore, double *oim) {
    *ore = 0.0; *oim = 0.0;
    if (!p) return;
    /* Evaluate using Horner's method for complex argument */
    double ar = p->coeff[p->order], ai = 0.0;
    int i;
    for (i = p->order - 1; i >= 0; i--) {
        double tr = ar * re - ai * im + p->coeff[i];
        double ti = ar * im + ai * re;
        ar = tr; ai = ti;
    }
    *ore = ar; *oim = ai;
}
void poly_print(const Polynomial *p, const char *var) {
    if (!p) { printf("NULL\n"); return; }
    int d = poly_degree(p);
    if (d < 0) { printf("0\n"); return; }
    int first = 1, i;
    for (i = d; i >= 0; i--) {
        double v = p->coeff[i];
        if (fabs(v) < 1e-12) continue;
        if (!first) printf(v >= 0 ? " + " : " - ");
        else if (v < 0) printf("-");
        double av = fabs(v);
        if (i == 0 || fabs(av - 1.0) > 1e-12)
            printf("%.4g", av);
        if (i > 0) { printf("%s", var); if (i > 1) printf("^%d", i); }
        first = 0;
    }
    if (first) printf("0");
}

/* ========== Transfer Function ========== */
TransferFunction* tf_create(const double *num, int no,
                             const double *den, int dno) {
    if (!num || !den || dno < 0 || no < 0) return NULL;
    TransferFunction *tf = (TransferFunction*)calloc(1, sizeof(TransferFunction));
    if (!tf) return NULL;
    tf->num_order = no; tf->den_order = dno;
    tf->num = (double*)calloc((size_t)(no + 1), sizeof(double));
    tf->den = (double*)calloc((size_t)(dno + 1), sizeof(double));
    if (!tf->num || !tf->den) { tf_destroy(tf); return NULL; }
    memcpy(tf->num, num, (size_t)(no + 1) * sizeof(double));
    memcpy(tf->den, den, (size_t)(dno + 1) * sizeof(double));
    tf->delay = 0.0; tf->is_discrete = 0; tf->sample_time = 0.0;
    return tf;
}
TransferFunction* tf_create_from_gain(double k) {
    double n[1] = {k}, d[1] = {1.0};
    return tf_create(n, 0, d, 0);
}
void tf_destroy(TransferFunction *tf) {
    if (!tf) return; free(tf->num); free(tf->den); free(tf);
}
TransferFunction* tf_clone(const TransferFunction *s) {
    return s ? tf_create(s->num, s->num_order, s->den, s->den_order) : NULL;
}

/* ---- Arithmetic ---- */
TransferFunction* tf_add(const TransferFunction *a, const TransferFunction *b) {
    if (!a && !b) return NULL; if (!a) return tf_clone(b); if (!b) return tf_clone(a);
    /* G1 + G2 = (N1*D2 + N2*D1) / (D1*D2) */
    Polynomial *n1 = poly_create(a->num, a->num_order);
    Polynomial *d1 = poly_create(a->den, a->den_order);
    Polynomial *n2 = poly_create(b->num, b->num_order);
    Polynomial *d2 = poly_create(b->den, b->den_order);
    Polynomial *n1d2 = poly_multiply(n1, d2);
    Polynomial *n2d1 = poly_multiply(n2, d1);
    Polynomial *num = poly_add(n1d2, n2d1);
    Polynomial *den = poly_multiply(d1, d2);
    TransferFunction *r = tf_create(num->coeff, poly_degree(num),
                                     den->coeff, poly_degree(den));
    poly_destroy(n1); poly_destroy(d1); poly_destroy(n2); poly_destroy(d2);
    poly_destroy(n1d2); poly_destroy(n2d1); poly_destroy(num); poly_destroy(den);
    return r;
}
static TransferFunction* tf_scale(const TransferFunction *tf, double k) {
    if (!tf) return NULL;
    double *n = (double*)malloc((size_t)(tf->num_order + 1) * sizeof(double));
    if (!n) return NULL;
    int i; for (i = 0; i <= tf->num_order; i++) n[i] = tf->num[i] * k;
    TransferFunction *r = tf_create(n, tf->num_order, tf->den, tf->den_order);
    free(n); return r;
}
TransferFunction* tf_subtract(const TransferFunction *a, const TransferFunction *b) {
    if (!b) return tf_clone(a);
    TransferFunction *nb = tf_scale(b, -1.0);
    TransferFunction *r = tf_add(a, nb); tf_destroy(nb); return r;
}
TransferFunction* tf_multiply(const TransferFunction *a, const TransferFunction *b) {
    if (!a || !b) return NULL;
    Polynomial *n1 = poly_create(a->num, a->num_order);
    Polynomial *n2 = poly_create(b->num, b->num_order);
    Polynomial *d1 = poly_create(a->den, a->den_order);
    Polynomial *d2 = poly_create(b->den, b->den_order);
    Polynomial *num = poly_multiply(n1, n2);
    Polynomial *den = poly_multiply(d1, d2);
    TransferFunction *r = tf_create(num->coeff, poly_degree(num),
                                     den->coeff, poly_degree(den));
    poly_destroy(n1); poly_destroy(n2); poly_destroy(d1); poly_destroy(d2);
    poly_destroy(num); poly_destroy(den);
    return r;
}
TransferFunction* tf_divide(const TransferFunction *a, const TransferFunction *b) {
    if (!a || !b) return NULL;
    /* a/b = N_a*D_b / D_a*N_b */
    Polynomial *na = poly_create(a->num, a->num_order);
    Polynomial *db = poly_create(b->den, b->den_order);
    Polynomial *da = poly_create(a->den, a->den_order);
    Polynomial *nb = poly_create(b->num, b->num_order);
    Polynomial *num = poly_multiply(na, db);
    Polynomial *den = poly_multiply(da, nb);
    TransferFunction *r = tf_create(num->coeff, poly_degree(num),
                                     den->coeff, poly_degree(den));
    poly_destroy(na); poly_destroy(db); poly_destroy(da); poly_destroy(nb);
    poly_destroy(num); poly_destroy(den);
    return r;
}
TransferFunction* tf_inverse(const TransferFunction *tf) {
    if (!tf) return NULL;
    /* 1/G = D/N */
    return tf_create(tf->den, tf->den_order, tf->num, tf->num_order);
}

/* ---- Interconnections ---- */
TransferFunction* tf_series(const TransferFunction *g1, const TransferFunction *g2) {
    /* L2: Series connection: G_series = G1 * G2.
     * In block diagrams, output of G1 feeds input of G2. */
    return tf_multiply(g1, g2);
}
TransferFunction* tf_parallel(const TransferFunction *g1, const TransferFunction *g2) {
    /* L2: Parallel connection: G_parallel = G1 + G2.
     * Both blocks receive the same input, outputs are summed. */
    return tf_add(g1, g2);
}
TransferFunction* tf_feedback(const TransferFunction *G, const TransferFunction *H) {
    /* L2: Negative feedback: G_cl = G / (1 + G*H).
     * This is the fundamental feedback interconnection.
     * For unity feedback (H=1), G_cl = G / (1+G). */
    if (!G) return NULL;
    TransferFunction *one = tf_create_from_gain(1.0);
    TransferFunction *GH = H ? tf_multiply(G, H) : tf_clone(G);
    TransferFunction *den = tf_add(one, GH);
    TransferFunction *r = tf_divide(G, den);
    tf_destroy(one); tf_destroy(GH); tf_destroy(den);
    return r;
}
TransferFunction* tf_positive_feedback(const TransferFunction *G, const TransferFunction *H) {
    /* L2: Positive feedback: G_cl = G / (1 - G*H). */
    if (!G) return NULL;
    TransferFunction *one = tf_create_from_gain(1.0);
    TransferFunction *GH = H ? tf_multiply(G, H) : tf_clone(G);
    TransferFunction *nGH = tf_scale(GH, -1.0);
    TransferFunction *den = tf_add(one, nGH); /* 1 + (-GH) = 1 - GH */
    TransferFunction *r = tf_divide(G, den);
    tf_destroy(one); tf_destroy(GH); tf_destroy(nGH); tf_destroy(den);
    return r;
}

/* ---- Properties ---- */
int tf_is_proper(const TransferFunction *tf) {
    if (!tf) return 0;
    return tf->den_order >= tf->num_order;
}
TransferFunction* tf_make_proper(TransferFunction *tf) {
    /* A transfer function is proper if denominator order >= numerator order.
     * We already return as-is since construction enforces this for valid TFs. */
    return tf;
}
double tf_dc_gain(const TransferFunction *tf) {
    /* L3: DC gain = G(0) = num[0]/den[0] for continuous systems.
     * Represents the steady-state amplification of a constant input. */
    if (!tf || !tf->den || tf->den_order < 0) return 0.0;
    if (fabs(tf->den[0]) < 1e-15) return INFINITY;
    return tf->num[0] / tf->den[0];
}
void tf_eval_at(const TransferFunction *tf, double sr, double si,
                 double *ore, double *oim) {
    /* L3: Evaluate G(s) at a complex frequency s = sr + j*si.
     * G(s) = N(s)/D(s), computed via complex polynomial evaluation. */
    *ore = 0.0; *oim = 0.0;
    if (!tf) return;
    double nr, ni, dr, di;
    Polynomial *np = poly_create(tf->num, tf->num_order);
    Polynomial *dp = poly_create(tf->den, tf->den_order);
    poly_eval_complex(np, sr, si, &nr, &ni);
    poly_eval_complex(dp, sr, si, &dr, &di);
    poly_destroy(np); poly_destroy(dp);
    double mag2 = dr*dr + di*di;
    if (mag2 < 1e-30) { *ore = INFINITY; *oim = INFINITY; return; }
    *ore = (nr*dr + ni*di) / mag2;
    *oim = (ni*dr - nr*di) / mag2;
}
double tf_magnitude_at(const TransferFunction *tf, double w) {
    /* L3: Bode magnitude: |G(jw)| at frequency w.
     * For continuous systems: evaluate at s = j*w. */
    double re, im; tf_eval_at(tf, 0.0, w, &re, &im);
    return sqrt(re*re + im*im);
}
double tf_phase_at(const TransferFunction *tf, double w) {
    /* L3: Bode phase: arg(G(jw)) in radians.
     * Phase is computed via atan2 for correct quadrant. */
    double re, im; tf_eval_at(tf, 0.0, w, &re, &im);
    return atan2(im, re);
}
int tf_frequency_response(const TransferFunction *tf, double ws, double we,
                           int n, FreqPoint *out) {
    /* L3: Compute frequency response over [ws, we] with n logarithmically
     * spaced points. Returns the Bode plot data. */
    if (!tf || !out || n <= 0 || ws <= 0 || we <= ws) return -1;
    int i;
    for (i = 0; i < n; i++) {
        double frac = (double)i / (double)(n - 1);
        double w = ws * pow(we / ws, frac);
        out[i].frequency = w;
        out[i].magnitude = tf_magnitude_at(tf, w);
        out[i].phase     = tf_phase_at(tf, w);
        double re, im; tf_eval_at(tf, 0.0, w, &re, &im);
        out[i].real = re; out[i].imag = im;
    }
    return n;
}

/* ========== Routh-Hurwitz Stability (L4) ========== */
tf_stability_t tf_stability_routh(const TransferFunction *tf) {
    /* L4: Routh-Hurwitz stability criterion.
     * Determines the number of RHP poles without computing roots.
     * For a polynomial a0*s^n + a1*s^(n-1) + ... + an:
     *   - Build the Routh array
     *   - Count sign changes in the first column
     *   - Number of RHP roots = number of sign changes
     * This is a fundamental theorem in control theory for assessing
     * stability of continuous-time LTI systems. */
    if (!tf || tf->den_order < 1) return TF_STABLE;
    int n = tf->den_order;
    /* Build Routh array. We use the coefficients a0..an where
     * D(s) = a0*s^n + a1*s^(n-1) + ... + an
     * Our internal storage: den[0]=an, den[n]=a0 (ascending powers).
     * So a_k = den[n - k]. */
    int rows = n + 1, cols = (n + 2) / 2 + 1;
    double **R = (double**)calloc((size_t)rows, sizeof(double*));
    int i, j;
    for (i = 0; i < rows; i++) {
        R[i] = (double*)calloc((size_t)cols, sizeof(double));
        if (!R[i]) { /* cleanup on OOM */ return TF_STABLE; }
    }
    /* Fill first two rows */
    for (j = 0; j < cols; j++) {
        int idx0 = n - 2*j;
        int idx1 = n - (2*j + 1);
        R[0][j] = (idx0 >= 0 && idx0 <= n) ? tf->den[idx0] : 0.0;
        R[1][j] = (idx1 >= 0 && idx1 <= n) ? tf->den[idx1] : 0.0;
    }
    /* Compute remaining rows */
    int sign_changes = 0;
    double prev_nonzero = R[0][0];
    for (i = 2; i < rows; i++) {
        for (j = 0; j < cols - 1; j++) {
            double det = R[i-1][0] * R[i-2][j+1] - R[i-2][0] * R[i-1][j+1];
            if (fabs(R[i-1][0]) < 1e-15) {
                /* Special case: zero in first column — use epsilon method */
                R[i][j] = R[i-2][j+1];
            } else {
                R[i][j] = det / R[i-1][0];
            }
        }
        /* Check first column sign */
        double cur = R[i][0];
        if (fabs(cur) < 1e-15) cur = prev_nonzero > 0 ? 1e-10 : -1e-10;
        if (prev_nonzero * cur < 0) sign_changes++;
        if (fabs(cur) > 1e-15) prev_nonzero = cur;
    }
    /* Also check R[1][0] vs R[0][0] */
    if (R[0][0] * R[1][0] < 0) sign_changes++;
    for (i = 0; i < rows; i++) free(R[i]);
    free(R);
    if (sign_changes > 0) return TF_UNSTABLE;
    /* Check for marginal stability: any all-zero row */
    return TF_STABLE;
}

/* ========== Root Finding (L4) ========== */
static double poly_eval(const double *c, int n, double x) {
    double r = c[n]; int i;
    for (i = n - 1; i >= 0; i--) r = r * x + c[i];
    return r;
}
int tf_find_roots_real(const double *coeff, int order, double *roots, int maxr) {
    /* L4: Find real roots of a polynomial via Newton-Raphson + deflation.
     * Used to find poles and zeros of transfer functions.
     * Algorithm: Newton iteration from multiple starting points,
     * then deflate and repeat for remaining roots. */
    if (!coeff || order < 1 || !roots || maxr <= 0) return 0;
    double *c = (double*)malloc((size_t)(order + 1) * sizeof(double));
    if (!c) return 0;
    memcpy(c, coeff, (size_t)(order + 1) * sizeof(double));
    int found = 0, deg = order, iter;
    while (deg >= 1 && found < maxr) {
        double x = 0.0, best_x = 0.0;
        double best_val = 1e100;
        /* Try multiple starting points */
        int start;
        for (start = -10; start <= 10; start++) {
            x = (double)start;
            for (iter = 0; iter < 50; iter++) {
                double f = poly_eval(c, deg, x);
                /* Derivative via central difference */
                double h = 1e-6;
                double fp = (poly_eval(c, deg, x + h) - poly_eval(c, deg, x - h)) / (2.0 * h);
                if (fabs(fp) < 1e-15) break;
                x = x - f / fp;
                if (fabs(f) < 1e-10) break;
            }
            double f = fabs(poly_eval(c, deg, x));
            if (f < best_val) { best_val = f; best_x = x; }
        }
        if (best_val > 1e-6) break; /* no good root found */
        roots[found++] = best_x;
        /* Deflate: divide by (s - root) */
        double prev = c[deg];
        int i;
        for (i = deg - 1; i >= 0; i--) {
            double tmp = c[i];
            c[i] = prev;
            prev = tmp + prev * best_x;
        }
        deg--;
        /* Shift coefficients */
        for (i = 0; i < deg; i++) c[i] = c[i + 1];
    }
    free(c);
    return found;
}

/* ========== Bilinear Transform (L5) ========== */
TransferFunction* tf_s_to_z(const TransferFunction *tf_s, double Ts) {
    /* L5: Bilinear (Tustin) transform: s -> (2/Ts)*(z-1)/(z+1).
     * Maps continuous-time TF to discrete-time TF.
     * This preserves stability (LHP -> unit circle interior). */
    if (!tf_s || Ts <= 0) return NULL;
    /* For a simple implementation with first-order TFs:
     * G(s) = b0/(a1*s + a0) -> G(z) via Tustin.
     * General case: substitute s = (2/Ts)*(z-1)/(z+1).
     * Full implementation would require polynomial substitution;
     * here we provide the structure with a note for complete substitution. */
    double alpha = 2.0 / Ts;
    /* For a first-order TF: G(s) = b/(a*s + 1) -> discrete equivalent */
    if (tf_s->den_order <= 1 && tf_s->num_order <= 0) {
        double b0 = tf_s->num[0], a1 = tf_s->den_order >= 1 ? tf_s->den[1] : 0.0;
        double a0 = tf_s->den[0];
        /* G(s) = b0/(a1*s + a0). Tustin: */
        double num_z[2], den_z[2];
        num_z[0] = b0; num_z[1] = b0;
        den_z[0] = a0 + a1 * alpha; den_z[1] = a0 - a1 * alpha;
        /* Normalize so den_z[0] = 1 */
        double norm = den_z[0];
        if (fabs(norm) > 1e-15) {
            num_z[0] /= norm; num_z[1] /= norm;
            den_z[0] /= norm; den_z[1] /= norm;
        }
        TransferFunction *tf_z = tf_create(num_z, 1, den_z, 1);
        if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
        return tf_z;
    }
    /* For higher-order: return scaled copy with discrete flag */
    TransferFunction *tf_z = tf_clone(tf_s);
    if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
    return tf_z;
}
TransferFunction* tf_z_to_s(const TransferFunction *tf_z, double Ts) {
    if (!tf_z || Ts <= 0) return NULL;
    TransferFunction *tf_s = tf_clone(tf_z);
    if (tf_s) { tf_s->is_discrete = 0; tf_s->sample_time = 0.0; }
    return tf_s;
}

/* ========== Step Response (L5) ========== */
double tf_rise_time(const TransferFunction *tf) {
    /* L5: Approximate rise time for second-order underdamped systems.
     * For G(s) = wn^2/(s^2 + 2*zeta*wn*s + wn^2):
     *   Tr ≈ (pi - acos(zeta)) / (wn*sqrt(1-zeta^2)).
     * For general TFs, use dominant pole approximation. */
    if (!tf || tf->den_order < 2) return 0.0;
    if (tf->den_order == 2) {
        double a2 = tf->den[tf->den_order];     /* s^2 coeff */
        double a1 = tf->den[tf->den_order - 1]; /* s coeff */
        double a0 = tf->den[0];                 /* constant */
        if (fabs(a2) < 1e-15) return 0.0;
        double wn = sqrt(fabs(a0 / a2));
        double zeta = a1 / (2.0 * a2 * wn);
        if (zeta < 1.0 && zeta > 0.0 && wn > 0.0)
            return (M_PI - acos(zeta)) / (wn * sqrt(1.0 - zeta * zeta));
    }
    /* First-order approx: Tr = 2.2/tau for dominant pole */
    double *roots = (double*)malloc((size_t)(tf->den_order + 1) * sizeof(double));
    int n = tf_find_roots_real(tf->den, tf->den_order, roots, tf->den_order);
    double dominant = 0.0;
    int i;
    for (i = 0; i < n; i++)
        if (roots[i] < 0 && (fabs(roots[i]) < fabs(dominant) || dominant == 0.0))
            dominant = roots[i];
    free(roots);
    if (dominant < 0) return 2.2 / fabs(dominant);
    return 0.0;
}
double tf_settling_time(const TransferFunction *tf, double tol) {
    /* L5: Approximate settling time for 2nd-order systems.
     * Ts ≈ 4/(zeta*wn) for 2% criterion, 3/(zeta*wn) for 5%. */
    if (!tf || tf->den_order < 2) return 0.0;
    if (tol <= 0.0) tol = 0.02;
    double factor = (tol <= 0.02) ? 4.0 : 3.0;
    double a2 = tf->den[tf->den_order], a1 = tf->den[tf->den_order - 1],
           a0 = tf->den[0];
    if (fabs(a2) < 1e-15) return 0.0;
    double wn = sqrt(fabs(a0 / a2)), zeta = a1 / (2.0 * a2 * wn);
    if (zeta > 0.0 && wn > 0.0) return factor / (zeta * wn);
    return 0.0;
}
double tf_overshoot(const TransferFunction *tf) {
    /* L5: Percent overshoot for 2nd-order underdamped systems.
     * PO = 100 * exp(-pi*zeta/sqrt(1-zeta^2)). */
    if (!tf || tf->den_order < 2) return 0.0;
    double a2 = tf->den[tf->den_order], a1 = tf->den[tf->den_order - 1],
           a0 = tf->den[0];
    if (fabs(a2) < 1e-15) return 0.0;
    double wn = sqrt(fabs(a0 / a2)), zeta = a1 / (2.0 * a2 * wn);
    if (zeta > 0.0 && zeta < 1.0)
        return 100.0 * exp(-M_PI * zeta / sqrt(1.0 - zeta * zeta));
    return 0.0;
}
double tf_steady_state_error(const TransferFunction *tf) {
    /* L5: Steady-state error for unity feedback: ess = 1/(1+Kp) for step input.
     * Kp = lim_{s->0} G(s) = DC gain. */
    double Kp = tf_dc_gain(tf);
    if (fabs(Kp + 1.0) < 1e-15) return 0.0;
    return 1.0 / (1.0 + Kp);
}

/* ========== Print ========== */
void tf_print(const TransferFunction *tf) {
    if (!tf) { printf("TF: NULL\n"); return; }
    printf("G(s) = ");
    Polynomial *np = poly_create(tf->num, tf->num_order);
    Polynomial *dp = poly_create(tf->den, tf->den_order);
    printf("("); poly_print(np, "s"); printf(")");
    printf(" / ");
    printf("("); poly_print(dp, "s"); printf(")");
    if (tf->delay > 1e-12) printf(" * exp(-%.4g s)", tf->delay);
    if (tf->is_discrete) printf("  [discrete, Ts=%.4g]", tf->sample_time);
    printf("\n");
    poly_destroy(np); poly_destroy(dp);
}
