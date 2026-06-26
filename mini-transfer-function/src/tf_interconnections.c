/**
 * @file tf_interconnections.c
 * @brief Transfer Function Interconnection Algebra
 *
 * Knowledge points:
 *   L1/L2: Series, parallel, feedback interconnections
 *   L2: Algebraic operations (add, subtract, multiply, divide, scale)
 *   L3/L4: Block diagram reduction, sensitivity, complementary sensitivity
 *   L4: Internal stability, pole-zero cancellation detection
 *   L5: Multiple interconnections (cascade, sum)
 *   L8: Kharitonov polynomials for robust stability
 *
 * Ref: Franklin Ch.3-4; Ogata Ch.3; Mason (1956); Kharitonov (1978)
 */

#include "tf_interconnections.h"
#include "tf_analysis.h"
#include "tf_polynomial.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TF_EPS 1e-15

/* ================================================================
 * L1/L2: Basic Interconnections
 * ================================================================ */

TransferFunction* tf_series(const TransferFunction *g1, const TransferFunction *g2) {
    /* Series: G_series(s) = G1(s) * G2(s).
     * Output of G1 feeds input of G2. Multiplication of TFs. */
    return tf_multiply(g1, g2);
}

TransferFunction* tf_parallel(const TransferFunction *g1, const TransferFunction *g2) {
    /* Parallel: G_parallel(s) = G1(s) + G2(s).
     * Both blocks share input, outputs summed. Addition of TFs. */
    return tf_add(g1, g2);
}

TransferFunction* tf_feedback(const TransferFunction *G, const TransferFunction *H) {
    /* Negative feedback: G_cl(s) = G(s) / (1 + G(s)*H(s)).
     * This is THE fundamental feedback formula in control theory.
     * For unity feedback (H=1): G_cl = G / (1+G).
     *
     * Derivation: Y = G*E, E = R - H*Y => Y = G*(R - H*Y)
     * => Y*(1 + G*H) = G*R => Y/R = G/(1+G*H) */
    if (!G) return NULL;

    TransferFunction *one = tf_create_from_gain(1.0);
    TransferFunction *GH = H ? tf_multiply(G, H) : tf_clone(G);
    TransferFunction *den = tf_add(one, GH);
    TransferFunction *result = tf_divide(G, den);

    tf_destroy(one); tf_destroy(GH); tf_destroy(den);
    return result;
}

TransferFunction* tf_positive_feedback(const TransferFunction *G, const TransferFunction *H) {
    /* Positive feedback: G_cl(s) = G(s) / (1 - G(s)*H(s)).
     * Used in oscillators and some specialized control structures. */
    if (!G) return NULL;

    TransferFunction *one = tf_create_from_gain(1.0);
    TransferFunction *GH = H ? tf_multiply(G, H) : tf_clone(G);
    TransferFunction *negGH = tf_scale(GH, -1.0);
    TransferFunction *den = tf_add(one, negGH); /* 1 + (-GH) = 1 - GH */
    TransferFunction *result = tf_divide(G, den);

    tf_destroy(one); tf_destroy(GH); tf_destroy(negGH); tf_destroy(den);
    return result;
}

TransferFunction* tf_unity_feedback(const TransferFunction *G) {
    /* Unity feedback: G_cl = G / (1 + G).
     * Most common feedback configuration - sensor gain is 1. */
    return tf_feedback(G, NULL);
}

/* ================================================================
 * L2: Algebraic Operations
 * ================================================================ */

TransferFunction* tf_add(const TransferFunction *a, const TransferFunction *b) {
    /* G1 + G2 = (N1*D2 + N2*D1) / (D1*D2). */
    if (!a && !b) return NULL;
    if (!a) return tf_clone(b);
    if (!b) return tf_clone(a);

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

TransferFunction* tf_subtract(const TransferFunction *a, const TransferFunction *b) {
    if (!b) return tf_clone(a);
    TransferFunction *nb = tf_scale(b, -1.0);
    TransferFunction *r = tf_add(a, nb);
    tf_destroy(nb);
    return r;
}

TransferFunction* tf_multiply(const TransferFunction *a, const TransferFunction *b) {
    /* G1 * G2 = (N1*N2) / (D1*D2). */
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
    /* G1 / G2 = (N1*D2) / (D1*N2). */
    if (!a || !b) return NULL;

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

TransferFunction* tf_scale(const TransferFunction *tf, double k) {
    /* k * G(s): multiply numerator by k. */
    if (!tf) return NULL;
    double *n = (double*)malloc((size_t)(tf->num_order + 1) * sizeof(double));
    if (!n) return NULL;
    int i;
    for (i = 0; i <= tf->num_order; i++) n[i] = tf->num[i] * k;
    TransferFunction *r = tf_create(n, tf->num_order, tf->den, tf->den_order);
    free(n);
    if (r) {
        r->delay = tf->delay;
        r->is_discrete = tf->is_discrete;
        r->sample_time = tf->sample_time;
    }
    return r;
}

/* ================================================================
 * L3/L4: Block Diagram Reduction and Sensitivity
 * ================================================================ */

TransferFunction* tf_disturbance_to_output(const TransferFunction *G,
                                            const TransferFunction *H,
                                            const TransferFunction *Gd) {
    /* Disturbance response: Y(s)/D(s) = Gd(s) / (1 + G(s)*H(s)).
     * With input disturbance entering before the plant. */
    if (!G || !Gd) return NULL;
    TransferFunction *den = tf_feedback(G, H); /* this returns G/(1+GH)... */
    /* Actually: Y/D = Gd / (1+G*H) = Gd * (1/(1+GH)) */
    TransferFunction *one = tf_create_from_gain(1.0);
    TransferFunction *GH = H ? tf_multiply(G, H) : tf_clone(G);
    TransferFunction *one_plus_GH = tf_add(one, GH);
    TransferFunction *result = tf_divide(Gd, one_plus_GH);

    tf_destroy(one); tf_destroy(GH); tf_destroy(one_plus_GH);
    tf_destroy(den);
    return result;
}

TransferFunction* tf_sensitivity(const TransferFunction *G, const TransferFunction *H) {
    /* Sensitivity function: S(s) = 1 / (1 + G(s)*H(s)).
     * S(s) = transfer function from reference to error.
     * |S(jw)| << 1 at low frequencies for good disturbance rejection. */
    if (!G) return tf_create_from_gain(1.0);
    TransferFunction *one = tf_create_from_gain(1.0);
    TransferFunction *GH = H ? tf_multiply(G, H) : tf_clone(G);
    TransferFunction *den = tf_add(one, GH);
    TransferFunction *result = tf_divide(one, den);

    tf_destroy(one); tf_destroy(GH); tf_destroy(den);
    return result;
}

TransferFunction* tf_complementary_sensitivity(const TransferFunction *G,
                                                 const TransferFunction *H) {
    /* Complementary sensitivity: T(s) = G*H / (1 + G*H).
     * T(s) = transfer function from reference to output.
     * Note: S(s) + T(s) = 1 (fundamental identity). */
    if (!G) return NULL;
    TransferFunction *one = tf_create_from_gain(1.0);
    TransferFunction *GH = H ? tf_multiply(G, H) : tf_clone(G);
    TransferFunction *den = tf_add(one, GH);
    TransferFunction *result = tf_divide(GH, den);

    tf_destroy(one); tf_destroy(GH); tf_destroy(den);
    return result;
}

TransferFunction* tf_loop_transfer(const TransferFunction *G, const TransferFunction *H) {
    /* Loop transfer: L(s) = G(s)*H(s).
     * Open-loop TF used for Nyquist/Bode analysis. */
    if (!G) return NULL;
    return H ? tf_multiply(G, H) : tf_clone(G);
}

TransferFunction* tf_nested_feedback(const TransferFunction *G1,
                                      const TransferFunction *G2,
                                      const TransferFunction *H1,
                                      const TransferFunction *H2) {
    /* Nested loops: inner loop G2 with H2, outer G1 with H1.
     * Inner: G_inner = G2/(1 + G2*H2)
     * Outer: G_total = G1*G_inner/(1 + G1*G_inner*H1) */
    if (!G1 || !G2) return NULL;

    TransferFunction *inner = tf_feedback(G2, H2);
    if (!inner) return NULL;

    TransferFunction *G1_inner = tf_multiply(G1, inner);
    TransferFunction *result = tf_feedback(G1_inner, H1);

    tf_destroy(inner); tf_destroy(G1_inner);
    return result;
}

TransferFunction* tf_feedforward(const TransferFunction *G_fb,
                                  const TransferFunction *G_ff,
                                  const TransferFunction *G) {
    /* Feedforward + feedback: Y = (G_ff + G_fb)*G*R.
     * 2-DOF control structure. */
    if (!G || !G_fb) return NULL;
    TransferFunction *combined = G_ff ? tf_add(G_fb, G_ff) : tf_clone(G_fb);
    TransferFunction *result = tf_multiply(combined, G);
    tf_destroy(combined);
    return result;
}

/* ================================================================
 * L4: Interconnection Theorems
 * ================================================================ */

int tf_has_pole_zero_cancellation(const TransferFunction *G, const TransferFunction *H) {
    /* Detect if G(s)*H(s) has common factors in num and den.
     * Pole-zero cancellation can hide internal instability.
     * Check if GCD(num, den) is non-trivial (degree > 0). */
    if (!G) return 0;

    TransferFunction *GH = H ? tf_multiply(G, H) : tf_clone(G);
    Polynomial *np = poly_create(GH->num, GH->num_order);
    Polynomial *dp = poly_create(GH->den, GH->den_order);
    Polynomial *gcd = poly_gcd(np, dp);

    int has_cancel = (poly_degree(gcd) > 0);

    poly_destroy(np); poly_destroy(dp); poly_destroy(gcd);
    tf_destroy(GH);
    return has_cancel;
}

tf_stability_t tf_internal_stability(const TransferFunction *G, const TransferFunction *H) {
    /* Internal stability: all four TFs must have LHP poles.
     * The four transfer functions are:
     *   G/(1+GH), H/(1+GH), 1/(1+GH), GH/(1+GH).
     *
     * Any RHP pole-zero cancellation between G and H causes
     * internal instability even if the closed-loop TF appears stable. */
    if (!G) return TF_STABLE;

    TransferFunction *S  = tf_sensitivity(G, H);
    TransferFunction *T  = tf_complementary_sensitivity(G, H);
    TransferFunction *GS = tf_multiply(G, S);

    tf_stability_t s_stable  = tf_stability_routh(S);
    tf_stability_t t_stable  = tf_stability_routh(T);
    tf_stability_t gs_stable = tf_stability_routh(GS);

    tf_destroy(S); tf_destroy(T); tf_destroy(GS);

    if (s_stable == TF_UNSTABLE || t_stable == TF_UNSTABLE ||
        gs_stable == TF_UNSTABLE)
        return TF_UNSTABLE;
    return TF_STABLE;
}

/* ================================================================
 * L5: Multiple Interconnections
 * ================================================================ */

TransferFunction* tf_cascade_n(const TransferFunction **tfs, int n) {
    /* Product of N TFs: G1 * G2 * ... * Gn. */
    if (!tfs || n <= 0) return NULL;
    TransferFunction *result = tf_clone(tfs[0]);
    if (!result) return NULL;

    int i;
    for (i = 1; i < n; i++) {
        TransferFunction *tmp = tf_multiply(result, tfs[i]);
        tf_destroy(result);
        result = tmp;
        if (!result) return NULL;
    }
    return result;
}

TransferFunction* tf_sum_n(const TransferFunction **tfs, int n) {
    /* Sum of N TFs: G1 + G2 + ... + Gn. */
    if (!tfs || n <= 0) return NULL;
    TransferFunction *result = tf_clone(tfs[0]);
    if (!result) return NULL;

    int i;
    for (i = 1; i < n; i++) {
        TransferFunction *tmp = tf_add(result, tfs[i]);
        tf_destroy(result);
        result = tmp;
        if (!result) return NULL;
    }
    return result;
}

int tf_kharitonov_polynomials(const double *num_low, const double *num_high,
                               const double *den_low, const double *den_high,
                               int order, TransferFunction **k_out) {
    /* Kharitonov's theorem (1978): An interval polynomial family
     * P(s) = sum [ak_min, ak_max] * s^k is robustly stable iff
     * four specific "Kharitonov polynomials" are stable.
     *
     * The four polynomials use alternating min/max coefficients:
     *   K1(s) = a0_min + a1_min*s + a2_max*s^2 + a3_max*s^3 + ...
     *   K2(s) = a0_max + a1_max*s + a2_min*s^2 + a3_min*s^3 + ...
     *   K3(s) = a0_min + a1_max*s + a2_max*s^2 + a3_min*s^3 + ...
     *   K4(s) = a0_max + a1_min*s + a2_min*s^2 + a3_max*s^3 + ...
     *
     * This is a cornerstone of robust parametric stability analysis.
     * L8: Advanced stability methods. */
    if (!den_low || !den_high || !k_out || order < 1) return -1;

    int i, k;
    for (k = 0; k < 4; k++) {
        double *den = (double*)malloc((size_t)(order + 1) * sizeof(double));
        double *num = (double*)malloc((size_t)(order + 1) * sizeof(double));
        if (!den || !num) { free(den); free(num); return -1; }

        for (i = 0; i <= order; i++) {
            int even_i = (i % 2 == 0);
            /* Select min or max based on Kharitonov pattern k and parity i */
            int use_max_den = 0;
            switch (k) {
                case 0: use_max_den = !even_i; break;
                case 1: use_max_den = even_i; break;
                case 2: use_max_den = !even_i; break;
                case 3: use_max_den = even_i; break;
            }
            /* K1: even=min, odd=min; K2: even=max, odd=max;
             * K3: even=min, odd=max; K4: even=max, odd=min */
            if (k == 0) use_max_den = 0;
            else if (k == 1) use_max_den = 1;
            else if (k == 2) use_max_den = even_i ? 1 : 0;
            else use_max_den = even_i ? 0 : 1;

            den[i] = use_max_den ? den_high[i] : den_low[i];
            num[i] = use_max_den ? num_high[i] : num_low[i];
        }

        k_out[k] = tf_create(num, order, den, order);
        free(num); free(den);
    }
    return 4;
}
