/**
 * @file transfer_function.c
 * @brief Core Transfer Function Implementation
 *
 * Implements the fundamental data structures and operations for
 * transfer functions in control systems. Covers L1 definitions,
 * L2 lifecycle management, L3 basic properties, and L5 transformations.
 *
 * Knowledge points implemented:
 *   L1: TransferFunction struct construction, data types, enums
 *   L2: Create/clone/destroy lifecycle, factory functions for standard forms
 *   L3: Properness, relative degree, system type, phase classification
 *   L3: DC gain, HF gain, pole/zero counting, equality check
 *   L4: Normalization (monic constant-term and monic leading-coeff)
 *   L5: Inverse, derivative, integral, frequency shift, time scaling
 *
 * References:
 *   G.F. Franklin, J.D. Powell, A. Emami-Naeini,
 *     "Feedback Control of Dynamic Systems", 8th ed., Ch. 3
 *   K. Ogata, "Modern Control Engineering", 5th ed., Ch. 2-3
 *   R.C. Dorf, R.H. Bishop, "Modern Control Systems", 13th ed., Ch. 2
 *   MIT 6.302 / Stanford ENGR105 / Berkeley ME132 / ETH 151-0591
 */

#include "transfer_function.h"
#include "tf_polynomial.h"
#include "tf_conversion.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Tolerance for floating-point comparisons */
#define TF_EPS 1e-15

/* ================================================================
 * L1: Lifecycle - Construction
 * ================================================================ */

TransferFunction* tf_create(const double *num, int num_order,
                             const double *den, int den_order) {
    /* Validate inputs: both arrays required, orders non-negative,
     * denominator must be non-trivial for proper systems */
    if (!num || !den || den_order < 0 || num_order < 0) return NULL;

    TransferFunction *tf = (TransferFunction*)calloc(1, sizeof(TransferFunction));
    if (!tf) return NULL;

    tf->num_order = num_order;
    tf->den_order = den_order;

    tf->num = (double*)calloc((size_t)(num_order + 1), sizeof(double));
    tf->den = (double*)calloc((size_t)(den_order + 1), sizeof(double));

    if (!tf->num || !tf->den) {
        tf_destroy(tf);
        return NULL;
    }

    memcpy(tf->num, num, (size_t)(num_order + 1) * sizeof(double));
    memcpy(tf->den, den, (size_t)(den_order + 1) * sizeof(double));

    tf->delay = 0.0;
    tf->is_discrete = 0;
    tf->sample_time = 0.0;

    return tf;
}

TransferFunction* tf_create_from_gain(double k) {
    /* G(s) = k = k/1.
     * A pure gain is both a transfer function and a constant.
     * This is the simplest non-trivial transfer function. */
    double num[1] = { k };
    double den[1] = { 1.0 };
    return tf_create(num, 0, den, 0);
}

TransferFunction* tf_create_integrator(void) {
    /* G(s) = 1/s.
     * The ideal integrator is the fundamental building block for
     * type-1 systems. It provides infinite DC gain and -90 degree
     * phase shift across all frequencies. */
    double num[1] = { 1.0 };
    double den[2] = { 0.0, 1.0 };  /* den[0]=0, den[1]=1 => 0 + 1*s = s */
    return tf_create(num, 0, den, 1);
}

TransferFunction* tf_create_first_order_lag(double K, double tau) {
    /* G(s) = K / (tau*s + 1).
     * First-order lag is the simplest dynamic system.
     * Time constant tau determines speed of response.
     * Pole at s = -1/tau (stable for tau > 0).
     * Step response: y(t) = K*(1 - exp(-t/tau)).
     * Settling time (2%): Ts = 4*tau. */
    if (tau <= 0.0) return NULL;
    double num[1] = { K };
    double den[2] = { 1.0, tau };
    return tf_create(num, 0, den, 1);
}

TransferFunction* tf_create_second_order(double K, double wn, double zeta) {
    /* G(s) = K*wn^2 / (s^2 + 2*zeta*wn*s + wn^2).
     *
     * The second-order system is the workhorse of control analysis:
     *  - zeta = 0: undamped (marginal stability, pure oscillation)
     *  - 0 < zeta < 1: underdamped (overshoot present)
     *  - zeta = 1: critically damped (fastest non-oscillatory)
     *  - zeta > 1: overdamped (slower, two real poles)
     *
     * Natural frequency wn determines speed.
     * Damping ratio zeta determines oscillation characteristics.
     *
     * Ref: Franklin §3.3, Ogata §5-3 */
    if (wn <= 0.0 || zeta < 0.0) return NULL;

    double num[1] = { K * wn * wn };
    double den[3] = { wn * wn, 2.0 * zeta * wn, 1.0 };
    return tf_create(num, 0, den, 2);
}

TransferFunction* tf_create_lead(double K, double T, double alpha) {
    /* G_lead(s) = K * (alpha*T*s + 1) / (T*s + 1), alpha > 1.
     *
     * Lead compensator adds positive phase (phase lead) to improve
     * stability margins. Maximum phase lead:
     *   phi_max = asin((alpha-1)/(alpha+1))
     * occurs at frequency w_max = 1/(T*sqrt(alpha)).
     *
     * Typical design: alpha = 5..20, providing 30..60 deg phase boost. */
    if (T <= 0.0 || alpha <= 1.0) return NULL;
    double num[2] = { K, K * alpha * T };
    double den[2] = { 1.0, T };
    return tf_create(num, 1, den, 1);
}

TransferFunction* tf_create_lag(double K, double T, double beta) {
    /* G_lag(s) = K * (T*s + 1) / (beta*T*s + 1), beta > 1.
     *
     * Lag compensator reduces high-frequency gain while preserving
     * low-frequency characteristics. Improves steady-state accuracy
     * without affecting transient response (if placed at low enough
     * frequency). Typical beta = 5..20. */
    if (T <= 0.0 || beta <= 1.0) return NULL;
    double num[2] = { K, K * T };
    double den[2] = { 1.0, beta * T };
    return tf_create(num, 1, den, 1);
}

TransferFunction* tf_create_pid(double Kp, double Ki, double Kd) {
    /* G_PID(s) = Kp + Ki/s + Kd*s = (Kd*s^2 + Kp*s + Ki) / s.
     *
     * PID is the most widely used controller in industry (>90% of loops).
     * - Proportional (Kp): affects present error
     * - Integral (Ki/s): eliminates steady-state error (type increase)
     * - Derivative (Kd*s): anticipates future error (phase lead)
     *
     * The derivative term is often filtered: Kd*s/(tau_f*s+1).
     * Here we provide the ideal (unfiltered) form.
     *
     * Ref: Astrom & Hagglund, "PID Controllers", 2nd ed. */
    double num[3] = { Ki, Kp, Kd };
    double den[2] = { 0.0, 1.0 };  /* s */
    return tf_create(num, 2, den, 1);
}

TransferFunction* tf_create_delay(double Td) {
    /* G(s) = exp(-Td*s).
     * Pure time delay (transport lag). Cannot be exactly represented as
     * a finite-order rational TF. Stored as metadata + unit TF.
     * For rational approximation, use Pade (see tf_advanced.c). */
    if (Td < 0.0) return NULL;
    double num[1] = { 1.0 };
    double den[1] = { 1.0 };
    TransferFunction *tf = tf_create(num, 0, den, 0);
    if (tf) tf->delay = Td;
    return tf;
}

TransferFunction* tf_clone(const TransferFunction *src) {
    if (!src) return NULL;
    TransferFunction *cpy = tf_create(src->num, src->num_order,
                                       src->den, src->den_order);
    if (cpy) {
        cpy->delay = src->delay;
        cpy->is_discrete = src->is_discrete;
        cpy->sample_time = src->sample_time;
    }
    return cpy;
}

void tf_destroy(TransferFunction *tf) {
    if (!tf) return;
    free(tf->num);
    free(tf->den);
    free(tf);
}

/* ================================================================
 * L3: Core Properties
 * ================================================================ */

tf_properness_t tf_properness(const TransferFunction *tf) {
    /* Properness classification determines whether a TF can be realized
     * as a causal state-space system.
     *   Strictly proper: more poles than zeros -> causal, no direct feedthrough
     *   Proper: equal poles and zeros -> causal, direct feedthrough = leading ratio
     *   Improper: more zeros than poles -> non-causal, differentiator present
     *
     * Ref: Franklin §3.1.4 */
    if (!tf) return TF_STRICTLY_PROPER;
    if (tf->den_order > tf->num_order) return TF_STRICTLY_PROPER;
    if (tf->den_order == tf->num_order) return TF_PROPER;
    return TF_IMPROPER;
}

int tf_relative_degree(const TransferFunction *tf) {
    /* Relative degree r = den_order - num_order.
     * Determines high-frequency roll-off: -20r dB/decade.
     * Also equals the number of integrators needed in the forward path
     * to make the closed-loop system proper. */
    if (!tf) return -1;
    return tf->den_order - tf->num_order;
}

tf_system_type_t tf_system_type(const TransferFunction *tf) {
    /* System type = number of integrators (poles at s=0).
     * Determines steady-state error characteristics:
     *   Type 0: finite step error, infinite ramp error
     *   Type 1: zero step error, finite ramp error, infinite parabolic error
     *   Type 2: zero step and ramp error, finite parabolic error
     *
     * Count by checking consecutive zero coefficients at the start
     * of the denominator polynomial. */
    if (!tf || !tf->den) return TF_TYPE_0;
    int count = 0, i;
    for (i = 0; i <= tf->den_order; i++) {
        if (fabs(tf->den[i]) < TF_EPS) count++;
        else break;
    }
    if (count >= 3) return TF_TYPE_3;
    if (count == 2) return TF_TYPE_2;
    if (count == 1) return TF_TYPE_1;
    return TF_TYPE_0;
}

tf_phase_type_t tf_phase_type(const TransferFunction *tf) {
    /* Minimum-phase systems have all zeros in LHP.
     * They have the smallest possible phase lag for a given magnitude
     * response (Bode's gain-phase relationship).
     *
     * Non-minimum phase systems have RHP zeros causing:
     *  - Initial undershoot in step response
     *  - Additional phase lag
     *  - Fundamental performance limitations */
    if (!tf) return TF_MINIMUM_PHASE;
    /* Use root finding to check for RHP zeros */
    double *zeros = (double*)malloc((size_t)(tf->num_order + 1) * sizeof(double));
    if (!zeros) return TF_MINIMUM_PHASE;
    int nz = poly_find_roots_real(tf->num, tf->num_order,
                                   zeros, tf->num_order);
    int i;
    for (i = 0; i < nz; i++) {
        if (zeros[i] > TF_EPS) {
            free(zeros);
            return TF_NONMINIMUM_PHASE;
        }
    }
    free(zeros);
    return TF_MINIMUM_PHASE;
}

int tf_has_rhp_zeros(const TransferFunction *tf) {
    return (tf_phase_type(tf) == TF_NONMINIMUM_PHASE);
}

double tf_dc_gain(const TransferFunction *tf) {
    /* DC gain = lim_{s->0} G(s) = num[0] / den[0].
     * For discrete systems evaluated at z=1:
     *   G(1) = sum(num[i]) / sum(den[i])
     *
     * DC gain represents steady-state amplification of a constant input.
     * It is the final value of the step response for stable systems.
     *
     * If den[0] = 0 (system has integrator), DC gain is infinite.
     * In that case, the system tracks step inputs with zero steady-state error. */
    if (!tf || !tf->den) return 0.0;

    if (tf->is_discrete) {
        /* Discrete DC gain: G(z=1) = sum(num)/sum(den) */
        double num_sum = 0.0, den_sum = 0.0;
        int i;
        for (i = 0; i <= tf->num_order; i++) num_sum += tf->num[i];
        for (i = 0; i <= tf->den_order; i++) den_sum += tf->den[i];
        if (fabs(den_sum) < TF_EPS) return INFINITY;
        return num_sum / den_sum;
    }

    /* Continuous DC gain */
    if (fabs(tf->den[0]) < TF_EPS) return INFINITY;
    return tf->num[0] / tf->den[0];
}

double tf_hf_gain(const TransferFunction *tf) {
    /* High-frequency gain: lim_{s->inf} G(s) = leading coefficient ratio.
     *
     * For s -> inf, G(s) ~ num[num_order] * s^{no} / (den[den_order] * s^{do})
     *                  = (num[num_order]/den[den_order]) * s^{no-do}
     *
     * For strictly proper systems (no < do): HF gain -> 0
     * For proper systems (no = do): HF gain = num_no/den_do
     * For improper systems: diverges
     *
     * This determines the initial value of the step response. */
    if (!tf) return 0.0;
    int no = tf->num_order, dno = tf->den_order;
    if (dno < 0) return 0.0;

    /* Trim trailing zeros */
    while (no >= 0 && fabs(tf->num[no]) < TF_EPS) no--;
    while (dno >= 0 && fabs(tf->den[dno]) < TF_EPS) dno--;

    if (dno < 0) return INFINITY; /* den = 0 (invalid) */
    if (no < dno) return 0.0;     /* strictly proper */
    if (no > dno) return INFINITY; /* improper, diverges */
    return tf->num[no] / tf->den[dno];
}

int tf_num_poles(const TransferFunction *tf) {
    /* Count true number of poles by trimming trailing near-zero coefficients.
     * A polynomial a_n s^n + ... + a_0 has degree n, but if a_n ~ 0,
     * the effective degree is lower. */
    if (!tf) return 0;
    int d = tf->den_order;
    while (d >= 0 && fabs(tf->den[d]) < TF_EPS) d--;
    return (d >= 0) ? d : 0;
}

int tf_num_zeros(const TransferFunction *tf) {
    /* Count true number of zeros. Same trimming logic as poles. */
    if (!tf) return 0;
    int n = tf->num_order;
    while (n >= 0 && fabs(tf->num[n]) < TF_EPS) n--;
    return (n >= 0) ? n : 0;
}

/* ================================================================
 * L3: Input/Output and Utility
 * ================================================================ */

void tf_print(const TransferFunction *tf) {
    if (!tf) { printf("TF: NULL\n"); return; }

    printf("G(s) = ");
    Polynomial *np = poly_create(tf->num, tf->num_order);
    Polynomial *dp = poly_create(tf->den, tf->den_order);

    printf("("); poly_print(np, "s"); printf(")");
    printf(" / ");
    printf("("); poly_print(dp, "s"); printf(")");

    if (tf->delay > TF_EPS) printf(" * exp(-%.4g s)", tf->delay);
    if (tf->is_discrete) printf("  [discrete, Ts=%.4g]", tf->sample_time);
    printf("\n");

    poly_destroy(np);
    poly_destroy(dp);
}

void tf_print_matlab(const TransferFunction *tf) {
    /* Print in MATLAB-compatible format for easy copy-paste.
     * Output: num = [b0 b1 ... bm]; den = [a0 a1 ... an]; G = tf(num, den) */
    if (!tf) { printf("%% TF: NULL\n"); return; }
    printf("%% MATLAB representation:\n");
    printf("num = [");
    int i;
    for (i = 0; i <= tf->num_order; i++) {
        printf("%.6g", tf->num[i]);
        if (i < tf->num_order) printf(", ");
    }
    printf("];\n");
    printf("den = [");
    for (i = 0; i <= tf->den_order; i++) {
        printf("%.6g", tf->den[i]);
        if (i < tf->den_order) printf(", ");
    }
    printf("];\n");
    printf("G = tf(num, den);\n");
    if (tf->delay > TF_EPS) printf("G = set(G, 'ioDelay', %.6g);\n", tf->delay);
}

void tf_print_zpk(const TransferFunction *tf) {
    /* Print transfer function in Zero-Pole-Gain form.
     * Requires root computation. */
    if (!tf) { printf("ZPK: NULL\n"); return; }
    ZeroPoleGain *zpk = tf_to_zpk(tf);
    if (!zpk) { printf("ZPK: unable to compute\n"); return; }
    zpk_print(zpk);
    zpk_destroy(zpk);
}

int tf_equals(const TransferFunction *a, const TransferFunction *b, double tol) {
    /* Two TFs are equal if their normalized forms have identical coefficients.
     * Cross-multiply to check: N_a*D_b == N_b*D_a (within tolerance).
     * This handles common-factor cancellation cases correctly. */
    if (!a && !b) return 1;
    if (!a || !b) return 0;

    /* Compare delays first */
    if (fabs(a->delay - b->delay) > tol) return 0;
    if (a->is_discrete != b->is_discrete) return 0;

    /* Cross-multiplication: N_a*D_b and N_b*D_a should match */
    Polynomial *na = poly_create(a->num, a->num_order);
    Polynomial *db = poly_create(b->den, b->den_order);
    Polynomial *nb = poly_create(b->num, b->num_order);
    Polynomial *da = poly_create(a->den, a->den_order);

    Polynomial *p1 = poly_multiply(na, db);
    Polynomial *p2 = poly_multiply(nb, da);

    int eq = poly_equals(p1, p2, tol);

    poly_destroy(na); poly_destroy(db);
    poly_destroy(nb); poly_destroy(da);
    poly_destroy(p1); poly_destroy(p2);

    return eq;
}

void tf_normalize(TransferFunction *tf) {
    /* Normalize so that den[0] = 1 (monic constant term).
     * This is convenient for type-0 systems where DC gain is clear.
     * For systems with integrators (den[0]=0), this is not applicable. */
    if (!tf) return;
    if (fabs(tf->den[0]) < TF_EPS) return; /* cannot normalize integrator */

    double factor = tf->den[0];
    int i;
    for (i = 0; i <= tf->num_order; i++) tf->num[i] /= factor;
    for (i = 0; i <= tf->den_order; i++) tf->den[i] /= factor;
}

void tf_normalize_monic(TransferFunction *tf) {
    /* Normalize so that den[den_order] = 1 (monic polynomial).
     * This is the standard form for control analysis: the highest power
     * of s in the denominator has coefficient 1.
     *
     * Example: (2s+1)/(3s^2+4s+5) -> (0.67s+0.33)/(s^2+1.33s+1.67)
     *
     * This form matches the characteristic polynomial convention
     * used in most textbooks: D(s) = s^n + a_{n-1}s^{n-1} + ... + a_0. */
    if (!tf) return;

    /* Find leading non-zero denominator coefficient */
    int d = tf->den_order;
    while (d >= 0 && fabs(tf->den[d]) < TF_EPS) d--;
    if (d < 0) return; /* zero denominator - invalid state */

    double factor = tf->den[d];
    if (fabs(factor) < TF_EPS) return;

    int i;
    for (i = 0; i <= tf->num_order; i++) tf->num[i] /= factor;
    for (i = 0; i <= d; i++) tf->den[i] /= factor;
}

void tf_trim(TransferFunction *tf) {
    /* Remove leading and trailing zero coefficients.
     * Trailing zeros: reduce num_order/den_order if highest coef ~ 0.
     * Leading zeros: shift and reduce if lowest coef ~ 0 for BOTH num and den
     * (common factor of s can be cancelled). */
    if (!tf) return;

    /* Trim trailing (highest power) near-zero coefficients */
    while (tf->num_order >= 0 && fabs(tf->num[tf->num_order]) < TF_EPS)
        tf->num_order--;
    while (tf->den_order >= 0 && fabs(tf->den[tf->den_order]) < TF_EPS)
        tf->den_order--;

    /* Trim leading (constant term) if BOTH have zero constant term */
    while (tf->num_order >= 0 && tf->den_order >= 0 &&
           fabs(tf->num[0]) < TF_EPS && fabs(tf->den[0]) < TF_EPS) {
        /* Shift arrays left by 1 (common factor s cancels) */
        int i;
        for (i = 0; i < tf->num_order; i++) tf->num[i] = tf->num[i + 1];
        for (i = 0; i < tf->den_order; i++) tf->den[i] = tf->den[i + 1];
        tf->num_order--;
        tf->den_order--;
    }

    /* Ensure orders are non-negative */
    if (tf->num_order < 0) tf->num_order = 0;
    if (tf->den_order < 0) tf->den_order = 0;

    /* Recompute stored orders based on actual non-zero content */
    while (tf->num_order > 0 && fabs(tf->num[tf->num_order]) < TF_EPS)
        tf->num_order--;
    while (tf->den_order > 0 && fabs(tf->den[tf->den_order]) < TF_EPS)
        tf->den_order--;
}

/* ================================================================
 * L5: Canonical Transformations
 * ================================================================ */

TransferFunction* tf_inverse(const TransferFunction *tf) {
    /* G_inv(s) = 1 / G(s) = D(s) / N(s).
     *
     * The inverse transfer function exchanges poles and zeros.
     * If original G(s) is strictly proper, G_inv(s) is improper
     * and represents a non-causal system (predictor/differentiator).
     *
     * Applications:
     *  - Feedforward control: G_inv * desired_output = required input
     *  - Model reference control: inverse plant model
     *  - Disturbance observer: nominal inverse */
    if (!tf) return NULL;

    TransferFunction *inv = tf_create(tf->den, tf->den_order,
                                       tf->num, tf->num_order);
    if (inv) {
        inv->delay = -tf->delay; /* time advance instead of delay */
        inv->is_discrete = tf->is_discrete;
        inv->sample_time = tf->sample_time;
    }
    return inv;
}

TransferFunction* tf_derivative(const TransferFunction *tf) {
    /* G_der(s) = s * G(s).
     *
     * Multiply numerator by s: s * N(s)/D(s) = (s*N(s)) / D(s).
     * Time domain: s * G(s) <-> d/dt g(t).
     *
     * This operation increases relative degree by -1 (more zeros than poles).
     * Applications: derivative kick compensation, lead filter design. */
    if (!tf) return NULL;

    /* Create s * N(s): shift numerator coefficients up by 1 */
    int new_no = tf->num_order + 1;
    double *new_num = (double*)calloc((size_t)(new_no + 1), sizeof(double));
    if (!new_num) return NULL;
    int i;
    for (i = 0; i <= tf->num_order; i++) new_num[i + 1] = tf->num[i];

    TransferFunction *der = tf_create(new_num, new_no, tf->den, tf->den_order);
    free(new_num);
    if (der) {
        der->delay = tf->delay;
        der->is_discrete = tf->is_discrete;
        der->sample_time = tf->sample_time;
    }
    return der;
}

TransferFunction* tf_integral(const TransferFunction *tf) {
    /* G_int(s) = G(s) / s.
     *
     * Multiply denominator by s: N(s) / (s * D(s)).
     * Time domain: G(s)/s <-> integral_0^t g(tau) dtau.
     *
     * This adds an integrator (pole at s=0), increasing system type by 1. */
    if (!tf) return NULL;

    int new_dno = tf->den_order + 1;
    double *new_den = (double*)calloc((size_t)(new_dno + 1), sizeof(double));
    if (!new_den) return NULL;
    int i;
    for (i = 0; i <= tf->den_order; i++) new_den[i + 1] = tf->den[i];

    TransferFunction *integ = tf_create(tf->num, tf->num_order,
                                         new_den, new_dno);
    free(new_den);
    if (integ) {
        integ->delay = tf->delay;
        integ->is_discrete = tf->is_discrete;
        integ->sample_time = tf->sample_time;
    }
    return integ;
}

TransferFunction* tf_frequency_shift(const TransferFunction *tf, double a) {
    /* G_shift(s) = G(s + a). Frequency/Laplace domain shift.
     *
     * This corresponds to modulation in time domain:
     *   L^{-1}[G(s+a)] = e^{-a*t} * g(t)
     *
     * Computed via polynomial substitution: replace s with (s+a) in both
     * numerator and denominator polynomials.
     *
     * Applications:
     *  - Shifting stability boundary (moving poles left)
     *  - Exponential weighting for convergence acceleration
     *  - Smith predictor analysis */
    if (!tf) return NULL;

    /* Create polynomial p(s) = s + a for substitution */
    double s_plus_a[2] = { a, 1.0 };
    Polynomial *shift = poly_create(s_plus_a, 1);
    Polynomial *np = poly_create(tf->num, tf->num_order);
    Polynomial *dp = poly_create(tf->den, tf->den_order);

    /* Compose: N(s+a) and D(s+a) */
    Polynomial *new_num = poly_compose(np, shift);
    Polynomial *new_den = poly_compose(dp, shift);

    TransferFunction *result = NULL;
    if (new_num && new_den) {
        result = tf_create(new_num->coeff, poly_degree(new_num),
                            new_den->coeff, poly_degree(new_den));
        if (result) {
            result->delay = tf->delay;
            result->is_discrete = tf->is_discrete;
            result->sample_time = tf->sample_time;
        }
    }

    poly_destroy(shift); poly_destroy(np); poly_destroy(dp);
    poly_destroy(new_num); poly_destroy(new_den);
    return result;
}

TransferFunction* tf_time_scale(const TransferFunction *tf, double alpha) {
    /* G_scaled(s) = G(alpha * s).
     *
     * Time domain: g(t/alpha)/alpha <-> G(alpha*s).
     *
     * Implementation: multiply each coefficient by alpha^power.
     *   If G(s) = N(s)/D(s) where N(s) = sum b_k s^k,
     *   then G(alpha*s) = sum b_k alpha^k s^k / sum a_k alpha^k s^k.
     *
     * For alpha > 1: speeds up time response (stretch frequency axis).
     * For alpha < 1: slows down time response (compress frequency axis). */
    if (!tf || alpha <= 0.0) return NULL;

    double *new_num = (double*)calloc((size_t)(tf->num_order + 1), sizeof(double));
    double *new_den = (double*)calloc((size_t)(tf->den_order + 1), sizeof(double));
    if (!new_num || !new_den) {
        free(new_num); free(new_den);
        return NULL;
    }

    int i;
    double pow_alpha = 1.0;
    for (i = 0; i <= tf->num_order; i++) {
        new_num[i] = tf->num[i] * pow_alpha;
        pow_alpha *= alpha;
    }
    pow_alpha = 1.0;
    for (i = 0; i <= tf->den_order; i++) {
        new_den[i] = tf->den[i] * pow_alpha;
        pow_alpha *= alpha;
    }

    TransferFunction *result = tf_create(new_num, tf->num_order,
                                          new_den, tf->den_order);
    free(new_num); free(new_den);
    if (result) {
        result->delay = tf->delay * alpha;
        result->is_discrete = tf->is_discrete;
        if (tf->is_discrete) result->sample_time = tf->sample_time / alpha;
    }
    return result;
}
