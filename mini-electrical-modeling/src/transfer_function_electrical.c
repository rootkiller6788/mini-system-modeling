/**
 * @file transfer_function_electrical.c
 * @brief Transfer functions of electrical networks — s-domain circuit analysis
 *
 * Independent knowledge points:
 *   - Transfer function evaluation via Horner method
 *   - DC gain computation
 *   - Pole/zero finding via Durand-Kerner method
 *   - Frequency response magnitude and phase computation
 *   - Bode asymptotic approximation
 *   - Cutoff frequency search via binary search
 *   - Bandwidth estimation from pole locations
 *   - Step response via partial fraction based inverse Laplace
 *   - Standard filter transfer functions (RC, RL, RLC, Sallen-Key, MFB)
 *   - Transfer function arithmetic (cascade, parallel, feedback)
 *   - Filter design with standard component values
 *   - Group delay computation
 */

#include "transfer_function_electrical.h"
#include <math.h>
#include <complex.h>
#include <string.h>
#include <float.h>

/* ===== Utility: polynomial evaluation via Horner method ===== */

static double complex poly_eval(const double *coeff, size_t order, double complex s)
{
    /* Horner method: (...(a_n*s + a_{n-1})*s + ... )*s + a_0 */
    /* More numerically stable than explicit power computation */
    if (order == 0 && coeff) return coeff[0];

    double complex result = 0.0;
    for (size_t i_int = order + 1; i_int > 0; i_int--) {
        size_t i = i_int - 1;
        result = result * s + coeff[i];
    }

    /* Handle case where order=0 properly */
    if (order == 0) result = coeff[0];
    return result;
}

/* ===== Evaluate TF at complex s ===== */

double complex tf_evaluate(const TransferFunction *tf, double complex s)
{
    if (!tf) return 0.0;
    double complex num_val = poly_eval(tf->num, tf->num_order, s);
    double complex den_val = poly_eval(tf->den, tf->den_order, s);

    double den_mag2 = creal(den_val) * creal(den_val) + cimag(den_val) * cimag(den_val);
    if (den_mag2 < 1e-30) {
        /* Pole at this s — return large value with correct sign */
        return DBL_MAX + I * 0.0;
    }
    return num_val / den_val;
}

/* ===== DC Gain ===== */

double tf_dc_gain(const TransferFunction *tf)
{
    /* H(0) = n0/d0 */
    if (!tf) return 0.0;
    if (fabs(tf->den[0]) < 1e-30) return DBL_MAX;
    return tf->num[0] / tf->den[0];
}

/* ===== Pole/Zero Finding: Durand-Kerner (Weierstrass) Method ===== */

/**
 * Durand-Kerner method for simultaneous polynomial root finding.
 * Uses Aberth correction for numerical stability.
 *
 * Reference: Durand (1960), Kerner (1966), Weierstrass (1891)
 * Modern implementation based on Bini (1996) for numerical stability.
 */
static int durand_kerner(const double *coeff, size_t order,
                          double complex *roots, int max_iter)
{
    if (order == 0 || !coeff || !roots) return 0;
    if (order > TF_MAX_ORDER) return -1;

    /* Leading coefficient */
    double a_n = coeff[order];
    if (fabs(a_n) < 1e-30) return -1;

    /* Initialize roots around unit circle with Chebyshev distribution */
    for (size_t k = 0; k < order; k++) {
        double angle = 2.0 * M_PI * (double)k / (double)order + 0.4;
        roots[k] = cos(angle) + I * sin(angle);
    }

    for (int iter = 0; iter < max_iter; iter++) {
        int converged = 1;
        for (size_t k = 0; k < order; k++) {
            /* Evaluate P(z_k) */
            double complex pz = coeff[order];
            for (size_t i_int = order; i_int > 0; i_int--) {
                size_t i = i_int - 1;
                pz = pz * roots[k] + coeff[i];
            }
            /* Also need leading coefficient contribution */
            double complex p_val = a_n;
            for (size_t i_int = order; i_int > 0; i_int--) {
                size_t i = i_int - 1;
                p_val = p_val * roots[k] + coeff[i];
            }

            /* Compute denominator: a_n * Π_{j≠k} (z_k - z_j) */
            double complex denom = a_n;
            for (size_t j = 0; j < order; j++) {
                if (j == k) continue;
                denom *= (roots[k] - roots[j]);
            }

            if (cabs(denom) < 1e-15) continue;

            /* Aberth correction */
            double complex delta = p_val / denom;
            roots[k] -= delta;

            if (cabs(delta) > 1e-12) converged = 0;
        }
        if (converged) break;
    }

    return (int)order;
}

int tf_find_poles(const TransferFunction *tf, double complex *poles)
{
    if (!tf || !poles) return 0;
    return durand_kerner(tf->den, tf->den_order, poles, 200);
}

int tf_find_zeros(const TransferFunction *tf, double complex *zeros)
{
    if (!tf || !zeros) return 0;
    return durand_kerner(tf->num, tf->num_order, zeros, 200);
}

/* ===== Frequency Response ===== */

double tf_magnitude_db(const TransferFunction *tf, double omega)
{
    /* 20 * log10(|H(jω)|) */
    if (!tf) return -200.0;
    double complex s = I * omega;
    double complex h = tf_evaluate(tf, s);
    double mag = cabs(h);
    if (mag < 1e-30) return -200.0;
    return 20.0 * log10(mag);
}

double tf_phase_rad(const TransferFunction *tf, double omega)
{
    /* ∠H(jω) in radians */
    if (!tf) return 0.0;
    double complex s = I * omega;
    double complex h = tf_evaluate(tf, s);
    return carg(h);
}

double tf_phase_deg(const TransferFunction *tf, double omega)
{
    return tf_phase_rad(tf, omega) * 180.0 / M_PI;
}

/* ===== Cutoff Frequency via Binary Search ===== */

double tf_find_cutoff_freq(const TransferFunction *tf, double f_low, double f_high)
{
    /* Find ω where |H(jω)| = |H(j0)|/√2 (=-3.01 dB from DC) */
    if (!tf) return -1.0;

    double dc_gain_val = tf_dc_gain(tf);
    if (fabs(dc_gain_val) < 1e-30) return -1.0;

    double target_mag = fabs(dc_gain_val) / sqrt(2.0);

    /* Check that crossing exists */
    double mag_low = cabs(tf_evaluate(tf, I * f_low));
    double mag_high = cabs(tf_evaluate(tf, I * f_high));

    if (mag_low < target_mag || mag_high > target_mag) return -1.0;

    double lo = f_low, hi = f_high;
    for (int iter = 0; iter < 80; iter++) {
        double mid = 0.5 * (lo + hi);
        double mag_mid = cabs(tf_evaluate(tf, I * mid));

        if (fabs(mag_mid - target_mag) / target_mag < 1e-9) {
            return mid;
        }

        if (mag_mid > target_mag) {
            lo = mid;
        } else {
            hi = mid;
        }

        if ((hi - lo) / fmax(hi, 1e-6) < 1e-12) {
            return 0.5 * (lo + hi);
        }
    }
    return 0.5 * (lo + hi);
}

/* ===== Bandwidth Estimate ===== */

double tf_bandwidth_estimate_from_poles(double natural_freq, double damping)
{
    /* For 2nd-order: BW = ω_n * sqrt(1-2ζ^2 + sqrt(2 - 4ζ^2 + 4ζ^4)) */
    /* Valid for ζ >= 0 */
    if (natural_freq <= 0.0) return 0.0;

    /* Clamp damping to meaningful range */
    if (damping < 0.0) damping = 0.0;

    double z2 = damping * damping;
    double inner = 1.0 - 2.0 * z2 + sqrt(fmax(0.0, 2.0 - 4.0 * z2 + 4.0 * z2 * z2));
    if (inner < 0.0) inner = 0.0;

    return natural_freq * sqrt(inner);
}

/* ===== Step Response ===== */

double tf_step_response(const TransferFunction *tf, double t)
{
    /* Step response y(t) = L^{-1}{H(s)/s} */
    /* Implementation: numerical inverse Laplace via series approximation */

    if (!tf || t < 0.0) return 0.0;

    /* For first-order: analytical solution */
    if (tf->den_order == 1) {
        /* H(s) = (b1*s + b0) / (a1*s + a0) */
        /* Step response: y(t) = b0/a0 + (b1*a0 - b0*a1)/(a0*a1) * e^{-(a0/a1)*t} */
        double a0 = tf->den[0], a1 = tf->den[1];
        double b0 = tf->num[0], b1 = (tf->num_order >= 1) ? tf->num[1] : 0.0;

        if (fabs(a1) < 1e-30) return 0.0;
        double tau = a1 / a0;
        if (tau < 0.0) return DBL_MAX; /* unstable */

        double dc = b0 / a0;
        double transient_coef = (b1 * a0 - b0 * a1) / (a0 * a1);

        return dc + transient_coef * exp(-t / tau);
    }

    /* For second-order: analytical solution */
    if (tf->den_order == 2) {
        double a0 = tf->den[0], a1 = tf->den[1], a2 = tf->den[2];
        double b0 = tf->num[0], b1 = (tf->num_order >= 1) ? tf->num[1] : 0.0;
        double b2 = (tf->num_order >= 2) ? tf->num[2] : 0.0;

        if (fabs(a2) < 1e-30) return 0.0;

        double omega_n = sqrt(fmax(0.0, a0 / a2));
        double zeta = a1 / (2.0 * a2 * omega_n);

        double dc = b0 / a0;

        if (zeta < 1.0) {
            /* Underdamped */
            double omega_d = omega_n * sqrt(1.0 - zeta * zeta);
            double sigma = zeta * omega_n;

            /* Residue computation at poles -σ ± jω_d */
            /* y(t) = DC + e^{-σt}[A·cos(ω_d·t) + B·sin(ω_d·t)] */
            double complex p1 = -sigma + I * omega_d;
            double complex r1 = (b2 * p1 * p1 + b1 * p1 + b0) / (a2 * 2.0 * I * omega_d);
            /* Conjugate residue */
            double A = 2.0 * creal(r1);
            double B = -2.0 * cimag(r1);

            double transient = exp(-sigma * t) * (A * cos(omega_d * t) + B * sin(omega_d * t));
            return dc + transient;
        } else if (fabs(zeta - 1.0) < 1e-10) {
            /* Critically damped */
            double sigma = omega_n;
            double numerator_at_pole = b2 * sigma * sigma - b1 * sigma + b0;
            double K1 = numerator_at_pole / (a2 * sigma * sigma);
            double K2 = (b1 - 2.0 * b2 * sigma) / (a2 * sigma);

            double transient = exp(-sigma * t) * (K1 + K2 * t);
            return dc + transient;
        } else {
            /* Overdamped */
            double s1 = omega_n * (-zeta + sqrt(zeta * zeta - 1.0));
            double s2 = omega_n * (-zeta - sqrt(zeta * zeta - 1.0));

            double K1 = (b2 * s1 * s1 + b1 * s1 + b0) / (a2 * (s1 - s2));
            double K2 = (b2 * s2 * s2 + b1 * s2 + b0) / (a2 * (s2 - s1));

            double transient = K1 * exp(s1 * t) + K2 * exp(s2 * t);
            return dc + transient;
        }
    }

    /* Higher order: numerical integration (simplified convolution) */
    /* Fall back to simple evaluation of inverse Laplace using dominant poles */
    double complex poles[TF_MAX_ORDER];
    int n_poles = tf_find_poles(tf, poles);
    if (n_poles <= 0) return 0.0;

    /* Sum of residues at poles: y(t) = Σ Res{H(s)/s·e^{st}} */
    double result = 0.0;
    for (int k = 0; k < n_poles; k++) {
        double complex pk = poles[k];
        /* Skip unstable poles */
        if (creal(pk) > 1e-9) return DBL_MAX;

        /* Compute residue at pole pk: H(pk)/pk * e^{pk·t} */
        double complex h_at_pk = tf_evaluate(tf, pk);
        if (cabs(pk) > 1e-12) {
            double complex residue = h_at_pk / pk;
            double complex term = residue * cexp(pk * t);
            result += creal(term);
        } else {
            /* Pole at origin: DC gain */
            result += tf_dc_gain(tf);
        }
    }
    return result;
}

/* ===== Standard Filter Transfer Functions ===== */

TransferFunction tf_rc_lowpass(double r, double c)
{
    /* H(s) = 1/(RCs + 1) = 1/(τs + 1) */
    TransferFunction tf;
    memset(&tf, 0, sizeof(tf));
    double tau = r * c;
    tf.num[0] = 1.0;
    tf.den[0] = 1.0;
    tf.den[1] = tau;
    tf.num_order = 0;
    tf.den_order = 1;
    return tf;
}

TransferFunction tf_rc_highpass(double r, double c)
{
    /* H(s) = RCs/(RCs + 1) = τs/(τs + 1) */
    TransferFunction tf;
    memset(&tf, 0, sizeof(tf));
    double tau = r * c;
    tf.num[1] = tau;
    tf.den[0] = 1.0;
    tf.den[1] = tau;
    tf.num_order = 1;
    tf.den_order = 1;
    return tf;
}

TransferFunction tf_rl_lowpass(double r, double l)
{
    /* H(s) = R/(Ls + R) = 1/((L/R)s + 1) */
    TransferFunction tf;
    memset(&tf, 0, sizeof(tf));
    double tau = (fabs(r) > 1e-30) ? l / r : DBL_MAX;
    tf.num[0] = 1.0;
    tf.den[0] = 1.0;
    tf.den[1] = tau;
    tf.num_order = 0;
    tf.den_order = 1;
    return tf;
}

TransferFunction tf_rl_highpass(double r, double l)
{
    /* H(s) = Ls/(Ls + R) */
    TransferFunction tf;
    memset(&tf, 0, sizeof(tf));
    double tau = (fabs(r) > 1e-30) ? l / r : DBL_MAX;
    tf.num[1] = tau;
    tf.den[0] = 1.0;
    tf.den[1] = tau;
    tf.num_order = 1;
    tf.den_order = 1;
    return tf;
}

TransferFunction tf_rlc_lowpass(double r, double l, double c)
{
    /* H(s) = 1/(LC·s^2 + RC·s + 1) */
    TransferFunction tf;
    memset(&tf, 0, sizeof(tf));
    tf.num[0] = 1.0;
    tf.den[0] = 1.0;
    tf.den[1] = r * c;
    tf.den[2] = l * c;
    tf.num_order = 0;
    tf.den_order = 2;
    return tf;
}

TransferFunction tf_rlc_bandpass(double r, double l, double c)
{
    /* H(s) = RCs/(LC·s^2 + RC·s + 1) */
    TransferFunction tf;
    memset(&tf, 0, sizeof(tf));
    tf.num[1] = r * c;
    tf.den[0] = 1.0;
    tf.den[1] = r * c;
    tf.den[2] = l * c;
    tf.num_order = 1;
    tf.den_order = 2;
    return tf;
}

TransferFunction tf_rlc_highpass(double r, double l, double c)
{
    /* H(s) = LC·s^2/(LC·s^2 + RC·s + 1) */
    TransferFunction tf;
    memset(&tf, 0, sizeof(tf));
    tf.num[2] = l * c;
    tf.den[0] = 1.0;
    tf.den[1] = r * c;
    tf.den[2] = l * c;
    tf.num_order = 2;
    tf.den_order = 2;
    return tf;
}

TransferFunction tf_sallen_key_lowpass(double r1, double r2, double c1, double c2)
{
    /* H(s) = 1/(R1·R2·C1·C2·s^2 + (R1+R2)·C2·s + 1) */
    TransferFunction tf;
    memset(&tf, 0, sizeof(tf));
    tf.num[0] = 1.0;
    tf.den[0] = 1.0;
    tf.den[1] = (r1 + r2) * c2;
    tf.den[2] = r1 * r2 * c1 * c2;
    tf.num_order = 0;
    tf.den_order = 2;
    return tf;
}

TransferFunction tf_mfb_lowpass(double r1, double r2, double r3, double c1, double c2)
{
    /* H(s) = -(R3/R1)/(R2·R3·C1·C2·s^2 + R3·(C1+C2)·s + 1) */
    TransferFunction tf;
    memset(&tf, 0, sizeof(tf));
    if (fabs(r1) < 1e-30) {
        tf.num[0] = 0.0;
        tf.den[0] = 1.0;
        tf.num_order = 0;
        tf.den_order = 0;
        return tf;
    }
    double gain = r3 / r1;
    tf.num[0] = -gain;
    tf.den[0] = 1.0;
    tf.den[1] = r3 * (c1 + c2);
    tf.den[2] = r2 * r3 * c1 * c2;
    tf.num_order = 0;
    tf.den_order = 2;
    return tf;
}

/* ===== Transfer Function Arithmetic ===== */

static void poly_multiply(const double *a, size_t na,
                           const double *b, size_t nb,
                           double *result, size_t *nr)
{
    /* Polynomial multiplication: convolution of coefficients */
    size_t n_result = na + nb;
    if (n_result > TF_MAX_ORDER) n_result = TF_MAX_ORDER;

    for (size_t k = 0; k <= n_result; k++) result[k] = 0.0;

    for (size_t i = 0; i <= na && i <= TF_MAX_ORDER; i++) {
        for (size_t j = 0; j <= nb && j <= TF_MAX_ORDER; j++) {
            if (i + j <= TF_MAX_ORDER) {
                result[i + j] += a[i] * b[j];
            }
        }
    }
    *nr = (n_result <= TF_MAX_ORDER) ? n_result : TF_MAX_ORDER;
}

static void poly_add(const double *a, size_t na,
                      const double *b, size_t nb,
                      double *result, size_t *nr)
{
    size_t max_order = (na > nb) ? na : nb;
    if (max_order > TF_MAX_ORDER) max_order = TF_MAX_ORDER;

    for (size_t k = 0; k <= max_order; k++) result[k] = 0.0;

    for (size_t i = 0; i <= na && i <= TF_MAX_ORDER; i++) result[i] += a[i];
    for (size_t i = 0; i <= nb && i <= TF_MAX_ORDER; i++) result[i] += b[i];

    /* Trim trailing zeros */
    *nr = max_order;
    while (*nr > 0 && fabs(result[*nr]) < 1e-30) (*nr)--;
}

TransferFunction tf_cascade(const TransferFunction *h1, const TransferFunction *h2)
{
    /* H(s) = H1(s) * H2(s) */
    /* N(s) = N1·N2, D(s) = D1·D2 */
    TransferFunction result;
    memset(&result, 0, sizeof(result));
    if (!h1 || !h2) return result;

    poly_multiply(h1->num, h1->num_order, h2->num, h2->num_order,
                  result.num, &result.num_order);
    poly_multiply(h1->den, h1->den_order, h2->den, h2->den_order,
                  result.den, &result.den_order);
    return result;
}

TransferFunction tf_parallel(const TransferFunction *h1, const TransferFunction *h2)
{
    /* H(s) = H1(s) + H2(s) = (N1·D2 + N2·D1)/(D1·D2) */
    TransferFunction result;
    memset(&result, 0, sizeof(result));
    if (!h1 || !h2) return result;

    /* N1 * D2 */
    double n1d2[2 * TF_MAX_ORDER + 1] = {0};
    size_t n1d2_order = 0;
    poly_multiply(h1->num, h1->num_order, h2->den, h2->den_order, n1d2, &n1d2_order);

    /* N2 * D1 */
    double n2d1[2 * TF_MAX_ORDER + 1] = {0};
    size_t n2d1_order = 0;
    poly_multiply(h2->num, h2->num_order, h1->den, h1->den_order, n2d1, &n2d1_order);

    /* Numerator = N1·D2 + N2·D1 */
    double num_tmp[2 * TF_MAX_ORDER + 1] = {0};
    size_t num_tmp_order = 0;
    poly_add(n1d2, n1d2_order, n2d1, n2d1_order, num_tmp, &num_tmp_order);
    for (size_t i = 0; i <= num_tmp_order && i <= TF_MAX_ORDER; i++)
        result.num[i] = num_tmp[i];
    result.num_order = num_tmp_order;

    /* Denominator = D1·D2 */
    poly_multiply(h1->den, h1->den_order, h2->den, h2->den_order,
                  result.den, &result.den_order);
    return result;
}

TransferFunction tf_negative_feedback(const TransferFunction *h1, const TransferFunction *h2)
{
    /* H(s) = H1(s) / (1 + H1(s)·H2(s)) */
    /* = N1·D2 / (D1·D2 + N1·N2) */
    TransferFunction result;
    memset(&result, 0, sizeof(result));
    if (!h1 || !h2) return result;

    /* Numerator: N1 * D2 */
    double num_tmp[2 * TF_MAX_ORDER + 1] = {0};
    size_t num_tmp_order = 0;
    poly_multiply(h1->num, h1->num_order, h2->den, h2->den_order, num_tmp, &num_tmp_order);
    for (size_t i = 0; i <= num_tmp_order && i <= TF_MAX_ORDER; i++)
        result.num[i] = num_tmp[i];
    result.num_order = num_tmp_order;

    /* D1 * D2 */
    double d1d2[2 * TF_MAX_ORDER + 1] = {0};
    size_t d1d2_order = 0;
    poly_multiply(h1->den, h1->den_order, h2->den, h2->den_order, d1d2, &d1d2_order);

    /* N1 * N2 */
    double n1n2[2 * TF_MAX_ORDER + 1] = {0};
    size_t n1n2_order = 0;
    poly_multiply(h1->num, h1->num_order, h2->num, h2->num_order, n1n2, &n1n2_order);

    /* Denominator: D1·D2 + N1·N2 */
    double den_tmp[2 * TF_MAX_ORDER + 1] = {0};
    size_t den_tmp_order = 0;
    poly_add(d1d2, d1d2_order, n1n2, n1n2_order, den_tmp, &den_tmp_order);
    for (size_t i = 0; i <= den_tmp_order && i <= TF_MAX_ORDER; i++)
        result.den[i] = den_tmp[i];
    result.den_order = den_tmp_order;
    return result;
}

/* ===== Asymptotic Bode Magnitude ===== */

double tf_bode_magnitude_asymptotic(const TransferFunction *tf, double omega)
{
    /* Piecewise linear approximation:
     *   - Constant at DC
     *   - Each pole adds -20dB/dec after its corner frequency
     *   - Each zero adds +20dB/dec after its corner frequency
     */
    if (!tf) return -200.0;

    double dc_db = 20.0 * log10(fmax(fabs(tf->num[0] / fmax(fabs(tf->den[0]), 1e-30)), 1e-30));
    if (omega <= 1e-30) return dc_db;

    /* Get corner frequencies from poles and zeros */
    double complex poles[TF_MAX_ORDER];
    double complex zeros[TF_MAX_ORDER];
    int n_poles = tf_find_poles(tf, poles);
    int n_zeros = tf_find_zeros(tf, zeros);

    double db = dc_db;

    /* Contribution from zeros: +20dB/dec per zero for ω > |z_i| */
    for (int i = 0; i < n_zeros; i++) {
        double corner = cabs(zeros[i]);
        if (omega > corner && corner > 1e-30) {
            db += 20.0 * log10(omega / corner);
        }
    }

    /* Contribution from poles: -20dB/dec per pole for ω > |p_i| */
    for (int i = 0; i < n_poles; i++) {
        double corner = cabs(poles[i]);
        if (omega > corner && corner > 1e-30) {
            db -= 20.0 * log10(omega / corner);
        }
    }

    return db;
}

/* ===== Filter Design with Standard Component Values ===== */

/* E12 standard resistor values (1Ω to 10MΩ) */
static const double e12_values[] = {
    1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2
};
static const double e12_decades[] = {
    1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7
};

/* E6 standard capacitor values (1pF to 1000µF) */
static const double e6_cap_values[] = {
    1.0, 1.5, 2.2, 3.3, 4.7, 6.8
};
static const double e6_cap_decades[] = {
    1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4
};

FilterDesign design_rc_lowpass_to_cutoff(double cutoff_hz)
{
    /* Choose C first (typically 10nF to 1µF for audio, 1nF for higher freq) */
    /* Then compute R = 1/(2π·fc·C) */
    FilterDesign result;
    memset(&result, 0, sizeof(result));

    if (cutoff_hz <= 0.0) return result;

    /* Select capacitor from E6 series */
    double target_rc = 1.0 / (2.0 * M_PI * cutoff_hz);

    /* Pick a reasonable capacitor */
    double best_c = 1e-9; /* Start with 1nF */
    double best_r = target_rc / best_c;

    /* Try different capacitor values */
    for (int di = 0; di < (int)(sizeof(e6_cap_decades)/sizeof(double)); di++) {
        for (int vi = 0; vi < (int)(sizeof(e6_cap_values)/sizeof(double)); vi++) {
            double c = e6_cap_values[vi] * e6_cap_decades[di];
            double r = target_rc / c;

            /* R should be in reasonable range (1k to 1M) */
            if (r >= 1e3 && r <= 1e6) {
                /* Find nearest standard R */
                double best_std_r = 1e3;
                double best_diff = DBL_MAX;
                for (int dd = 0; dd < (int)(sizeof(e12_decades)/sizeof(double)); dd++) {
                    for (int dv = 0; dv < (int)(sizeof(e12_values)/sizeof(double)); dv++) {
                        double std_r = e12_values[dv] * e12_decades[dd];
                        double diff = fabs(std_r - r) / r;
                        if (diff < best_diff) {
                            best_diff = diff;
                            best_std_r = std_r;
                        }
                    }
                }
                best_c = c;
                best_r = best_std_r;
                /* Accept if the resulting cutoff is close */
                break;
            }
        }
    }

    result.c_standard = best_c;
    result.r_standard = best_r;
    result.actual_cutoff = 1.0 / (2.0 * M_PI * best_r * best_c);
    if (cutoff_hz > 1e-30)
        result.error_percent = 100.0 * (result.actual_cutoff - cutoff_hz) / cutoff_hz;
    return result;
}

/* ===== Group Delay ===== */

double tf_group_delay(const TransferFunction *tf, double omega, double delta)
{
    /* τ_g(ω) = -dφ/dω ≈ -(φ(ω+δ) - φ(ω-δ))/(2δ) */
    if (!tf) return 0.0;
    if (delta <= 0.0) delta = omega * 1e-6 + 1e-9; /* small perturbation */

    double phi_plus  = tf_phase_rad(tf, omega + delta);
    double phi_minus = tf_phase_rad(tf, omega - delta);

    /* Handle phase wrapping */
    double dphi = phi_plus - phi_minus;
    /* Unwrap if discontinuity > π */
    if (dphi > M_PI)  dphi -= 2.0 * M_PI;
    if (dphi < -M_PI) dphi += 2.0 * M_PI;

    return -dphi / (2.0 * delta);
}
