/**
 * @file mechanical_transfer_function.c
 * @brief Implementation of mechanical transfer functions
 * Each function implements one independent knowledge point.
 */

#include "mechanical_transfer_function.h"
#include <math.h>
#include <complex.h>
#include <float.h>
#include <string.h>

/* Helper: set TF coefficients from arrays */
static void tf_set(MechanicalTF *tf, int n_num, int n_den,
                   const double complex *num, const double complex *den)
{
    tf->num_order = n_num;
    tf->den_order = n_den;
    tf->is_strictly_proper = (n_num < n_den) ? 1 : 0;
    memset(tf->num, 0, sizeof(tf->num));
    memset(tf->den, 0, sizeof(tf->den));
    for (int i = 0; i <= n_num; i++) tf->num[i] = num[i];
    for (int i = 0; i <= n_den; i++) tf->den[i] = den[i];
    /* DC gain */
    if (fabs(cabs(tf->den[n_den])) > 1e-30)
        tf->dc_gain = creal(tf->num[n_num] / tf->den[n_den]);
    else
        tf->dc_gain = 0.0;
}

/* ===== SDOF Transfer Functions ===== */

MechanicalTF tf_sdof_displacement(double mass, double damping, double stiffness)
{
    /* H(s) = X(s)/F(s) = 1 / (m*s^2 + b*s + k)
     * Polynomial convention: coeffs[0]*s^n + ... + coeffs[n] */
    MechanicalTF tf = {0};
    double complex num[] = {1.0};            /* constant 1 */
    double complex den[] = {mass, damping, stiffness}; /* m*s^2 + b*s + k */
    tf_set(&tf, 0, 2, num, den);
    return tf;
}

MechanicalTF tf_sdof_velocity(double mass, double damping, double stiffness)
{
    /* H(s) = s / (m*s^2 + b*s + k) */
    MechanicalTF tf = {0};
    double complex num[] = {1.0, 0.0};       /* 1*s + 0 */
    double complex den[] = {mass, damping, stiffness};
    tf_set(&tf, 1, 2, num, den);
    return tf;
}

MechanicalTF tf_sdof_acceleration(double mass, double damping, double stiffness)
{
    /* H(s) = s^2 / (m*s^2 + b*s + k) */
    MechanicalTF tf = {0};
    double complex num[] = {1.0, 0.0, 0.0};  /* 1*s^2 + 0*s + 0 */
    double complex den[] = {mass, damping, stiffness};
    tf_set(&tf, 2, 2, num, den);
    return tf;
}

MechanicalTF tf_base_excitation_displacement(double mass, double damping, double stiffness)
{
    /* H(s) = (b*s + k) / (m*s^2 + b*s + k) */
    MechanicalTF tf = {0};
    double complex num[] = {damping, stiffness};        /* b*s + k */
    double complex den[] = {mass, damping, stiffness};
    tf_set(&tf, 1, 2, num, den);
    return tf;
}

MechanicalTF tf_base_excitation_force(double mass, double damping, double stiffness)
{
    /* H(s) = F_t(s)/Y(s) = s*(b*s + k) / (m*s^2 + b*s + k)
     * Transmitted force to base per unit base displacement.
     */
    MechanicalTF tf = {0};
    double complex num[] = {damping, stiffness, 0.0};
    double complex den[] = {mass, damping, stiffness};
    tf_set(&tf, 2, 2, num, den);
    return tf;
}

MechanicalTF tf_rotating_unbalance_displacement(double total_mass, double damping, double stiffness)
{
    /* H(s) = X(s) / F_unbal = 1 / (M_total*s^2 + b*s + k)
     * Same form as displacement TF with total mass.
     */
    return tf_sdof_displacement(total_mass, damping, stiffness);
}

/* ===== 2-DOF Transfer Functions ===== */

MechanicalTF tf_2dof_direct(double m1, double m2, double k1, double k2, double kc)
{
    /* H11(s) = X1(s)/F1(s) for undamped 2-DOF system.
     * H11 = (m2*s^2 + k2 + kc) / [(m1*s^2+k1+kc)*(m2*s^2+k2+kc) - kc^2]
     * Denominator: m1*m2*s^4 + (m1*k2+m2*k1+m1*kc+m2*kc)*s^2 + (k1*k2+k1*kc+k2*kc)
     */
    MechanicalTF tf = {0};
    double complex num[] = {m2, 0.0, k2 + kc};
    double complex den[] = {m1*m2, 0.0,
        m1*k2 + m2*k1 + m1*kc + m2*kc, 0.0,
        k1*k2 + k1*kc + k2*kc};
    tf_set(&tf, 2, 4, num, den);
    return tf;
}

MechanicalTF tf_2dof_cross(double m1, double m2, double k1, double k2, double kc)
{
    /* H12(s) = X2(s)/F1(s) = kc / [denominator]
     * Cross transfer function: response of mass 2 due to force on mass 1.
     */
    MechanicalTF tf = {0};
    double complex num[] = {kc}; /* constant */
    double complex den[] = {m1*m2, 0.0,
        m1*k2 + m2*k1 + m1*kc + m2*kc, 0.0,
        k1*k2 + k1*kc + k2*kc};
    tf_set(&tf, 0, 4, num, den);
    return tf;
}

MechanicalTF tf_2dof_damped_direct(double m1, double m2, double k1, double k2,
                                    double kc, double b1, double b2, double bc)
{
    /* Full damped 2-DOF: H11(s).
     * Den: (m1*s^2+(b1+bc)*s+(k1+kc))*(m2*s^2+(b2+bc)*s+(k2+kc)) - (bc*s+kc)^2
     * = m1*m2*s^4 + (m1*(b2+bc)+m2*(b1+bc))*s^3
     *   + (m1*(k2+kc)+m2*(k1+kc)+(b1+bc)*(b2+bc)-bc^2)*s^2
     *   + ((b1+bc)*(k2+kc)+(b2+bc)*(k1+kc)-2*bc*kc)*s
     *   + ((k1+kc)*(k2+kc)-kc^2)
     * Num: m2*s^2 + (b2+bc)*s + (k2+kc)
     */
    MechanicalTF tf = {0};
    double complex num[] = {m2, b2+bc, k2+kc};
    double complex den[] = {
        m1*m2,
        m1*(b2+bc) + m2*(b1+bc),
        m1*(k2+kc) + m2*(k1+kc) + (b1+bc)*(b2+bc) - bc*bc,
        (b1+bc)*(k2+kc) + (b2+bc)*(k1+kc) - 2.0*bc*kc,
        (k1+kc)*(k2+kc) - kc*kc
    };
    tf_set(&tf, 2, 4, num, den);
    return tf;
}

MechanicalTF tf_vibration_absorber(double M, double ma, double K, double ka,
                                    double B, double ba)
{
    /* Dynamic vibration absorber: main mass M + absorber mass ma.
     * Tuned when ka/ma = K/M (absorber natural freq = main natural freq).
     */
    MechanicalTF tf = {0};
    double complex num[] = {ma, ba, ka};
    double complex den[] = {
        M*ma,
        M*ba + B*ma + B*ba,
        M*ka + K*ma + B*ba + ka*ma + B*B,
        B*ka + K*ba,
        K*ka
    };
    tf_set(&tf, 2, 4, num, den);
    return tf;
}

/* ===== Special Transfer Functions ===== */

MechanicalTF tf_seismic_accelerometer(double mass, double damping, double stiffness)
{
    /* H(s) = X(s)/A_base(s) = 1 / (s^2 + (b/m)*s + k/m)
     * Normalized accelerometer: X*m/(m*a_base) = 1/(s^2 + 2*zeta*wn*s + wn^2).
     * For large wn, H ~ 1/wn^2 (displacement proportional to base acceleration).
     */
    MechanicalTF tf = {0};
    double complex num[] = {1.0}; /* constant */
    double complex den[] = {1.0, damping/mass, stiffness/mass};
    tf_set(&tf, 0, 2, num, den);
    return tf;
}

MechanicalTF tf_seismic_velocimeter(double mass, double damping, double stiffness)
{
    /* H(s) = V(s)/V_base(s) — relative velocity output = s/(s^2+(b/m)*s+k/m) */
    MechanicalTF tf = {0};
    double complex num[] = {1.0, 0.0}; /* 1*s + 0 */
    double complex den[] = {1.0, damping/mass, stiffness/mass};
    tf_set(&tf, 1, 2, num, den);
    return tf;
}

MechanicalTF tf_vibration_isolator(double mass, double damping, double stiffness)
{
    /* Identical to base excitation displacement TF.
     * H(s) = (b*s + k) / (m*s^2 + b*s + k)
     */
    return tf_base_excitation_displacement(mass, damping, stiffness);
}

/* ===== TF Operations ===== */

double complex tf_eval(const MechanicalTF *tf, double complex s)
{
    if (!tf) return 0.0;
    double complex num_val = poly_eval(tf->num, tf->num_order, s);
    double complex den_val = poly_eval(tf->den, tf->den_order, s);
    double mag2 = creal(den_val)*creal(den_val) + cimag(den_val)*cimag(den_val);
    if (mag2 < 1e-60) return DBL_MAX + I * 0.0;
    return num_val / den_val;
}

double tf_mag(const MechanicalTF *tf, double angular_freq)
{
    return cabs(tf_eval(tf, I * angular_freq));
}

double tf_mag_db(const MechanicalTF *tf, double angular_freq)
{
    double mag = tf_mag(tf, angular_freq);
    if (mag < 1e-30) return -200.0;
    return 20.0 * log10(mag);
}

double tf_phase(const MechanicalTF *tf, double angular_freq)
{
    double complex val = tf_eval(tf, I * angular_freq);
    return atan2(cimag(val), creal(val));
}

double complex poly_eval(const double complex *coeffs, int order, double complex s)
{
    /* Horner's method: a0*s^n + a1*s^{n-1} + ... + an = ((a0*s + a1)*s + ... )*s + an
     * Complexity: O(n).
     */
    if (!coeffs || order < 0) return 0.0;
    double complex result = coeffs[0];
    for (int i = 1; i <= order; i++) {
        result = result * s + coeffs[i];
    }
    return result;
}

/* Polynomial convolution for TF cascade/feedback */
static void poly_conv(const double complex *a, int na, const double complex *b, int nb,
                       double complex *c)
{
    /* c = a * b, c[k] = sum_{i+j=k} a[i]*b[j] */
    for (int k = 0; k <= na + nb; k++) c[k] = 0.0;
    for (int i = 0; i <= na; i++) {
        for (int j = 0; j <= nb; j++) {
            c[i+j] += a[i] * b[j];
        }
    }
}

MechanicalTF tf_cascade(const MechanicalTF *tf1, const MechanicalTF *tf2)
{
    /* H(s) = H1(s) * H2(s) */
    MechanicalTF result = {0};
    if (!tf1 || !tf2) return result;

    double complex num[TF_MAX_ORDER*2+1], den[TF_MAX_ORDER*2+1];
    poly_conv(tf1->num, tf1->num_order, tf2->num, tf2->num_order, num);
    poly_conv(tf1->den, tf1->den_order, tf2->den, tf2->den_order, den);

    int n_num = tf1->num_order + tf2->num_order;
    int n_den = tf1->den_order + tf2->den_order;

    tf_set(&result, n_num, n_den, num, den);
    return result;
}

MechanicalTF tf_parallel(const MechanicalTF *tf1, const MechanicalTF *tf2)
{
    /* H(s) = H1(s) + H2(s) = (N1*D2 + N2*D1) / (D1*D2) */
    MechanicalTF result = {0};
    if (!tf1 || !tf2) return result;

    double complex N1D2[TF_MAX_ORDER*2+1], N2D1[TF_MAX_ORDER*2+1];
    double complex num[TF_MAX_ORDER*2+1], den[TF_MAX_ORDER*2+1];

    poly_conv(tf1->num, tf1->num_order, tf2->den, tf2->den_order, N1D2);
    poly_conv(tf2->num, tf2->num_order, tf1->den, tf1->den_order, N2D1);
    poly_conv(tf1->den, tf1->den_order, tf2->den, tf2->den_order, den);

    int n_num = (tf1->num_order+tf2->den_order > tf2->num_order+tf1->den_order)
                ? tf1->num_order+tf2->den_order : tf2->num_order+tf1->den_order;
    int n_den = tf1->den_order + tf2->den_order;

    /* Add polynomials N1D2 + N2D1 */
    for (int i = n_num; i >= 0; i--) {
        double complex n1 = (i <= tf1->num_order+tf2->den_order) ? N1D2[i] : 0.0;
        double complex n2 = (i <= tf2->num_order+tf1->den_order) ? N2D1[i] : 0.0;
        num[i] = n1 + n2;
    }

    tf_set(&result, n_num, n_den, num, den);
    return result;
}

MechanicalTF tf_neg_feedback(const MechanicalTF *tf1, const MechanicalTF *tf2)
{
    /* H(s) = H1 / (1 + H1*H2) = N1*D2 / (D1*D2 + N1*N2) */
    MechanicalTF result = {0};
    if (!tf1 || !tf2) return result;

    double complex N1D2[TF_MAX_ORDER*2+1];
    double complex N1N2[TF_MAX_ORDER*2+1], D1D2[TF_MAX_ORDER*2+1];
    double complex num[TF_MAX_ORDER*2+1], den[TF_MAX_ORDER*2+1];

    poly_conv(tf1->num, tf1->num_order, tf2->den, tf2->den_order, N1D2);
    poly_conv(tf1->num, tf1->num_order, tf2->num, tf2->num_order, N1N2);
    poly_conv(tf1->den, tf1->den_order, tf2->den, tf2->den_order, D1D2);

    int n_num = tf1->num_order + tf2->den_order;
    int n_den_max = (tf1->den_order+tf2->den_order > tf1->num_order+tf2->num_order)
                     ? tf1->den_order+tf2->den_order : tf1->num_order+tf2->num_order;

    for (int i = 0; i <= n_num; i++) num[i] = N1D2[i];
    for (int i = 0; i <= n_den_max; i++) {
        double complex d1d2_v = D1D2[i];
        double complex n1n2_v = N1N2[i];
        den[i] = d1d2_v + n1n2_v;
    }

    tf_set(&result, n_num, n_den_max, num, den);
    return result;
}

MechanicalTF tf_pos_feedback(const MechanicalTF *tf1, const MechanicalTF *tf2)
{
    /* H(s) = H1 / (1 - H1*H2) = N1*D2 / (D1*D2 - N1*N2) */
    MechanicalTF result = {0};
    if (!tf1 || !tf2) return result;

    double complex N1D2[TF_MAX_ORDER*2+1];
    double complex N1N2[TF_MAX_ORDER*2+1], D1D2[TF_MAX_ORDER*2+1];
    double complex num_pf[TF_MAX_ORDER*2+1], den_pf[TF_MAX_ORDER*2+1];

    poly_conv(tf1->num, tf1->num_order, tf2->den, tf2->den_order, N1D2);
    poly_conv(tf1->num, tf1->num_order, tf2->num, tf2->num_order, N1N2);
    poly_conv(tf1->den, tf1->den_order, tf2->den, tf2->den_order, D1D2);

    int n_num = tf1->num_order + tf2->den_order;
    int n_den_max = (tf1->den_order+tf2->den_order > tf1->num_order+tf2->num_order)
                     ? tf1->den_order+tf2->den_order : tf1->num_order+tf2->num_order;

    for (int i = 0; i <= n_num; i++) num_pf[i] = N1D2[i];
    for (int i = 0; i <= n_den_max; i++) den_pf[i] = D1D2[i] - N1N2[i];

    tf_set(&result, n_num, n_den_max, num_pf, den_pf);
    return result;
}

/* ===== Frequency Response Analysis ===== */

void tf_sdof_poles(double m, double b, double k, double complex *p1, double complex *p2)
{
    /* s^2 + (b/m)*s + (k/m) = 0
     * s = -b/(2*m) +/- sqrt((b/(2*m))^2 - k/m)
     */
    if (m <= 0.0) { if(p1)*p1=0; if(p2)*p2=0; return; }
    double alpha = b / (2.0 * m);
    double omega_n2 = k / m;
    double disc = alpha*alpha - omega_n2;
    if (disc >= 0) {
        double sqrt_disc = sqrt(disc);
        if(p1)*p1 = -alpha + sqrt_disc;
        if(p2)*p2 = -alpha - sqrt_disc;
    } else {
        double omega_d = sqrt(-disc);
        if(p1)*p1 = -alpha + I * omega_d;
        if(p2)*p2 = -alpha - I * omega_d;
    }
}

int tf_sdof_is_overdamped(double m, double b, double k)
{
    if (m <= 0.0 || k < 0.0) return 0;
    double zeta = b / (2.0 * sqrt(k*m));
    return (zeta > 1.0) ? 1 : 0;
}

int tf_sdof_is_critically_damped(double m, double b, double k)
{
    if (m <= 0.0 || k < 0.0) return 0;
    double zeta = b / (2.0 * sqrt(k*m));
    return (fabs(zeta - 1.0) < 1e-9) ? 1 : 0;
}

int tf_sdof_is_underdamped(double m, double b, double k)
{
    if (m <= 0.0 || k < 0.0) return 0;
    double zeta = b / (2.0 * sqrt(k*m));
    return (zeta < 1.0 && zeta >= 0.0) ? 1 : 0;
}

double tf_steady_state_amp(const MechanicalTF *tf, double force_amp, double w)
{
    return force_amp * tf_mag(tf, w);
}

double tf_steady_state_phase(const MechanicalTF *tf, double w)
{
    return tf_phase(tf, w);
}

int tf_is_proper(const MechanicalTF *tf)
{
    if (!tf) return 0;
    return (tf->den_order >= tf->num_order) ? 1 : 0;
}

double tf_estimate_damping(const MechanicalTF *tf)
{
    /* Extract damping from dominant complex pole pair */
    if (!tf || tf->den_order < 2) return 0.0;
    /* For SDOF case: evaluate denominator at s=jw form */
    if (tf->den_order == 2) {
        double m2 = creal(tf->den[0]);
        double b_coeff = creal(tf->den[1]);
        double k_coeff = creal(tf->den[2]);
        if (m2 <= 0.0) return 0.0;
        return b_coeff / (2.0 * sqrt(k_coeff * m2));
    }
    return 0.0;
}

double tf_estimate_natural_freq(const MechanicalTF *tf)
{
    if (!tf || tf->den_order < 2) return 0.0;
    if (tf->den_order == 2) {
        double m2 = creal(tf->den[0]);
        double k_coeff = creal(tf->den[2]);
        if (m2 <= 0.0) return 0.0;
        return sqrt(k_coeff / m2);
    }
    return 0.0;
}

double tf_dc_gain(const MechanicalTF *tf)
{
    /* H(0) = num[den_order] / den[den_order] (constant term ratio)
     * With coeffs organized as [a_n*s^n + ... + a_0],
     * the constant term is the last coefficient.
     */
    if (!tf) return 0.0;
    double complex den0 = tf->den[tf->den_order];
    if (cabs(den0) < 1e-30) return DBL_MAX;
    /* Transfer function written in s-domain convention:
     * num[0]*s^{n} + ... + num[n_order]
     * DC gain = eval at s=0 => last coefficient ratio */
    double complex num_dc = (tf->num_order <= tf->den_order)
        ? tf->num[tf->num_order] : DBL_MAX;
    return creal(num_dc / den0);
}

double tf_step_response(const MechanicalTF *tf, double t)
{
    /* Approximate step response using SDOF model parameters */
    if (!tf) return 0.0;
    double wn = tf_estimate_natural_freq(tf);
    double zeta = tf_estimate_damping(tf);
    double dc = tf_dc_gain(tf);
    double delta_st = dc; /* DC gain = static deflection for displacement TF */

    double wn_val = wn;
    double zeta_val = zeta;
    return delta_st * (1.0 - exp(-zeta_val*wn_val*t) *
           (cos(wn_val*sqrt(1.0-zeta_val*zeta_val)*t) +
            (zeta_val/sqrt(1.0-zeta_val*zeta_val))*sin(wn_val*sqrt(1.0-zeta_val*zeta_val)*t)));
}

/* ===== Modal TF Decomposition ===== */

void tf_2dof_modal_decomposition(double m1, double m2, double k1, double k2, double kc,
                                  double phi[4], double omega[2], double zeta[2],
                                  MechanicalTF modal_tf[2])
{
    /* 2-DOF undamped modal decomposition.
     * M = [[m1,0],[0,m2]], K = [[k1+kc,-kc],[-kc,k2+kc]]
     * Solve K*phi = omega^2 * M * phi.
     */
    double a = m1*m2;
    double b = -(m1*(k2+kc) + m2*(k1+kc));
    double c = (k1+kc)*(k2+kc) - kc*kc;

    double disc = b*b - 4.0*a*c;
    if (disc < 0) disc = 0.0;

    double lambda1 = (-b - sqrt(disc)) / (2.0*a);
    double lambda2 = (-b + sqrt(disc)) / (2.0*a);

    omega[0] = sqrt(lambda1);
    omega[1] = sqrt(lambda2);

    /* Mode shape 1: (K - lambda1*M)*phi1 = 0 */
    double k11 = k1 + kc;
    double k22 = k2 + kc;
    double k12 = -kc;

    if (fabs(k11 - lambda1*m1) > 1e-9) {
        phi[0] = 1.0;
        phi[1] = -k12 / (k11 - lambda1*m1);
    } else {
        phi[0] = -k12 / (k22 - lambda1*m2);
        phi[1] = 1.0;
    }
    /* Normalize */
    double norm1 = sqrt(phi[0]*phi[0] + phi[1]*phi[1]);
    phi[0] /= norm1;
    phi[1] /= norm1;

    if (fabs(k11 - lambda2*m1) > 1e-9) {
        phi[2] = 1.0;
        phi[3] = -k12 / (k11 - lambda2*m1);
    } else {
        phi[2] = -k12 / (k22 - lambda2*m2);
        phi[3] = 1.0;
    }
    double norm2 = sqrt(phi[2]*phi[2] + phi[3]*phi[3]);
    phi[2] /= norm2;
    phi[3] /= norm2;

    /* Undamped => zeta = 0 */
    zeta[0] = 0.0;
    zeta[1] = 0.0;

    /* Build modal TFs as SDOF systems */
    modal_tf[0] = tf_sdof_displacement(1.0, 0.0, omega[0]*omega[0]);
    modal_tf[1] = tf_sdof_displacement(1.0, 0.0, omega[1]*omega[1]);
}

/* ===== Frequency Response Data ===== */

void tf_freq_response(const MechanicalTF *tf, double w_start, double w_end,
                      int n_pts, FreqResp *fr)
{
    if (!tf || !fr || n_pts <= 0) return;
    fr->n_points = n_pts;
    fr->freq = (double*)malloc(n_pts * sizeof(double));
    fr->mag_db = (double*)malloc(n_pts * sizeof(double));
    fr->phase_deg = (double*)malloc(n_pts * sizeof(double));
    if (!fr->freq || !fr->mag_db || !fr->phase_deg) { fr->n_points = 0; return; }

    double log_start = log10(w_start);
    double log_end = log10(w_end);
    for (int i = 0; i < n_pts; i++) {
        double log_w = log_start + (log_end - log_start) * (double)i / (double)(n_pts - 1);
        fr->freq[i] = pow(10.0, log_w);
        fr->mag_db[i] = tf_mag_db(tf, fr->freq[i]);
        fr->phase_deg[i] = tf_phase(tf, fr->freq[i]) * 180.0 / M_PI;
    }
}

void freq_resp_free(FreqResp *fr)
{
    if (!fr) return;
    free(fr->freq);
    free(fr->mag_db);
    free(fr->phase_deg);
    fr->n_points = 0;
}
