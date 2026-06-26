/**
 * @file tf_analysis.c
 * @brief Transfer Function Analysis: Stability, Frequency, Time Response
 *
 * Knowledge points: L4 Routh-Hurwitz, FVT, IVT, Gain/Phase margins
 * L5: Frequency response, time response (RK4), step characteristics
 * L6: Partial fraction expansion, dominant pole
 *
 * Ref: Ogata Ch.5,8; Franklin Ch.3,6; Bode (1945); Nyquist (1932)
 */

#include "tf_analysis.h"
#include "tf_conversion.h"
#include "tf_polynomial.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TF_EPS 1e-15

/* ================================================================
 * L4: Routh-Hurwitz Stability Criterion
 * ================================================================ */

tf_stability_t tf_stability_routh(const TransferFunction *tf) {
    if (!tf || tf->den_order < 1) return TF_STABLE;
    int n = tf_num_poles(tf);
    if (n < 1) return TF_STABLE;
    int rows = n + 1, cols = (n + 2) / 2 + 1;
    double **R = (double**)calloc((size_t)rows, sizeof(double*));
    if (!R) return TF_STABLE;
    int i, j;
    for (i = 0; i < rows; i++) {
        R[i] = (double*)calloc((size_t)cols, sizeof(double));
        if (!R[i]) { int k; for (k = 0; k < i; k++) free(R[k]); free(R); return TF_STABLE; }
    }
    int dno = tf->den_order;
    for (j = 0; j < cols; j++) {
        int idx0 = dno - 2 * j, idx1 = dno - (2 * j + 1);
        R[0][j] = (idx0 >= 0 && idx0 <= dno) ? tf->den[idx0] : 0.0;
        R[1][j] = (idx1 >= 0 && idx1 <= dno) ? tf->den[idx1] : 0.0;
    }
    int sign_changes = 0;
    double prev = (fabs(R[0][0]) > TF_EPS) ? R[0][0] : 1e-10;
    for (i = 2; i < rows; i++) {
        for (j = 0; j < cols - 1; j++) {
            double det = R[i-1][0] * R[i-2][j+1] - R[i-2][0] * R[i-1][j+1];
            R[i][j] = (fabs(R[i-1][0]) < TF_EPS) ? R[i-2][j+1] : det / R[i-1][0];
        }
        double cur = R[i][0];
        if (fabs(cur) < TF_EPS) cur = (prev > 0) ? 1e-10 : -1e-10;
        if (prev * cur < 0) sign_changes++;
        if (fabs(cur) > TF_EPS) prev = cur;
    }
    if (R[0][0] * R[1][0] < 0) sign_changes++;
    for (i = 0; i < rows; i++) free(R[i]);
    free(R);
    return (sign_changes > 0) ? TF_UNSTABLE : TF_STABLE;
}

tf_stability_t tf_stability_poles(const TransferFunction *tf) {
    if (!tf) return TF_STABLE;
    int np = tf_num_poles(tf);
    if (np < 1) return TF_STABLE;
    double *rr = (double*)malloc((size_t)(np + 1) * sizeof(double));
    double *ri = (double*)malloc((size_t)(np + 1) * sizeof(double));
    if (!rr || !ri) { free(rr); free(ri); return TF_STABLE; }
    int n = poly_find_roots_complex(tf->den, tf->den_order, rr, ri);
    int rhp = 0, jw = 0, i;
    for (i = 0; i < n; i++) {
        if (rr[i] > TF_EPS) rhp++;
        else if (fabs(rr[i]) < TF_EPS) jw++;
    }
    free(rr); free(ri);
    if (rhp > 0) return TF_UNSTABLE;
    if (jw > 0) return TF_MARGINALLY_STABLE;
    return TF_STABLE;
}

int tf_stability_routh_array(const double *den, int order, double **array_out) {
    if (!den || order < 1 || !array_out) return 0;
    int rows = order + 1, cols = (order + 2) / 2 + 1;
    double *arr = (double*)calloc((size_t)(rows * cols), sizeof(double));
    if (!arr) return -1;
    int j;
    for (j = 0; j < cols; j++) {
        int i0 = order - 2 * j, i1 = order - (2 * j + 1);
        arr[j] = (i0 >= 0 && i0 <= order) ? den[i0] : 0.0;
        arr[cols + j] = (i1 >= 0 && i1 <= order) ? den[i1] : 0.0;
    }
    int sc = 0, i;
    double prev = (fabs(arr[0]) > TF_EPS) ? arr[0] : 1e-10;
    for (i = 2; i < rows; i++) {
        for (j = 0; j < cols - 1; j++) {
            double det = arr[(i-1)*cols] * arr[(i-2)*cols+j+1] - arr[(i-2)*cols] * arr[(i-1)*cols+j+1];
            arr[i*cols+j] = (fabs(arr[(i-1)*cols]) < TF_EPS) ? arr[(i-2)*cols+j+1] : det / arr[(i-1)*cols];
        }
        double cur = arr[i*cols];
        if (fabs(cur) < TF_EPS) cur = (prev > 0) ? 1e-10 : -1e-10;
        if (prev * cur < 0) sc++;
        if (fabs(cur) > TF_EPS) prev = cur;
    }
    if (arr[0] * arr[cols] < 0) sc++;
    *array_out = arr;
    return sc;
}

int tf_routh_hurwitz_stable(const double *den, int order) {
    if (!den || order < 1) return 1;
    double *arr = NULL;
    int rhp = tf_stability_routh_array(den, order, &arr);
    free(arr);
    return (rhp == 0);
}

int tf_count_rhp_poles(const TransferFunction *tf) {
    if (!tf) return 0;
    int np = tf_num_poles(tf);
    if (np < 1) return 0;
    double *rr = (double*)malloc((size_t)(np + 1) * sizeof(double));
    double *ri = (double*)malloc((size_t)(np + 1) * sizeof(double));
    if (!rr || !ri) { free(rr); free(ri); return 0; }
    int n = poly_find_roots_complex(tf->den, tf->den_order, rr, ri);
    int c = 0, i;
    for (i = 0; i < n; i++) if (rr[i] > TF_EPS) c++;
    free(rr); free(ri);
    return c;
}

double tf_gain_margin(const TransferFunction *tf) {
    if (!tf) return INFINITY;
    double bw = 0.0, bd = INFINITY;
    int i;
    for (i = 0; i < 1000; i++) {
        double w = 0.001 * pow(1e6, (double)i / 999.0);
        double d = fabs(fabs(tf_phase_at(tf, w)) - M_PI);
        if (d < bd) { bd = d; bw = w; }
    }
    if (bw <= 0) return INFINITY;
    double m = tf_magnitude_at(tf, bw);
    return (m < TF_EPS) ? INFINITY : -20.0 * log10(m);
}

double tf_phase_margin(const TransferFunction *tf) {
    if (!tf) return 0.0;
    double bw = 0.0, bd = INFINITY;
    int i;
    for (i = 0; i < 1000; i++) {
        double w = 0.001 * pow(1e6, (double)i / 999.0);
        double d = fabs(tf_magnitude_at(tf, w) - 1.0);
        if (d < bd) { bd = d; bw = w; }
    }
    if (bw <= 0) return 0.0;
    return (M_PI + tf_phase_at(tf, bw)) * 180.0 / M_PI;
}

double tf_final_value_step(const TransferFunction *tf) {
    if (!tf) return 0.0;
    if (tf_stability_routh(tf) == TF_UNSTABLE) return NAN;
    return tf_dc_gain(tf);
}

double tf_initial_value(const TransferFunction *tf) {
    if (!tf) return 0.0;
    return (tf->den_order > tf->num_order) ? 0.0 : tf_hf_gain(tf);
}

int tf_final_value_theorem(const TransferFunction *tf, double *result) {
    if (!tf || !result) return 0;
    if (tf_stability_routh(tf) == TF_UNSTABLE) { *result = NAN; return 0; }
    *result = tf_dc_gain(tf);
    return 1;
}

void tf_eval_at(const TransferFunction *tf, double sr, double si, double *ore, double *oim) {
    *ore = 0.0; *oim = 0.0;
    if (!tf) return;
    double nr, ni, dr, di;
    Polynomial *np = poly_create(tf->num, tf->num_order);
    Polynomial *dp = poly_create(tf->den, tf->den_order);
    poly_eval_complex(np, sr, si, &nr, &ni);
    poly_eval_complex(dp, sr, si, &dr, &di);
    poly_destroy(np); poly_destroy(dp);
    double m2 = dr*dr + di*di;
    if (m2 < 1e-30) { *ore = INFINITY; *oim = INFINITY; return; }
    *ore = (nr*dr + ni*di) / m2;
    *oim = (ni*dr - nr*di) / m2;
}

double tf_magnitude_at(const TransferFunction *tf, double w) {
    double re, im; tf_eval_at(tf, 0.0, w, &re, &im);
    return sqrt(re*re + im*im);
}

double tf_phase_at(const TransferFunction *tf, double w) {
    double re, im; tf_eval_at(tf, 0.0, w, &re, &im);
    return atan2(im, re);
}

int tf_frequency_response(const TransferFunction *tf, double ws, double we, int n, FreqPoint *out) {
    if (!tf || !out || n <= 0 || ws <= 0 || we <= ws) return -1;
    int i;
    for (i = 0; i < n; i++) {
        double w = ws * pow(we / ws, (double)i / (double)(n - 1));
        out[i].frequency = w;
        out[i].magnitude = tf_magnitude_at(tf, w);
        out[i].phase = tf_phase_at(tf, w);
        tf_eval_at(tf, 0.0, w, &out[i].real, &out[i].imag);
    }
    return n;
}

int tf_nyquist_eval(const TransferFunction *tf, double ws, double we, int n, FreqPoint *out) {
    return tf_frequency_response(tf, ws, we, n, out);
}

double tf_bandwidth(const TransferFunction *tf) {
    if (!tf) return -1.0;
    double dc = fabs(tf_dc_gain(tf));
    if (dc < TF_EPS || isinf(dc)) return -1.0;
    double tgt = dc / sqrt(2.0);
    double lo = 1e-3, hi = 1e4;
    if (tf_magnitude_at(tf, hi) > tgt) return hi;
    int i;
    for (i = 0; i < 50; i++) {
        double mid = (lo + hi) / 2.0;
        if (tf_magnitude_at(tf, mid) > tgt) lo = mid;
        else hi = mid;
    }
    return (lo + hi) / 2.0;
}

double tf_resonant_peak(const TransferFunction *tf) {
    if (!tf) return 1.0;
    double dc = fabs(tf_dc_gain(tf));
    if (dc < TF_EPS) return INFINITY;
    double mx = 0.0;
    int i;
    for (i = 0; i < 500; i++) {
        double m = tf_magnitude_at(tf, 0.01 * pow(1e4, (double)i / 499.0));
        if (m > mx) mx = m;
    }
    return mx / dc;
}

double tf_resonant_frequency(const TransferFunction *tf) {
    if (!tf) return -1.0;
    double mx = 0.0, bw = -1.0;
    int i;
    for (i = 0; i < 500; i++) {
        double w = 0.01 * pow(1e4, (double)i / 499.0);
        double m = tf_magnitude_at(tf, w);
        if (m > mx) { mx = m; bw = w; }
    }
    return bw;
}

static void rk4_step(const double *A, const double *B, const double *C,
                      const double *D, int n, double *x, double u,
                      double dt, double *y) {
    double *k = (double*)malloc((size_t)(4 * n) * sizeof(double));
    if (!k) { *y = 0.0; return; }
    double *k1=k, *k2=k+n, *k3=k+2*n, *k4=k+3*n;
    int i, j;
    for (i = 0; i < n; i++) {
        double s = 0.0;
        for (j = 0; j < n; j++) s += A[i*n+j] * x[j];
        k1[i] = s + B[i] * u;
    }
    for (i = 0; i < n; i++) {
        double s = 0.0;
        for (j = 0; j < n; j++) s += A[i*n+j] * (x[j] + 0.5*dt*k1[j]);
        k2[i] = s + B[i] * u;
    }
    for (i = 0; i < n; i++) {
        double s = 0.0;
        for (j = 0; j < n; j++) s += A[i*n+j] * (x[j] + 0.5*dt*k2[j]);
        k3[i] = s + B[i] * u;
    }
    for (i = 0; i < n; i++) {
        double s = 0.0;
        for (j = 0; j < n; j++) s += A[i*n+j] * (x[j] + dt*k3[j]);
        k4[i] = s + B[i] * u;
    }
    for (i = 0; i < n; i++)
        x[i] += dt * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
    double ys = D[0] * u;
    for (j = 0; j < n; j++) ys += C[j] * x[j];
    *y = ys;
    free(k);
}

int tf_step_response(const TransferFunction *tf, double t_end, int n, TimePoint *out) {
    if (!tf || !out || n < 2 || t_end <= 0.0) return -1;
    StateSpace *ss = tf_to_ss_controllable(tf);
    if (!ss) return -1;
    double dt = t_end / (double)(n - 1);
    double *x = (double*)calloc((size_t)ss->n, sizeof(double));
    if (!x) { ss_destroy(ss); return -1; }
    int i;
    for (i = 0; i < n; i++) {
        out[i].t = (double)i * dt;
        rk4_step(ss->A, ss->B, ss->C, ss->D, ss->n, x, 1.0, dt, &out[i].y);
    }
    free(x); ss_destroy(ss);
    return n;
}

int tf_impulse_response(const TransferFunction *tf, double t_end, int n, TimePoint *out) {
    if (!tf || !out || n < 2 || t_end <= 0.0) return -1;
    StateSpace *ss = tf_to_ss_controllable(tf);
    if (!ss) return -1;
    double dt = t_end / (double)(n - 1);
    double *x = (double*)calloc((size_t)ss->n, sizeof(double));
    if (!x) { ss_destroy(ss); return -1; }
    int i;
    for (i = 0; i < ss->n; i++) x[i] = ss->B[i];
    for (i = 0; i < n; i++) {
        out[i].t = (double)i * dt;
        rk4_step(ss->A, ss->B, ss->C, ss->D, ss->n, x, 0.0, dt, &out[i].y);
    }
    free(x); ss_destroy(ss);
    return n;
}

double tf_rise_time(const TransferFunction *tf) {
    if (!tf || tf->den_order < 2) return 0.0;
    if (tf->den_order == 2 && tf->num_order <= 1) {
        int dno = tf->den_order;
        double a2 = tf->den[dno], a1 = tf->den[dno-1], a0 = tf->den[0];
        if (fabs(a2) < TF_EPS) a2 = 1.0;
        double wn = sqrt(fabs(a0/a2)), zeta = a1/(2.0*a2*wn);
        if (zeta > 0.0 && zeta < 1.0 && wn > 0.0)
            return (M_PI - acos(zeta)) / (wn * sqrt(1.0 - zeta*zeta));
    }
    return 0.0;
}

double tf_settling_time(const TransferFunction *tf, double tol) {
    if (!tf || tf->den_order < 2) return 0.0;
    if (tol <= 0.0) tol = 0.02;
    if (tf->den_order == 2 && tf->num_order <= 1) {
        int dno = tf->den_order;
        double a2 = tf->den[dno], a1 = tf->den[dno-1], a0 = tf->den[0];
        if (fabs(a2) < TF_EPS) a2 = 1.0;
        double wn = sqrt(fabs(a0/a2)), zeta = a1/(2.0*a2*wn);
        if (zeta > 0.0 && wn > 0.0)
            return (tol <= 0.02 ? 4.0 : 3.0) / (zeta * wn);
    }
    return 0.0;
}

double tf_overshoot(const TransferFunction *tf) {
    if (!tf || tf->den_order < 2) return 0.0;
    if (tf->den_order == 2 && tf->num_order <= 1) {
        int dno = tf->den_order;
        double a2 = tf->den[dno], a1 = tf->den[dno-1], a0 = tf->den[0];
        if (fabs(a2) < TF_EPS) a2 = 1.0;
        double wn = sqrt(fabs(a0/a2)), zeta = a1/(2.0*a2*wn);
        if (zeta > 0.0 && zeta < 1.0)
            return 100.0 * exp(-M_PI * zeta / sqrt(1.0 - zeta*zeta));
    }
    return 0.0;
}

double tf_peak_time(const TransferFunction *tf) {
    if (!tf || tf->den_order != 2) return 0.0;
    int dno = tf->den_order;
    double a2 = tf->den[dno], a1 = tf->den[dno-1], a0 = tf->den[0];
    if (fabs(a2) < TF_EPS) a2 = 1.0;
    double wn = sqrt(fabs(a0/a2)), zeta = a1/(2.0*a2*wn);
    if (zeta > 0.0 && zeta < 1.0 && wn > 0.0)
        return M_PI / (wn * sqrt(1.0 - zeta*zeta));
    return 0.0;
}

double tf_steady_state_error(const TransferFunction *tf) {
    if (!tf) return 1.0;
    if (tf_system_type(tf) >= TF_TYPE_1) return 0.0;
    double Kp = tf_dc_gain(tf);
    return isinf(Kp) ? 0.0 : 1.0 / (1.0 + Kp);
}

PartialFractionExpansion* tf_partial_fraction(const TransferFunction *tf) {
    if (!tf) return NULL;
    PartialFractionExpansion *pfe = (PartialFractionExpansion*)calloc(1, sizeof(PartialFractionExpansion));
    if (!pfe) return NULL;
    pfe->direct = (tf->num_order >= tf->den_order) ? tf_hf_gain(tf) : 0.0;
    int np = tf_num_poles(tf);
    if (np < 1) { return pfe; }
    double *pr = (double*)malloc((size_t)(np+1)*sizeof(double));
    double *pi = (double*)malloc((size_t)(np+1)*sizeof(double));
    if (!pr || !pi) { free(pr); free(pi); free(pfe); return NULL; }
    int n = poly_find_roots_complex(tf->den, tf->den_order, pr, pi);
    pfe->num_terms = n;
    pfe->terms = (PFE_Term*)calloc((size_t)n, sizeof(PFE_Term));
    if (!pfe->terms) { free(pr); free(pi); free(pfe); return NULL; }
    Polynomial *dp = poly_create(tf->den, tf->den_order);
    Polynomial *dpd = poly_derivative(dp);
    int i;
    for (i = 0; i < n; i++) {
        pfe->terms[i].pole = pr[i];
        pfe->terms[i].multiplicity = 1;
        pfe->terms[i].is_complex = (fabs(pi[i]) > TF_EPS);
        double nr, ni, dr, di;
        Polynomial *npol = poly_create(tf->num, tf->num_order);
        poly_eval_complex(npol, pr[i], pi[i], &nr, &ni);
        poly_eval_complex(dpd, pr[i], pi[i], &dr, &di);
        poly_destroy(npol);
        double m2 = dr*dr + di*di;
        pfe->terms[i].residue = (m2 > TF_EPS) ? (nr*dr + ni*di)/m2 : 0.0;
    }
    poly_destroy(dp); poly_destroy(dpd);
    free(pr); free(pi);
    pfe->delay = tf->delay;
    return pfe;
}

void pfe_destroy(PartialFractionExpansion *pfe) {
    if (pfe) { free(pfe->terms); free(pfe); }
}

void pfe_print(const PartialFractionExpansion *pfe) {
    if (!pfe) { printf("PFE: NULL\n"); return; }
    printf("G(s) = ");
    if (fabs(pfe->direct) > TF_EPS) printf("%.4g + ", pfe->direct);
    int i;
    for (i = 0; i < pfe->num_terms; i++) {
        printf("%.4g/(s - %.4g)", pfe->terms[i].residue, pfe->terms[i].pole);
        if (pfe->terms[i].multiplicity > 1) printf("^%d", pfe->terms[i].multiplicity);
        if (i < pfe->num_terms - 1) printf(" + ");
    }
    if (pfe->delay > TF_EPS) printf(" * exp(-%.4g s)", pfe->delay);
    printf("\n");
}

int tf_dominant_pole(const TransferFunction *tf, double *pole) {
    if (!tf || !pole) return 0;
    int np = tf_num_poles(tf);
    if (np < 1) return 0;
    double *rr = (double*)malloc((size_t)(np+1)*sizeof(double));
    double *ri = (double*)malloc((size_t)(np+1)*sizeof(double));
    if (!rr || !ri) { free(rr); free(ri); return 0; }
    int n = poly_find_roots_complex(tf->den, tf->den_order, rr, ri);
    double best = INFINITY;
    int best_i = -1, i;
    for (i = 0; i < n; i++) {
        if (rr[i] < 0 && fabs(rr[i]) < best) { best = fabs(rr[i]); best_i = i; }
    }
    if (best_i >= 0) { *pole = rr[best_i]; free(rr); free(ri); return 1; }
    free(rr); free(ri);
    return 0;
}

double tf_bode_sensitivity_integral(const TransferFunction *tf) {
    if (!tf) return 0.0;
    int np = tf_num_poles(tf);
    if (np < 1) return 0.0;
    double *rr = (double*)malloc((size_t)(np+1)*sizeof(double));
    double *ri = (double*)malloc((size_t)(np+1)*sizeof(double));
    if (!rr || !ri) { free(rr); free(ri); return 0.0; }
    int n = poly_find_roots_complex(tf->den, tf->den_order, rr, ri);
    double sum = 0.0;
    int i;
    for (i = 0; i < n; i++) if (rr[i] > TF_EPS) sum += rr[i];
    free(rr); free(ri);
    return M_PI * sum;
}

tf_stability_t tf_nyquist_stability(const TransferFunction *L, int *num_enc) {
    if (!L) { if (num_enc) *num_enc = 0; return TF_STABLE; }
    int P = tf_count_rhp_poles(L), N = 0;
    double prev_im = 0.0;
    int i;
    for (i = 0; i < 1000; i++) {
        double w = 0.001 * pow(1e6, (double)i/999.0);
        double re, im; tf_eval_at(L, 0.0, w, &re, &im);
        if (i > 0 && prev_im * im < 0 && re < -1.0) N += (im > prev_im) ? 1 : -1;
        prev_im = im;
    }
    if (num_enc) *num_enc = N;
    return (N + P == 0) ? TF_STABLE : TF_UNSTABLE;
}
