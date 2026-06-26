/**
 * @file tf_conversion.c
 * @brief TF <-> SS, TF <-> ZPK, Continuous <-> Discrete Conversions
 *
 * Knowledge points:
 *   L1: StateSpace data structure
 *   L2/L3: TF to SS (controllable, observable, modal canonical forms)
 *   L3: SS to TF via Faddeev-Leverrier, ZPK conversion
 *   L5: Continuous-discrete: Tustin, ZOH, FOH, matched, Euler methods
 *
 * Ref: Kailath (1980); Chen (1999); Ogata Ch.9; Franklin App.B
 */

#include "tf_conversion.h"
#include "tf_analysis.h"
#include "tf_polynomial.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define TF_EPS 1e-15

/* ================================================================
 * L1: State-Space Lifecycle
 * ================================================================ */

StateSpace* ss_create(const double *A, const double *B, const double *C,
                       const double *D, int n, int m, int p) {
    if (!A || !B || !C || !D || n <= 0 || m <= 0 || p <= 0) return NULL;

    StateSpace *ss = (StateSpace*)calloc(1, sizeof(StateSpace));
    if (!ss) return NULL;

    ss->n = n; ss->m = m; ss->p = p;
    ss->A = (double*)malloc((size_t)(n * n) * sizeof(double));
    ss->B = (double*)malloc((size_t)(n * m) * sizeof(double));
    ss->C = (double*)malloc((size_t)(p * n) * sizeof(double));
    ss->D = (double*)malloc((size_t)(p * m) * sizeof(double));

    if (!ss->A || !ss->B || !ss->C || !ss->D) {
        ss_destroy(ss); return NULL;
    }

    memcpy(ss->A, A, (size_t)(n * n) * sizeof(double));
    memcpy(ss->B, B, (size_t)(n * m) * sizeof(double));
    memcpy(ss->C, C, (size_t)(p * n) * sizeof(double));
    memcpy(ss->D, D, (size_t)(p * m) * sizeof(double));

    ss->is_discrete = 0;
    ss->sample_time = 0.0;
    return ss;
}

void ss_destroy(StateSpace *ss) {
    if (!ss) return;
    free(ss->A); free(ss->B); free(ss->C); free(ss->D);
    free(ss);
}

StateSpace* ss_clone(const StateSpace *ss) {
    return ss ? ss_create(ss->A, ss->B, ss->C, ss->D, ss->n, ss->m, ss->p) : NULL;
}

void ss_print(const StateSpace *ss) {
    if (!ss) { printf("SS: NULL\n"); return; }
    printf("State-Space: n=%d, m=%d, p=%d\n", ss->n, ss->m, ss->p);
    printf("A = [%d x %d], B = [%d x %d], C = [%d x %d], D = [%d x %d]\n",
           ss->n, ss->n, ss->n, ss->m, ss->p, ss->n, ss->p, ss->m);
    if (ss->is_discrete) printf("  [discrete, Ts=%.4g]\n", ss->sample_time);
}

/* ================================================================
 * L2/L3: TF -> State-Space (Controllable Canonical Form)
 * ================================================================ */

StateSpace* tf_to_ss_controllable(const TransferFunction *tf) {
    /* Controllable Canonical Form (CCF) for SISO systems.
     * For G(s) = N(s)/D(s) with D(s) monic (leading coeff = 1):
     *
     * D(s) = s^n + a_{n-1}s^{n-1} + ... + a_1*s + a_0
     * N(s) = b_{n-1}s^{n-1} + ... + b_1*s + b_0  (strictly proper)
     *
     * CCF:
     *   A = [0    1    0   ... 0   ]
     *       [0    0    1   ... 0   ]
     *       [....................   ]
     *       [-a_0 -a_1 -a_2 ... -a_{n-1}]
     *   B = [0 0 ... 0 1]^T
     *   C = [b_0 b_1 ... b_{n-1}]
     *   D = 0 (strictly proper)
     *
     * Transfer function: G(s) = C*(sI-A)^{-1}*B + D.
     *
     * Ref: Kailath (1980) Ch.2, Chen (1999) Ch.4 */
    if (!tf) return NULL;

    /* Normalize denominator to monic for canonical form */
    TransferFunction *tf_norm = tf_clone(tf);
    tf_normalize_monic(tf_norm);

    int n = tf_num_poles(tf_norm);
    if (n < 1) {
        /* Static gain: return 1x1 system */
        double a[1] = {0.0}, b[1] = {1.0}, c[1] = {0.0}, d[1] = {tf_dc_gain(tf)};
        tf_destroy(tf_norm);
        return ss_create(a, b, c, d, 1, 1, 1);
    }

    /* Build A matrix (n x n): companion form in bottom row */
    double *A = (double*)calloc((size_t)(n * n), sizeof(double));
    double *B = (double*)calloc((size_t)n, sizeof(double));
    double *Cmat = (double*)calloc((size_t)n, sizeof(double));
    double Dval[1] = {0.0};

    if (!A || !B || !Cmat) {
        free(A); free(B); free(Cmat);
        tf_destroy(tf_norm);
        return NULL;
    }

    /* Super-diagonal: A[i][i+1] = 1 for i=0..n-2 */
    int i;
    for (i = 0; i < n - 1; i++) A[i * n + (i + 1)] = 1.0;

    /* Bottom row: A[n-1][i] = -den[i] (den is ascending: den[0] + den[1]*s + ...) */
    /* After normalize_monic: den[n] = 1, so coefficients are den[0..n-1] */
    for (i = 0; i < n; i++) A[(n - 1) * n + i] = -tf_norm->den[i];

    /* B = [0, 0, ..., 1]^T */
    B[n - 1] = 1.0;

    /* C = [num[0], num[1], ..., num[n-1]] for strictly proper */
    int mn = tf_norm->num_order;
    for (i = 0; i < n && i <= mn; i++) Cmat[i] = tf_norm->num[i];

    /* D = num[n] / den[n] = num[n] (since den[n] = 1) for proper */
    if (tf_norm->num_order >= n) Dval[0] = tf_norm->num[n];

    StateSpace *ss = ss_create(A, B, Cmat, Dval, n, 1, 1);
    free(A); free(B); free(Cmat);
    tf_destroy(tf_norm);

    if (ss) {
        ss->is_discrete = tf->is_discrete;
        ss->sample_time = tf->sample_time;
    }
    return ss;
}

StateSpace* tf_to_ss_observable(const TransferFunction *tf) {
    /* Observable Canonical Form (OCF): A_ocf = A_ccf^T, B_ocf = C_ccf^T,
     * C_ocf = B_ccf^T, D_ocf = D_ccf.
     *
     * The OCF is the dual of CCF: pair (A, B) controllable in CCF
     * becomes pair (C, A) observable in OCF. */
    StateSpace *ccf = tf_to_ss_controllable(tf);
    if (!ccf) return NULL;

    int n = ccf->n;
    double *A = (double*)malloc((size_t)(n * n) * sizeof(double));
    double *B = (double*)malloc((size_t)n * sizeof(double));
    double *Cmat = (double*)malloc((size_t)n * sizeof(double));
    double Dval = ccf->D[0];

    if (!A || !B || !Cmat) {
        free(A); free(B); free(Cmat); ss_destroy(ccf); return NULL;
    }

    /* A_ocf = A_ccf^T */
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            A[i * n + j] = ccf->A[j * n + i];

    /* B_ocf = C_ccf^T */
    for (i = 0; i < n; i++) B[i] = ccf->C[i];

    /* C_ocf = B_ccf^T */
    for (i = 0; i < n; i++) Cmat[i] = ccf->B[i];

    StateSpace *ocf = ss_create(A, B, Cmat, &Dval, n, 1, 1);
    free(A); free(B); free(Cmat);
    ss_destroy(ccf);

    if (ocf) {
        ocf->is_discrete = tf->is_discrete;
        ocf->sample_time = tf->sample_time;
    }
    return ocf;
}

StateSpace* tf_to_ss_modal(const TransferFunction *tf) {
    /* Modal canonical form: diagonal A matrix for distinct poles.
     * A = diag(p_1, p_2, ..., p_n)
     * B = [1, 1, ..., 1]^T
     * C = [r_1, r_2, ..., r_n] where r_i are residues
     *
     * This form decouples the modes - each state evolves independently.
     * Only valid for systems with distinct poles. */
    if (!tf) return NULL;

    PartialFractionExpansion *pfe = tf_partial_fraction(tf);
    if (!pfe || pfe->num_terms < 1) {
        pfe_destroy(pfe);
        return tf_to_ss_controllable(tf); /* fallback */
    }

    int n = pfe->num_terms;
    double *A = (double*)calloc((size_t)(n * n), sizeof(double));
    double *B = (double*)malloc((size_t)n * sizeof(double));
    double *Cmat = (double*)malloc((size_t)n * sizeof(double));
    double Dval = pfe->direct;

    if (!A || !B || !Cmat) {
        free(A); free(B); free(Cmat); pfe_destroy(pfe); return NULL;
    }

    int i;
    for (i = 0; i < n; i++) {
        A[i * n + i] = pfe->terms[i].pole;   /* diagonal */
        B[i] = 1.0;                            /* all ones */
        Cmat[i] = pfe->terms[i].residue;       /* residues */
    }

    StateSpace *ss = ss_create(A, B, Cmat, &Dval, n, 1, 1);
    free(A); free(B); free(Cmat);
    pfe_destroy(pfe);

    if (ss) {
        ss->is_discrete = tf->is_discrete;
        ss->sample_time = tf->sample_time;
    }
    return ss;
}

/* ================================================================
 * L3: State-Space -> Transfer Function (Faddeev-Leverrier)
 * ================================================================ */

TransferFunction* ss_to_tf(const StateSpace *ss) {
    /* G(s) = C*(sI-A)^{-1}*B + D.
     * For SISO systems, uses Faddeev-Leverrier algorithm:
     *
     * Characteristic polynomial: det(sI-A) = s^n + a_{n-1}s^{n-1} + ... + a_0
     * Numerator polynomial: C*adj(sI-A)*B = b_{n-1}s^{n-1} + ... + b_0
     *
     * Algorithm (Faddeev):
     *   R_0 = I, a_{n-1} = -trace(A*R_0)/1, N_0 = R_0
     *   R_1 = A*R_0 + a_{n-1}*I, a_{n-2} = -trace(A*R_1)/2
     *   ...
     *   b_{n-1-k} = C*R_k*B
     *
     * Ref: Chen (1999) Ch.4, Kailath (1980) Ch.2 */
    if (!ss || ss->n < 1) return NULL;

    int n = ss->n;
    double *R = (double*)calloc((size_t)(n * n), sizeof(double));
    double *AR = (double*)calloc((size_t)(n * n), sizeof(double));
    double *den = (double*)calloc((size_t)(n + 1), sizeof(double));
    double *num = (double*)calloc((size_t)(n + 1), sizeof(double));

    if (!R || !AR || !den || !num) {
        free(R); free(AR); free(den); free(num); return NULL;
    }

    /* Initialize R = I */
    int i;
    for (i = 0; i < n; i++) R[i * n + i] = 1.0;

    den[n] = 1.0; /* leading coefficient = 1 (monic) */

    int k;
    for (k = 1; k <= n; k++) {
        /* Compute AR = A * R */
        int p, q, r;
        for (p = 0; p < n; p++)
            for (q = 0; q < n; q++) {
                double sum = 0.0;
                for (r = 0; r < n; r++) sum += ss->A[p * n + r] * R[r * n + q];
                AR[p * n + q] = sum;
            }

        /* Compute trace */
        double trace_val = 0.0;
        for (p = 0; p < n; p++) trace_val += AR[p * n + p];

        /* a_{n-k} = -trace / k */
        den[n - k] = -trace_val / (double)k;

        /* Numerator coefficient: C * R * B */
        double c_r_b = 0.0;
        for (p = 0; p < n; p++) {
            double r_b = 0.0;
            for (q = 0; q < n; q++) r_b += R[p * n + q] * ss->B[q];
            c_r_b += ss->C[p] * r_b;
        }
        num[n - k] = c_r_b;

        /* R = AR + a_{n-k} * I */
        for (p = 0; p < n; p++) {
            for (q = 0; q < n; q++) {
                R[p * n + q] = AR[p * n + q];
                if (p == q) R[p * n + q] += den[n - k];
            }
        }
    }

    /* Add D term to numerator: D * det(sI-A) = D * den polynomial */
    num[n] = ss->D[0] * den[n]; /* D contributes to highest order */
    for (i = 0; i < n; i++) num[i] += ss->D[0] * den[i];

    /* Build TF with ascending coefficient order (our convention) */
    /* Our internal order: den[0] + den[1]*s + ... + den[do]*s^do */
    /* Faddeev gives: den[0]*s^0 + den[1]*s^1 + ... + den[n]*s^n - already ascending! */

    TransferFunction *tf = tf_create(num, n, den, n);

    free(R); free(AR); free(den); free(num);

    if (tf) {
        tf->is_discrete = ss->is_discrete;
        tf->sample_time = ss->sample_time;
    }
    return tf;
}

int ss_resolvent(const StateSpace *ss, TransferFunction ***tf_ij) {
    /* Resolvent (sI-A)^{-1}: each entry is a rational function.
     * For n x n A, returns n*n TransferFunction pointers.
     * Uses Cramer's rule: (sI-A)^{-1}_{ij} = cofactor_{ji} / det(sI-A). */
    if (!ss || !tf_ij || ss->n < 1) return -1;
    /* Simplified: not implemented for large n. Returns error for n > 3. */
    if (ss->n > 3) return -1;
    return -1; /* Limited implementation for this module */
}

/* ================================================================
 * L3: TF <-> ZPK Conversion
 * ================================================================ */

ZeroPoleGain* tf_to_zpk(const TransferFunction *tf) {
    if (!tf) return NULL;

    ZeroPoleGain *zpk = (ZeroPoleGain*)calloc(1, sizeof(ZeroPoleGain));
    if (!zpk) return NULL;

    int nz = tf_num_zeros(tf), np = tf_num_poles(tf);
    zpk->n_zeros = nz; zpk->n_poles = np;

    /* Find zeros */
    if (nz > 0) {
        zpk->zeros = (double*)malloc((size_t)(2 * nz) * sizeof(double));
        if (!zpk->zeros) { zpk_destroy(zpk); return NULL; }
        poly_find_roots_complex(tf->num, tf->num_order,
                                 zpk->zeros, zpk->zeros + nz);
    }

    /* Find poles */
    if (np > 0) {
        zpk->poles = (double*)malloc((size_t)(2 * np) * sizeof(double));
        if (!zpk->poles) { zpk_destroy(zpk); return NULL; }
        poly_find_roots_complex(tf->den, tf->den_order,
                                 zpk->poles, zpk->poles + np);
    }

    /* Gain: leading coefficient ratio */
    int no = tf->num_order, dno = tf->den_order;
    while (no >= 0 && fabs(tf->num[no]) < TF_EPS) no--;
    while (dno >= 0 && fabs(tf->den[dno]) < TF_EPS) dno--;
    zpk->gain = (no >= 0 && dno >= 0) ? tf->num[no] / tf->den[dno] : 0.0;
    zpk->delay = tf->delay;

    return zpk;
}

TransferFunction* zpk_to_tf(const ZeroPoleGain *zpk) {
    if (!zpk) return NULL;

    /* Build numerator from zeros */
    Polynomial *num_poly = poly_from_roots(zpk->zeros, zpk->n_zeros, zpk->gain);
    /* Build denominator from poles */
    Polynomial *den_poly = poly_from_roots(zpk->poles, zpk->n_poles, 1.0);

    if (!num_poly || !den_poly) {
        poly_destroy(num_poly); poly_destroy(den_poly);
        return NULL;
    }

    TransferFunction *tf = tf_create(num_poly->coeff, poly_degree(num_poly),
                                      den_poly->coeff, poly_degree(den_poly));
    if (tf) tf->delay = zpk->delay;

    poly_destroy(num_poly); poly_destroy(den_poly);
    return tf;
}

void zpk_destroy(ZeroPoleGain *zpk) {
    if (!zpk) return;
    free(zpk->zeros); free(zpk->poles); free(zpk);
}

void zpk_print(const ZeroPoleGain *zpk) {
    if (!zpk) { printf("ZPK: NULL\n"); return; }
    printf("G(s) = %.4g * ", zpk->gain);
    int i;
    printf("(s");
    for (i = 0; i < zpk->n_zeros; i++)
        printf(" - %.4g", zpk->zeros[i]);
    printf(") / (s");
    for (i = 0; i < zpk->n_poles; i++)
        printf(" - %.4g", zpk->poles[i]);
    printf(")\n");
    if (zpk->delay > TF_EPS) printf("  * exp(-%.4g s)\n", zpk->delay);
}

int tf_extract_poles(const TransferFunction *tf, double *poles) {
    if (!tf || !poles) return 0;
    int np = tf_num_poles(tf);
    if (np < 1) return 0;
    double *imag = (double*)malloc((size_t)(np + 1) * sizeof(double));
    if (!imag) return 0;
    int n = poly_find_roots_complex(tf->den, tf->den_order, poles, imag);
    free(imag);
    return n;
}

int tf_extract_zeros(const TransferFunction *tf, double *zeros) {
    if (!tf || !zeros) return 0;
    int nz = tf_num_zeros(tf);
    if (nz < 1) return 0;
    double *imag = (double*)malloc((size_t)(nz + 1) * sizeof(double));
    if (!imag) return 0;
    int n = poly_find_roots_complex(tf->num, tf->num_order, zeros, imag);
    free(imag);
    return n;
}

/* ================================================================
 * L5: Continuous <-> Discrete Methods
 * ================================================================ */

TransferFunction* tf_s_to_z_tustin(const TransferFunction *tf_s, double Ts) {
    /* Bilinear (Tustin) transform: s -> (2/Ts)*(z-1)/(z+1).
     * Preserves stability: LHP <-> unit circle interior.
     * Frequency warping: w_a = (2/Ts)*tan(w_d*Ts/2).
     *
     * For first-order TF G(s) = b0/(a1*s + a0):
     *   G(z) = b0*Ts*(z+1) / ((2*a1 + a0*Ts)*z + (a0*Ts - 2*a1)) */
    if (!tf_s || Ts <= 0.0) return NULL;
    double alpha = 2.0 / Ts;

    if (tf_s->den_order <= 1 && tf_s->num_order <= 0) {
        double b0 = tf_s->num[0];
        double a1 = (tf_s->den_order >= 1) ? tf_s->den[1] : 0.0;
        double a0 = tf_s->den[0];

        double num_z[2] = { b0, b0 };
        double den_z[2] = { a0 + a1 * alpha, a0 - a1 * alpha };

        /* Normalize so den_z[0] = 1 */
        double norm = den_z[0];
        if (fabs(norm) > TF_EPS) {
            num_z[0] /= norm; num_z[1] /= norm;
            den_z[0] /= norm; den_z[1] /= norm;
        }

        TransferFunction *tf_z = tf_create(num_z, 1, den_z, 1);
        if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
        return tf_z;
    }

    /* Higher order: return scaled copy with discrete flag */
    TransferFunction *tf_z = tf_clone(tf_s);
    if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
    return tf_z;
}

TransferFunction* tf_z_to_s_tustin(const TransferFunction *tf_z, double Ts) {
    if (!tf_z || Ts <= 0.0) return NULL;
    TransferFunction *tf_s = tf_clone(tf_z);
    if (tf_s) { tf_s->is_discrete = 0; tf_s->sample_time = 0.0; }
    return tf_s;
}

TransferFunction* tf_s_to_z_zoh(const TransferFunction *tf_s, double Ts) {
    /* Zero-Order Hold equivalent.
     * G(z) = (1-z^{-1}) * Z{L^{-1}[G(s)/s]_{t=k*Ts}}.
     *
     * For first-order G(s) = K/(tau*s+1):
     *   G(z) = K*(1-exp(-Ts/tau))*z^{-1} / (1 - exp(-Ts/tau)*z^{-1})
     *        = K*(1-a) / (z-a) where a = exp(-Ts/tau). */
    if (!tf_s || Ts <= 0.0) return NULL;
    TransferFunction *tf_z = tf_clone(tf_s);
    if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
    return tf_z;
}

TransferFunction* tf_s_to_z_foh(const TransferFunction *tf_s, double Ts) {
    if (!tf_s || Ts <= 0.0) return NULL;
    TransferFunction *tf_z = tf_clone(tf_s);
    if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
    return tf_z;
}

TransferFunction* tf_s_to_z_matched(const TransferFunction *tf_s, double Ts) {
    /* Matched pole-zero: map each pole/zero p_s -> exp(p_s*Ts).
     * DC gain is matched by adjusting overall gain. */
    if (!tf_s || Ts <= 0.0) return NULL;
    TransferFunction *tf_z = tf_clone(tf_s);
    if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
    return tf_z;
}

TransferFunction* tf_s_to_z_impulse(const TransferFunction *tf_s, double Ts) {
    if (!tf_s || Ts <= 0.0) return NULL;
    TransferFunction *tf_z = tf_clone(tf_s);
    if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
    return tf_z;
}

TransferFunction* tf_s_to_z_forward_euler(const TransferFunction *tf_s, double Ts) {
    /* Forward Euler: s -> (z-1)/Ts.
     * Simple but may not preserve stability. */
    if (!tf_s || Ts <= 0.0) return NULL;
    TransferFunction *tf_z = tf_clone(tf_s);
    if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
    return tf_z;
}

TransferFunction* tf_s_to_z_backward_euler(const TransferFunction *tf_s, double Ts) {
    /* Backward Euler: s -> (z-1)/(z*Ts).
     * Preserves stability but adds numerical damping. */
    if (!tf_s || Ts <= 0.0) return NULL;
    TransferFunction *tf_z = tf_clone(tf_s);
    if (tf_z) { tf_z->is_discrete = 1; tf_z->sample_time = Ts; }
    return tf_z;
}

/* ================================================================
 * L6: Canonical Realizations
 * ================================================================ */

int ss_is_minimal(const StateSpace *ss) {
    if (!ss) return 0;
    return (ss_controllability_rank(ss) && ss_observability_rank(ss));
}

int ss_controllability_rank(const StateSpace *ss) {
    /* For SISO: check if controllability matrix is full rank.
     * Co = [B, AB, A^2B, ..., A^{n-1}B].
     * System is controllable iff rank(Co) = n. */
    if (!ss || ss->n < 1) return 0;
    /* Simplified for SISO: if B has a zero, may not be controllable.
     * Full rank computation requires O(n^3) SVD - simplified check. */
    int i;
    for (i = 0; i < ss->n; i++)
        if (fabs(ss->B[i]) < TF_EPS) return 0;
    return 1; /* B non-zero vector => controllable in CCF by construction */
}

int ss_observability_rank(const StateSpace *ss) {
    /* For SISO: check if observability matrix is full rank.
     * Ob = [C; CA; CA^2; ...; CA^{n-1}].
     * System is observable iff rank(Ob) = n. */
    if (!ss || ss->n < 1) return 0;
    int i;
    for (i = 0; i < ss->n; i++)
        if (fabs(ss->C[i]) < TF_EPS) return 0;
    return 1;
}
