/**
 * @file electrical_state_space.c
 * @brief State-space models for electrical circuits
 *
 * Knowledge points — each function implements one independent concept:
 *   L1: State-space initialization and data access patterns
 *   L2: Eigenvalue computation via QR algorithm, Hurwitz stability
 *   L3: Controllability and observability Gramians and rank
 *   L4: Matrix exponential via Pade with scaling-and-squaring
 *   L5: RK4 simulation, RLC-to-state-space, canonical forms, ZOH discretization
 *   L6: Modal parameter extraction, dominant time constant
 *
 * Reference:
 *   MIT 6.241 Dynamic Systems and Control
 *   Higham (2009) The Scaling and Squaring Method for Matrix Exponential
 *   Van Loan (1978) Computing Integrals Involving the Matrix Exponential
 */

#include "electrical_state_space.h"
#include <math.h>
#include <complex.h>
#include <string.h>
#include <float.h>

/* ===== State-Space Initialization (L1) ===== */

void ss_init(StateSpace *ss, size_t n, size_t m, size_t p)
{
    if (!ss) return;
    memset(ss, 0, sizeof(StateSpace));
    ss->n = (n <= SS_MAX_ORDER) ? n : SS_MAX_ORDER;
    ss->m = (m <= SS_MAX_INPUTS) ? m : SS_MAX_INPUTS;
    ss->p = (p <= SS_MAX_OUTPUTS) ? p : SS_MAX_OUTPUTS;
}

/* ===== Hessenberg Reduction via Householder Reflections (L2) ===== */
/* Preprocesses matrix A for QR iteration. Complexity: O(n^3) */

static void hessenberg_reduce(double *H, size_t n)
{
    for (size_t k = 0; k + 1 < n; k++) {
        double sigma = 0.0;
        for (size_t i = k + 1; i < n; i++)
            sigma += H[i*n + k] * H[i*n + k];
        if (sigma < 1e-30) continue;
        double x1 = H[(k+1)*n + k];
        double alpha = (x1 > 0) ? -sqrt(sigma) : sqrt(sigma);
        double v1 = x1 - alpha;
        double beta = 1.0 / (sigma - x1 * alpha);
        /* Left transformation: H = (I - beta*v*v^T) * H */
        for (size_t j = k; j < n; j++) {
            double dot = v1 * H[(k+1)*n + j];
            for (size_t i = k + 2; i < n; i++)
                dot += H[i*n + k] * H[i*n + j];
            dot *= beta;
            H[(k+1)*n + j] -= dot * v1;
            for (size_t i = k + 2; i < n; i++)
                H[i*n + j] -= dot * H[i*n + k];
        }
        /* Right transformation: H = H * (I - beta*v*v^T) */
        for (size_t i = 0; i < n; i++) {
            double dot = H[i*n + (k+1)] * v1;
            for (size_t j = k + 2; j < n; j++)
                dot += H[i*n + j] * H[j*n + k];
            dot *= beta;
            H[i*n + (k+1)] -= dot * v1;
            for (size_t j = k + 2; j < n; j++)
                H[i*n + j] -= dot * H[j*n + k];
        }
    }
}

/* ===== QR Double-Step with Wilkinson Shift (L2) ===== */
/* Implicit double-shift bulge-chase on Hessenberg form */

static void qr_double_step(double *H, size_t n)
{
    /* Wilkinson shift from trailing 2x2 submatrix */
    double tr = H[(n-2)*n+(n-2)] + H[(n-1)*n+(n-1)];
    double det = H[(n-2)*n+(n-2)]*H[(n-1)*n+(n-1)]
               - H[(n-2)*n+(n-1)]*H[(n-1)*n+(n-2)];
    /* Shift = eigenvalue of trailing 2x2 closest to a_{nn}.
     * The implicit QR double-shift uses (A - mu1*I)(A - mu2*I)
     * which is real even when mu1,mu2 are complex conjugates.
     * The real matrix M = A^2 - tr*A + det*I encodes both shifts.
     * disc = tr^2 - 4*det is the discriminant for eigenvalue calculation. */
    (void)(tr*tr - 4.0*det); /* Suppress unused warning — used conceptually */

    /* Implicit double-shift: create bulge at top of Hessenberg
     * using M = (H - mu1*I)(H - mu2*I) = H^2 - tr*H + det*I
     * The first column of M is [x; y; z; 0; ...] */
    double x = H[0]*H[0] - tr*H[0] + det;
    double y = H[1*n+0]*(H[0] + H[1*n+1] - tr);
    double z = H[1*n+0]*H[2*n+1];

    for (size_t k = 0; k + 2 < n; k++) {
        /* Compute 3x3 Householder reflector for column [x; y; z] */
        double r = sqrt(x*x + y*y + z*z);
        if (r < 1e-15) { x = H[(k+1)*n+k]; y = H[(k+2)*n+k]; z = 0; continue; }
        double s = (x < 0) ? -r : r;
        double vx = x + s, vy = y, vz = z;
        double b = 1.0 / (r*(r + fabs(x)));

        /* Apply from left to columns k..n-1 */
        for (size_t j = k; j < n; j++) {
            double dot = vx*H[k*n+j] + vy*H[(k+1)*n+j] + vz*H[(k+2)*n+j];
            dot *= b;
            H[k*n+j]     -= dot * vx;
            H[(k+1)*n+j] -= dot * vy;
            H[(k+2)*n+j] -= dot * vz;
        }
        /* Apply from right to rows 0..min(k+3,n-1) */
        size_t rmax = (k + 3 < n) ? k + 3 : n;
        for (size_t i = 0; i < rmax; i++) {
            double dot = H[i*n+k]*vx + H[i*n+(k+1)]*vy + H[i*n+(k+2)]*vz;
            dot *= b;
            H[i*n+k]     -= dot * vx;
            H[i*n+(k+1)] -= dot * vy;
            H[i*n+(k+2)] -= dot * vz;
        }
        /* Next bulge location */
        if (k + 3 < n) {
            x = H[(k+2)*n+(k+1)];
            y = H[(k+3)*n+(k+1)];
            z = 0;
        }
    }
}

/* ===== Eigenvalue Extraction (L2) ===== */

int ss_eigenvalues(const StateSpace *ss, double complex *eigenvalues)
{
    if (!ss || !eigenvalues || ss->n == 0) return 0;
    size_t n = ss->n;
    double *H = (double*)malloc(n * n * sizeof(double));
    if (!H) return 0;
    for (size_t i = 0; i < n*n; i++) H[i] = ss->A[i];

    hessenberg_reduce(H, n);

    for (int iter = 0; iter < 200; iter++) {
        int done = 1;
        for (size_t i = 0; i + 1 < n; i++) {
            if (fabs(H[(i+1)*n + i]) > 1e-12 *
                (fabs(H[i*n+i]) + fabs(H[(i+1)*n + i + 1]))) {
                done = 0;
                break;
            }
        }
        if (done) break;
        qr_double_step(H, n);
    }

    /* Extract eigenvalues from diagonal and 2x2 blocks */
    size_t cnt = 0;
    size_t i = 0;
    while (i < n && cnt < SS_MAX_ORDER) {
        if (i + 1 < n && fabs(H[(i+1)*n + i]) > 1e-10) {
            /* 2x2 block -> complex conjugate pair */
            double a = H[i*n+i], b = H[i*n+i+1];
            double c = H[(i+1)*n+i], d = H[(i+1)*n+i+1];
            double tr2 = a + d;
            double det2 = a*d - b*c;
            double disc2 = tr2*tr2 - 4.0*det2;
            if (disc2 < 0) {
                double re = 0.5 * tr2;
                double im = 0.5 * sqrt(-disc2);
                eigenvalues[cnt++] = re + I*im;
                eigenvalues[cnt++] = re - I*im;
            } else {
                double sd = sqrt(disc2);
                eigenvalues[cnt++] = 0.5*(tr2 + sd) + I*0.0;
                eigenvalues[cnt++] = 0.5*(tr2 - sd) + I*0.0;
            }
            i += 2;
        } else {
            eigenvalues[cnt++] = H[i*n+i] + I*0.0;
            i++;
        }
    }
    free(H);
    return (int)cnt;
}

/* ===== Hurwitz Stability Check (L2) ===== */

int ss_is_stable(const StateSpace *ss)
{
    /* All eigenvalues must have Re(lambda) < 0 for asymptotic stability */
    if (!ss || ss->n == 0) return 0;
    double complex evals[SS_MAX_ORDER];
    int nev = ss_eigenvalues(ss, evals);
    for (int i = 0; i < nev; i++)
        if (creal(evals[i]) >= -1e-12) return 0;
    return 1;
}

/* ===== Controllability Gramian (L3) ===== */
/* W_c = sum_{k=0}^{n-1} A^k B B^T (A^T)^k */

int ss_controllability_gramian(const StateSpace *ss, double *Gramian)
{
    if (!ss || !Gramian || ss->n == 0) return 0;
    size_t n = ss->n, m = ss->m;
    for (size_t i = 0; i < n*n; i++) Gramian[i] = 0.0;

    for (size_t p = 0; p < m && p < SS_MAX_INPUTS; p++) {
        double x[SS_MAX_ORDER] = {0};
        for (size_t i = 0; i < n; i++) x[i] = ss->B[i*n + p];

        for (size_t k = 0; k < n; k++) {
            /* Accumulate x*x^T into Gramian */
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                    Gramian[i*n + j] += x[i] * x[j];

            /* x = A * x */
            double xn[SS_MAX_ORDER] = {0};
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                    xn[i] += SS_A(ss, i, j) * x[j];
            for (size_t i = 0; i < n; i++) x[i] = xn[i];
        }
    }
    return (int)n;
}

/* ===== Controllability Rank via QR with Column Pivoting (L3) ===== */

int ss_controllability_rank(const StateSpace *ss)
{
    if (!ss || ss->n == 0) return 0;
    size_t n = ss->n, m = ss->m;
    size_t ncols = n * m;
    if (ncols > SS_MAX_ORDER * SS_MAX_INPUTS) ncols = SS_MAX_ORDER * SS_MAX_INPUTS;

    double *C = (double*)calloc(n * ncols, sizeof(double));
    if (!C) return 0;

    /* Build controllability matrix: C = [B, AB, A^2B, ..., A^{n-1}B] */
    for (size_t p = 0; p < m && p < SS_MAX_INPUTS; p++) {
        double x[SS_MAX_ORDER] = {0};
        for (size_t i = 0; i < n; i++) x[i] = ss->B[i*n + p];
        for (size_t k = 0; k < n; k++) {
            for (size_t i = 0; i < n; i++)
                C[i*ncols + (k*m + p)] = x[i];
            double xn[SS_MAX_ORDER] = {0};
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                    xn[i] += SS_A(ss, i, j) * x[j];
            for (size_t i = 0; i < n; i++) x[i] = xn[i];
        }
    }

    /* Rank via modified Gram-Schmidt with column pivoting */
    double *col_norms = (double*)malloc(ncols * sizeof(double));
    int *pivoted = (int*)calloc(ncols, sizeof(int));
    if (!col_norms || !pivoted) {
        free(C); free(col_norms); free(pivoted);
        return 0;
    }
    for (size_t j = 0; j < ncols; j++) {
        col_norms[j] = 0.0;
        for (size_t i = 0; i < n; i++)
            col_norms[j] += C[i*ncols + j] * C[i*ncols + j];
        col_norms[j] = sqrt(col_norms[j]);
    }

    int rank = 0;
    for (size_t step = 0; step < n && step < ncols; step++) {
        size_t max_j = 0;
        double max_norm = -1.0;
        for (size_t j = 0; j < ncols; j++) {
            if (!pivoted[j] && col_norms[j] > max_norm) {
                max_norm = col_norms[j];
                max_j = j;
            }
        }
        if (max_norm < 1e-12) break;
        pivoted[max_j] = 1;
        rank++;

        /* Orthogonalize remaining columns against chosen column */
        for (size_t j = 0; j < ncols; j++) {
            if (pivoted[j]) continue;
            double dot = 0.0, qn2 = 0.0;
            for (size_t i = 0; i < n; i++) {
                dot += C[i*ncols + max_j] * C[i*ncols + j];
                qn2 += C[i*ncols + max_j] * C[i*ncols + max_j];
            }
            double proj = dot / qn2;
            for (size_t i = 0; i < n; i++)
                C[i*ncols + j] -= proj * C[i*ncols + max_j];
            col_norms[j] = 0.0;
            for (size_t i = 0; i < n; i++)
                col_norms[j] += C[i*ncols + j] * C[i*ncols + j];
            col_norms[j] = sqrt(col_norms[j]);
        }
    }
    free(pivoted);
    free(col_norms);
    free(C);
    return rank;
}

/* ===== Observability Gramian (L3) ===== */
/* W_o = sum_{k=0}^{n-1} (A^T)^k C^T C A^k */

int ss_observability_gramian(const StateSpace *ss, double *Gramian)
{
    if (!ss || !Gramian || ss->n == 0) return 0;
    size_t n = ss->n;
    for (size_t i = 0; i < n*n; i++) Gramian[i] = 0.0;

    for (size_t q = 0; q < ss->p && q < SS_MAX_OUTPUTS; q++) {
        double y[SS_MAX_ORDER] = {0};
        for (size_t j = 0; j < n; j++) y[j] = ss->C[q*n + j];

        for (size_t k = 0; k < n; k++) {
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                    Gramian[i*n + j] += y[i] * y[j];

            /* y = A^T * y */
            double yn[SS_MAX_ORDER] = {0};
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                    yn[i] += SS_A(ss, j, i) * y[j];
            for (size_t i = 0; i < n; i++) y[i] = yn[i];
        }
    }
    return (int)n;
}

/* ===== Observability Rank (L3) ===== */

int ss_observability_rank(const StateSpace *ss)
{
    if (!ss || ss->n == 0) return 0;
    size_t n = ss->n, p = ss->p;
    size_t nrows = n * p;
    if (nrows > SS_MAX_ORDER * SS_MAX_OUTPUTS) nrows = SS_MAX_ORDER * SS_MAX_OUTPUTS;

    double *O = (double*)calloc(nrows * n, sizeof(double));
    if (!O) return 0;

    /* Build observability matrix: O = [C; CA; CA^2; ...; CA^{n-1}] */
    for (size_t q = 0; q < p && q < SS_MAX_OUTPUTS; q++)
        for (size_t j = 0; j < n; j++)
            O[(q*n) + j] = ss->C[q*n + j];

    double *cur = (double*)malloc(p * n * sizeof(double));
    if (!cur) { free(O); return 0; }
    memcpy(cur, ss->C, p * n * sizeof(double));

    for (size_t k = 1; k < n; k++) {
        double *nxt = (double*)calloc(p * n, sizeof(double));
        if (!nxt) { free(cur); free(O); return 0; }
        for (size_t q = 0; q < p; q++)
            for (size_t j = 0; j < n; j++)
                for (size_t i = 0; i < n; i++)
                    nxt[q*n + j] += cur[q*n + i] * SS_A(ss, i, j);
        for (size_t q = 0; q < p && q < SS_MAX_OUTPUTS; q++)
            for (size_t j = 0; j < n; j++)
                O[(k*p + q)*n + j] = nxt[q*n + j];
        memcpy(cur, nxt, p * n * sizeof(double));
        free(nxt);
    }
    free(cur);

    /* Rank via column-pivoted QR on transposed matrix (columns = states) */
    double *cn2 = (double*)malloc(n * sizeof(double));
    int *pv2 = (int*)calloc(n, sizeof(int));
    if (!cn2 || !pv2) { free(cn2); free(pv2); free(O); return 0; }
    for (size_t j = 0; j < n; j++) {
        cn2[j] = 0.0;
        for (size_t i = 0; i < nrows; i++)
            cn2[j] += O[i*n + j] * O[i*n + j];
        cn2[j] = sqrt(cn2[j]);
    }

    int rank = 0;
    for (size_t step = 0; step < n; step++) {
        size_t mj = 0;
        double mn = -1.0;
        for (size_t j = 0; j < n; j++)
            if (!pv2[j] && cn2[j] > mn) { mn = cn2[j]; mj = j; }
        if (mn < 1e-12) break;
        pv2[mj] = 1;
        rank++;

        for (size_t j = 0; j < n; j++) {
            if (pv2[j]) continue;
            double dot = 0.0, qn2 = 0.0;
            for (size_t i = 0; i < nrows; i++) {
                dot += O[i*n + mj] * O[i*n + j];
                qn2 += O[i*n + mj] * O[i*n + mj];
            }
            double proj = dot / qn2;
            for (size_t i = 0; i < nrows; i++)
                O[i*n + j] -= proj * O[i*n + mj];
            cn2[j] = 0.0;
            for (size_t i = 0; i < nrows; i++)
                cn2[j] += O[i*n + j] * O[i*n + j];
            cn2[j] = sqrt(cn2[j]);
        }
    }
    free(pv2); free(cn2); free(O);
    return rank;
}

/* ===== Matrix Exponential via Pade [8/8] with Scaling-and-Squaring (L4) ===== */
/* Reference: Higham (2009) SIAM Review 51(4), Algorithm 10.20 */

int ss_state_transition(const StateSpace *ss, double t, double *Phi)
{
    if (!ss || !Phi || ss->n == 0) return -1;
    size_t n = ss->n;

    /* Determine scaling factor: ||A*t||_inf <= 0.5 after scaling */
    double normA = 0.0;
    for (size_t i = 0; i < n; i++) {
        double rs = 0.0;
        for (size_t j = 0; j < n; j++) rs += fabs(SS_A(ss, i, j) * t);
        if (rs > normA) normA = rs;
    }
    int k = 0;
    double sc = 1.0;
    while (normA / sc > 0.5) { k++; sc *= 2.0; }

    double *As = (double*)malloc(n * n * sizeof(double));
    double *F  = (double*)calloc(n * n, sizeof(double));
    double *Fn = (double*)calloc(n * n, sizeof(double));
    double *Fd = (double*)calloc(n * n, sizeof(double));
    double *T  = (double*)calloc(n * n, sizeof(double));
    if (!As || !F || !Fn || !Fd || !T) {
        free(As); free(F); free(Fn); free(Fd); free(T);
        return -1;
    }

    for (size_t i = 0; i < n*n; i++) As[i] = ss->A[i] * t / sc;
    for (size_t i = 0; i < n; i++) F[i*n+i] = Fn[i*n+i] = Fd[i*n+i] = 1.0;

    /* Pade [8/8] coefficients: e^x = P_8(x)/Q_8(x) */
    const double pc[9] = {
        1.0, 0.5, 3.0/28.0, 1.0/84.0, 5.0/3024.0,
        1.0/15120.0, 1.0/332640.0, 1.0/8648640.0, 1.0/518918400.0
    };
    const double qc[9] = {
        1.0, -0.5, 3.0/28.0, -1.0/84.0, 5.0/3024.0,
        -1.0/15120.0, 1.0/332640.0, -1.0/8648640.0, 1.0/518918400.0
    };

    /* Initialize P_8 and Q_8 with highest-degree coefficient * I */
    for (size_t ii = 0; ii < n; ii++) {
        for (size_t jj = 0; jj < n; jj++) {
            Fn[ii*n+jj] = (ii == jj) ? pc[8] : 0.0;
            Fd[ii*n+jj] = (ii == jj) ? qc[8] : 0.0;
        }
    }

    /* Horner evaluation: P_k = P_{k+1}*A_s + p_k*I */
    for (int d = 7; d >= 0; d--) {
        /* Fn = Fn * As + pc[d] * I */
        for (size_t i = 0; i < n*n; i++) T[i] = 0.0;
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                for (size_t r = 0; r < n; r++)
                    T[i*n+j] += Fn[i*n+r] * As[r*n+j];
        for (size_t i = 0; i < n*n; i++) Fn[i] = T[i];
        for (size_t i = 0; i < n; i++) Fn[i*n+i] += pc[d];

        /* Fd = Fd * As + qc[d] * I */
        for (size_t i = 0; i < n*n; i++) T[i] = 0.0;
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                for (size_t r = 0; r < n; r++)
                    T[i*n+j] += Fd[i*n+r] * As[r*n+j];
        for (size_t i = 0; i < n*n; i++) Fd[i] = T[i];
        for (size_t i = 0; i < n; i++) Fd[i*n+i] += qc[d];
    }

    /* Solve Fd * F = Fn column by column (Gaussian elimination) */
    for (size_t col = 0; col < n; col++) {
        double rhs2[SS_MAX_ORDER], lhs2[SS_MAX_ORDER * SS_MAX_ORDER];
        for (size_t i = 0; i < n; i++) {
            rhs2[i] = Fn[i*n + col];
            for (size_t j = 0; j < n; j++)
                lhs2[i*n + j] = Fd[i*n + j];
        }
        /* Forward elimination with partial pivoting */
        for (size_t kk = 0; kk < n; kk++) {
            size_t mr = kk;
            double mv = fabs(lhs2[kk*n + kk]);
            for (size_t i2 = kk + 1; i2 < n; i2++)
                if (fabs(lhs2[i2*n + kk]) > mv)
                    { mv = fabs(lhs2[i2*n + kk]); mr = i2; }
            if (mv < 1e-15) continue;
            if (mr != kk) {
                for (size_t j2 = kk; j2 < n; j2++) {
                    double tp = lhs2[kk*n + j2];
                    lhs2[kk*n + j2] = lhs2[mr*n + j2];
                    lhs2[mr*n + j2] = tp;
                }
                double tp = rhs2[kk]; rhs2[kk] = rhs2[mr]; rhs2[mr] = tp;
            }
            double pv = lhs2[kk*n + kk];
            for (size_t i2 = kk + 1; i2 < n; i2++) {
                double fc = lhs2[i2*n + kk] / pv;
                for (size_t j2 = kk; j2 < n; j2++)
                    lhs2[i2*n + j2] -= fc * lhs2[kk*n + j2];
                rhs2[i2] -= fc * rhs2[kk];
            }
        }
        /* Back substitution */
        double sol2[SS_MAX_ORDER] = {0};
        for (size_t ii2 = n; ii2 > 0; ii2--) {
            size_t i2 = ii2 - 1;
            double sm = rhs2[i2];
            for (size_t j2 = i2 + 1; j2 < n; j2++)
                sm -= lhs2[i2*n + j2] * sol2[j2];
            sol2[i2] = sm / lhs2[i2*n + i2];
        }
        for (size_t i2 = 0; i2 < n; i2++) F[i2*n + col] = sol2[i2];
    }

    /* Repeated squaring: F = F^(2^k) to undo scaling */
    for (int s = 0; s < k; s++) {
        for (size_t i = 0; i < n*n; i++) T[i] = 0.0;
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                for (size_t r = 0; r < n; r++)
                    T[i*n+j] += F[i*n+r] * F[r*n+j];
        for (size_t i = 0; i < n*n; i++) F[i] = T[i];
    }

    for (size_t i = 0; i < n*n; i++) Phi[i] = F[i];
    free(As); free(F); free(Fn); free(Fd); free(T);
    return 0;
}

/* ===== RK4 State Simulation (L5) ===== */

int ss_simulate_rk4(StateSpace *ss, double *x, const double *u,
                     double t0, double tf, size_t steps)
{
    if (!ss || !x || !u || steps == 0) return -1;
    double dt = (tf - t0) / (double)steps;
    size_t n = ss->n, m = ss->m;

    double *k1 = (double*)malloc(n * sizeof(double));
    double *k2 = (double*)malloc(n * sizeof(double));
    double *k3 = (double*)malloc(n * sizeof(double));
    double *k4 = (double*)malloc(n * sizeof(double));
    double *xt = (double*)malloc(n * sizeof(double));
    if (!k1 || !k2 || !k3 || !k4 || !xt) {
        free(k1); free(k2); free(k3); free(k4); free(xt);
        return -1;
    }

    for (size_t step = 0; step < steps; step++) {
        /* k1 = A*x + B*u */
        for (size_t i = 0; i < n; i++) {
            k1[i] = 0.0;
            for (size_t j = 0; j < n; j++) k1[i] += SS_A(ss, i, j) * x[j];
            for (size_t j = 0; j < m; j++) k1[i] += SS_B(ss, i, j) * u[j];
        }
        /* k2: evaluate at midpoint using k1 */
        for (size_t i = 0; i < n; i++) xt[i] = x[i] + 0.5 * dt * k1[i];
        for (size_t i = 0; i < n; i++) {
            k2[i] = 0.0;
            for (size_t j = 0; j < n; j++) k2[i] += SS_A(ss, i, j) * xt[j];
            for (size_t j = 0; j < m; j++) k2[i] += SS_B(ss, i, j) * u[j];
        }
        /* k3: evaluate at midpoint using k2 */
        for (size_t i = 0; i < n; i++) xt[i] = x[i] + 0.5 * dt * k2[i];
        for (size_t i = 0; i < n; i++) {
            k3[i] = 0.0;
            for (size_t j = 0; j < n; j++) k3[i] += SS_A(ss, i, j) * xt[j];
            for (size_t j = 0; j < m; j++) k3[i] += SS_B(ss, i, j) * u[j];
        }
        /* k4: evaluate at endpoint using k3 */
        for (size_t i = 0; i < n; i++) xt[i] = x[i] + dt * k3[i];
        for (size_t i = 0; i < n; i++) {
            k4[i] = 0.0;
            for (size_t j = 0; j < n; j++) k4[i] += SS_A(ss, i, j) * xt[j];
            for (size_t j = 0; j < m; j++) k4[i] += SS_B(ss, i, j) * u[j];
        }
        /* Update: x += (dt/6)*(k1 + 2k2 + 2k3 + k4) */
        for (size_t i = 0; i < n; i++)
            x[i] += (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    free(k1); free(k2); free(k3); free(k4); free(xt);
    return 0;
}

/* ===== RLC Circuit to State-Space (L5) ===== */

StateSpace ss_from_series_rlc(double R, double L, double C)
{
    /* State: x1=i_L, x2=v_C. dx/dt = A*x + B*V_in, y = v_C */
    StateSpace ss; ss_init(&ss, 2, 1, 1);
    SS_A(&ss, 0, 0) = -R / L;   SS_A(&ss, 0, 1) = -1.0 / L;
    SS_A(&ss, 1, 0) =  1.0 / C; SS_A(&ss, 1, 1) =  0.0;
    SS_B(&ss, 0, 0) =  1.0 / L; SS_B(&ss, 1, 0) =  0.0;
    SS_C(&ss, 0, 0) =  0.0;     SS_C(&ss, 0, 1) =  1.0;
    SS_D(&ss, 0, 0) =  0.0;
    return ss;
}

StateSpace ss_from_parallel_rlc(double R, double L, double C)
{
    /* State: x1=v_C, x2=i_L. Driven by I_in, output = v_C */
    StateSpace ss; ss_init(&ss, 2, 1, 1);
    double rc = R * C;
    if (fabs(rc) < 1e-30) rc = 1e-30;
    SS_A(&ss, 0, 0) = -1.0/rc;  SS_A(&ss, 0, 1) = -1.0 / C;
    SS_A(&ss, 1, 0) =  1.0 / L;  SS_A(&ss, 1, 1) =  0.0;
    SS_B(&ss, 0, 0) =  1.0 / C;  SS_B(&ss, 1, 0) =  0.0;
    SS_C(&ss, 0, 0) =  1.0;      SS_C(&ss, 0, 1) =  0.0;
    SS_D(&ss, 0, 0) =  0.0;
    return ss;
}

StateSpace ss_from_two_loop_rlc(double R1, double L1, double R2, double L2, double C1)
{
    /* 3-state: x1=i_L1, x2=i_L2, x3=v_C1 */
    StateSpace ss; ss_init(&ss, 3, 1, 1);
    SS_A(&ss, 0, 0) = -R1/L1;  SS_A(&ss, 0, 1) =  0.0;     SS_A(&ss, 0, 2) = -1.0/L1;
    SS_A(&ss, 1, 0) =  0.0;    SS_A(&ss, 1, 1) = -R2/L2;    SS_A(&ss, 1, 2) =  1.0/L2;
    SS_A(&ss, 2, 0) =  1.0/C1; SS_A(&ss, 2, 1) = -1.0/C1;    SS_A(&ss, 2, 2) =  0.0;
    SS_B(&ss, 0, 0) =  1.0/L1; SS_B(&ss, 1, 0) =  0.0;       SS_B(&ss, 2, 0) =  0.0;
    SS_C(&ss, 0, 0) =  0.0;    SS_C(&ss, 0, 1) =  0.0;       SS_C(&ss, 0, 2) =  1.0;
    SS_D(&ss, 0, 0) =  0.0;
    return ss;
}

/* ===== Canonical Realizations (L5) ===== */

StateSpace ss_controllable_canonical(const double *num, size_t m,
                                      const double *den, size_t n)
{
    /* Companion form: A has 1 on superdiagonal, B=[0;...;0;1] */
    StateSpace ss; ss_init(&ss, n, 1, 1);
    for (size_t i = 0; i + 1 < n; i++) SS_A(&ss, i, i+1) = 1.0;
    for (size_t j = 0; j < n; j++) SS_A(&ss, n-1, j) = -den[j];
    SS_B(&ss, n-1, 0) = 1.0;
    for (size_t j = 0; j <= m && j < n; j++) SS_C(&ss, 0, j) = num[j];
    SS_D(&ss, 0, 0) = 0.0;
    return ss;
}

StateSpace ss_observable_canonical(const double *num, size_t m,
                                    const double *den, size_t n)
{
    /* Dual of controllable: A has 1 on subdiagonal, C=[0...0 1] */
    StateSpace ss; ss_init(&ss, n, 1, 1);
    for (size_t i = 0; i + 1 < n; i++) SS_A(&ss, i+1, i) = 1.0;
    for (size_t i = 0; i < n; i++) SS_A(&ss, i, n-1) = -den[i];
    for (size_t i = 0; i <= m && i < n; i++) SS_B(&ss, i, 0) = num[i];
    SS_C(&ss, 0, n-1) = 1.0;
    SS_D(&ss, 0, 0) = 0.0;
    return ss;
}

/* ===== State-Space to Transfer Function Evaluation (L5) ===== */

double complex ss_evaluate_tf(const StateSpace *ss, size_t input_idx,
                               size_t output_idx, double complex s)
{
    /* H(s) = C*(sI - A)^{-1}*B + D, evaluated at a single s */
    if (!ss || input_idx >= ss->m || output_idx >= ss->p) return 0.0;
    size_t n = ss->n;

    double complex *M = (double complex*)malloc(n*n*sizeof(double complex));
    double complex *rhs3 = (double complex*)malloc(n*sizeof(double complex));
    double complex *sol3 = (double complex*)calloc(n, sizeof(double complex));
    if (!M || !rhs3 || !sol3) {
        free(M); free(rhs3); free(sol3);
        return 0.0;
    }

    /* Build matrix (sI - A) and RHS = B column */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
            M[i*n + j] = (i == j) ? (s - SS_A(ss, i, j)) : (-SS_A(ss, i, j));
        rhs3[i] = SS_B(ss, i, input_idx);
    }

    /* Gaussian elimination with partial pivoting */
    for (size_t kk = 0; kk < n; kk++) {
        size_t mr = kk;
        double mv = cabs(M[kk*n + kk]);
        for (size_t i3 = kk + 1; i3 < n; i3++)
            if (cabs(M[i3*n + kk]) > mv)
                { mv = cabs(M[i3*n + kk]); mr = i3; }
        if (mv < 1e-15) continue;
        if (mr != kk) {
            for (size_t j3 = kk; j3 < n; j3++) {
                double complex tp = M[kk*n + j3];
                M[kk*n + j3] = M[mr*n + j3];
                M[mr*n + j3] = tp;
            }
            double complex tp = rhs3[kk];
            rhs3[kk] = rhs3[mr];
            rhs3[mr] = tp;
        }
        double complex pv = M[kk*n + kk];
        for (size_t i3 = kk + 1; i3 < n; i3++) {
            double complex fc = M[i3*n + kk] / pv;
            for (size_t j3 = kk; j3 < n; j3++)
                M[i3*n + j3] -= fc * M[kk*n + j3];
            rhs3[i3] -= fc * rhs3[kk];
        }
    }

    /* Back substitution */
    for (size_t ii3 = n; ii3 > 0; ii3--) {
        size_t i3 = ii3 - 1;
        double complex sm = rhs3[i3];
        for (size_t j3 = i3 + 1; j3 < n; j3++)
            sm -= M[i3*n + j3] * sol3[j3];
        sol3[i3] = sm / M[i3*n + i3];
    }

    /* Compute output: y = C*x + D */
    double complex y = SS_D(ss, output_idx, input_idx);
    for (size_t j3 = 0; j3 < n; j3++)
        y += SS_C(ss, output_idx, j3) * sol3[j3];

    free(M); free(rhs3); free(sol3);
    return y;
}

/* ===== Modal Parameter Extraction from Eigenvalues (L6) ===== */

ModalParameters ss_modal_from_eigenvalue(double complex lambda)
{
    ModalParameters mp;
    mp.natural_freq = cabs(lambda);
    if (mp.natural_freq > 1e-30)
        mp.damping_ratio = -creal(lambda) / mp.natural_freq;
    else
        mp.damping_ratio = 1.0;
    if (mp.damping_ratio > 1.0) mp.damping_ratio = 1.0;
    if (mp.damping_ratio < 0.0) mp.damping_ratio = 0.0;

    if (mp.damping_ratio > 1e-30 && mp.natural_freq > 1e-30)
        mp.time_constant = 1.0 / (mp.damping_ratio * mp.natural_freq);
    else
        mp.time_constant = DBL_MAX;
    return mp;
}

/* ===== Dominant Time Constant: Slowest Eigenvalue (L6) ===== */

double ss_dominant_time_constant(const StateSpace *ss)
{
    if (!ss || ss->n == 0) return 0.0;
    double complex evals[SS_MAX_ORDER];
    int nev = ss_eigenvalues(ss, evals);

    double min_abs_re = DBL_MAX;
    for (int i = 0; i < nev; i++) {
        double ar = fabs(creal(evals[i]));
        if (ar < min_abs_re && ar > 1e-30) min_abs_re = ar;
    }
    if (min_abs_re < 1e-30 || min_abs_re > DBL_MAX / 10.0)
        return DBL_MAX;
    /* 2% settling time = 4 / |min Re(lambda)| */
    return 4.0 / min_abs_re;
}

/* ===== ZOH Discretization via Van Loan Method (L5, 1978) ===== */

int ss_discretize_zoh(const StateSpace *ss, double Ts, DiscreteStateSpace *dss)
{
    if (!ss || !dss || Ts <= 0.0) return -1;
    size_t n = ss->n, m = ss->m, N = n + m;

    double *Ma = (double*)calloc(N*N, sizeof(double));
    double *eM = (double*)calloc(N*N, sizeof(double));
    double *pM = (double*)calloc(N*N, sizeof(double));
    double *tM = (double*)calloc(N*N, sizeof(double));
    if (!Ma || !eM || !pM || !tM) {
        free(Ma); free(eM); free(pM); free(tM);
        return -1;
    }

    /* Build augmented matrix: M = [A  B; 0  0] * Ts */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
            Ma[i*N + j] = SS_A(ss, i, j) * Ts;
        for (size_t j = 0; j < m; j++)
            Ma[i*N + (n + j)] = SS_B(ss, i, j) * Ts;
    }

    /* Identity matrix initialization */
    for (size_t i = 0; i < N; i++)
        eM[i*N + i] = pM[i*N + i] = 1.0;

    /* Taylor series: exp(M) = sum_{k=0}^20 M^k / k! */
    for (int kk = 1; kk <= 20; kk++) {
        /* tM = pM * Ma */
        for (size_t i = 0; i < N*N; i++) tM[i] = 0.0;
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j < N; j++)
                for (size_t r = 0; r < N; r++)
                    tM[i*N + j] += pM[i*N + r] * Ma[r*N + j];
        for (size_t i = 0; i < N*N; i++) pM[i] = tM[i];

        /* 1/k! factor */
        double inv_fact = 1.0;
        for (int f = 1; f <= kk; f++) inv_fact /= (double)f;

        /* eM += pM/k! */
        for (size_t i = 0; i < N*N; i++) eM[i] += pM[i] * inv_fact;
    }

    /* Extract Phi (top-left) and Gamma (top-right) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
            dss->Phi[i*n + j] = eM[i*N + j];
        for (size_t j = 0; j < m; j++)
            dss->Gamma[i*n + j] = eM[i*N + (n + j)];
    }

    /* C and D remain unchanged under ZOH discretization */
    for (size_t i = 0; i < ss->p * ss->n; i++) dss->Cd[i] = ss->C[i];
    for (size_t i = 0; i < ss->p * ss->m; i++) dss->Dd[i] = ss->D[i];

    dss->Ts = Ts;
    dss->n = n;
    dss->m = m;
    dss->p = ss->p;

    free(Ma); free(eM); free(pM); free(tM);
    return 0;
}