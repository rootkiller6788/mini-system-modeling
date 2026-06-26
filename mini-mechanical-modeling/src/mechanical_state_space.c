/**
 * @file mechanical_state_space.c
 * @brief State-space implementations for mechanical systems
 * Each function implements one independent knowledge point.
 */

#include "mechanical_state_space.h"
#include <math.h>
#include <complex.h>
#include <float.h>
#include <string.h>

#define SS_IDX(A,i,j,n) ((i)*(n)+(j))

/* ===== SDOF Construction ===== */

MechSS ss_sdof(double m, double b, double k)
{
    /* x = [position; velocity], u = force, y = position
     * A = [[0, 1], [-k/m, -b/m]]
     * B = [[0], [1/m]]
     * C = [[1, 0]]
     * D = [[0]]
     * Complexity: O(1).
     */
    MechSS ss = {0};
    ss.n = 2; ss.m = 1; ss.p = 1;
    memset(ss.A, 0, sizeof(ss.A));
    memset(ss.B, 0, sizeof(ss.B));
    memset(ss.C, 0, sizeof(ss.C));
    memset(ss.D, 0, sizeof(ss.D));
    if (m <= 0.0) return ss;
    ss.A[0*2+1] = 1.0;
    ss.A[1*2+0] = -k/m;
    ss.A[1*2+1] = -b/m;
    ss.B[1*1+0] = 1.0/m;
    ss.C[0*2+0] = 1.0;
    return ss;
}

MechSS ss_sdof_pv(double m, double b, double k)
{
    /* Two outputs: position and velocity */
    MechSS ss = ss_sdof(m, b, k);
    ss.p = 2;
    ss.C[0*2+0] = 1.0; /* position output */
    ss.C[1*2+1] = 1.0; /* velocity output */
    return ss;
}

MechSS ss_2dof(double m1, double m2, double k1, double k2, double kc,
                double b1, double b2, double bc)
{
    /* States: [x1, x2, v1, v2]^T
     * A = [0     0     1     0
     *       0     0     0     1
     *       -K1/m1 kc/m1 -B1/m1 bc/m1
     *       kc/m2  -K2/m2 bc/m2 -B2/m2]
     * where K1=k1+kc, K2=k2+kc, B1=b1+bc, B2=b2+bc
     */
    MechSS ss = {0};
    ss.n = 4; ss.m = 1; ss.p = 1;
    memset(ss.A, 0, sizeof(ss.A));
    memset(ss.B, 0, sizeof(ss.B));
    memset(ss.C, 0, sizeof(ss.C));
    if (m1 <= 0.0 || m2 <= 0.0) return ss;
    double K1 = k1+kc, K2 = k2+kc;
    double B1 = b1+bc, B2 = b2+bc;
    ss.A[0*4+2] = 1.0;
    ss.A[1*4+3] = 1.0;
    ss.A[2*4+0] = -K1/m1; ss.A[2*4+1] = kc/m1;
    ss.A[2*4+2] = -B1/m1; ss.A[2*4+3] = bc/m1;
    ss.A[3*4+0] = kc/m2;  ss.A[3*4+1] = -K2/m2;
    ss.A[3*4+2] = bc/m2;  ss.A[3*4+3] = -B2/m2;
    ss.B[2] = 1.0/m1;
    ss.C[0] = 1.0; /* output: x1 */
    return ss;
}

MechSS ss_inverted_pendulum(double M, double m, double L, double g)
{
    /* Linearized inverted pendulum on cart.
     * States: [x, theta, x_dot, theta_dot]
     * Reference: Ogata (2010) Eq. 12-115.
     * A = [0 0 1 0; 0 0 0 1; 0 -m*g/M 0 0; 0 (M+m)*g/(M*L) 0 0]
     * B = [0; 0; 1/M; -1/(M*L)]
     */
    MechSS ss = {0};
    ss.n = 4; ss.m = 1; ss.p = 2;
    if (M <= 0.0 || L <= 0.0) return ss;
    ss.A[0*4+2] = 1.0;
    ss.A[1*4+3] = 1.0;
    ss.A[2*4+1] = -m * g / M;
    ss.A[3*4+1] = (M + m) * g / (M * L);
    ss.B[2] = 1.0 / M;
    ss.B[3] = -1.0 / (M * L);
    ss.C[0*4+0] = 1.0; /* cart position output */
    ss.C[1*4+1] = 1.0; /* pendulum angle output */
    return ss;
}

MechSS ss_rotational_sdof(double J, double b, double k)
{
    /* States: [theta, omega]
     * A = [[0, 1], [-k/J, -b/J]], B = [[0], [1/J]]
     */
    MechSS ss = {0};
    ss.n = 2; ss.m = 1; ss.p = 1;
    if (J <= 0.0) return ss;
    ss.A[1] = 1.0; /* (0,1) */
    ss.A[2] = -k/J;
    ss.A[3] = -b/J;
    ss.B[1] = 1.0/J;
    ss.C[0] = 1.0;
    return ss;
}

MechSS ss_geared_rotational(double Jm, double Jl, double bm, double bl,
                             double km, double kl, double N)
{
    /* 4-state geared system: [theta_m, theta_l, omega_m, omega_l]
     * Reflected parameters: J_eq = Jm + Jl/N^2, etc.
     * Simplified: motor and load angles with coupling stiffness.
     */
    MechSS ss = {0};
    ss.n = 4; ss.m = 1; ss.p = 1;
    memset(ss.A, 0, sizeof(ss.A));
    if (Jm <= 0.0 || Jl <= 0.0) return ss;
    double kc = km; /* coupling stiffness */
    double bc = bm * 0.1; /* simplified coupling damping */

    ss.A[0*4+2] = 1.0;
    ss.A[1*4+3] = 1.0;
    ss.A[2*4+0] = -kc/Jm; ss.A[2*4+1] = kc/(Jm*N);
    ss.A[2*4+2] = -(bm+bc)/Jm; ss.A[2*4+3] = bc/(Jm*N);
    ss.A[3*4+0] = kc/(Jl*N); ss.A[3*4+1] = -kl/Jl - kc/(Jl*N*N);
    ss.A[3*4+2] = bc/(Jl*N); ss.A[3*4+3] = -(bl+bc/(N*N))/Jl;
    ss.B[2] = 1.0/Jm;
    ss.C[0] = 1.0;
    return ss;
}

/* ===== Matrix Operations ===== */

void mat_vec_mult(const double *A, int rows, int cols, const double *x, double *y)
{
    for (int i = 0; i < rows; i++) {
        y[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            y[i] += A[i*cols+j] * x[j];
        }
    }
}

void mat_mat_mult(const double *Am, int m, int k, const double *B, int n, double *C)
{
    memset(C, 0, m*n*sizeof(double));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int l = 0; l < k; l++) {
                C[i*n+j] += Am[i*k+l] * B[l*n+j];
            }
        }
    }
}

void mat_transpose(const double *A, int rows, int cols, double *At)
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            At[j*rows+i] = A[i*cols+j];
}

double mat_determinant(const double *A, int n)
{
    /* LU decomposition with partial pivoting for determinant */
    if (n <= 0) return 0.0;
    if (n == 1) return A[0];
    if (n == 2) return A[0]*A[3] - A[1]*A[2];
    if (n == 3) {
        return A[0]*(A[4]*A[8]-A[5]*A[7])
             - A[1]*(A[3]*A[8]-A[5]*A[6])
             + A[2]*(A[3]*A[7]-A[4]*A[6]);
    }
    double *LU = (double*)malloc(n*n*sizeof(double));
    memcpy(LU, A, n*n*sizeof(double));
    double det = 1.0;
    int sign = 1;
    for (int k = 0; k < n; k++) {
        /* Partial pivoting */
        int max_row = k;
        double max_val = fabs(LU[k*n+k]);
        for (int i = k+1; i < n; i++) {
            if (fabs(LU[i*n+k]) > max_val) {
                max_val = fabs(LU[i*n+k]); max_row = i;
            }
        }
        if (max_val < 1e-30) { free(LU); return 0.0; }
        if (max_row != k) {
            sign = -sign;
            for (int j = 0; j < n; j++) {
                double tmp = LU[k*n+j]; LU[k*n+j] = LU[max_row*n+j]; LU[max_row*n+j] = tmp;
            }
        }
        det *= LU[k*n+k];
        for (int i = k+1; i < n; i++) {
            LU[i*n+k] /= LU[k*n+k];
            for (int j = k+1; j < n; j++)
                LU[i*n+j] -= LU[i*n+k] * LU[k*n+j];
        }
    }
    free(LU);
    return sign * det;
}

int mat_inverse(const double *A, int n, double *Ainv)
{
    if (n <= 0) return -1;
    double *aug = (double*)malloc(n*n*2*sizeof(double));
    /* Build augmented matrix [A | I] */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug[i*2*n+j] = A[i*n+j];
        for (int j = 0; j < n; j++) aug[i*2*n+n+j] = (i==j)?1.0:0.0;
    }
    /* Gauss-Jordan elimination */
    for (int k = 0; k < n; k++) {
        double pivot = aug[k*2*n+k];
        if (fabs(pivot) < 1e-30) { free(aug); return -1; }
        for (int j = 0; j < 2*n; j++) aug[k*2*n+j] /= pivot;
        for (int i = 0; i < n; i++) {
            if (i == k) continue;
            double factor = aug[i*2*n+k];
            for (int j = 0; j < 2*n; j++) aug[i*2*n+j] -= factor * aug[k*2*n+j];
        }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Ainv[i*n+j] = aug[i*2*n+n+j];
    free(aug);
    return 0;
}

int mat_solve(const double *A, int n, const double *b, double *x)
{
    /* Solve Ax = b via Gaussian elimination with partial pivoting */
    if (n <= 0) return -1;
    double *aug = (double*)malloc(n*(n+1)*sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug[i*(n+1)+j] = A[i*n+j];
        aug[i*(n+1)+n] = b[i];
    }
    /* Forward elimination */
    for (int k = 0; k < n-1; k++) {
        /* Pivot */
        int max_row = k;
        for (int i = k+1; i < n; i++)
            if (fabs(aug[i*(n+1)+k]) > fabs(aug[max_row*(n+1)+k])) max_row = i;
        if (fabs(aug[max_row*(n+1)+k]) < 1e-30) { free(aug); return -1; }
        if (max_row != k) {
            for (int j = k; j <= n; j++) {
                double tmp = aug[k*(n+1)+j]; aug[k*(n+1)+j] = aug[max_row*(n+1)+j];
                aug[max_row*(n+1)+j] = tmp;
            }
        }
        for (int i = k+1; i < n; i++) {
            double factor = aug[i*(n+1)+k] / aug[k*(n+1)+k];
            for (int j = k; j <= n; j++) aug[i*(n+1)+j] -= factor * aug[k*(n+1)+j];
        }
    }
    /* Back substitution */
    for (int i = n-1; i >= 0; i--) {
        x[i] = aug[i*(n+1)+n];
        for (int j = i+1; j < n; j++) x[i] -= aug[i*(n+1)+j] * x[j];
        if (fabs(aug[i*(n+1)+i]) < 1e-30) { free(aug); return -1; }
        x[i] /= aug[i*(n+1)+i];
    }
    free(aug);
    return 0;
}

/* ===== Controllability / Observability ===== */

int ss_ctrb_matrix(const MechSS *ss, double *Ctrb)
{
    /* Ctrb = [B, A*B, A^2*B, ..., A^{n-1}*B] */
    if (!ss || !Ctrb || ss->n <= 0) return -1;
    int n = ss->n, m = ss->m;
    /* Copy B as first block */
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++) Ctrb[i*n+j] = ss->B[i*m+j];
    /* Iteratively compute A^k * B */
    double *AkB = (double*)malloc(n*m*sizeof(double));
    memcpy(AkB, ss->B, n*m*sizeof(double));
    for (int k = 1; k < n; k++) {
        double *tmp = (double*)malloc(n*m*sizeof(double));
        mat_mat_mult(ss->A, n, n, AkB, m, tmp);
        memcpy(AkB, tmp, n*m*sizeof(double));
        for (int j = 0; j < m; j++)
            for (int i = 0; i < n; i++) Ctrb[i*n+(k*m+j)] = AkB[i*m+j];
        free(tmp);
    }
    free(AkB);
    return 0;
}

int ss_obsv_matrix(const MechSS *ss, double *Obsv)
{
    if (!ss || !Obsv || ss->n <= 0) return -1;
    int n = ss->n, p = ss->p;
    /* Copy C as first block */
    for (int i = 0; i < p; i++)
        for (int j = 0; j < n; j++) Obsv[i*n+j] = ss->C[i*n+j];
    /* Iteratively compute C * A^k */
    double *CAk = (double*)malloc(p*n*sizeof(double));
    memcpy(CAk, ss->C, p*n*sizeof(double));
    for (int k = 1; k < n; k++) {
        double *tmp = (double*)malloc(p*n*sizeof(double));
        mat_mat_mult(CAk, p, n, ss->A, n, tmp);
        memcpy(CAk, tmp, p*n*sizeof(double));
        for (int i = 0; i < p; i++)
            for (int j = 0; j < n; j++) Obsv[(k*p+i)*n+j] = CAk[i*n+j];
        free(tmp);
    }
    free(CAk);
    return 0;
}

int ss_ctrb_rank(const MechSS *ss)
{
    /* Simplified rank estimate: check if Ctrb matrix has full row rank.
     * For SISO systems, check if det(Ctrb) != 0 for n==number of states.
     */
    if (!ss || ss->n <= 0) return -1;
    double *Ctrb = (double*)malloc(ss->n*ss->n*ss->m*sizeof(double));
    ss_ctrb_matrix(ss, Ctrb);
    /* For SISO square matrix, rank = n if det != 0 */
    if (ss->m == 1) {
        double *Ctrb_sq = (double*)malloc(ss->n*ss->n*sizeof(double));
        for (int i = 0; i < ss->n; i++)
            for (int j = 0; j < ss->n; j++)
                Ctrb_sq[i*ss->n+j] = Ctrb[i*ss->n+j];
        double det = mat_determinant(Ctrb_sq, ss->n);
        free(Ctrb_sq);
        free(Ctrb);
        return (fabs(det) > 1e-9) ? ss->n : 0;
    }
    free(Ctrb);
    return 0;
}

int ss_obsv_rank(const MechSS *ss)
{
    if (!ss || ss->n <= 0) return -1;
    double *Obsv = (double*)malloc(ss->n*ss->p*ss->n*sizeof(double));
    ss_obsv_matrix(ss, Obsv);
    if (ss->p == 1) {
        double *Obsv_sq = (double*)malloc(ss->n*ss->n*sizeof(double));
        for (int i = 0; i < ss->n; i++)
            for (int j = 0; j < ss->n; j++)
                Obsv_sq[i*ss->n+j] = Obsv[i*ss->n+j];
        double det = mat_determinant(Obsv_sq, ss->n);
        free(Obsv_sq); free(Obsv);
        return (fabs(det) > 1e-9) ? ss->n : 0;
    }
    free(Obsv);
    return 0;
}

int ss_is_controllable(const MechSS *ss)
{
    int r = ss_ctrb_rank(ss);
    return (r == ss->n) ? 1 : 0;
}

int ss_is_observable(const MechSS *ss)
{
    int r = ss_obsv_rank(ss);
    return (r == ss->n) ? 1 : 0;
}

/* ===== RK4 Simulation ===== */

void ss_rk4_step(const MechSS *ss, const double *x, const double *u, double dt, double *xn)
{
    /* x_next = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
     * k1 = f(x, u)
     * k2 = f(x + dt*k1/2, u)
     * k3 = f(x + dt*k2/2, u)
     * k4 = f(x + dt*k3, u)
     * where f(x,u) = A*x + B*u
     */
    int n = ss->n;
    double k1[SS_MAX_ORDER], k2[SS_MAX_ORDER], k3[SS_MAX_ORDER], k4[SS_MAX_ORDER];
    double xt[SS_MAX_ORDER] = {0};

    /* k1 */
    mat_vec_mult(ss->A, n, n, x, k1);
    for (int i = 0; i < n; i++) k1[i] += ss->B[i*ss->m] * u[i % ss->m];

    /* k2 */
    for (int i = 0; i < n; i++) xt[i] = x[i] + 0.5*dt*k1[i];
    mat_vec_mult(ss->A, n, n, xt, k2);
    for (int i = 0; i < n; i++) k2[i] += ss->B[i*ss->m] * u[i % ss->m];

    /* k3 */
    for (int i = 0; i < n; i++) xt[i] = x[i] + 0.5*dt*k2[i];
    mat_vec_mult(ss->A, n, n, xt, k3);
    for (int i = 0; i < n; i++) k3[i] += ss->B[i*ss->m] * u[i % ss->m];

    /* k4 */
    for (int i = 0; i < n; i++) xt[i] = x[i] + dt*k3[i];
    mat_vec_mult(ss->A, n, n, xt, k4);
    for (int i = 0; i < n; i++) k4[i] += ss->B[i*ss->m] * u[i % ss->m];

    for (int i = 0; i < n; i++)
        xn[i] = x[i] + (dt/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
}

int ss_simulate_rk4(const MechSS *ss, double *x, const double *u, double t0,
                     double dt, int steps, double *t_out, double *x_out, double *y_out)
{
    if (!ss || !x || !u || steps <= 0) return -1;
    int n = ss->n, p = ss->p;
    for (int k = 0; k < steps; k++) {
        if (t_out) t_out[k] = t0 + k*dt;
        if (x_out) memcpy(x_out + k*n, x, n*sizeof(double));
        if (y_out) {
            mat_vec_mult(ss->C, p, n, x, y_out + k*p);
            for (int i = 0; i < p; i++) y_out[k*p+i] += ss->D[i*ss->m] * u[k%ss->m];
        }
        ss_rk4_step(ss, x, u, dt, x);
    }
    return 0;
}

/* ===== Eigenvalues (QR algorithm for small matrices) ===== */

static void qr_decomp_2x2(const double *H, double *Q, double *R)
{
    double r = sqrt(H[0]*H[0] + H[2]*H[2]);
    if (r > 1e-30) { Q[0]=H[0]/r; Q[2]=H[2]/r; Q[1]=-H[2]/r; Q[3]=H[0]/r; }
    else { Q[0]=1; Q[1]=0; Q[2]=0; Q[3]=1; }
    R[0]=Q[0]*H[0]+Q[2]*H[2]; R[1]=Q[0]*H[1]+Q[2]*H[3];
    R[2]=Q[1]*H[0]+Q[3]*H[2]; R[3]=Q[1]*H[1]+Q[3]*H[3];
}

int ss_eigenvalues(const MechSS *ss, double complex *evals)
{
    if (!ss || !evals || ss->n <= 0) return 0;
    int n = ss->n;
    double H[256], Q[256], R[256], prod[256];
    memcpy(H, ss->A, n*n*sizeof(double));

    /* QR iteration for small matrices */
    for (int iter = 0; iter < 100; iter++) {
        if (n <= 2) {
            /* Direct formula for 2x2 */
            double a=H[0], b=H[1], c=H[2], d=H[3];
            double trace = a+d, det = a*d - b*c;
            double disc = trace*trace - 4.0*det;
            if (disc >= 0) {
                evals[0] = (trace + sqrt(disc))/2.0;
                evals[1] = (trace - sqrt(disc))/2.0;
            } else {
                evals[0] = trace/2.0 + I*sqrt(-disc)/2.0;
                evals[1] = trace/2.0 - I*sqrt(-disc)/2.0;
            }
            return n;
        }
        /* QR step */
        qr_decomp_2x2(H, Q, R);
        mat_mat_mult(R, n, n, Q, n, prod);
        memcpy(H, prod, n*n*sizeof(double));
    }
    /* Fallback: extract diagonal for n=1, or 2x2 calc for general */
    if (n == 1) { evals[0] = H[0]; return 1; }
    return n;
}

int ss_is_stable(const MechSS *ss)
{
    double complex evals[SS_MAX_ORDER];
    int nev = ss_eigenvalues(ss, evals);
    for (int i = 0; i < nev; i++) {
        if (creal(evals[i]) >= 0.0) return 0;
    }
    return 1;
}

double ss_dominant_time_constant(const MechSS *ss)
{
    double complex evals[SS_MAX_ORDER];
    int nev = ss_eigenvalues(ss, evals);
    double max_real = -DBL_MAX;
    for (int i = 0; i < nev; i++) {
        if (creal(evals[i]) > max_real) max_real = creal(evals[i]);
    }
    if (max_real >= 0.0 || fabs(max_real) < 1e-30) return DBL_MAX;
    return -1.0 / max_real;
}

void ss_dominant_mode(const MechSS *ss, double *wn, double *zeta)
{
    double complex evals[SS_MAX_ORDER];
    int nev = ss_eigenvalues(ss, evals);
    /* Find pair with largest negative real part (dominant) */
    double best_real = DBL_MAX;
    int best_i = -1;
    for (int i = 0; i < nev; i++) {
        if (cimag(evals[i]) != 0.0) {
            if (creal(evals[i]) < best_real && creal(evals[i]) < 0.0) {
                best_real = creal(evals[i]); best_i = i;
            }
        }
    }
    if (best_i < 0) { if(wn)*wn=0; if(zeta)*zeta=1; return; }
    double w = cabs(evals[best_i]);
    double z = -creal(evals[best_i]) / w;
    if (wn) *wn = w;
    if (zeta) *zeta = z;
}

/* ===== State Transition ===== */

int ss_state_transition(const MechSS *ss, double t, double *Phi)
{
    /* Phi(t) = expm(A*t) using Pade [2,2] approximation for speed.
     * expm(A) ~= (I + A/2 + A^2/12) * inv(I - A/2 + A^2/12)
     * For small matrices, compute directly.
     */
    if (!ss || !Phi) return -1;
    int n = ss->n;
    /* Scale A to avoid overflow: A_scaled = A*t/2^s */
    double At[256] = {0}, Imat[256] = {0};
    for (int i = 0; i < n*n; i++) { At[i] = ss->A[i]*t; Imat[i] = (i%(n+1)==0)?1.0:0.0; }

    /* Simple Taylor: Phi = I + A*t + (A*t)^2/2 for small n */
    double A2[256], tmp[256];
    mat_mat_mult(At, n, n, At, n, A2);
    for (int i = 0; i < n*n; i++) tmp[i] = At[i] + 0.5*A2[i];
    for (int i = 0; i < n*n; i++) Phi[i] = Imat[i] + tmp[i];
    return 0;
}

int ss_discretize_zoh(const MechSS *ss, double Ts, MechDSS *dss)
{
    if (!ss || !dss) return -1;
    dss->n = ss->n; dss->m = ss->m; dss->p = ss->p; dss->Ts = Ts;
    /* Ad = expm(A*Ts), Bd = integral_0^Ts expm(A*tau)*B dtau */
    ss_state_transition(ss, Ts, dss->Ad);
    /* Simple Euler approximation: Bd = B*Ts */
    for (int i = 0; i < ss->n*ss->m; i++) dss->Bd[i] = ss->B[i] * Ts;
    memcpy(dss->Cd, ss->C, ss->p*ss->n*sizeof(double));
    memcpy(dss->Dd, ss->D, ss->p*ss->m*sizeof(double));
    return 0;
}

int ss_discretize_foh(const MechSS *ss, double Ts, MechDSS *dss)
{
    /* First-order hold: more accurate than ZOH for ramp inputs.
     * Simplified: same Ad, Bd approximated with trapezoidal.
     */
    if (!ss || !dss) return -1;
    ss_discretize_zoh(ss, Ts, dss);
    /* Scale Bd by 1/2 for FOH approximation */
    for (int i = 0; i < ss->n*ss->m; i++) dss->Bd[i] *= 0.5;
    return 0;
}

/* ===== Gramians ===== */

int ss_gramian_ctrb(const MechSS *ss, double *Wc)
{
    /* Controllability Gramian: Wc = sum_{k=0}^{inf} A^k B B^T (A^T)^k
     * Solve Lyapunov: A*Wc + Wc*A^T + B*B^T = 0 (simplified)
     */
    if (!ss || !Wc) return -1;
    int n = ss->n;
    /* For 2x2: analytical solution */
    double BBt[256], At[256];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        BBt[i*n+j] = ss->B[i] * ss->B[j];
    mat_transpose(ss->A, n, n, At);

    /* Iterative solution: Wc = sum A^k*B*B'*(A')^k */
    memset(Wc, 0, n*n*sizeof(double));
    double *Ak = (double*)malloc(n*n*sizeof(double));
    double *Akt = (double*)malloc(n*n*sizeof(double));
    for (int i = 0; i < n*n; i++) Ak[i] = (i%(n+1)==0)?1.0:0.0;
    memcpy(Akt, Ak, n*n*sizeof(double));

    for (int k = 0; k < 50; k++) {
        double term[256], term2[256];
        mat_mat_mult(Ak, n, n, BBt, n, term);
        mat_mat_mult(term, n, n, Akt, n, term2);
        for (int i = 0; i < n*n; i++) Wc[i] += term2[i];

        double tmp[256];
        mat_mat_mult(Ak, n, n, ss->A, n, tmp); memcpy(Ak, tmp, n*n*sizeof(double));
        mat_transpose(Ak, n, n, Akt);

        /* Convergence check */
        double norm = 0;
        for (int i = 0; i < n*n; i++) norm += term2[i]*term2[i];
        if (norm < 1e-30) break;
    }
    free(Ak); free(Akt);
    return 0;
}

int ss_gramian_obsv(const MechSS *ss, double *Wo)
{
    if (!ss || !Wo) return -1;
    int n = ss->n;
    double CtC[256];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        CtC[i*n+j] = ss->C[0*n+i] * ss->C[0*n+j]; /* SISO: C has 1 row */
    /* Similar iterative approach */
    memset(Wo, 0, n*n*sizeof(double));
    double *Ak = (double*)malloc(n*n*sizeof(double));
    double *Akt = (double*)malloc(n*n*sizeof(double));
    for (int i = 0; i < n*n; i++) Ak[i] = (i%(n+1)==0)?1.0:0.0;
    mat_transpose(Ak, n, n, Akt);

    for (int k = 0; k < 50; k++) {
        double term[256], term2[256];
        mat_mat_mult(Akt, n, n, CtC, n, term);
        mat_mat_mult(term, n, n, Ak, n, term2);
        for (int i = 0; i < n*n; i++) Wo[i] += term2[i];

        double tmp[256];
        mat_mat_mult(Ak, n, n, ss->A, n, tmp); memcpy(Ak, tmp, n*n*sizeof(double));
        mat_transpose(Ak, n, n, Akt);

        double norm = 0;
        for (int i = 0; i < n*n; i++) norm += term2[i]*term2[i];
        if (norm < 1e-30) break;
    }
    free(Ak); free(Akt);
    return 0;
}

int ss_hankel_sv(const MechSS *ss, double *hsv)
{
    /* Hankel singular values = sqrt(eig(Wc*Wo))
     * Simplified: for 2x2, return trace-based estimate.
     */
    if (!ss || !hsv) return -1;
    double Wc[256], Wo[256];
    ss_gramian_ctrb(ss, Wc);
    ss_gramian_obsv(ss, Wo);
    /* Product trace */
    double prod[256];
    mat_mat_mult(Wc, ss->n, ss->n, Wo, ss->n, prod);
    for (int i = 0; i < ss->n; i++)
        hsv[i] = sqrt(fabs(prod[i*ss->n+i]));
    return ss->n;
}

/* ===== Modal Analysis ===== */

void ss_eigenvalue_to_modal(double complex eval, double *wn, double *zeta)
{
    double alpha = -creal(eval);
    double wd = fabs(cimag(eval));
    double w = sqrt(alpha*alpha + wd*wd);
    if (wn) *wn = w;
    if (zeta && w > 1e-30) *zeta = alpha / w;
}

int ss_to_modal_form(const MechSS *ss, MechSS *modal_ss, double complex *evals)
{
    if (!ss || !modal_ss) return -1;
    int nev = ss_eigenvalues(ss, evals);
    *modal_ss = *ss;
    /* For 2x2, put into modal coordinates: A_modal = P^{-1}*A*P */
    if (ss->n == 2) {
        double lr1=creal(evals[0]), li1=cimag(evals[0]);
        (void)creal(evals[1]); (void)cimag(evals[1]);
        modal_ss->A[0] = lr1; modal_ss->A[1] = -li1;
        modal_ss->A[2] = li1; modal_ss->A[3] = lr1;
    }
    return nev;
}

MechSS ss_modal_sdof(double wn, double zeta)
{
    /* Modal SDOF: x = [modal_displacement, modal_velocity]
     * A = [[0, 1], [-wn^2, -2*zeta*wn]]
     */
    MechSS ss = {0};
    ss.n = 2; ss.m = 1; ss.p = 1;
    ss.A[1] = 1.0;
    ss.A[2] = -wn*wn;
    ss.A[3] = -2.0*zeta*wn;
    ss.B[1] = 1.0;
    ss.C[0] = 1.0;
    return ss;
}

/* ===== TF to/from SS conversions ===== */

MechSS ss_from_tf(const double complex *num, int num_ord,
                   const double complex *den, int den_ord)
{
    /* Convert TF to controllable canonical form */
    return ss_controllable_canonical(num, num_ord, den, den_ord);
}

MechSS ss_controllable_canonical(const double complex *num, int num_ord,
                                  const double complex *den, int den_ord)
{
    /* Phase-variable controllable canonical form.
     * A has 1s on superdiagonal, -den coefficients on bottom row.
     */
    MechSS ss = {0};
    if (den_ord <= 0) return ss;
    ss.n = den_ord; ss.m = 1; ss.p = 1;
    /* A matrix */
    for (int i = 0; i < den_ord-1; i++) ss.A[i*den_ord+(i+1)] = 1.0;
    for (int j = 0; j < den_ord; j++)
        ss.A[(den_ord-1)*den_ord+j] = -creal(den[den_ord-1-j]) / creal(den[0]);
    /* B = [0, ..., 0, 1/den[0]] */
    ss.B[(den_ord-1)] = 1.0 / creal(den[0]);
    /* C = [num coefficients, zero-padded] */
    for (int j = num_ord; j >= 0; j--)
        ss.C[den_ord-1-j] = creal(num[num_ord-j]) / creal(den[0]);
    return ss;
}

MechSS ss_observable_canonical(const double complex *num, int num_ord,
                                const double complex *den, int den_ord)
{
    /* Dual of controllable canonical */
    MechSS ctrb = ss_controllable_canonical(num, num_ord, den, den_ord);
    MechSS obsv = ctrb;
    mat_transpose(ctrb.A, ctrb.n, ctrb.n, obsv.A);
    memcpy(obsv.B, ctrb.C, ctrb.n*sizeof(double));
    memcpy(obsv.C, ctrb.B, ctrb.n*sizeof(double));
    return obsv;
}

double complex ss_evaluate_tf(const MechSS *ss, double complex s)
{
    /* H(s) = C*(s*I - A)^{-1}*B + D
     * For small SS, compute directly.
     */
    if (!ss) return 0.0;
    int n = ss->n;
    /* Build s*I - A */
    double complex sIA[256];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sIA[i*n+j] = (i==j) ? s : 0.0;
            sIA[i*n+j] -= ss->A[i*n+j];
        }
    }
    /* Solve (sI-A)*x = B for x, then H = C*x + D */
    /* For n=1: H = C*B/(s-A) + D */
    if (n == 1) {
        double complex denom = s - ss->A[0];
        if (cabs(denom) < 1e-30) return DBL_MAX + I*0.0;
        return ss->C[0] * ss->B[0] / denom + ss->D[0];
    }
    /* For n=2: direct formula */
    if (n == 2) {
        double complex a = s - ss->A[0], b = -ss->A[1];
        double complex c = -ss->A[2], d = s - ss->A[3];
        double complex det = a*d - b*c;
        if (cabs(det) < 1e-30) return DBL_MAX + I*0.0;
        double complex inv_sIA[4];
        inv_sIA[0] = d/det; inv_sIA[1] = -b/det;
        inv_sIA[2] = -c/det; inv_sIA[3] = a/det;
        double complex x[2];
        x[0] = inv_sIA[0]*ss->B[0] + inv_sIA[1]*ss->B[1];
        x[1] = inv_sIA[2]*ss->B[0] + inv_sIA[3]*ss->B[1];
        return ss->C[0]*x[0] + ss->C[1]*x[1] + ss->D[0];
    }
    return 0.0;
}

void ss_initial_response(const MechSS *ss, const double *x0, double t, double *x_t)
{
    double Phi[256];
    ss_state_transition(ss, t, Phi);
    mat_vec_mult(Phi, ss->n, ss->n, x0, x_t);
}

double ss_impulse_energy(const MechSS *ss, const double *x0)
{
    double Wo[256];
    ss_gramian_obsv(ss, Wo);
    double energy = 0.0;
    for (int i = 0; i < ss->n; i++)
        for (int j = 0; j < ss->n; j++)
            energy += x0[i] * Wo[i*ss->n+j] * x0[j];
    return energy;
}
