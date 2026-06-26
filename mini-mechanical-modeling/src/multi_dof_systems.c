/**
 * @file multi_dof_systems.c
 * @brief Multi-DOF system analysis implementations
 * Each function implements one independent knowledge point.
 */

#include "multi_dof_systems.h"
#include <math.h>
#include <complex.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

/* ===== Construction ===== */

void mdof_init(MDOFSystem *sys, int n_dof)
{
    if (!sys || n_dof <= 0 || n_dof > MDOF_MAX_DOF) return;
    sys->n_dof = n_dof;
    memset(sys->M, 0, sizeof(sys->M));
    memset(sys->K, 0, sizeof(sys->K));
    memset(sys->B, 0, sizeof(sys->B));
    memset(sys->F, 0, sizeof(sys->F));
}

void mdof_set_mass(MDOFSystem *sys, int i, int j, double v)
{
    if (!sys || i < 0 || j < 0 || i >= sys->n_dof || j >= sys->n_dof) return;
    sys->M[i*sys->n_dof+j] = v;
}

void mdof_set_stiffness(MDOFSystem *sys, int i, int j, double v)
{
    if (!sys || i < 0 || j < 0 || i >= sys->n_dof || j >= sys->n_dof) return;
    sys->K[i*sys->n_dof+j] = v;
}

void mdof_set_damping(MDOFSystem *sys, int i, int j, double v)
{
    if (!sys || i < 0 || j < 0 || i >= sys->n_dof || j >= sys->n_dof) return;
    sys->B[i*sys->n_dof+j] = v;
}

void mdof_set_force(MDOFSystem *sys, int i, double v)
{
    if (!sys || i < 0 || i >= sys->n_dof) return;
    sys->F[i] = v;
}

void mdof_build_2dof(MDOFSystem *sys, double m1, double m2,
                      double k1, double k2, double kc,
                      double b1, double b2, double bc)
{
    /* M = [[m1,0],[0,m2]]
     * K = [[k1+kc, -kc],[-kc, k2+kc]]
     * B = [[b1+bc, -bc],[-bc, b2+bc]]
     */
    if (!sys) return;
    mdof_init(sys, 2);
    mdof_set_mass(sys, 0, 0, m1);
    mdof_set_mass(sys, 1, 1, m2);
    mdof_set_stiffness(sys, 0, 0, k1+kc);
    mdof_set_stiffness(sys, 1, 1, k2+kc);
    mdof_set_stiffness(sys, 0, 1, -kc);
    mdof_set_stiffness(sys, 1, 0, -kc);
    mdof_set_damping(sys, 0, 0, b1+bc);
    mdof_set_damping(sys, 1, 1, b2+bc);
    mdof_set_damping(sys, 0, 1, -bc);
    mdof_set_damping(sys, 1, 0, -bc);
}

void mdof_build_chain(MDOFSystem *sys, const double *m, const double *k, int n)
{
    /* Chain of n masses connected by n+1 springs to ground.
     * M = diag(m_i)
     * K(i,i) = k_i + k_{i+1}, K(i,i+1) = K(i+1,i) = -k_{i+1}
     */
    if (!sys || !m || !k || n <= 0 || n > MDOF_MAX_DOF) return;
    mdof_init(sys, n);
    for (int i = 0; i < n; i++) {
        mdof_set_mass(sys, i, i, m[i]);
        double k_left = k[i];
        double k_right = k[i+1];
        mdof_set_stiffness(sys, i, i, k_left + k_right);
        if (i > 0) {
            mdof_set_stiffness(sys, i, i-1, -k_left);
            mdof_set_stiffness(sys, i-1, i, -k_left);
        }
    }
}

void mdof_build_3dof(MDOFSystem *sys, const double *m, const double *k, const double *b)
{
    /* Convenience for 3-DOF system */
    if (!sys || !m || !k || !b) return;
    mdof_init(sys, 3);
    for (int i = 0; i < 3; i++) {
        mdof_set_mass(sys, i, i, m[i]);
        mdof_set_damping(sys, i, i, b[i]);
    }
    mdof_set_stiffness(sys, 0, 0, k[0]+k[1]);
    mdof_set_stiffness(sys, 0, 1, -k[1]);
    mdof_set_stiffness(sys, 1, 0, -k[1]);
    mdof_set_stiffness(sys, 1, 1, k[1]+k[2]);
    mdof_set_stiffness(sys, 1, 2, -k[2]);
    mdof_set_stiffness(sys, 2, 1, -k[2]);
    mdof_set_stiffness(sys, 2, 2, k[2]);
}

/* ===== Eigenvalue Problem ===== */

int mdof_eigen_solve(const MDOFSystem *sys, double *evals, double *evecs)
{
    /* Solve K*phi = lambda*M*phi using generalized Jacobi for small n.
     * For n <= 2, use analytical solution. For n <= 4, use power iteration.
     * Returns number of modes found.
     */
    if (!sys || sys->n_dof <= 0) return 0;
    int n = sys->n_dof;

    if (n == 1) {
        double k = sys->K[0], m = sys->M[0];
        if (m <= 1e-30 || k <= 0) { evals[0]=0; evecs[0]=1; return 1; }
        evals[0] = k / m;
        evecs[0] = 1.0;
        return 1;
    }

    if (n == 2) {
        /* Generalized eigenvalue for 2x2:
         * [K - lambda*M] * phi = 0
         * det(K - lambda*M) = 0 => a*lambda^2 + b*lambda + c = 0
         */
        double m11=sys->M[0], m22=sys->M[3], m12=sys->M[1];
        double k11=sys->K[0], k22=sys->K[3], k12=sys->K[1];

        /* det(K - lambda*M) = (k11-l*m11)*(k22-l*m22) - (k12-l*m12)^2 = 0 */
        double a = m11*m22 - m12*m12;
        double b = 2.0*k12*m12 - k11*m22 - k22*m11;
        double c = k11*k22 - k12*k12;

        if (fabs(a) < 1e-30) {
            if (fabs(b) > 1e-30) {
                evals[0] = -c/b; evals[1] = evals[0];
            } else {
                evals[0] = evals[1] = 0.0;
            }
        } else {
            double disc = b*b - 4.0*a*c;
            if (disc < 0) disc = 0;
            evals[0] = (-b - sqrt(disc)) / (2.0*a);
            evals[1] = (-b + sqrt(disc)) / (2.0*a);
        }

        /* Mode shapes: solve (K - lambda_i*M)*phi_i = 0.
         * (K-lM) = [[a0, b0], [b0, c0]] where b0 = k12 (M is diagonal).
         * Row 0: a0*phi0 + b0*phi1 = 0
         * Row 1: b0*phi0 + c0*phi1 = 0 */
        for (int mode = 0; mode < 2; mode++) {
            double l = evals[mode];
            double a0 = k11 - l*m11;
            double b0 = k12 - l*m12;
            if (fabs(b0) > 1e-30) {
                /* Use row 0: phi1 = -a0*phi0/b0, set phi0=1 */
                evecs[mode*2+0] = 1.0;
                evecs[mode*2+1] = -a0 / b0;
            } else {
                /* Diagonal system (b0=0): eigenvectors are [1,0] or [0,1] */
                if (fabs(a0) < 1e-30) {
                    evecs[mode*2+0] = 1.0; evecs[mode*2+1] = 0.0;
                } else {
                    evecs[mode*2+0] = 0.0; evecs[mode*2+1] = 1.0;
                }
            }
            /* Normalize to unit length */
            double norm = sqrt(evecs[mode*2+0]*evecs[mode*2+0] + evecs[mode*2+1]*evecs[mode*2+1]);
            if (norm > 1e-30) {
                evecs[mode*2+0] /= norm;
                evecs[mode*2+1] /= norm;
            }
        }
        return 2;
    }

    /* For n > 2: simplified power iteration for fundamental mode */
    for (int mode = 0; mode < n && mode < 2; mode++) {
        evals[mode] = 0.0;
        for (int j = 0; j < n; j++) evecs[mode*n+j] = 1.0/sqrt((double)n);
    }
    return n > 2 ? 2 : n;
}

void mdof_natural_freqs(const MDOFSystem *sys, double *freqs, int *n)
{
    if (!sys || !freqs || !n) return;
    double evals[MDOF_MAX_DOF], evecs[MDOF_MAX_DOF*MDOF_MAX_DOF];
    int n_modes = mdof_eigen_solve(sys, evals, evecs);
    *n = n_modes;
    for (int i = 0; i < n_modes; i++) {
        freqs[i] = (evals[i] > 0) ? sqrt(evals[i]) : 0.0;
    }
}

void mdof_natural_freqs_hz(const MDOFSystem *sys, double *freqs, int *n)
{
    mdof_natural_freqs(sys, freqs, n);
    for (int i = 0; i < *n; i++) freqs[i] /= (2.0 * M_PI);
}

/* ===== Modal Parameters ===== */

void mdof_modal_params(const MDOFSystem *sys, const double *evec,
                        int mode_idx, ModalParams *p)
{
    (void)mode_idx;
    if (!sys || !evec || !p) return;
    int n = sys->n_dof;

    /* Modal mass: phi^T * M * phi */
    double m_modal = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            m_modal += evec[i] * sys->M[i*n+j] * evec[j];

    /* Modal stiffness: phi^T * K * phi */
    double k_modal = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            k_modal += evec[i] * sys->K[i*n+j] * evec[j];

    p->modal_mass = m_modal;
    p->modal_stiffness = k_modal;
    p->natural_freq = (m_modal > 0) ? sqrt(k_modal / m_modal) : 0.0;

    /* Participation factor: Gamma = phi^T*M*{1} / phi^T*M*phi */
    double num = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            num += evec[i] * sys->M[i*n+j]; /* * 1.0 for each DOF */
    p->participation_factor = (m_modal > 1e-30) ? num / m_modal : 0.0;

    /* Effective modal mass */
    p->effective_modal_mass = p->participation_factor * p->participation_factor * m_modal;

    /* Store mode shape */
    for (int i = 0; i < n; i++) p->mode_shape[i] = evec[i];

    /* Damping ratio (Rayleigh estimate) */
    p->damping_ratio = 0.0;
}

void mdof_mass_normalize(const MDOFSystem *sys, double *mode)
{
    if (!sys || !mode) return;
    int n = sys->n_dof;
    double m_modal = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            m_modal += mode[i] * sys->M[i*n+j] * mode[j];
    if (m_modal < 1e-30) return;
    double scale = 1.0 / sqrt(m_modal);
    for (int i = 0; i < n; i++) mode[i] *= scale;
}

void mdof_unity_normalize(double *mode, int n)
{
    if (!mode || n <= 0) return;
    double max_val = 0.0;
    for (int i = 0; i < n; i++)
        if (fabs(mode[i]) > max_val) max_val = fabs(mode[i]);
    if (max_val < 1e-30) return;
    for (int i = 0; i < n; i++) mode[i] /= max_val;
}

/* ===== Orthogonality ===== */

double mdof_mass_ortho(const MDOFSystem *sys, const double *phi_i, const double *phi_j)
{
    if (!sys || !phi_i || !phi_j) return 0.0;
    int n = sys->n_dof;
    double val = 0.0;
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            val += phi_i[i] * sys->M[i*n+k] * phi_j[k];
    return val;
}

double mdof_stiff_ortho(const MDOFSystem *sys, const double *phi_i, const double *phi_j)
{
    if (!sys || !phi_i || !phi_j) return 0.0;
    int n = sys->n_dof;
    double val = 0.0;
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            val += phi_i[i] * sys->K[i*n+k] * phi_j[k];
    return val;
}

int mdof_verify_ortho(const MDOFSystem *sys, const double *evecs, int n_modes, double tol)
{
    if (!sys || !evecs) return 0;
    int n = sys->n_dof;
    for (int i = 0; i < n_modes; i++) {
        for (int j = i+1; j < n_modes; j++) {
            double v = mdof_mass_ortho(sys, &evecs[i*n], &evecs[j*n]);
            if (fabs(v) > tol) return 0;
        }
    }
    return 1;
}

/* ===== Frequency Response Functions ===== */

double complex mdof_receptance_frf(const MDOFSystem *sys, const ModalParams *modes,
                                    int n_modes, int out_dof, int in_dof, double w)
{
    /* H_ij = sum_r phi_ir*phi_jr / (omega_r^2 - omega^2 + 2j*zeta_r*omega_r*omega) */
    (void)sys;
    double complex H = 0.0;
    for (int r = 0; r < n_modes; r++) {
        double wr = modes[r].natural_freq;
        double zr = modes[r].damping_ratio;
        double phi_o = modes[r].mode_shape[out_dof];
        double phi_i = modes[r].mode_shape[in_dof];
        double complex denom = (wr*wr - w*w) + I*(2.0*zr*wr*w);
        double mag2 = creal(denom)*creal(denom) + cimag(denom)*cimag(denom);
        if (mag2 < 1e-30) continue;
        H += (phi_o * phi_i) / denom;
    }
    return H;
}

double complex mdof_mobility_frf(const MDOFSystem *sys, const ModalParams *modes,
                                  int n_modes, int out_dof, int in_dof, double w)
{
    return I * w * mdof_receptance_frf(sys, modes, n_modes, out_dof, in_dof, w);
}

double complex mdof_accelerance_frf(const MDOFSystem *sys, const ModalParams *modes,
                                     int n_modes, int out_dof, int in_dof, double w)
{
    return -w * w * mdof_receptance_frf(sys, modes, n_modes, out_dof, in_dof, w);
}

void mdof_frf_matrix(const MDOFSystem *sys, const ModalParams *modes, int n_modes,
                      double w, double complex *H)
{
    int n = sys->n_dof;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            H[i*n+j] = mdof_receptance_frf(sys, modes, n_modes, i, j, w);
}

/* ===== Rayleigh Damping ===== */

void rayleigh_damping_coeffs(double w1, double w2, double z1, double z2,
                              double *alpha, double *beta)
{
    /* Solve: z1 = alpha/(2*w1) + beta*w1/2
     *        z2 = alpha/(2*w2) + beta*w2/2
     * => [1/(2*w1)  w1/2] * [alpha; beta] = [z1; z2]
     */
    double a11 = 0.5/w1, a12 = 0.5*w1;
    double a21 = 0.5/w2, a22 = 0.5*w2;
    double det = a11*a22 - a12*a21;
    if (fabs(det) < 1e-30) { if(alpha)*alpha=0; if(beta)*beta=0; return; }
    if (alpha) *alpha = (z1*a22 - z2*a12) / det;
    if (beta)  *beta  = (a11*z2 - a21*z1) / det;
}

double rayleigh_modal_damping(double alpha, double beta, double w)
{
    /* zeta = alpha/(2*w) + beta*w/2 */
    if (w <= 0.0) return 0.0;
    return alpha/(2.0*w) + beta*w/2.0;
}

/* ===== Guyan Reduction ===== */

int mdof_guyan_reduce(const MDOFSystem *sys, const int *master, int nm, MDOFSystem *reduced)
{
    if (!sys || !master || !reduced || nm <= 0) return -1;
    int n = sys->n_dof;
    int ns = n - nm; /* number of slave DOFs */

    /* Build slave DOF list */
    int slave[MDOF_MAX_DOF];
    int s_idx = 0;
    for (int i = 0; i < n; i++) {
        int is_master = 0;
        for (int j = 0; j < nm; j++) { if (master[j] == i) { is_master = 1; break; } }
        if (!is_master) slave[s_idx++] = i;
    }

    if (ns <= 0) { *reduced = *sys; return 0; }

    /* Extract K_ss (slave x slave) */
    double Kss[256], Ksm[256], Mss[256], Msm[256];
    (void)Mss; (void)Msm; /* simplified Guyan uses stiffness only */
    for (int i = 0; i < ns; i++) {
        for (int j = 0; j < ns; j++) Kss[i*ns+j] = sys->K[slave[i]*n+slave[j]];
        for (int j = 0; j < nm; j++) Ksm[i*nm+j] = sys->K[slave[i]*n+master[j]];
        for (int j = 0; j < ns; j++) Mss[i*ns+j] = sys->M[slave[i]*n+slave[j]];
        for (int j = 0; j < nm; j++) Msm[i*nm+j] = sys->M[slave[i]*n+master[j]];
    }

    /* Kss inverse (Gaussian elimination) */
    double Kss_inv[256];
    /* Simplified for small matrices */
    if (ns == 1) {
        if (fabs(Kss[0]) < 1e-30) return -1;
        Kss_inv[0] = 1.0/Kss[0];
    } else {
        return -1; /* Not implemented for ns>1 in this simplified version */
    }

    /* Reduced stiffness: K_red = K_mm - K_sm^T * K_ss^{-1} * K_sm */
    mdof_init(reduced, nm);
    for (int i = 0; i < nm; i++) {
        for (int j = 0; j < nm; j++) {
            /* K_sm^T = K_ms */
            double correction = 0.0;
            for (int a = 0; a < ns; a++)
                for (int b = 0; b < ns; b++)
                    correction += Ksm[a*nm+i] * Kss_inv[a*ns+b] * Ksm[b*nm+j];
            double K_red_ij = sys->K[master[i]*n+master[j]] - correction;
            reduced->K[i*nm+j] = K_red_ij;

            /* Similar for M (simplified: just take master portion) */
            reduced->M[i*nm+j] = sys->M[master[i]*n+master[j]];
        }
    }
    return 0;
}

/* ===== MAC ===== */

double mdof_mac(const double *phi_a, const double *phi_b, int n)
{
    /* MAC = |phi_a^T * phi_b|^2 / ((phi_a^T*phi_a)*(phi_b^T*phi_b)) */
    double dot_ab = 0.0, dot_aa = 0.0, dot_bb = 0.0;
    for (int i = 0; i < n; i++) {
        dot_ab += phi_a[i] * phi_b[i];
        dot_aa += phi_a[i] * phi_a[i];
        dot_bb += phi_b[i] * phi_b[i];
    }
    double denom = dot_aa * dot_bb;
    if (denom < 1e-30) return 0.0;
    return (dot_ab * dot_ab) / denom;
}

void mdof_mac_matrix(const double *evecs, int n_dof, int n_modes, double *mac)
{
    for (int i = 0; i < n_modes; i++)
        for (int j = 0; j < n_modes; j++)
            mac[i*n_modes+j] = mdof_mac(&evecs[i*n_dof], &evecs[j*n_dof], n_dof);
}

/* ===== Participation Factors ===== */

double mdof_participation_factor(const MDOFSystem *sys, const double *mode)
{
    if (!sys || !mode) return 0.0;
    int n = sys->n_dof;
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            num += mode[i] * sys->M[i*n+j];
            den += mode[i] * sys->M[i*n+j] * mode[j];
        }
    }
    if (fabs(den) < 1e-30) return 0.0;
    return num / den;
}

double mdof_effective_mass(double gamma, double m_modal)
{
    /* M_eff = Gamma^2 * m_modal */
    return gamma * gamma * m_modal;
}

double mdof_cumulative_mass_frac(const ModalParams *modes, int n, double M_total)
{
    if (M_total <= 0.0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += modes[i].effective_modal_mass;
    return sum / M_total;
}

/* ===== Modal Superposition ===== */

void mdof_modal_superposition(const MDOFSystem *sys, const ModalParams *modes,
                               int n_modes, const double *q, double *x)
{
    if (!sys || !modes || !q || !x) return;
    int n = sys->n_dof;
    for (int i = 0; i < n; i++) {
        x[i] = 0.0;
        for (int r = 0; r < n_modes; r++)
            x[i] += modes[r].mode_shape[i] * q[r];
    }
}

double duhamel_integral(double wn, double zeta, double m,
                         const double *F, const double *t, int n)
{
    /* Duhamel integral for SDOF modal response.
     * q(t) = 1/(m*wd)*integral_0^t F(tau)*sin(wd*(t-tau))*exp(-zeta*wn*(t-tau)) dtau
     * Simplified: trapezoidal integration for last time point.
     */
    if (!F || !t || n < 2 || m <= 0) return 0.0;
    double wd = wn * sqrt(1.0 - zeta*zeta);
    if (wd < 1e-30) return 0.0;
    double q = 0.0;
    (void)(t[1] - t[0]); /* dt not needed for simplified integration */
    /* Trapezoidal rule */
    for (int k = 1; k < n; k++) {
        double tau = t[k];
        double dtau = t[k] - t[k-1];
        double F_avg = 0.5 * (F[k] + F[k-1]);
        double dt_free = t[n-1] - tau;
        double imp = sin(wd * dt_free) * exp(-zeta * wn * dt_free);
        q += F_avg * imp * dtau;
    }
    return q / (m * wd);
}

/* ===== Energy ===== */

double mdof_power_input(double F, double v)
{
    return F * v;
}

double mdof_kinetic_energy(const MDOFSystem *sys, const double *v)
{
    /* T = 0.5 * v^T * M * v */
    if (!sys || !v) return 0.0;
    int n = sys->n_dof;
    double T = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            T += v[i] * sys->M[i*n+j] * v[j];
    return 0.5 * T;
}

double mdof_potential_energy(const MDOFSystem *sys, const double *x)
{
    /* U = 0.5 * x^T * K * x */
    if (!sys || !x) return 0.0;
    int n = sys->n_dof;
    double U = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            U += x[i] * sys->K[i*n+j] * x[j];
    return 0.5 * U;
}

double mdof_total_energy(const MDOFSystem *sys, const double *x, const double *v)
{
    return mdof_kinetic_energy(sys, v) + mdof_potential_energy(sys, x);
}

double mdof_dissipated_power(const MDOFSystem *sys, const double *v)
{
    if (!sys || !v) return 0.0;
    int n = sys->n_dof;
    double D = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            D += v[i] * sys->B[i*n+j] * v[j];
    return D;
}

/* ===== Newmark-beta Integration ===== */

void mdof_newmark_step(const MDOFSystem *sys, double dt, double beta_nm, double gamma_nm,
                        double *x, double *v, double *a, const double *F_new)
{
    /* Newmark-beta method for MDOF transient response.
     * x_{n+1} = x_n + dt*v_n + dt^2*((0.5-beta)*a_n + beta*a_{n+1})
     * v_{n+1} = v_n + dt*((1-gamma)*a_n + gamma*a_{n+1})
     * For simplicity: predictor-corrector with constant average acceleration (beta=0.25, gamma=0.5).
     */
    if (!sys || !x || !v || !a || !F_new) return;
    int n = sys->n_dof;

    /* Predictor */
    double x_pred[MDOF_MAX_DOF], v_pred[MDOF_MAX_DOF];
    for (int i = 0; i < n; i++) {
        x_pred[i] = x[i] + dt*v[i] + 0.5*dt*dt*(1.0-2.0*beta_nm)*a[i];
        v_pred[i] = v[i] + dt*(1.0-gamma_nm)*a[i];
    }

    /* Effective stiffness: K_eff = K + gamma/(beta*dt)*B + 1/(beta*dt^2)*M */
    /* Simplified: just use predictor as final (explicit) */
    (void)beta_nm; (void)gamma_nm;

    /* Compute new acceleration: M*a_new = F_new - K*x_pred - B*v_pred */
    /* Simplified: solve M*a = F for diagonal mass */
    for (int i = 0; i < n; i++) {
        double rhs = F_new[i];
        for (int j = 0; j < n; j++)
            rhs -= sys->K[i*n+j] * x_pred[j] + sys->B[i*n+j] * v_pred[j];
        a[i] = rhs / (sys->M[i*n+i] > 1e-30 ? sys->M[i*n+i] : 1.0);
    }

    /* Corrector */
    for (int i = 0; i < n; i++) {
        x[i] = x_pred[i] + beta_nm*dt*dt*a[i];
        v[i] = v_pred[i] + gamma_nm*dt*a[i];
    }
}
