/**
 * @file multi_dof_systems.h
 * @brief Multi-degree-of-freedom mechanical systems - modal analysis, eigenvalue problems
 *
 * Knowledge Coverage:
 *   L1 - Definitions: MDOF, mass/stiffness/damping matrices, mode shapes,
 *        natural frequencies, modal mass/stiffness, FRF
 *   L2 - Core Concepts: modal superposition, orthogonality, Rayleigh damping,
 *        dynamic stiffness, frequency response functions
 *   L3 - Mathematical Structures: generalized eigenvalue problem,
 *        matrix symmetry, spectral decomposition, modal coordinates
 *   L4 - Fundamental Laws: D"Alembert principle, Lagrange equations for MDOF
 *   L5 - Engineering Methods: modal analysis, static condensation (Guyan),
 *        modal participation factors, effective modal mass
 *   L8 - Advanced: Guyan reduction, MAC criterion, Newmark-beta integration
 *
 * Reference:
 *   MIT 2.060 Structural Dynamics
 *   Stanford ME 351A
 *   Ewins, Modal Testing: Theory and Application (2000)
 *   Craig & Kurdila, Fundamentals of Structural Dynamics (2006)
 */

#ifndef MULTI_DOF_SYSTEMS_H
#define MULTI_DOF_SYSTEMS_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define MDOF_MAX_DOF 24

/** MDOF system: M*x_ddot + B*x_dot + K*x = F(t) */
typedef struct {
    int n_dof;
    double M[MDOF_MAX_DOF * MDOF_MAX_DOF];
    double K[MDOF_MAX_DOF * MDOF_MAX_DOF];
    double B[MDOF_MAX_DOF * MDOF_MAX_DOF];
    double F[MDOF_MAX_DOF];
} MDOFSystem;

/** Modal parameters */
typedef struct {
    double natural_freq;
    double damping_ratio;
    double modal_mass;
    double modal_stiffness;
    double participation_factor;
    double effective_modal_mass;
    double mode_shape[MDOF_MAX_DOF];
} ModalParams;

/* Construction */
void mdof_init(MDOFSystem *sys, int n_dof);
void mdof_set_mass(MDOFSystem *sys, int i, int j, double v);
void mdof_set_stiffness(MDOFSystem *sys, int i, int j, double v);
void mdof_set_damping(MDOFSystem *sys, int i, int j, double v);
void mdof_set_force(MDOFSystem *sys, int i, double v);
void mdof_build_2dof(MDOFSystem *sys, double m1, double m2,
                      double k1, double k2, double kc,
                      double b1, double b2, double bc);
void mdof_build_chain(MDOFSystem *sys, const double *m, const double *k, int n);
void mdof_build_3dof(MDOFSystem *sys, const double *m, const double *k, const double *b);

/* Eigenvalue problem */
int mdof_eigen_solve(const MDOFSystem *sys, double *evals, double *evecs);
void mdof_natural_freqs(const MDOFSystem *sys, double *freqs, int *n);
void mdof_natural_freqs_hz(const MDOFSystem *sys, double *freqs, int *n);

/* Modal parameters */
void mdof_modal_params(const MDOFSystem *sys, const double *evec,
                        int mode_idx, ModalParams *p);
void mdof_mass_normalize(const MDOFSystem *sys, double *mode);
void mdof_unity_normalize(double *mode, int n);

/* Orthogonality */
double mdof_mass_ortho(const MDOFSystem *sys, const double *phi_i, const double *phi_j);
double mdof_stiff_ortho(const MDOFSystem *sys, const double *phi_i, const double *phi_j);
int mdof_verify_ortho(const MDOFSystem *sys, const double *evecs, int n_modes, double tol);

/* Frequency Response Functions */
double complex mdof_receptance_frf(const MDOFSystem *sys, const ModalParams *modes,
                                    int n_modes, int out_dof, int in_dof, double w);
double complex mdof_mobility_frf(const MDOFSystem *sys, const ModalParams *modes,
                                  int n_modes, int out_dof, int in_dof, double w);
double complex mdof_accelerance_frf(const MDOFSystem *sys, const ModalParams *modes,
                                     int n_modes, int out_dof, int in_dof, double w);
void mdof_frf_matrix(const MDOFSystem *sys, const ModalParams *modes, int n_modes,
                      double w, double complex *H);

/* Rayleigh damping */
void rayleigh_damping_coeffs(double w1, double w2, double z1, double z2,
                              double *alpha, double *beta);
double rayleigh_modal_damping(double alpha, double beta, double w);

/* Guyan reduction */
int mdof_guyan_reduce(const MDOFSystem *sys, const int *master, int nm, MDOFSystem *reduced);

/* MAC */
double mdof_mac(const double *phi_a, const double *phi_b, int n);
void mdof_mac_matrix(const double *evecs, int n_dof, int n_modes, double *mac);

/* Participation factors */
double mdof_participation_factor(const MDOFSystem *sys, const double *mode);
double mdof_effective_mass(double gamma, double m_modal);
double mdof_cumulative_mass_frac(const ModalParams *modes, int n, double M_total);

/* Modal superposition */
void mdof_modal_superposition(const MDOFSystem *sys, const ModalParams *modes,
                               int n_modes, const double *q, double *x);
double duhamel_integral(double wn, double zeta, double m,
                         const double *F, const double *t, int n);

/* Energy */
double mdof_power_input(double F, double v);
double mdof_kinetic_energy(const MDOFSystem *sys, const double *v);
double mdof_potential_energy(const MDOFSystem *sys, const double *x);
double mdof_total_energy(const MDOFSystem *sys, const double *x, const double *v);
double mdof_dissipated_power(const MDOFSystem *sys, const double *v);

/* Newmark-beta integration */
void mdof_newmark_step(const MDOFSystem *sys, double dt, double beta_nm, double gamma_nm,
                        double *x, double *v, double *a, const double *F_new);

#endif /* MULTI_DOF_SYSTEMS_H */
