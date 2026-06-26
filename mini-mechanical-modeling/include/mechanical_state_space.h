/**
 * @file mechanical_state_space.h
 * @brief State-space representations of mechanical systems
 *
 * Knowledge Coverage:
 *   L1 - Definitions: state variables, ABCD matrices, canonical forms
 *   L2 - Core Concepts: controllability, observability, stability, eigenvalues
 *   L3 - Mathematical Structures: matrix algebra, eigen decomposition
 *   L4 - Fundamental Laws: state equations from Newton-Euler
 *   L5 - Engineering Methods: RK4 simulation, ZOH discretization, Pade matrix exp
 *
 * Reference:
 *   MIT 6.241 Dynamic Systems & Control
 *   Stanford ENGR207B Linear Control Systems
 *   Ogata, Modern Control Engineering (2010)
 *   Brogan, Modern Control Theory (1991)
 */

#ifndef MECHANICAL_STATE_SPACE_H
#define MECHANICAL_STATE_SPACE_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define SS_MAX_ORDER 16

/** Continuous-time state-space: dx/dt = A*x + B*u, y = C*x + D*u */
typedef struct {
    int n, m, p;
    double A[SS_MAX_ORDER * SS_MAX_ORDER];
    double B[SS_MAX_ORDER * SS_MAX_ORDER];
    double C[SS_MAX_ORDER * SS_MAX_ORDER];
    double D[SS_MAX_ORDER * SS_MAX_ORDER];
} MechSS;

/** Discrete-time: x[k+1] = Ad*x[k] + Bd*u[k] */
typedef struct {
    int n, m, p;
    double Ts;
    double Ad[SS_MAX_ORDER * SS_MAX_ORDER];
    double Bd[SS_MAX_ORDER * SS_MAX_ORDER];
    double Cd[SS_MAX_ORDER * SS_MAX_ORDER];
    double Dd[SS_MAX_ORDER * SS_MAX_ORDER];
} MechDSS;

/* Construction */
MechSS ss_sdof(double m, double b, double k);
MechSS ss_sdof_pv(double m, double b, double k);
MechSS ss_2dof(double m1, double m2, double k1, double k2, double kc,
                double b1, double b2, double bc);
MechSS ss_inverted_pendulum(double M, double m, double L, double g);
MechSS ss_rotational_sdof(double J, double b, double k);
MechSS ss_geared_rotational(double Jm, double Jl, double bm, double bl,
                             double km, double kl, double N);

/* Matrix operations */
void mat_vec_mult(const double *A, int rows, int cols, const double *x, double *y);
void mat_mat_mult(const double *Am, int m, int k, const double *B, int n, double *C);
void mat_transpose(const double *A, int rows, int cols, double *At);
double mat_determinant(const double *A, int n);
int mat_inverse(const double *A, int n, double *Ainv);
int mat_solve(const double *A, int n, const double *b, double *x);

/* Controllability / Observability */
int ss_ctrb_matrix(const MechSS *ss, double *Ctrb);
int ss_obsv_matrix(const MechSS *ss, double *Obsv);
int ss_ctrb_rank(const MechSS *ss);
int ss_obsv_rank(const MechSS *ss);
int ss_is_controllable(const MechSS *ss);
int ss_is_observable(const MechSS *ss);

/* Simulation */
void ss_rk4_step(const MechSS *ss, const double *x, const double *u, double dt, double *xn);
int ss_simulate_rk4(const MechSS *ss, double *x, const double *u, double t0,
                     double dt, int steps, double *t_out, double *x_out, double *y_out);

/* Eigenvalues */
int ss_eigenvalues(const MechSS *ss, double complex *evals);
int ss_is_stable(const MechSS *ss);
double ss_dominant_time_constant(const MechSS *ss);
void ss_dominant_mode(const MechSS *ss, double *wn, double *zeta);

/* State transition and discretization */
int ss_state_transition(const MechSS *ss, double t, double *Phi);
int ss_discretize_zoh(const MechSS *ss, double Ts, MechDSS *dss);
int ss_discretize_foh(const MechSS *ss, double Ts, MechDSS *dss);

/* Gramians */
int ss_gramian_ctrb(const MechSS *ss, double *Wc);
int ss_gramian_obsv(const MechSS *ss, double *Wo);
int ss_hankel_sv(const MechSS *ss, double *hsv);

/* Modal analysis */
void ss_eigenvalue_to_modal(double complex eval, double *wn, double *zeta);
int ss_to_modal_form(const MechSS *ss, MechSS *modal_ss, double complex *evals);
MechSS ss_modal_sdof(double wn, double zeta);

/* TF from/to SS */
MechSS ss_from_tf(const double complex *num, int num_ord,
                   const double complex *den, int den_ord);
MechSS ss_controllable_canonical(const double complex *num, int num_ord,
                                  const double complex *den, int den_ord);
MechSS ss_observable_canonical(const double complex *num, int num_ord,
                                const double complex *den, int den_ord);
double complex ss_evaluate_tf(const MechSS *ss, double complex s);

/** Initial condition response: x(t) = expm(A*t)*x0. O(n^3). */
void ss_initial_response(const MechSS *ss, const double *x0, double t, double *x_t);

/** Impulse response energy: integral of y^2 dt = x0"*Wo*x0. O(n^2). */
double ss_impulse_energy(const MechSS *ss, const double *x0);

#endif /* MECHANICAL_STATE_SPACE_H */
