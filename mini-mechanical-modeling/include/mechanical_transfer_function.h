/**
 * @file mechanical_transfer_function.h
 * @brief Transfer functions for mechanical systems - SDOF, 2-DOF, base excitation,
 *        rotating unbalance, vibration absorber, seismic sensor.
 *
 * Knowledge Coverage:
 *   L1 - Definitions: mechanical TF, resonance, anti-resonance, dynamic stiffness
 *   L2 - Core Concepts: frequency response, Bode representation, vibration modes
 *   L3 - Mathematical Structures: rational polynomials, complex s-domain, pole-zero
 *   L4 - Fundamental Laws: D"Alembert principle, Newton-Euler in s-domain
 *   L5 - Engineering Methods: TF derivation from ODE, cascading, feedback
 *
 * Reference:
 *   MIT 2.004 Dynamics and Control II
 *   Stanford ENGR105
 *   Ogata, Modern Control Engineering (2010)
 *   Inman, Engineering Vibration (2014)
 */

#ifndef MECHANICAL_TRANSFER_FUNCTION_H
#define MECHANICAL_TRANSFER_FUNCTION_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define TF_MAX_ORDER 8

/** Transfer function H(s) = N(s)/D(s) as rational polynomial */
typedef struct {
    int num_order;
    int den_order;
    double complex num[TF_MAX_ORDER + 1];
    double complex den[TF_MAX_ORDER + 1];
    double dc_gain;
    int is_strictly_proper;
} MechanicalTF;

/* ===== SDOF Transfer Functions ===== */

/** Displacement/Force: H(s)=X(s)/F(s)=1/(m*s^2+b*s+k). O(1). */
MechanicalTF tf_sdof_displacement(double mass, double damping, double stiffness);

/** Velocity/Force: H(s)=s/(m*s^2+b*s+k). O(1). */
MechanicalTF tf_sdof_velocity(double mass, double damping, double stiffness);

/** Acceleration/Force: H(s)=s^2/(m*s^2+b*s+k). O(1). */
MechanicalTF tf_sdof_acceleration(double mass, double damping, double stiffness);

/** Base excitation: X(s)/Y(s)=(b*s+k)/(m*s^2+b*s+k). O(1). */
MechanicalTF tf_base_excitation_displacement(double mass, double damping, double stiffness);

/** Base excitation transmitted force. O(1). */
MechanicalTF tf_base_excitation_force(double mass, double damping, double stiffness);

/** Rotating unbalance. O(1). */
MechanicalTF tf_rotating_unbalance_displacement(double total_mass, double damping, double stiffness);

/* ===== 2-DOF Transfer Functions ===== */

MechanicalTF tf_2dof_direct(double m1, double m2, double k1, double k2, double kc);
MechanicalTF tf_2dof_cross(double m1, double m2, double k1, double k2, double kc);
MechanicalTF tf_2dof_damped_direct(double m1, double m2, double k1, double k2,
                                    double kc, double b1, double b2, double bc);
MechanicalTF tf_vibration_absorber(double M, double ma, double K, double ka,
                                    double B, double ba);

/* ===== Special Transfer Functions ===== */

MechanicalTF tf_seismic_accelerometer(double mass, double damping, double stiffness);
MechanicalTF tf_seismic_velocimeter(double mass, double damping, double stiffness);
MechanicalTF tf_vibration_isolator(double mass, double damping, double stiffness);

/* ===== TF Operations ===== */

double complex tf_eval(const MechanicalTF *tf, double complex s);
double tf_mag(const MechanicalTF *tf, double angular_freq);
double tf_mag_db(const MechanicalTF *tf, double angular_freq);
double tf_phase(const MechanicalTF *tf, double angular_freq);
double complex poly_eval(const double complex *coeffs, int order, double complex s);

/** Cascade: H(s)=H1(s)*H2(s) via polynomial convolution. O(n1*n2). */
MechanicalTF tf_cascade(const MechanicalTF *tf1, const MechanicalTF *tf2);

/** Parallel: H(s)=H1(s)+H2(s). O(max(n1,n2)). */
MechanicalTF tf_parallel(const MechanicalTF *tf1, const MechanicalTF *tf2);

/** Negative feedback: H1/(1+H1*H2). O(n1*n2). */
MechanicalTF tf_neg_feedback(const MechanicalTF *tf1, const MechanicalTF *tf2);

/** Positive feedback: H1/(1-H1*H2). */
MechanicalTF tf_pos_feedback(const MechanicalTF *tf1, const MechanicalTF *tf2);

/* ===== Frequency Response Analysis ===== */

void tf_sdof_poles(double m, double b, double k, double complex *p1, double complex *p2);
int tf_sdof_is_overdamped(double m, double b, double k);
int tf_sdof_is_critically_damped(double m, double b, double k);
int tf_sdof_is_underdamped(double m, double b, double k);
double tf_steady_state_amp(const MechanicalTF *tf, double force_amp, double w);
double tf_steady_state_phase(const MechanicalTF *tf, double w);
int tf_is_proper(const MechanicalTF *tf);
double tf_estimate_damping(const MechanicalTF *tf);
double tf_estimate_natural_freq(const MechanicalTF *tf);
double tf_dc_gain(const MechanicalTF *tf);
double tf_step_response(const MechanicalTF *tf, double t);

/* ===== Modal TF Decomposition ===== */

void tf_2dof_modal_decomposition(double m1, double m2, double k1, double k2, double kc,
                                  double phi[4], double omega[2], double zeta[2],
                                  MechanicalTF modal_tf[2]);

/** Frequency response data structure */
typedef struct {
    int n_points;
    double *freq;
    double *mag_db;
    double *phase_deg;
} FreqResp;

void tf_freq_response(const MechanicalTF *tf, double w_start, double w_end,
                      int n_pts, FreqResp *fr);
void freq_resp_free(FreqResp *fr);

#endif /* MECHANICAL_TRANSFER_FUNCTION_H */
