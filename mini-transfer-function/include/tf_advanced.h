/**
 * @file tf_advanced.h
 * @brief Advanced Transfer Function Methods
 *
 * Knowledge Coverage:
 *   L5: Model order reduction (balanced truncation, Pade approximation,
 *       moment matching, modal truncation)
 *   L7: Application-specific TFs (DC motor, mass-spring-damper,
 *       quarter-car suspension, quadrotor, RLC, inverted pendulum, hydraulic)
 *   L8: Coprime factorization, Youla parameterization, robust stability
 *   L9: Fractional-order transfer functions (Oustaloup approximation)
 *
 * References:
 *   A.C. Antoulas, Approximation of Large-Scale Dynamical Systems, 2005
 *   S. Skogestad, I. Postlethwaite, Multivariable Feedback Control, 2nd ed.
 *   K. Zhou, J.C. Doyle, K. Glover, Robust and Optimal Control, 1996
 *   I. Podlubny, Fractional Differential Equations, 1999
 *   MIT 6.245 / Stanford ENGR210B / ETH 151-0563 / Caltech CDS 212
 */
#ifndef TF_ADVANCED_H
#define TF_ADVANCED_H

#include "transfer_function.h"
#include "tf_conversion.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L5: Model Order Reduction
 * ================================================================ */

/**
 * Balanced truncation model reduction.
 * Reduces a high-order TF to lower order while preserving dominant
 * input-output behavior based on Hankel singular values.
 *
 * Algorithm (Moore, 1981):
 *   1. Compute controllability Gramian P and observability Gramian Q
 *   2. Compute Cholesky factors: P = R*R^T, Q = L*L^T
 *   3. SVD of L^T*R = U*Sigma*V^T, Hankel singular values = diag(Sigma)
 *   4. Form balancing transformation, truncate smallest HSV states
 *
 * @param tf      Original high-order transfer function
 * @param order   Desired reduced model order
 * @return Reduced-order TF approximation, or NULL on error
 * Complexity: O(n^3)
 * Ref: Antoulas (2005) Ch. 7
 */
TransferFunction* tf_balanced_truncation(const TransferFunction *tf, int order);

/**
 * Pade approximation of time delay: exp(-Td*s) = N_p(s)/D_q(s).
 *
 * The [p/q] Pade approximant matches the Taylor series to order p+q:
 *   N_p(x) = sum_{k=0}^p C(p+q-k,q)*p!*(-x)^k/((p+q)!*k!)
 *   D_q(x) = sum_{k=0}^q C(p+q-k,p)*q!*x^k/((p+q)!*k!)
 * where x = Td*s.
 *
 * @param Td   Time delay (seconds, >0)
 * @param p    Numerator order
 * @param q    Denominator order
 * @return Rational approximation of exp(-Td*s)
 * Ref: Baker & Graves-Morris, Pade Approximants, 2nd ed.
 */
TransferFunction* tf_pade_approximation(double Td, int p, int q);

/**
 * Moment matching model reduction (Krylov subspace method).
 * Matches the first k Markov parameters: G(s) = sum m_i/s^{i+1}.
 * where m_i = C * A^i * B (impulse response moments).
 *
 * @param tf  Original transfer function
 * @param k   Number of moments to match
 * @return Reduced-order TF
 * Complex: O(k * n^2) via Arnoldi iteration
 */
TransferFunction* tf_moment_matching(const TransferFunction *tf, int k);

/**
 * Modal truncation: retain only the k most dominant poles.
 * Dominance = |residue / Re(pole)| - larger = more dominant.
 *
 * @param tf  Original transfer function
 * @param k   Number of modes to retain
 * @return Reduced TF with k dominant modes
 */
TransferFunction* tf_modal_truncation(const TransferFunction *tf, int k);

/* ================================================================
 * L7: Application-Specific Transfer Functions
 * ================================================================ */

/**
 * DC motor speed TF: Omega(s)/V_a(s).
 *
 * Electrical: L_a*di/dt + R_a*i = v_a - K_e*omega
 * Mechanical: J*domega/dt + b*omega = K_t*i
 *
 * G(s) = K_t / ((J*s+b)*(L_a*s+R_a) + K_t*K_e)
 *
 * L7 Application: DC motor control (Detroit, Toyota), robotics
 * Ref: Ogata (2010) Ch. 4, Franklin et al. (2019) Ch. 2
 */
TransferFunction* tf_dc_motor_speed(double J, double b, double K_t,
                                     double K_e, double R_a, double L_a);
TransferFunction* tf_dc_motor_position(double J, double b, double K_t,
                                        double K_e, double R_a, double L_a);

/**
 * Mass-spring-damper: G(s) = X(s)/F(s) = 1/(m*s^2 + b*s + k).
 *
 * L7 Application: Boeing 747 landing gear, automotive suspension
 */
TransferFunction* tf_mass_spring_damper(double m, double b, double k);

/**
 * Quarter-car suspension: sprung mass displacement / road input.
 *
 * G(s) = (Cs*Kt*s + Ks*Kt) / (Ms*Mus*s^4 + ...)
 *
 * L7 Application: Tesla Model S suspension, ISO 2631 ride comfort
 * Ref: Gillespie, Fundamentals of Vehicle Dynamics (1992)
 */
TransferFunction* tf_quarter_car_suspension(double Ms, double Mus,
                                              double Ks, double Cs, double Kt);

/**
 * Quadrotor attitude: G(s) = 1/(I*s^2). Double integrator.
 * L7 Application: Quadrotor UAV, SpaceX landing stability
 */
TransferFunction* tf_quadrotor_attitude(double I);

/**
 * RLC lowpass: G(s) = 1/(LC*s^2 + RC*s + 1).
 * L7 Application: Tesla inverter EMI filter design
 */
TransferFunction* tf_rlc_lowpass(double R, double L, double C);

/**
 * RLC bandpass: G(s) = RC*s/(LC*s^2 + RC*s + 1).
 */
TransferFunction* tf_rlc_bandpass(double R, double L, double C);

/**
 * Inverted pendulum (cart-pole) linearization:
 * G(s) = -1/(M*l) / (s^2 - (M+m)*g/(M*l))
 * Unstable pole at s = sqrt((M+m)*g/(M*l))
 *
 * L7 Application: Segway, Falcon 9 landing, bipedal robots
 */
TransferFunction* tf_inverted_pendulum(double M, double m, double l);

/**
 * Hydraulic actuator: G(s) = K/(s*(tau*s + 1)).
 * L7 Application: Boeing 747 flight controls, construction equipment
 */
TransferFunction* tf_hydraulic_actuator(double K, double tau);

/* ================================================================
 * L8: Coprime Factorization and Youla Parameterization
 * ================================================================ */

/**
 * Left coprime factorization: G = M^{-1} * N with M, N stable and proper.
 * Bezout identity: N*X + M*Y = 1 for some stable X, Y.
 */
int tf_left_coprime_factors(const TransferFunction *G,
                             TransferFunction **M, TransferFunction **N);
int tf_right_coprime_factors(const TransferFunction *G,
                              TransferFunction **N, TransferFunction **M);

/**
 * Youla parameterization of all stabilizing controllers:
 * K(Q) = (X + M*Q) / (Y - N*Q) where Q is any stable proper TF.
 * Ref: Youla, Jabr, Bongiorno (1976)
 */
TransferFunction* tf_youla_controller(const TransferFunction *G,
                                       const TransferFunction *Q);

/**
 * Verify controller K internally stabilizes plant G.
 */
tf_stability_t tf_is_stabilizing(const TransferFunction *G,
                                  const TransferFunction *K);

/* ================================================================
 * L9: Fractional-Order Transfer Functions (Research Frontier)
 * ================================================================ */

/**
 * Fractional integrator 1/s^alpha via Oustaloup recursive filter.
 * s^alpha = K * prod (s + w'_k)/(s + w_k), geometrically spaced.
 * @param alpha  Fractional order (0, 2)
 * @param N      Approximation order
 * @param wL,wH  Frequency bounds (rad/s)
 * L9 Frontier: Fractional calculus in control (CRONE, fractional PID)
 */
TransferFunction* tf_fractional_integrator(double alpha, int N,
                                             double wL, double wH);
TransferFunction* tf_fractional_differentiator(double alpha, int N,
                                                 double wL, double wH);

/**
 * Fractional PID: Kp + Ki/s^lambda + Kd*s^mu.
 * L9: 5 tuning parameters vs 3 for classical PID.
 */
TransferFunction* tf_fractional_pid(double Kp, double Ki, double lambda,
                                     double Kd, double mu,
                                     int N, double wL, double wH);

#ifdef __cplusplus
}
#endif
#endif
