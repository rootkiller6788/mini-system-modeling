#ifndef MECHATRONIC_STATE_SPACE_H
#define MECHATRONIC_STATE_SPACE_H

/**
 * @file mechatronic_state_space.h
 * @brief L3-L5 State-Space Modeling for Mechatronic Systems
 *
 * Builds state-space representations for multi-domain mechatronic systems
 * combining motor dynamics, transmission mechanics, and structural modes.
 *
 * Key References:
 *  - Ogata "Modern Control Engineering" (2010)
 *  - Isermann, R. "Mechatronic Systems: Fundamentals" (2005)
 *  - Spong et al. "Robot Modeling and Control" (2006)
 *
 * 9-Layer Coverage:
 *   L3 — State-space matrices, canonical forms, eigenvalue analysis
 *   L4 — Conservation laws embedded in A,B,C,D matrices
 *   L5 — Controllability, observability, discretization, simulation
 */

#include "mechatronic_definitions.h"
#include "electromechanical_coupling.h"
#include <stddef.h>

/*===========================================================================*/
/* L3: State-Space Linear System Definition                                  */
/*===========================================================================*/

/** Maximum state-space dimension supported */
#define SS_MAX_DIM 12

/** Continuous-time state-space system: dx/dt = A x + B u, y = C x + D u */
typedef struct {
    int    n_states;        /**< Number of states                         */
    int    n_inputs;        /**< Number of inputs                         */
    int    n_outputs;       /**< Number of outputs                        */
    double A[SS_MAX_DIM * SS_MAX_DIM]; /**< System matrix (row-major)    */
    double B[SS_MAX_DIM * SS_MAX_DIM]; /**< Input matrix                 */
    double C[SS_MAX_DIM * SS_MAX_DIM]; /**< Output matrix                */
    double D[SS_MAX_DIM * SS_MAX_DIM]; /**< Feedthrough matrix           */
} LinearStateSpace;

/*===========================================================================*/
/* L3: DC Motor State-Space Models                                           */
/*===========================================================================*/

/** Build 2nd-order DC motor state-space: [I_a, ω]^T
 *
 *  A = [-R_a/L_a   -K_e/L_a ]
 *      [  K_t/J     -B/J    ]
 *  B = [1/L_a   0   ]
 *      [ 0     -1/J  ]
 *  C = [0  1]  (speed output)
 *  D = [0  0]
 *
 *  @param motor  Motor parameters
 *  @param J_total Total reflected inertia at motor shaft [kg·m²]
 *  @param B_total Total damping at motor shaft [N·m·s/rad]
 *  @return LinearStateSpace with 2 states, 2 inputs (V_a, T_L), 1 output (ω)
 *
 *  Ref: Ogata (2010) Example 3.6 */
LinearStateSpace ss_dc_motor_second_order(const DCMotorParams *motor,
                                           double J_total, double B_total);

/** Build 3rd-order DC motor state-space: [I_a, ω, θ]^T
 *  Adds position (rotor angle) as third state.
 *
 *  A = [-R_a/L_a  -K_e/L_a   0]
 *      [  K_t/J     -B/J     0]
 *      [    0         1      0]
 *
 *  @param motor  Motor parameters
 *  @param J_total Total inertia [kg·m²]
 *  @param B_total Total damping [N·m·s/rad]
 *  @return LinearStateSpace with 3 states (I, ω, θ), 2 inputs, 2 outputs */
LinearStateSpace ss_dc_motor_third_order(const DCMotorParams *motor,
                                          double J_total, double B_total);

/*===========================================================================*/
/* L3: Motor + Flexible Transmission State-Space                             */
/*===========================================================================*/

/** Build 4th-order state-space: motor + flexible gear + load.
 *
 *  States: [I_a, ω_m, ω_L, Δθ]^T
 *    where Δθ = θ_m/N - θ_L  (torsional deflection across gear)
 *
 *  Dynamics:
 *    L_a·dI_a/dt  = V_a - R_a·I_a - K_e·ω_m
 *    J_m·dω_m/dt  = K_t·I_a - B_m·ω_m - (K_s/N)·Δθ - (D_s/N)·(ω_m/N - ω_L)
 *    J_L·dω_L/dt  = K_s·Δθ + D_s·(ω_m/N - ω_L) - B_L·ω_L - T_L
 *    dΔθ/dt       = ω_m/N - ω_L
 *
 *  @param motor   Motor parameters
 *  @param J_m     Motor-side inertia [kg·m²]
 *  @param B_m     Motor-side damping [N·m·s/rad]
 *  @param gear    Gear train parameters
 *  @param J_L     Load inertia [kg·m²]
 *  @param B_L     Load damping [N·m·s/rad]
 *  @return 4-state system [I_a, ω_m, ω_L, Δθ], inputs [V_a, T_L]
 *
 *  Ref: Spong et al. "Robot Modeling and Control" Ch.6 */
LinearStateSpace ss_motor_flexible_gear(const DCMotorParams *motor,
                                         double J_m, double B_m,
                                         const GearTrain *gear,
                                         double J_L, double B_L);

/*===========================================================================*/
/* L3: Ball Screw Linear Axis State-Space                                    */
/*===========================================================================*/

/** Build state-space for ball screw axis (motor + screw + linear load)
 *
 *  States: [I_a, ω, x, v]^T
 *    ω = motor angular velocity
 *    x = linear position
 *    v = linear velocity
 *
 *  Screw converts: v = ω·pitch/(2π),  F = T·2π·η/pitch
 *
 *  @param motor   Motor parameters
 *  @param screw   Lead screw parameters
 *  @param M_load  Mass of linear load [kg]
 *  @param B_linear Linear guide damping [N·s/m]
 *  @return 4-state system, inputs [V_a, F_ext]
 *
 *  Ref: Altintas "Manufacturing Automation" (2012) Ch.4 */
LinearStateSpace ss_ball_screw_axis(const DCMotorParams *motor,
                                     const LeadScrew *screw,
                                     double M_load, double B_linear);

/*===========================================================================*/
/* L3: Belt Drive Axis State-Space                                           */
/*===========================================================================*/

/** Build state-space for belt-driven linear axis
 *
 *  States: [I_a, ω_m, ω_L]^T  (includes belt compliance as spring-damper)
 *
 *  @param motor   Motor parameters
 *  @param belt    Belt drive parameters
 *  @param M_carriage Carriage mass [kg]
 *  @param B_guide Linear guide damping [N·s/m]
 *  @return 3-state system */
LinearStateSpace ss_belt_drive_axis(const DCMotorParams *motor,
                                     const BeltDrive *belt,
                                     double M_carriage, double B_guide);

/*===========================================================================*/
/* L4: Controllability and Observability Analysis                            */
/*===========================================================================*/

/** Build controllability matrix: C = [B AB A²B ... A^(n-1)B]
 *  Stored as n × (n·m) matrix in row-major order.
 *  @param ss     State-space system
 *  @param C_mat  Output buffer of size n_states × (n_states*n_inputs)
 *  @return Rank of controllability matrix (via non-zero singular values)
 *
 *  Complexity: O(n³) for rank via Gaussian elimination
 *
 *  Ref: Kalman, R.E. "On the general theory of control systems" (1960) */
int ss_controllability_matrix(const LinearStateSpace *ss, double *C_mat);

/** Check if system is controllable (returns 1 if fully controllable)
 *  @param ss  State-space system
 *  @return 1 if controllable, 0 otherwise */
int ss_is_controllable(const LinearStateSpace *ss);

/** Build observability matrix: O = [C; CA; CA²; ...; CA^(n-1)]
 *  @param ss     State-space system
 *  @param O_mat  Output buffer of size (n_states*n_outputs) × n_states
 *  @return Rank of observability matrix */
int ss_observability_matrix(const LinearStateSpace *ss, double *O_mat);

/** Check if system is observable
 *  @param ss  State-space system
 *  @return 1 if observable, 0 otherwise */
int ss_is_observable(const LinearStateSpace *ss);

/*===========================================================================*/
/* L5: Transfer Function from State-Space                                    */
/*===========================================================================*/

/** Compute transfer function G(s) = C·(sI - A)^(-1)·B + D
 *  for a SISO system at a given frequency point.
 *
 *  Uses Cramer's rule: G(s) = [det(sI - A + B·C) - det(sI - A)] / det(sI - A)
 *  plus D for the full expression.
 *
 *  @param ss    SISO state-space system
 *  @param s_re  Real part of complex frequency
 *  @param s_im  Imaginary part of complex frequency
 *  @param G_re  Output: real part of G(s)
 *  @param G_im  Output: imaginary part of G(s)
 *  @return 0 on success, -1 if singular sI-A
 *
 *  Ref: Kailath "Linear Systems" (1980) */
int ss_frequency_response(const LinearStateSpace *ss,
                           double s_re, double s_im,
                           double *G_re, double *G_im);

/** Compute DC gain matrix: G_dc = -C·A^(-1)·B + D
 *  @param ss     State-space system
 *  @param G_dc   Output matrix of size n_outputs × n_inputs (row-major)
 *  @return 0 on success, -1 if A is singular */
int ss_dc_gain(const LinearStateSpace *ss, double *G_dc);

/** Compute eigenvalues of A via QR algorithm
 *  Real eigenvalues stored in eval_re[i], imaginary in eval_im[i]
 *  @param A        n×n matrix (row-major)
 *  @param n        dimension
 *  @param eval_re  Output: real parts of eigenvalues
 *  @param eval_im  Output: imaginary parts of eigenvalues
 *  @return 0 on convergence, -1 if max iterations exceeded
 *
 *  Complexity: O(n³)
 *  Ref: Golub & Van Loan "Matrix Computations" (2013) Ch.7 */
int ss_eigenvalues(const double *A, int n, double *eval_re, double *eval_im);

/** Determine stability from eigenvalues
 *  @param eval_re  Real parts
 *  @param eval_im  Imaginary parts
 *  @param n        Dimension
 *  @return 1 if all eigenvalues have negative real part, 0 otherwise */
int ss_is_stable(const double *eval_re, const double *eval_im, int n);

/** Natural frequency and damping ratio for a pair of complex eigenvalues
 *  For λ = σ ± jω_d:
 *    ω_n = √(σ² + ω_d²)
 *    ζ = -σ / ω_n
 *  @param sigma    Real part of eigenvalue
 *  @param omega_d  Imaginary part (damped natural frequency)
 *  @param omega_n  Output: natural frequency [rad/s]
 *  @param zeta     Output: damping ratio */
void ss_mode_properties(double sigma, double omega_d,
                         double *omega_n, double *zeta);

/*===========================================================================*/
/* L5: Time-Domain Simulation (Runge-Kutta 4)                                */
/*===========================================================================*/

/** Simulate state-space system forward using RK4
 *  dx/dt = A x + B u, y = C x + D u
 *
 *  @param ss      State-space system
 *  @param x0      Initial state [n_states]
 *  @param u       Input array [n_inputs], constant over step
 *  @param dt      Time step [s]
 *  @param x_next  Output: state at t+dt [n_states]
 *  @param y_next  Output: output at t+dt [n_outputs] (may be NULL)
 *
 *  Complexity: O(n² + n·m) per step
 *  Ref: Isermann (2005) "Mechatronic Systems" */
void ss_rk4_step(const LinearStateSpace *ss,
                  const double *x, const double *u, double dt,
                  double *x_next, double *y_next);

/** Simulate trajectory for N steps, returns state and output history
 *  @param ss       State-space system
 *  @param x0       Initial state [n_states]
 *  @param u_history Input history [N_steps * n_inputs] row-major (may be NULL → zero)
 *  @param N_steps  Number of time steps
 *  @param dt       Step size [s]
 *  @param x_history Output state trajectory [N_steps * n_states]
 *  @param y_history Output output trajectory [N_steps * n_outputs]
 *
 *  Complexity: O(N·n²) */
void ss_simulate(const LinearStateSpace *ss,
                  const double *x0, const double *u_history,
                  int N_steps, double dt,
                  double *x_history, double *y_history);

/*===========================================================================*/
/* L5: Discretization (ZOH)                                                  */
/*===========================================================================*/

/** Zero-Order Hold discretization:
 *  x_{k+1} = Φ·x_k + Γ·u_k
 *  y_k = C·x_k + D·u_k
 *
 *  where Φ = exp(A·Ts), Γ = A^(-1)·(Φ - I)·B
 *
 *  Uses Pade [8/8] approximation for matrix exponential.
 *
 *  @param ss_cont   Continuous-time system
 *  @param Ts        Sampling period [s]
 *  @param ss_disc   Output: discrete-time system (same dimensions)
 *  @return 0 on success
 *
 *  Complexity: O(n³)
 *  Ref: Van Loan (1978), Higham "Functions of Matrices" (2008) */
int ss_discretize_zoh(const LinearStateSpace *ss_cont, double Ts,
                       LinearStateSpace *ss_disc);

/*===========================================================================*/
/* L5: Cascaded Controller Design                                            */
/*===========================================================================*/

/** Design current PI controller for DC motor
 *
 *  Desired closed-loop bandwidth ω_c_rad_s.
 *  Plant: G_i(s) = 1 / (L_a·s + R_a)
 *
 *  PI: Kp_cur + Ki_cur/s
 *  Pole placement: Kp_cur = L_a·ω_c,  Ki_cur = R_a·ω_c
 *
 *  @param motor   Motor parameters
 *  @param omega_c Desired current loop bandwidth [rad/s]
 *  @param Kp      Output: proportional gain [V/A]
 *  @param Ki      Output: integral gain [V/(A·s)] */
void design_current_pi(const DCMotorParams *motor, double omega_c,
                        double *Kp, double *Ki);

/** Design velocity PI controller (cascaded outside current loop)
 *
 *  Approximate current loop as first-order: 1/(τ_cur·s + 1)
 *  Plant (speed): K_t / (J_total·s + B_total)
 *
 *  PI poles placed via symmetrical optimum or bandwidth specification.
 *
 *  @param motor   Motor parameters
 *  @param J_total Total reflected inertia [kg·m²]
 *  @param B_total Total reflected damping [N·m·s/rad]
 *  @param omega_c Desired velocity loop bandwidth [rad/s]
 *  @param tau_cur Current loop time constant [s]
 *  @param Kp      Output: proportional gain [(N·m)/(rad/s)]
 *  @param Ki      Output: integral gain [N·m/(rad)] */
void design_velocity_pi(const DCMotorParams *motor,
                         double J_total, double B_total,
                         double omega_c, double tau_cur,
                         double *Kp, double *Ki);

/** Design position P controller (cascaded outside velocity loop)
 *
 *  Approximate velocity loop as 1/(τ_vel·s + 1)
 *  Plant (position): 1/s
 *
 *  P-gain: Kp_pos = ω_c_pos (proportional only for servo)
 *
 *  @param omega_c_pos Desired position loop bandwidth [rad/s]
 *  @param Kp          Output: position proportional gain [1/s] */
void design_position_p(double omega_c_pos, double *Kp);

#endif /* MECHATRONIC_STATE_SPACE_H */
