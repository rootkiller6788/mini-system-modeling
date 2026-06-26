/**
 * @file electrical_state_space.h
 * @brief State-space models of electrical circuits
 *
 * Knowledge Coverage:
 *   L1 - Definitions: state variables, state vector, state matrix,
 *        input matrix, output matrix, feedthrough matrix, eigenvalues
 *   L2 - Core Concepts: state-space vs transfer function duality,
 *        controllability, observability, minimal realization
 *   L3 - Mathematical Structures: linear ODE systems, matrix exponential,
 *        eigenvalue decomposition, controllability Gramian
 *   L4 - Fundamental Laws: state equations from KVL/KCL, energy-based
 *        state selection for electrical networks
 *   L5 - Engineering Methods: systematic state equation derivation for RLC,
 *        normal tree method, bond graph to state-space
 *
 * Reference:
 *   MIT 6.241 Dynamic Systems and Control — state-space fundamentals
 *   Berkeley ME 232 — State-space for electromechanical systems
 *   Ogata, Modern Control Engineering (2010) §9
 */

#ifndef ELECTRICAL_STATE_SPACE_H
#define ELECTRICAL_STATE_SPACE_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define SS_MAX_ORDER 32
#define SS_MAX_INPUTS 8
#define SS_MAX_OUTPUTS 8

/* L1: Continuous-time state-space model */
/* dx/dt = A*x + B*u, y = C*x + D*u */
typedef struct {
    double A[SS_MAX_ORDER * SS_MAX_ORDER]; /* State matrix (n x n) */
    double B[SS_MAX_ORDER * SS_MAX_INPUTS]; /* Input matrix (n x m) */
    double C[SS_MAX_OUTPUTS * SS_MAX_ORDER]; /* Output matrix (p x n) */
    double D[SS_MAX_OUTPUTS * SS_MAX_INPUTS]; /* Feedthrough matrix (p x m) */
    size_t n; /* Number of states */
    size_t m; /* Number of inputs */
    size_t p; /* Number of outputs */
} StateSpace;

/* L1: Initialize state-space model to zero */
void ss_init(StateSpace *ss, size_t n, size_t m, size_t p);

/* L1: Accessor macros for flat arrays (row-major) */
#define SS_A(ss, i, j) ((ss)->A[(i)*(ss)->n + (j)])
#define SS_B(ss, i, j) ((ss)->B[(i)*(ss)->n + (j)])
#define SS_C(ss, i, j) ((ss)->C[(i)*(ss)->p + (j)])
#define SS_D(ss, i, j) ((ss)->D[(i)*(ss)->p + (j)])

/* L2: Eigenvalues of state matrix A */

/**
 * @brief Compute eigenvalues of A matrix using QR algorithm
 *
 * Theorem: Eigenvalues λ_i of A satisfy det(λI - A) = 0.
 *   For stability: all Re(λ_i) < 0 (continuous-time).
 *   Natural frequencies: ω_i = |Im(λ_i)|
 *   Damping ratios: ζ_i = -Re(λ_i)/|λ_i|
 *
 * For electrical networks with only RLC (no controlled sources),
 * all eigenvalues have Re(λ) <= 0 (passivity ensures stability).
 *
 * @param eigenvalues Output array (caller allocates n * sizeof(double complex))
 * @return Number of eigenvalues found
 * Complexity: O(n^3) for QR iteration
 */
int ss_eigenvalues(const StateSpace *ss, double complex *eigenvalues);

/**
 * @brief Check if system is Hurwitz stable (all eigenvalues strictly in LHP)
 *
 * Theorem (Routh-Hurwitz 1874/1895): A continuous-time LTI system is
 *   asymptotically stable iff all eigenvalues of A have negative real parts.
 *   For physically realizable passive RLC networks, this holds automatically.
 *
 * @return 1 if stable, 0 if marginally stable or unstable
 * Complexity: O(n^3)
 */
int ss_is_stable(const StateSpace *ss);

/* L2: Controllability matrix */

/**
 * @brief Compute controllability Gramian W_c = ∫_0^∞ e^{At}BB'e^{A't} dt
 *
 * Theorem (Kalman 1960): The system is controllable iff the controllability
 *   Gramian is positive definite (for continuous-time LTI systems).
 *   Alternatively, the controllability matrix
 *   C = [B, AB, A^2B, ..., A^{n-1}B] must have full rank n.
 *
 * @param Gramian Output: n x n matrix (caller allocates n*n doubles)
 * @return rank of controllability matrix
 * Complexity: O(n^3)
 */
int ss_controllability_gramian(const StateSpace *ss, double *Gramian);

/**
 * @brief Compute rank of controllability matrix
 *
 * Theorem: rank(C) = n iff system is completely controllable.
 *   Uses QR decomposition with column pivoting for numerical rank.
 *
 * @return numerical rank of controllability matrix
 * Complexity: O(n^3)
 */
int ss_controllability_rank(const StateSpace *ss);

/* L2: Observability matrix */

/**
 * @brief Compute observability Gramian W_o = ∫_0^∞ e^{A't}C'Ce^{At} dt
 *
 * Theorem (Kalman 1960): The system is observable iff the observability
 *   Gramian is positive definite. Alternatively, the observability matrix
 *   O = [C; CA; CA^2; ...; CA^{n-1}] must have full column rank n.
 *
 * @param Gramian Output: n x n matrix (caller allocates n*n doubles)
 * @return rank of observability matrix
 * Complexity: O(n^3)
 */
int ss_observability_gramian(const StateSpace *ss, double *Gramian);

/**
 * @brief Compute rank of observability matrix
 * Complexity: O(n^3)
 */
int ss_observability_rank(const StateSpace *ss);

/* L3: Matrix exponential — state transition */

/**
 * @brief Compute state transition matrix Φ(t) = e^{At} using Padé approximation
 *
 * Theorem: For dx/dt = Ax, the solution is x(t) = e^{At} x(0).
 *   Matrix exponential defined via power series: e^{At} = Σ_{k=0}^∞ (At)^k/k!
 *   Implementation uses [8/8] Padé with scaling-and-squaring
 *   (Higham 2009, "The Scaling and Squaring Method for the Matrix Exponential Revisited").
 *
 * @param Phi Output: n x n matrix (caller allocates n*n doubles)
 * @param t Time [s]
 * Complexity: O(n^3)
 */
int ss_state_transition(const StateSpace *ss, double t, double *Phi);

/**
 * @brief Simulate state-space system forward in time using Runge-Kutta 4
 *
 * Solves dx/dt = Ax + Bu from t0 to tf given initial state x0 and
 * constant input u over the interval.
 *
 * @param x Current state (input: x(t0), output: x(tf))
 * @param u Input vector [m]
 * @param t0 Start time [s]
 * @param tf End time [s]
 * @param steps Number of RK4 steps
 * Complexity: O(n^3 * steps) if using full state transition, O(n^2 * steps) for RK4
 */
int ss_simulate_rk4(StateSpace *ss, double *x, const double *u,
                     double t0, double tf, size_t steps);

/* L5: RLC circuit to state-space (systematic derivation) */

/**
 * @brief Build state-space model for series RLC circuit
 *
 * Circuit: V_in → R → L → C → GND
 * Output: voltage across C
 *
 * State variables: x1 = i_L (inductor current), x2 = v_C (capacitor voltage)
 * State equations:
 *   di_L/dt = -(R/L)·i_L - (1/L)·v_C + (1/L)·V_in
 *   dv_C/dt = (1/C)·i_L
 *   y = v_C
 *
 * Matrices:
 *   A = [-R/L  -1/L;  1/C  0]
 *   B = [1/L;  0]
 *   C = [0  1]
 *   D = [0]
 *
 * Complexity: O(1)
 */
StateSpace ss_from_series_rlc(double R, double L, double C);

/**
 * @brief Build state-space model for parallel RLC circuit
 *
 * Circuit: I_in → R || L || C → GND
 * Output: voltage across all elements
 *
 * State variables: x1 = v_C (voltage), x2 = i_L (inductor current)
 * State equations:
 *   dv_C/dt = -(1/(RC))·v_C - (1/C)·i_L + (1/C)·I_in
 *   di_L/dt = (1/L)·v_C
 *   y = v_C
 *
 * Complexity: O(1)
 */
StateSpace ss_from_parallel_rlc(double R, double L, double C);

/**
 * @brief Build state-space model for two-loop RLC circuit
 *
 * Circuit topology: V_in → R1 → L1 → C1 → GND, with R2→L2 from node between L1 and C1
 *
 * State variables: x1 = i_L1, x2 = i_L2, x3 = v_C1
 * Useful for modeling coupled electrical dynamics.
 *
 * Complexity: O(1)
 */
StateSpace ss_from_two_loop_rlc(double R1, double L1, double R2, double L2, double C1);

/* L5: Transfer function to state-space conversion */

/**
 * @brief Convert SISO transfer function to controllable canonical form
 *
 * Given H(s) = (b0 + b1*s + ... + bm*s^m) / (a0 + a1*s + ... + an*s^n)
 * with an = 1 (monic denominator).
 *
 * Controllable canonical form (phase-variable form):
 *   A = [0  1  0 ... 0; 0 0 1 ... 0; ...; -a0 -a1 -a2 ... -a_{n-1}]
 *   B = [0; 0; ...; 1]
 *   C = [b0 b1 ... bm 0 ... 0]
 *   D = [0] (if strictly proper, m < n)
 *
 * Theorem: Any proper rational transfer function has infinitely many
 *   state-space realizations. The controllable canonical form is
 *   minimal iff numerator and denominator are coprime.
 *
 * Reference: Kailath, Linear Systems (1980) §2.4
 *
 * Complexity: O(n)
 */
StateSpace ss_controllable_canonical(const double *num, size_t m,
                                      const double *den, size_t n);

/**
 * @brief Convert SISO transfer function to observable canonical form
 *
 * Observable canonical form:
 *   A = [0 0 ... 0 -a0; 1 0 ... 0 -a1; 0 1 ... 0 -a2; ...; 0 0 ... 1 -a_{n-1}]
 *   B = [b0; b1; ...; b_{n-1}]
 *   C = [0 0 ... 0 1]
 *   D = [0]
 *
 * This is the dual of the controllable canonical form.
 * Complexity: O(n)
 */
StateSpace ss_observable_canonical(const double *num, size_t m,
                                    const double *den, size_t n);

/* L3: State-space to transfer function */

/**
 * @brief Compute SISO transfer function H(s) = C(sI - A)^{-1}B + D
 *
 * Uses Leverrier-Faddeev algorithm O(n^4) for symbolic polynomial computation,
 * or evaluate at specific s using Cramer rule.
 *
 * @param s Complex frequency point
 * @return H(s) evaluated at s
 * Complexity: O(n^3) for matrix inverse
 */
double complex ss_evaluate_tf(const StateSpace *ss, size_t input_idx,
                               size_t output_idx, double complex s);

/* L6: Natural frequency and damping from eigenvalues */

/**
 * @brief Extract modal parameters from a complex eigenvalue pair
 *
 * For a complex conjugate pair λ = -ζω_n ± jω_n√(1-ζ^2):
 *   ω_n = |λ| (natural frequency [rad/s])
 *   ζ = -Re(λ)/|λ| (damping ratio)
 *
 * For real eigenvalues: ζ = 1 (critically damped/overdamped)
 *
 * Complexity: O(1)
 */
typedef struct {
    double natural_freq;   /* ω_n [rad/s] */
    double damping_ratio;  /* ζ [dimensionless] */
    double time_constant;  /* τ = 1/(ζ·ω_n) [s] */
} ModalParameters;

ModalParameters ss_modal_from_eigenvalue(double complex lambda);

/**
 * @brief Dominant time constant of the system
 *
 * Theorem: The slowest mode (eigenvalue closest to imaginary axis) determines
 *   the settling time: t_s ≈ 4/|min Re(λ_i)| for 2% criterion.
 *
 * Complexity: O(n)
 */
double ss_dominant_time_constant(const StateSpace *ss);

/* L7: Discretization for digital implementation */

/**
 * @brief Exact discretization using zero-order hold (ZOH)
 *
 * Computes discrete-time state-space:
 *   x[k+1] = Φ·x[k] + Γ·u[k]
 *   y[k] = C·x[k] + D·u[k]
 *
 * where Φ = e^{A·Ts}, Γ = ∫_0^{Ts} e^{A·τ} dτ · B
 *
 * Uses Van Loan method (1978) for computing Φ and Γ simultaneously.
 * Reference: Van Loan, IEEE TAC 23(3), 1978.
 *
 * Complexity: O(n^3)
 */
typedef struct {
    double Phi[SS_MAX_ORDER * SS_MAX_ORDER];
    double Gamma[SS_MAX_ORDER * SS_MAX_INPUTS];
    double Cd[SS_MAX_OUTPUTS * SS_MAX_ORDER];
    double Dd[SS_MAX_OUTPUTS * SS_MAX_INPUTS];
    double Ts;
    size_t n, m, p;
} DiscreteStateSpace;

int ss_discretize_zoh(const StateSpace *ss, double Ts, DiscreteStateSpace *dss);

#endif /* ELECTRICAL_STATE_SPACE_H */
