/**
 * @file tf_conversion.h
 * @brief Transfer Function Conversion and Transformation Methods
 *
 * Knowledge Coverage:
 *   L1: State-space representation (A,B,C,D), ZPK form
 *   L2: TF <-> State-Space conversion (controllable/observable canonical)
 *   L3: TF <-> ZPK conversion, pole/zero extraction
 *   L4: Similarity transformations, canonical forms
 *   L5: Continuous <-> Discrete conversion (Tustin, ZOH, FOH, matched)
 *   L6: Canonical form realizations
 *   L8: Model order reduction via balanced truncation
 *
 * Conversions between different system representations are essential
 * for control system analysis and design. Each representation exposes
 * different properties.
 *
 * References:
 *   T. Kailath, Linear Systems, 1980
 *   C.T. Chen, Linear System Theory and Design, 4th ed.
 *   K. Ogata, Modern Control Engineering, Ch. 9
 *   MIT 6.241 (Dynamic Systems) / Stanford ENGR207B / Berkeley ME232
 */
#ifndef TF_CONVERSION_H
#define TF_CONVERSION_H

#include "transfer_function.h"
#include <stddef.h>

#ifdef __cplusplus
extern C {
#endif

/* ---- L1: State-Space Data Structure ---- */

/** State-space representation: dx/dt = A*x + B*u, y = C*x + D*u */
typedef struct {
    double *A;     /**< System matrix [n*n], row-major */
    double *B;     /**< Input matrix [n*m], row-major */
    double *C;     /**< Output matrix [p*n], row-major */
    double *D;     /**< Feedthrough matrix [p*m], row-major */
    int     n;     /**< Number of states */
    int     m;     /**< Number of inputs */
    int     p;     /**< Number of outputs */
    int     is_discrete;
    double  sample_time;
} StateSpace;

/* ---- L1: Lifecycle ---- */

/** Create a state-space model. All matrices are copied internally. */
StateSpace* ss_create(const double *A, const double *B, const double *C,
                       const double *D, int n, int m, int p);
void ss_destroy(StateSpace *ss);
StateSpace* ss_clone(const StateSpace *ss);
void ss_print(const StateSpace *ss);

/* ---- L2/L3: TF <-> SS Conversion ---- */

/**
 * Convert SISO transfer function to state-space in controllable canonical form.
 *
 * For G(s) = N(s)/D(s) = (b0 + b1*s + ... + bm*s^m)/(s^n + a_{n-1}*s^{n-1} + ... + a0)
 * with m <= n (proper), the controllable canonical form (CCF) is:
 *
 *   A = [[0, 1, 0, ..., 0],
 *        [0, 0, 1, ..., 0],
 *        ...
 *        [-a0, -a1, -a2, ..., -a_{n-1}]]
 *   B = [0, 0, ..., 0, 1]^T
 *   C = [b0-a0*bn, b1-a1*bn, ..., b_{n-1}-a_{n-1}*bn]
 *   D = bn  (where bn = b_m if m=n, else 0)
 *
 * Complexity: O(n^2)
 * Ref: Kailath ¡ì2.4, Chen ¡ì4.3
 */
StateSpace* tf_to_ss_controllable(const TransferFunction *tf);

/**
 * Convert to observable canonical form.
 * A = CCF_A^T, B = CCF_C^T, C = CCF_B^T, D = CCF_D.
 * Complexity: O(n^2)
 */
StateSpace* tf_to_ss_observable(const TransferFunction *tf);

/**
 * Convert to modal (diagonal/Jordan) canonical form.
 * Requires pole computation. Only works when all poles are distinct.
 * Complexity: O(n^3) due to eigenvalue decomposition
 */
StateSpace* tf_to_ss_modal(const TransferFunction *tf);

/**
 * Convert state-space to transfer function: G(s) = C*(sI-A)^{-1}*B + D.
 * Uses Faddeev-Leverrier algorithm for SISO systems.
 * Complexity: O(n^3)
 *
 * Theorem: G(s) = det(sI-A+B*C) / det(sI-A) for SISO (Faddeev formula)
 * Ref: Chen ¡ì4.3, Kailath ¡ì2.4
 */
TransferFunction* ss_to_tf(const StateSpace *ss);

/**
 * Compute the resolvent (sI-A)^{-1} as a matrix of rational functions.
 * Each entry (i,j) is returned as a TransferFunction.
 * @param ss    State-space model
 * @param tf_ij Output: flat array of n*n TransferFunction*
 * @return 0 on success, -1 on error
 *
 * For small n, uses Cramer's rule; for larger, Leverrier's algorithm.
 * Complexity: O(n^4)
 */
int ss_resolvent(const StateSpace *ss, TransferFunction ***tf_ij);

/* ---- L3: TF <-> ZPK Conversion ---- */

/** Convert transfer function to Zero-Pole-Gain form.
 *  Requires root-finding for num and den polynomials.
 *  Complexity: O((no+do)^3)
 */
ZeroPoleGain* tf_to_zpk(const TransferFunction *tf);

/** Convert ZPK to transfer function by polynomial expansion.
 *  Complexity: O((nz+np)^2)
 */
TransferFunction* zpk_to_tf(const ZeroPoleGain *zpk);

/** Free a ZPK structure. */
void zpk_destroy(ZeroPoleGain *zpk);

/** Print ZPK in human-readable format. */
void zpk_print(const ZeroPoleGain *zpk);

/** Extract poles from a transfer function.
 *  @param poles  Output array [tf->den_order]
 *  @return Number of real poles found (complex pairs counted separately)
 */
int tf_extract_poles(const TransferFunction *tf, double *poles);

/** Extract zeros from a transfer function.
 *  @param zeros  Output array [tf->num_order]
 *  @return Number of real zeros found
 */
int tf_extract_zeros(const TransferFunction *tf, double *zeros);

/* ---- L5: Continuous <-> Discrete Conversion ---- */

/**
 * Bilinear (Tustin) transform: s -> (2/Ts)*(z-1)/(z+1).
 * Preserves stability: LHP -> unit circle interior.
 * May cause frequency warping.
 * Complexity: O((no+do)^2) for full substitution
 */
TransferFunction* tf_s_to_z_tustin(const TransferFunction *tf_s, double Ts);

/**
 * Inverse bilinear transform: z -> (1 + s*Ts/2)/(1 - s*Ts/2).
 * Complexity: O((no+do)^2)
 */
TransferFunction* tf_z_to_s_tustin(const TransferFunction *tf_z, double Ts);

/**
 * Zero-Order Hold (ZOH) discretization:
 * G(z) = (1 - z^{-1}) * Z{L^{-1}[G(s)/s]_{t=k*Ts}}
 *
 * For simple TFs (up to 2nd order), uses analytical formulas.
 * For higher order, uses numerical approximation.
 * Complexity: O(2^n) for exact, O(1) for approximations
 */
TransferFunction* tf_s_to_z_zoh(const TransferFunction *tf_s, double Ts);

/**
 * First-Order Hold (FOH) discretization.
 * Complexity: similar to ZOH
 */
TransferFunction* tf_s_to_z_foh(const TransferFunction *tf_s, double Ts);

/**
 * Matched pole-zero discretization: map each pole/zero via z=e^{s*Ts}.
 * DC gain is matched.
 * Complexity: O((no+do)^3) due to root-finding
 */
TransferFunction* tf_s_to_z_matched(const TransferFunction *tf_s, double Ts);

/**
 * Impulse-invariant discretization.
 * Complexity: O((no+do)^3)
 */
TransferFunction* tf_s_to_z_impulse(const TransferFunction *tf_s, double Ts);

/**
 * Forward Euler discretization: s -> (z-1)/Ts.
 * Simple but may not preserve stability.
 */
TransferFunction* tf_s_to_z_forward_euler(const TransferFunction *tf_s, double Ts);

/**
 * Backward Euler discretization: s -> (z-1)/(z*Ts).
 * Preserves stability but introduces damping.
 */
TransferFunction* tf_s_to_z_backward_euler(const TransferFunction *tf_s, double Ts);

/* ---- L6: Canonical Realizations ---- */

/**
 * Check if a state-space realization is minimal (controllable + observable).
 * Uses Kalman rank tests.
 * @return 1 if minimal, 0 otherwise
 */
int ss_is_minimal(const StateSpace *ss);

/**
 * Compute the controllability matrix: Co = [B, AB, A^2B, ..., A^{n-1}B].
 * @return Rank (for SISO, returns 1 if controllable, 0 if not)
 */
int ss_controllability_rank(const StateSpace *ss);

/**
 * Compute the observability matrix: Ob = [C; CA; CA^2; ...; CA^{n-1}].
 * @return Rank (for SISO, returns 1 if observable, 0 if not)
 */
int ss_observability_rank(const StateSpace *ss);

#ifdef __cplusplus
}
#endif
#endif
