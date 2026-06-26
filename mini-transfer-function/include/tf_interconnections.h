/**
 * @file tf_interconnections.h
 * @brief Transfer Function Interconnection Algebra
 *
 * Knowledge Coverage:
 *   L1: Block interconnections (series, parallel, feedback)
 *   L2: Block diagram algebra rules (summing junction, take-off point)
 *   L3: Signal flow graph representation and Mason's rule
 *   L4: Feedback equivalence theorems
 *   L5: Systematic block diagram reduction
 *
 * The interconnection of transfer functions forms the foundation of
 * block diagram modeling in control systems.
 *
 * References:
 *   Franklin et al., Feedback Control of Dynamic Systems, Ch. 3
 *   Ogata, Modern Control Engineering, Ch. 3
 *   Mason, S.J. (1956), Feedback Theory - Further Properties of SFGs
 *   MIT 6.302 / Stanford ENGR105 / ETH 151-0591
 */
#ifndef TF_INTERCONNECTIONS_H
#define TF_INTERCONNECTIONS_H

#include "transfer_function.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- L1/L2: Basic Interconnections ---- */

/**
 * Series connection: G_series(s) = G1(s) * G2(s).
 * In block diagrams, the output of G1 feeds the input of G2.
 * Theorem: For LTI systems, series connection = multiplication of TFs.
 * Complexity: O(deg1*deg2)
 */
TransferFunction* tf_series(const TransferFunction *g1, const TransferFunction *g2);

/**
 * Parallel connection: G_parallel(s) = G1(s) + G2(s).
 * Both blocks receive the same input, outputs are summed.
 * Complexity: O(max(deg1, deg2))
 */
TransferFunction* tf_parallel(const TransferFunction *g1, const TransferFunction *g2);

/**
 * Negative feedback: G_cl(s) = G(s) / (1 + G(s)*H(s)).
 * This is the fundamental single-loop feedback interconnection.
 * For unity feedback (H=1), simplified to G_cl = G/(1+G).
 *
 * Theorem (Feedback Equivalence):
 *   Y(s)/R(s) = G(s)/(1+G(s)H(s)) for negative feedback
 *   Y(s)/R(s) = G(s)/(1-G(s)H(s)) for positive feedback
 *
 * Complexity: O(deg(G)*deg(H))
 * Ref: Franklin ¡ì4.2, Ogata ¡ì3-3
 */
TransferFunction* tf_feedback(const TransferFunction *G, const TransferFunction *H);

/**
 * Positive feedback: G_cl(s) = G(s) / (1 - G(s)*H(s)).
 * Complexity: O(deg(G)*deg(H))
 */
TransferFunction* tf_positive_feedback(const TransferFunction *G, const TransferFunction *H);

/**
 * Unity feedback (H=1): G_cl(s) = G(s) / (1 + G(s)).
 * This is the most common feedback configuration in practice.
 * Complexity: O(deg(G))
 */
TransferFunction* tf_unity_feedback(const TransferFunction *G);

/* ---- L2: Algebraic Operations ---- */

/** Add two transfer functions: G_add = G1 + G2 */
TransferFunction* tf_add(const TransferFunction *a, const TransferFunction *b);

/** Subtract: G_sub = G1 - G2 */
TransferFunction* tf_subtract(const TransferFunction *a, const TransferFunction *b);

/** Multiply: G_mult = G1 * G2 */
TransferFunction* tf_multiply(const TransferFunction *a, const TransferFunction *b);

/** Divide: G_div = G1 / G2 = N1*D2 / D1*N2 */
TransferFunction* tf_divide(const TransferFunction *a, const TransferFunction *b);

/** Scale: k * G(s) */
TransferFunction* tf_scale(const TransferFunction *tf, double k);

/* ---- L3: Block Diagram Reduction ---- */

/**
 * Reduce a feedback loop with disturbance input.
 * @param G    Plant transfer function
 * @param H    Sensor transfer function
 * @param Gd   Disturbance transfer function
 * @return Closed-loop TF from disturbance to output
 *
 * Formula: Y(s)/D(s) = G_d(s) / (1 + G(s)*H(s))
 */
TransferFunction* tf_disturbance_to_output(const TransferFunction *G,
                                            const TransferFunction *H,
                                            const TransferFunction *Gd);

/**
 * Compute the sensitivity function: S(s) = 1/(1 + G(s)*H(s)).
 * S(s) = complementary sensitivity of the feedback loop.
 * From reference to error: E(s)/R(s) = S(s).
 * Complexity: O(deg(G)*deg(H))
 */
TransferFunction* tf_sensitivity(const TransferFunction *G, const TransferFunction *H);

/**
 * Compute the complementary sensitivity: T(s) = G*H/(1 + G*H).
 * T(s) = closed-loop TF from reference to output.
 * Note: S(s) + T(s) = 1 (fundamental identity).
 * Complexity: O(deg(G)*deg(H))
 */
TransferFunction* tf_complementary_sensitivity(const TransferFunction *G,
                                                 const TransferFunction *H);

/**
 * Compute the loop transfer function: L(s) = G(s)*H(s).
 * L(s) is the open-loop transfer function used for Nyquist/Bode analysis.
 */
TransferFunction* tf_loop_transfer(const TransferFunction *G, const TransferFunction *H);

/**
 * Reduce two nested feedback loops.
 * Inner loop: G2 with H2 feedback, outer: G1 in series, with H1 feedback.
 * @return Equivalent single TF
 */
TransferFunction* tf_nested_feedback(const TransferFunction *G1,
                                      const TransferFunction *G2,
                                      const TransferFunction *H1,
                                      const TransferFunction *H2);

/**
 * Feedforward interconnection: output = (G_ff + G_fb*G)*input.
 * Used in two-degree-of-freedom (2DOF) control structures.
 * @param G_fb  Feedback path TF
 * @param G_ff  Feedforward path TF
 * @param G     Plant TF
 * @return Equivalent TF from reference to output
 */
TransferFunction* tf_feedforward(const TransferFunction *G_fb,
                                  const TransferFunction *G_ff,
                                  const TransferFunction *G);

/* ---- L4: Interconnection Theorems ---- */

/**
 * Check if pole-zero cancellation occurs in a feedback interconnection.
 * @return 1 if cancellation detected, 0 otherwise
 *
 * Cancellation can cause internal instability even if the closed-loop
 * TF appears stable (Franklin ¡ì7.4.3).
 */
int tf_has_pole_zero_cancellation(const TransferFunction *G, const TransferFunction *H);

/**
 * Internal stability check for feedback interconnection.
 * Uses the four transfer functions: G/(1+GH), H/(1+GH), 1/(1+GH), GH/(1+GH).
 * All four must have all poles in LHP for internal stability.
 * @return TF_STABLE if internally stable, TF_UNSTABLE otherwise
 */
tf_stability_t tf_internal_stability(const TransferFunction *G, const TransferFunction *H);

/* ---- L5: Multiple Interconnections ---- */

/**
 * Cascade N transfer functions in series.
 * @param tfs   Array of N transfer functions
 * @param n     Number of TFs
 * @return Product G1*G2*...*Gn
 * Complexity: O(n * max_degree^2)
 */
TransferFunction* tf_cascade_n(const TransferFunction **tfs, int n);

/**
 * Sum N transfer functions in parallel.
 * @param tfs   Array of N transfer functions
 * @param n     Number of TFs
 * @return Sum G1+G2+...+Gn
 */
TransferFunction* tf_sum_n(const TransferFunction **tfs, int n);

/**
 * Kharitonov's rectangle: build the four Kharitonov polynomials for
 * robust stability analysis of interval polynomials.
 * (L8: Advanced stability methods)
 */
int tf_kharitonov_polynomials(const double *num_low, const double *num_high,
                               const double *den_low, const double *den_high,
                               int order, TransferFunction **k_out);

#ifdef __cplusplus
}
#endif
#endif
