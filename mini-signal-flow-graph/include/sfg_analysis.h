/**
 * sfg_analysis.h — Signal Flow Graph Analysis: Sensitivity, Stability,
 *                   Frequency Response, and Monte Carlo Methods
 *
 * Once the SFG transfer function is obtained via Mason's formula or
 * reduction, this module provides engineering analysis tools that
 * leverage the SFG structure:
 *
 *   - Sensitivity analysis: dT/dg for any branch gain g
 *   - Stability analysis: denominator (Δ) root analysis
 *   - Frequency response: s → jω substitution in SFG
 *   - Monte Carlo tolerance analysis: statistical gain variation
 *   - Worst-case analysis: extreme value gain variation
 *
 * References:
 *  - Truxal, J.G. "Automatic Feedback Control System Synthesis"
 *    McGraw-Hill, 1955. Chapter 6 (Sensitivity).
 *  - Bode, H.W. "Network Analysis and Feedback Amplifier Design"
 *    Van Nostrand, 1945. (Bode sensitivity integral).
 *  - Spanos, C.J. "Statistical Process Control in Semiconductor
 *    Manufacturing" Proc. IEEE, 1992. (Monte Carlo methods).
 */

#ifndef SFG_ANALYSIS_H
#define SFG_ANALYSIS_H

#include "sfg_core.h"
#include "sfg_mason.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L5: Sensitivity Analysis
 * ================================================================ */

/**
 * sfg_sensitivity — Compute sensitivity of transfer function T
 *                   with respect to a specific branch gain g
 *
 * L5 Definition: The sensitivity S_g^T = (dT/T) / (dg/g) = (g/T) · (dT/dg)
 *
 * For SFG analysis, Bode's sensitivity formula states:
 *   S_g^T = (forward-path gain through g) · (cofactor without g) / T
 *
 * This uses Mason's formula structure: the sensitivity can be computed
 * by evaluating the transfer function contribution of the specific
 * branch, treating it as a separate forward path.
 *
 * L6 Application: Sensitivity analysis is critical for:
 *   - Identifying which component tolerances matter most
 *   - Robust design: minimize sensitivity to parameter variations
 *   - Fault detection: monitor sensitive parameters for drift
 *
 * Complexity: O(P·2^L) per sensitivity evaluation
 *
 * @param g       The signal flow graph
 * @param source  Source node
 * @param sink    Sink node
 * @param branch_src  Branch source node
 * @param branch_dst  Branch destination node
 * @param sens    Output: sensitivity value S_g^T
 * @return        0 on success, -1 on error
 */
int sfg_sensitivity(const sfg_graph_t* g,
                    sfg_node_id_t source, sfg_node_id_t sink,
                    sfg_node_id_t branch_src, sfg_node_id_t branch_dst,
                    sfg_real_t* sens);

/**
 * sfg_sensitivity_all — Compute sensitivity for all branches
 *
 * Returns an array sens[b] = dT/dg_b normalized by T
 *
 * Complexity: O(B · P · 2^L)
 */
int sfg_sensitivity_all(const sfg_graph_t* g,
                        sfg_node_id_t source, sfg_node_id_t sink,
                        sfg_real_t* sens_array, int max_branches);

/**
 * sfg_sensitivity_sort — Sort branches by sensitivity magnitude
 *
 * Returns indices of branches sorted by |sensitivity| (descending).
 * Most sensitive branch is index 0.
 */
void sfg_sensitivity_sort(const sfg_real_t* sens, int n, int* ranked);

/**
 * sfg_bode_sensitivity_integral — Verify Bode sensitivity integral
 *
 * L4 Theorem: For any stable feedback system, the integral of
 * ln|S(jω)| over all frequencies is zero:
 *   ∫₀^∞ ln|S(jω)| dω = 0
 * (for systems with ≥2 more poles than zeros)
 *
 * This function verifies the integral numerically using the SFG.
 *
 * Returns: computed integral value (should be ~0 for valid systems)
 */
double sfg_bode_sensitivity_integral(const sfg_graph_t* g,
                                      sfg_node_id_t source,
                                      sfg_node_id_t sink);

/* ================================================================
 * L5: Stability Analysis via SFG
 * ================================================================ */

/**
 * sfg_stability_routh — Build Routh array from SFG denominator
 *
 * L4 Method: The denominator Δ(s) from Mason's formula gives the
 * characteristic polynomial. The Routh-Hurwitz criterion determines
 * stability without finding roots explicitly.
 *
 * Note: This applies to the SFG with symbolic s-dependent branches
 * (integrators represented as 1/s). The denominator Δ(s) is the
 * characteristic polynomial.
 *
 * @param coef    Characteristic polynomial coefficients [a_n ... a_0]
 * @param order   Polynomial order
 * @param routh   Output: Routh array (order+1 × (order+2)/2, row-major)
 * @param stable  Output: 1 if stable, 0 if unstable
 * @return        0 on success, -1 on error
 */
int sfg_stability_routh(const sfg_real_t* coef, int order,
                         sfg_real_t* routh, int* stable);

/**
 * sfg_stability_nyquist — Analyze SFG stability via Nyquist criterion
 *
 * L4 Theorem: N = Z - P where N = net CCW encirclements of -1.
 * For the SFG, the loop transfer function L(s) extracted from
 * Mason's formula is used.
 *
 * Complexity: O(K·B) for K frequency points
 */
int sfg_stability_nyquist(const sfg_graph_t* g,
                          sfg_node_id_t source, sfg_node_id_t sink,
                          int* stable, double* gain_margin,
                          double* phase_margin);

/* ================================================================
 * L5: Frequency Response Analysis
 * ================================================================ */

/**
 * sfg_frequency_response — Compute frequency response from SFG
 *
 * L6: For frequency-domain analysis, replace each integrator gain
 * (1/s) with 1/(jω) and evaluate the transfer function at each
 * frequency ω.
 *
 * The SFG naturally supports frequency response analysis because
 * each branch gain can be a complex function of frequency.
 *
 * @param g          SFG with frequency-dependent branches
 * @param source     Source node
 * @param sink       Sink node
 * @param omega      Array of frequency points (rad/s)
 * @param n_freqs    Number of frequency points
 * @param magnitude  Output: |T(jω)| for each frequency
 * @param phase      Output: ∠T(jω) in radians for each frequency
 * @return           0 on success, -1 on error
 */
int sfg_frequency_response(const sfg_graph_t* g,
                           sfg_node_id_t source, sfg_node_id_t sink,
                           const double* omega, int n_freqs,
                           double* magnitude, double* phase);

/**
 * sfg_bode_plot_data — Generate data for Bode plot from SFG
 *
 * Convenience wrapper: generates logarithmically spaced frequencies
 * and computes magnitude (dB) and phase (degrees).
 */
int sfg_bode_plot_data(const sfg_graph_t* g,
                        sfg_node_id_t source, sfg_node_id_t sink,
                        double w_min, double w_max, int n_points,
                        double* omega, double* mag_db, double* phase_deg);

/**
 * sfg_bandwidth — Compute -3dB bandwidth from SFG
 *
 * Returns: frequency (rad/s) where |T(jω)| = |T(0)|/√2
 */
double sfg_bandwidth(const sfg_graph_t* g,
                     sfg_node_id_t source, sfg_node_id_t sink);

/**
 * sfg_resonant_peak — Find resonant peak of frequency response
 *
 * Returns: peak magnitude ratio M_p = max_ω |T(jω)|/|T(0)|
 */
double sfg_resonant_peak(const sfg_graph_t* g,
                          sfg_node_id_t source, sfg_node_id_t sink,
                          double* resonant_freq);

/* ================================================================
 * L8: Monte Carlo Tolerance Analysis
 * ================================================================ */

/**
 * sfg_monte_carlo — Monte Carlo simulation of gain variations
 *
 * L8 Advanced: Randomly perturb branch gains according to specified
 * distributions, compute the resulting transfer function for each
 * sample, and collect statistics.
 *
 * This is an agent-based approach to understanding how component
 * tolerances propagate through a system.
 *
 * @param g           Nominal SFG
 * @param source      Source node
 * @param sink        Sink node
 * @param tolerance   Relative tolerance (e.g., 0.05 = ±5%)
 * @param n_samples   Number of Monte Carlo samples
 * @param mean        Output: mean of |T|
 * @param stddev      Output: standard deviation of |T|
 * @param worst_case  Output: maximum |T| observed
 * @param best_case   Output: minimum |T| observed
 * @return            0 on success, -1 on error
 */
int sfg_monte_carlo(const sfg_graph_t* g,
                    sfg_node_id_t source, sfg_node_id_t sink,
                    double tolerance, int n_samples,
                    double* mean, double* stddev,
                    double* worst_case, double* best_case);

/**
 * sfg_worst_case_analysis — Extreme value (worst-case) analysis
 *
 * L6: Determine the maximum possible deviation of T given bounded
 * gain variations. Uses interval arithmetic on Mason's formula.
 *
 * Each gain g is replaced by interval [g(1-ε), g(1+ε)].
 * The worst-case T is computed by extremizing each contribution
 * in Mason's formula.
 */
int sfg_worst_case_analysis(const sfg_graph_t* g,
                             sfg_node_id_t source, sfg_node_id_t sink,
                             double tolerance,
                             double* T_nominal,
                             double* T_worst_high,
                             double* T_worst_low);

/**
 * sfg_yield_analysis — Estimate production yield given spec limits
 *
 * L7 Application: Given upper and lower specification limits on T,
 * estimate the fraction of systems that meet spec under gain
 * variations (using Monte Carlo).
 *
 * Returns: estimated yield (0.0 to 1.0)
 */
double sfg_yield_analysis(const sfg_graph_t* g,
                           sfg_node_id_t source, sfg_node_id_t sink,
                           double tolerance,
                           double T_spec_low, double T_spec_high,
                           int n_samples);

/* ================================================================
 * L8: Advanced Sensitivity Methods
 * ================================================================ */

/**
 * sfg_sensitivity_matrix — Compute sensitivity matrix for MIMO SFG
 *
 * L8: For MIMO systems, S_ij,kl = dT_ij / dg_kl
 * The sensitivity matrix captures how each element of the transfer
 * function matrix depends on each branch gain.
 */
int sfg_sensitivity_matrix(const sfg_graph_t* g,
                            const sfg_node_id_t* inputs, int M,
                            const sfg_node_id_t* outputs, int N,
                            sfg_real_t* S_matrix);

/**
 * sfg_balanced_sensitivity — Check if SFG is balanced (equal sensitivity
 *                              in parallel paths)
 *
 * L8: A balanced SFG has equal sensitivity magnitudes across parallel
 * forward paths, providing inherent robustness to gain variations.
 * This property is used in analog circuit design.
 *
 * Returns: balance metric (0 = perfectly balanced, higher = more imbalance)
 */
double sfg_balanced_sensitivity(const sfg_graph_t* g,
                                 sfg_node_id_t source, sfg_node_id_t sink);

#ifdef __cplusplus
}
#endif

#endif /* SFG_ANALYSIS_H */
