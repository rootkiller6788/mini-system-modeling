/**
 * @file tf_analysis.h
 * @brief Transfer Function Analysis: Stability, Response, Modal Methods
 *
 * Knowledge Coverage:
 *   L4: Routh-Hurwitz stability criterion, Final/Initial Value Theorems,
 *       Nyquist stability criterion, Bode sensitivity integrals
 *   L5: Frequency response (Bode, Nyquist, Nichols), time response
 *       simulation, step response characteristics, bandwidth computation
 *   L6: Modal decomposition via partial fraction expansion,
 *       dominant pole approximation, residue analysis
 *
 * This module provides the complete analytical toolbox for transfer
 * function evaluation, stability assessment, and response characterization.
 *
 * References:
 *   H. Nyquist, "Regeneration Theory", Bell System Tech. J., 1932
 *   E.J. Routh, "A Treatise on the Stability of a Given State of Motion", 1877
 *   A. Hurwitz, "Ueber die Bedingungen...", Math. Ann., 1895
 *   H.W. Bode, "Network Analysis and Feedback Amplifier Design", 1945
 *   MIT 6.302 / Stanford ENGR105 / ETH 151-0591 / Cambridge 3F2
 */
#ifndef TF_ANALYSIS_H
#define TF_ANALYSIS_H

#include "transfer_function.h"
#include "tf_polynomial.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L4: Stability Theorems
 * ================================================================ */

/**
 * Routh-Hurwitz stability criterion.
 *
 * Determines the number of RHP poles without computing roots.
 * Builds the Routh array from denominator coefficients and counts
 * sign changes in the first column.
 *
 * Theorem (Routh-Hurwitz, 1877/1895):
 *   Number of RHP roots = number of sign changes in first column
 *   of the Routh array. System is stable iff all first-column
 *   entries are positive (for monic polynomial with a_n > 0).
 *
 * Complexity: O(den_order^2)
 * Ref: Ogata §5-6, Franklin §3.5
 */
tf_stability_t tf_stability_routh(const TransferFunction *tf);

/**
 * Stability check via explicit pole computation.
 * More computationally expensive than Routh-Hurwitz but gives
 * exact pole locations.
 * Complexity: O(den_order^3)
 */
tf_stability_t tf_stability_poles(const TransferFunction *tf);

/**
 * Build Routh array for a given polynomial and analyze stability.
 *
 * @param den        Denominator coefficients [order+1], ascending powers
 * @param order      Polynomial degree
 * @param array_out  Output: Routh array as flat double* (caller frees)
 *                   Size: (order+1) * ((order+2)/2) elements
 * @return Number of RHP poles (0 = stable)
 * Complexity: O(order^2)
 */
int tf_stability_routh_array(const double *den, int order, double **array_out);

/**
 * Simple Routh-Hurwitz stable check: returns 1 if all poles in LHP, 0 otherwise.
 * Complexity: O(order^2)
 */
int tf_routh_hurwitz_stable(const double *den, int order);

/**
 * Count RHP poles using explicit computation.
 * @return Number of poles with positive real part
 * Complexity: O(den_order^3)
 */
int tf_count_rhp_poles(const TransferFunction *tf);

/**
 * Gain margin: amount of gain increase (in dB) before closed-loop instability.
 * Computed from open-loop frequency response where phase = -180 deg.
 *
 * GM = -20*log10(|L(j*w_pc)|) where w_pc is the phase crossover frequency.
 * Positive GM (in dB) indicates stability margin.
 *
 * Complexity: O(freq_points) for frequency sweep
 * Ref: Ogata §8-8, Franklin §6.4
 */
double tf_gain_margin(const TransferFunction *tf);

/**
 * Phase margin: additional phase lag (in degrees) at gain crossover frequency.
 * PM = 180 + arg(L(j*w_gc)) where w_gc is where |L(j*w_gc)| = 1.
 *
 * Typical design: PM > 45 deg for good damping, PM > 60 deg preferred.
 * Complexity: O(freq_points)
 * Ref: Ogata §8-8, Franklin §6.4
 */
double tf_phase_margin(const TransferFunction *tf);

/* ================================================================
 * L4: Laplace Transform Theorems
 * ================================================================ */

/**
 * Final Value Theorem (FVT): lim_{t->inf} y(t) = lim_{s->0} s*Y(s).
 *
 * For step response: y(inf) = lim_{s->0} G(s) = DC gain.
 * For ramp response: tracks if type >= 1.
 *
 * Precondition: All poles of s*Y(s) must be in LHP (stability).
 * If any pole is in RHP or on jw-axis (except single at origin),
 * the theorem does not apply.
 *
 * @return Final value, or NAN if theorem conditions not met
 * Complexity: O(den_order) for stability check
 * Ref: Ogata §2-5
 */
double tf_final_value_step(const TransferFunction *tf);

/**
 * Initial Value Theorem (IVT): lim_{t->0+} y(t) = lim_{s->inf} s*Y(s).
 *
 * For step response: y(0+) = lim_{s->inf} G(s) = high-frequency gain.
 * For strictly proper TF, y(0+) = 0.
 * For proper TF, y(0+) = num[num_order]/den[den_order].
 *
 * Complexity: O(1)
 * Ref: Ogata §2-5, Franklin §3.1
 */
double tf_initial_value(const TransferFunction *tf);

/**
 * Full Final Value Theorem: computes lim_{t->inf} step response
 * with rigorous condition checking.
 * @param result  Output: final value
 * @return 1 if FVT applies, 0 if conditions violated
 */
int tf_final_value_theorem(const TransferFunction *tf, double *result);

/* ================================================================
 * L3/L5: Frequency Domain Analysis
 * ================================================================ */

/**
 * Evaluate G(s) at complex frequency s = sr + j*si.
 * @param ore  Output: real part G(s)
 * @param oim  Output: imaginary part G(s)
 * Complexity: O(num_order + den_order)
 */
void tf_eval_at(const TransferFunction *tf, double sr, double si,
                 double *ore, double *oim);

/**
 * Bode magnitude at frequency w: |G(jw)|.
 * Complexity: O(num_order + den_order)
 */
double tf_magnitude_at(const TransferFunction *tf, double w);

/**
 * Bode phase at frequency w: arg(G(jw)) in radians [-pi, pi].
 * Complexity: O(num_order + den_order)
 */
double tf_phase_at(const TransferFunction *tf, double w);

/**
 * Compute frequency response over [ws, we] with n log-spaced points.
 *
 * @param ws   Start frequency (rad/s, >0)
 * @param we   End frequency (rad/s, > ws)
 * @param n    Number of points (>=2)
 * @param out  Output array [n] of FreqPoint
 * @return n on success, -1 on error
 * Complexity: O(n * (num_order + den_order))
 */
int tf_frequency_response(const TransferFunction *tf, double ws, double we,
                           int n, FreqPoint *out);

/**
 * Compute Nyquist contour evaluation (s = j*w for w in [ws, we]).
 * Returns complex values G(jw) for polar plot stability analysis.
 *
 * Nyquist stability criterion: Z = N + P where
 *   Z = number of closed-loop RHP poles
 *   N = number of CCW encirclements of (-1, 0)
 *   P = number of open-loop RHP poles
 *
 * Complexity: O(n * degree)
 * Ref: Nyquist (1932), Ogata §8-6
 */
int tf_nyquist_eval(const TransferFunction *tf, double ws, double we,
                     int n, FreqPoint *out);

/**
 * Bandwidth: frequency where |G(jw_bw)| = |G(0)|/sqrt(2) (-3 dB point).
 * For open-loop TF, this is the gain crossover frequency.
 * @return Bandwidth in rad/s, or -1 if not found
 * Complexity: O(100 * degree) for iterative search
 */
double tf_bandwidth(const TransferFunction *tf);

/**
 * Resonant peak: maximum magnitude |G(jw)|.
 * @return Maximum magnitude ratio > 1
 * Complexity: O(200 * degree) for search
 */
double tf_resonant_peak(const TransferFunction *tf);

/**
 * Resonant frequency: frequency w_r where |G(jw_r)| is maximum.
 * For second-order: w_r = wn*sqrt(1-2*zeta^2) for zeta < 1/sqrt(2).
 * @return Resonant frequency in rad/s, or -1 if none
 */
double tf_resonant_frequency(const TransferFunction *tf);

/* ================================================================
 * L5: Time Domain Response Simulation
 * ================================================================ */

/**
 * Simulate step response using numerical integration (RK4).
 *
 * The state-space realization is used for simulation:
 *   dx/dt = A*x + B*u, y = C*x + D*u, u(t) = 1 for t >= 0.
 *
 * @param t_end  Simulation duration (seconds)
 * @param n      Number of time points (>=10)
 * @param out    Output array [n] of TimePoint
 * @return n on success, -1 on error
 * Complexity: O(n * ss_order^2) for RK4 integration
 */
int tf_step_response(const TransferFunction *tf, double t_end, int n,
                      TimePoint *out);

/**
 * Simulate impulse response.
 * u(t) approximated as large initial condition in integrator states.
 * Complexity: O(n * ss_order^2)
 */
int tf_impulse_response(const TransferFunction *tf, double t_end, int n,
                         TimePoint *out);

/**
 * Rise time: time for step response to go from 10% to 90% of final value.
 * For second-order underdamped: Tr = (pi - acos(zeta))/(wn*sqrt(1-zeta^2)).
 * @return Rise time in seconds, or -1 if cannot be determined
 * Complexity: O(100) for iterative computation
 */
double tf_rise_time(const TransferFunction *tf);

/**
 * Settling time: time for step response to stay within +/-tol of final value.
 * For second-order: Ts = 4/(zeta*wn) for 2% criterion.
 * @param tol  Tolerance band (e.g., 0.02 for 2%, 0.05 for 5%)
 * @return Settling time in seconds, or -1
 * Complexity: O(100 * ss_order^2) for simulation
 */
double tf_settling_time(const TransferFunction *tf, double tol);

/**
 * Percent overshoot: PO = 100*(y_max - y_final)/y_final.
 * For second-order: PO = 100*exp(-pi*zeta/sqrt(1-zeta^2)).
 * @return Percent overshoot (0-100), or -1
 * Complexity: O(100 * ss_order^2)
 */
double tf_overshoot(const TransferFunction *tf);

/**
 * Peak time: time at which maximum overshoot occurs.
 * For second-order: Tp = pi/(wn*sqrt(1-zeta^2)).
 * @return Peak time in seconds, or -1
 */
double tf_peak_time(const TransferFunction *tf);

/**
 * Steady-state error for unity feedback step input.
 * System type 0: ess = 1/(1+Kp), type 1+: ess = 0.
 * System type 1: ess_ramp = 1/Kv, type 2+: ess_ramp = 0.
 * @return Steady-state error (0 to 1), or -1
 */
double tf_steady_state_error(const TransferFunction *tf);

/* ================================================================
 * L6: Modal Decomposition (Partial Fraction Expansion)
 * ================================================================ */

/**
 * Compute partial fraction expansion of a transfer function.
 *
 * G(s) = sum_{i} r_i/(s - p_i)^{m_i} + D (direct feedthrough)
 *
 * For distinct real poles: r_i = N(p_i)/D'(p_i) (residue formula).
 * For complex conjugate pairs: combined into 2nd-order terms.
 * For repeated poles: generalized residue formula with derivatives.
 *
 * @return PartialFractionExpansion, or NULL on error
 * Complexity: O(den_order^3) due to pole computation
 *
 * Ref: Ogata §2-4, Churchill & Brown "Complex Variables" Ch. 6
 */
PartialFractionExpansion* tf_partial_fraction(const TransferFunction *tf);

/** Free partial fraction expansion memory. */
void pfe_destroy(PartialFractionExpansion *pfe);

/** Print partial fraction expansion in human-readable form. */
void pfe_print(const PartialFractionExpansion *pfe);

/**
 * Find the dominant pole (pole with largest time constant = 1/|Re(p)|).
 * The dominant pole primarily determines transient response speed.
 *
 * @param pole  Output: dominant pole location (real part stored)
 * @return 1 if dominant pole found, 0 if none
 * Complexity: O(den_order^3)
 */
int tf_dominant_pole(const TransferFunction *tf, double *pole);

/* ================================================================
 * L4: Bode Sensitivity Integrals (Advanced Stability)
 * ================================================================ */

/**
 * Bode sensitivity integral:
 *   integral_0^inf log|S(jw)| dw = pi * sum(Re(p_i)) for stable S(s)
 *
 * This is a fundamental limitation theorem: reducing sensitivity at some
 * frequencies necessarily increases it at others (waterbed effect).
 *
 * For continuous-time with relative degree >= 2:
 *   integral_0^inf log|S(jw)| dw = pi * sum(Re(p_i of L(s)))
 *
 * @return Value of the sensitivity integral
 * Complexity: O(den_order) for pole extraction
 *
 * Ref: Bode (1945), Freudenberg & Looze (1985)
 * L8: Advanced stability limitations
 */
double tf_bode_sensitivity_integral(const TransferFunction *tf);

/**
 * Compute the full Nyquist stability assessment.
 *
 * @param L         Open-loop transfer function G(s)*H(s)
 * @param num_enc   Output: number of CCW encirclements of (-1,0j)
 * @return TF_STABLE if closed-loop stable, TF_UNSTABLE otherwise
 *
 * Complexity: O(1000 * degree) for Nyquist contour evaluation
 */
tf_stability_t tf_nyquist_stability(const TransferFunction *L, int *num_enc);

#ifdef __cplusplus
}
#endif
#endif /* TF_ANALYSIS_H */
