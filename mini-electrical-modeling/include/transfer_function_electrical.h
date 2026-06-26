/**
 * @file transfer_function_electrical.h
 * @brief Transfer functions of electrical networks s-domain analysis
 *
 * Knowledge Coverage:
 *   L1 - Definitions: transfer function, pole, zero, gain, bandwidth, cutoff frequency,
 *        Bode plot, frequency response, phase margin, gain margin
 *   L2 - Core Concepts: s-domain circuit representation, voltage divider in s-domain,
 *        impedance transfer functions, cascade/parallel interconnections
 *   L3 - Mathematical Structures: rational functions, partial fraction expansion,
 *        Laplace transform pairs for RLC circuits
 *   L4 - Fundamental Laws: s-domain Kirchhoff laws, pole location
 *   L5 - Engineering Methods: filter topology synthesis (Sallen-Key, MFB),
 *        frequency response analysis, asymptotic Bode approximation
 *
 * Reference:
 *   Stanford ENGR105 - Transfer functions from circuit equations
 *   MIT 6.302 Feedback Systems - Loop shaping with electrical networks
 *   Van Valkenburg, Analog Filter Design (1982)
 *   Zverev, Handbook of Filter Synthesis (1967)
 */

#ifndef TRANSFER_FUNCTION_ELECTRICAL_H
#define TRANSFER_FUNCTION_ELECTRICAL_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* Maximum order of transfer function (numerator and denominator) */
#define TF_MAX_ORDER 16

/* L1: Transfer function representation H(s) = N(s)/D(s) */
typedef struct {
    double num[TF_MAX_ORDER + 1]; /* n0 + n1*s + n2*s^2 + ... + nm*s^m */
    double den[TF_MAX_ORDER + 1]; /* d0 + d1*s + d2*s^2 + ... + dn*s^n */
    size_t num_order;              /* m = highest power in numerator */
    size_t den_order;              /* n = highest power in denominator */
} TransferFunction;

/* L1: Calculation at a specific s */

/**
 * @brief Evaluate transfer function H(s) = N(s)/D(s) at complex frequency s
 *
 * Uses Horner method for numerical stability.
 * Complexity: O(n) where n = max(num_order, den_order)
 */
double complex tf_evaluate(const TransferFunction *tf, double complex s);

/* L1: DC gain H(0) = n0/d0 */

/**
 * @brief DC gain = H(0) = n0/d0
 *
 * For physical circuits, DC gain is the steady-state amplification at zero frequency.
 * Complexity: O(1)
 */
double tf_dc_gain(const TransferFunction *tf);

/* L1: Poles and zeros via root-finding */

/**
 * @brief Find all real poles of transfer function (roots of denominator)
 *
 * Uses Durand-Kerner (Weierstrass) method for polynomial root-finding.
 * For denominator D(s) = d0 + d1*s + ... + dn*s^n, finds all n roots.
 *
 * Theorem: Poles determine stability. All poles must have Re(p) < 0
 *   for BIBO stability. Poles of electrical networks from passive RLC
 *   always satisfy Re(p) <= 0 (no active gain).
 *
 * @param poles Output array of complex roots (caller allocates den_order * sizeof(double complex))
 * @return Number of poles found
 * Complexity: O(n^2) per iteration for Durand-Kerner
 */
int tf_find_poles(const TransferFunction *tf, double complex *poles);

/**
 * @brief Find all real zeros of transfer function (roots of numerator)
 * Complexity: O(n^2) per iteration
 */
int tf_find_zeros(const TransferFunction *tf, double complex *zeros);

/* L2: Frequency response — magnitude and phase */

/**
 * @brief Magnitude response |H(jω)| at a frequency
 *
 * Useful for Bode magnitude plots.
 * Complexity: O(n) where n = max(num_order, den_order)
 */
double tf_magnitude_db(const TransferFunction *tf, double omega);

/**
 * @brief Phase response ∠H(jω) [rad] at a frequency
 * Complexity: O(n)
 */
double tf_phase_rad(const TransferFunction *tf, double omega);

/**
 * @brief Phase response ∠H(jω) [degrees] at a frequency
 * Complexity: O(n)
 */
double tf_phase_deg(const TransferFunction *tf, double omega);

/* L3: Bandwidth and cutoff frequency */

/**
 * @brief Find -3dB cutoff frequency via binary search
 *
 * Theorem: -3dB frequency is where |H(jω_c)| = |H(j0)|/√2.
 * For low-pass: ω_c = 1/RC (first-order) or more complex for higher orders.
 *
 * @param f_low Lower bound for search [rad/s]
 * @param f_high Upper bound for search [rad/s]
 * @return Cutoff frequency ω_c or -1 if not found
 * Complexity: O(n * log(range/eps))
 */
double tf_find_cutoff_freq(const TransferFunction *tf, double f_low, double f_high);

/**
 * @brief Estimate -3dB bandwidth from pole locations
 *
 * For a second-order system with complex poles at -ζω_n ± jω_n√(1-ζ^2),
 * bandwidth ≈ ω_n * sqrt(1 - 2ζ^2 + sqrt(2 - 4ζ^2 + 4ζ^4))
 *
 * Complexity: O(1)
 */
double tf_bandwidth_estimate_from_poles(double natural_freq, double damping);

/* L2: Step response from transfer function */

/**
 * @brief Compute step response at time t using inverse Laplace via partial fractions
 *
 * Theorem: Step response y(t) = L^{-1}{H(s)/s} computed by residue method.
 * For first-order H(s) = K/(τs + 1): y(t) = K*(1 - e^{-t/τ})
 * For second-order: superposition of exponentials or damped sinusoids.
 *
 * @param t Time [s]
 * @return Step response value y(t)
 * Complexity: O(n)
 */
double tf_step_response(const TransferFunction *tf, double t);

/* L5: Filter topologies */

/**
 * @brief First-order RC low-pass: H(s) = 1/(RCs + 1) = 1/(τs + 1)
 *
 * Circuit: R in series, C to ground. Output across C.
 * Cutoff: ω_c = 1/τ = 1/(RC)
 *
 * Reference: Van Valkenburg (1982) §3.2
 */
TransferFunction tf_rc_lowpass(double r, double c);

/**
 * @brief First-order RC high-pass: H(s) = RCs/(RCs + 1) = τs/(τs + 1)
 *
 * Circuit: C in series, R to ground. Output across R.
 * Cutoff: ω_c = 1/τ = 1/(RC)
 */
TransferFunction tf_rc_highpass(double r, double c);

/**
 * @brief First-order RL low-pass: H(s) = R/(Ls + R) = 1/((L/R)s + 1)
 *
 * Circuit: L in series, R to ground. Output across R.
 * Cutoff: ω_c = R/L
 */
TransferFunction tf_rl_lowpass(double r, double l);

/**
 * @brief First-order RL high-pass: H(s) = Ls/(Ls + R)
 *
 * Circuit: R in series, L to ground. Output across L.
 * Cutoff: ω_c = R/L
 */
TransferFunction tf_rl_highpass(double r, double l);

/**
 * @brief Second-order RLC low-pass: H(s) = 1/(LC·s^2 + RC·s + 1)
 *
 * Circuit: R and L in series, C to ground. Output across C.
 * Natural frequency: ω_n = 1/√(LC)
 * Damping: ζ = (R/2)·√(C/L)
 * This is the standard form for series RLC with output across capacitor.
 */
TransferFunction tf_rlc_lowpass(double r, double l, double c);

/**
 * @brief Second-order RLC band-pass: H(s) = RCs/(LC·s^2 + RC·s + 1)
 *
 * Circuit: Series RLC, output across R.
 * Center frequency: ω_0 = 1/√(LC)
 * Quality factor: Q = ω_0·L/R = 1/(ω_0·R·C)
 */
TransferFunction tf_rlc_bandpass(double r, double l, double c);

/**
 * @brief Second-order RLC high-pass: H(s) = LC·s^2/(LC·s^2 + RC·s + 1)
 *
 * Circuit: Series RLC, output across L.
 */
TransferFunction tf_rlc_highpass(double r, double l, double c);

/**
 * @brief Sallen-Key second-order low-pass (unity-gain topology)
 *
 * H(s) = 1/(R1·R2·C1·C2·s^2 + (R1+R2)·C2·s + 1)
 *
 * Theorem (Sallen and Key 1955): Active filter topology using a single
 *   op-amp, two resistors, and two capacitors. Provides second-order
 *   response without inductors.
 *
 * Reference: Sallen and Key, IRE Trans. Circuit Theory (1955)
 */
TransferFunction tf_sallen_key_lowpass(double r1, double r2, double c1, double c2);

/**
 * @brief Multiple-feedback (MFB) second-order low-pass
 *
 * H(s) = -(R3/R1)/(R2·R3·C1·C2·s^2 + R3·(C1+C2)·s + 1)
 *
 * Reference: Huelsman and Allen, Active Filter Design (1980)
 */
TransferFunction tf_mfb_lowpass(double r1, double r2, double r3, double c1, double c2);

/* L3: Transfer function arithmetic — series and feedback */

/**
 * @brief Series (cascade) connection: H(s) = H1(s) * H2(s)
 *
 * Theorem: For two systems in series, transfer functions multiply (in s-domain).
 *   Assumes no loading effect (ideal isolation between stages).
 *   Convolution in time domain → multiplication in s-domain.
 *
 * Complexity: O(n1 * n2) where ni = order of Hi
 */
TransferFunction tf_cascade(const TransferFunction *h1, const TransferFunction *h2);

/**
 * @brief Parallel connection: H(s) = H1(s) + H2(s)
 *
 * Theorem: For two systems with outputs summed, transfer functions add.
 * Complexity: O(n1 + n2)
 */
TransferFunction tf_parallel(const TransferFunction *h1, const TransferFunction *h2);

/**
 * @brief Negative feedback: H(s) = H1(s) / (1 + H1(s)*H2(s))
 *
 * Theorem (Black 1927): With forward path H1(s) and feedback path H2(s),
 *   closed-loop transfer function is H1/(1 + H1*H2) for negative feedback.
 *   Characteristic equation: 1 + H1(s)·H2(s) = 0 → determines stability.
 *
 * Complexity: O((n1+n2)^2)
 */
TransferFunction tf_negative_feedback(const TransferFunction *h1, const TransferFunction *h2);

/* L5: Asymptotic Bode approximation */

/**
 * @brief Compute asymptotic Bode magnitude at frequency ω (in dB)
 *
 * Approximates |H(jω)| using piece-wise linear segments:
 *   - DC segment (ω << first corner): constant
 *   - At each pole: -20 dB/dec
 *   - At each zero: +20 dB/dec
 *
 * Complexity: O(num_poles + num_zeros)
 */
double tf_bode_magnitude_asymptotic(const TransferFunction *tf, double omega);

/* L7: Real-world filter design */

/**
 * @brief Design RC filter for a given cutoff frequency
 *
 * Chooses standard capacitor value (E12 series) and computes matching
 * resistor value. Returns the nearest standard resistor value.
 *
 * Reference: IEC 60063 standard values; Horowitz and Hill (2015) §1.7.2
 */
typedef struct {
    double r_standard;     /* Nearest standard resistor value [ohm] */
    double c_standard;     /* Nearest standard capacitor value [F] */
    double actual_cutoff;  /* Actual achieved cutoff [Hz] */
    double error_percent;  /* Error from target [%] */
} FilterDesign;

FilterDesign design_rc_lowpass_to_cutoff(double cutoff_hz);

/**
 * @brief Compute group delay τ_g = -dφ/dω of transfer function
 *
 * Theorem: Group delay measures the time delay of the amplitude envelope
 *   of a signal passing through the filter.
 *   For a linear-phase filter, τ_g is constant over frequency.
 *   Implemented via finite-difference approximation.
 *
 * @param omega Frequency [rad/s]
 * @param delta Small perturbation for numerical derivative
 * @return Group delay [s]
 * Complexity: O(n) per evaluation
 */
double tf_group_delay(const TransferFunction *tf, double omega, double delta);

#endif /* TRANSFER_FUNCTION_ELECTRICAL_H */
