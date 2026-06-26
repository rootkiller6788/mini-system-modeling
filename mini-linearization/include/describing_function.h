/**
 * @file describing_function.h
 * @brief Describing function (DF) analysis -- harmonic linearization of
 *        nonlinear elements for limit cycle prediction.
 *
 * Knowledge Coverage:
 *   L1 - Describing function, sinusoidal-input describing function (SIDF),
 *        dual-input describing function (DIDF)
 *   L2 - Limit cycle prediction, harmonic balance, DF-based Nyquist criterion
 *   L3 - Fourier series representation, complex gain N(A, omega)
 *   L5 - DF for relay, saturation, dead zone, backlash, hysteresis
 *   L6 - Limit cycle existence, amplitude/frequency, stability of limit cycles
 *
 * Reference:
 *   Slotine & Li, Applied Nonlinear Control (1991), Ch.5
 *   Khalil, Nonlinear Systems (3rd ed, 2002), Ch.7
 *   Gelb & Vander Velde, Multiple-Input Describing Functions (1968)
 *   Atherton, Nonlinear Control Engineering (1982)
 */

#ifndef DESCRIBING_FUNCTION_H
#define DESCRIBING_FUNCTION_H

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>

/**
 * L1: Describing function N(A, omega) = complex gain of a nonlinear element
 * under sinusoidal input x = A*sin(omega*t).
 *
 * N(A) = (1/A) * (b_1(A) + j*a_1(A))
 * where a_1, b_1 are first Fourier coefficients of the output.
 *
 * For static nonlinearities (no frequency dependence), N = N(A) only.
 */
typedef struct {
    double complex gain;   /**< N(A, omega) */
    double amplitude;      /**< Input amplitude A */
    double frequency;      /**< Input frequency omega [rad/s] */
    double phase_deg;      /**< Phase shift in degrees */
    double magnitude_db;   /**< |N(A)| in dB */
} DescribingFunction;

/**
 * L6: Limit cycle prediction from describing function analysis.
 *
 * A limit cycle exists at amplitude A if:
 *   G(j*omega) = -1 / N(A)
 *
 * i.e., the Nyquist plot of G(jw) intersects the negative reciprocal
 * of the describing function -1/N(A).
 */
typedef struct {
    bool exists;           /**< true if limit cycle predicted */
    double amplitude;      /**< Predicted amplitude */
    double frequency_rad_s; /**< Predicted frequency [rad/s] */
    double frequency_hz;   /**< Predicted frequency [Hz] */
    bool is_stable;        /**< true if limit cycle is stable */
} LimitCyclePrediction;

/**
 * L5: Compute describing function for an ideal relay.
 * f(x) = M * sign(x),  N(A) = 4*M / (pi*A)
 *
 * Reference: Slotine & Li, section 5.2.1
 */
DescribingFunction df_ideal_relay(double M, double A);

/**
 * L5: Compute describing function for relay with hysteresis.
 * f(x) = M*sign(x) with switching at +/- Delta.
 * N(A) = 4*M/(pi*A) * (sqrt(1-(Delta/A)^2) - j*(Delta/A)) for A >= Delta
 *
 * Reference: Atherton (1982), Ch.4
 */
DescribingFunction df_relay_hysteresis(double M, double Delta, double A);

/**
 * L5: Compute describing function for saturation.
 * f(x) = k*x for |x|<=a, k*a*sign(x) for |x|>a
 * N(A) = k for A<=a, (2k/pi)*(asin(a/A) + (a/A)*sqrt(1-(a/A)^2)) for A>a
 *
 * Reference: Slotine & Li, section 5.2.2
 */
DescribingFunction df_saturation(double k, double a, double A);

/**
 * L5: Compute describing function for dead zone.
 * f(x) = 0 for |x|<=d, k*(x - d*sign(x)) for |x|>d
 * N(A) = 0 for A<=d, k*(1 - (2/pi)*(asin(d/A) + (d/A)*sqrt(1-(d/A)^2))) for A>d
 *
 * Reference: Slotine & Li, section 5.2.3
 */
DescribingFunction df_dead_zone(double k, double d, double A);

/**
 * L5: Compute describing function for backlash.
 * f(x) tracks x with gap b.
 * N(A) = (k/pi) * [pi/2 + asin(1-2b/A) + 2*(1-2b/A)*sqrt((b/A)*(1-b/A))]
 *        + j * (-4k*b*(1-b/A) / (pi*A))
 * for A >= b (zero otherwise).
 *
 * Reference: Gelb & Vander Velde (1968), section 3.5
 */
DescribingFunction df_backlash(double k, double b, double A);

/**
 * L5: Compute describing function for cubic nonlinearity.
 * f(x) = c * x^3,  N(A) = (3/4) * c * A^2
 *
 * Reference: Khalil (2002), Example 7.1
 */
DescribingFunction df_cubic(double c, double A);

/**
 * L5: Compute describing function for quantizer (uniform step q).
 * N(A) approximates quantization as a linear gain.
 *
 * Reference: Slotine & Li, section 5.2.4
 */
DescribingFunction df_quantizer(double q, double A);

/**
 * L6: Predict limit cycle using describing function and linear transfer function.
 *
 * Solve for (A, omega) such that G(j*omega) * N(A) = -1.
 * Uses numerical search with harmonic balance equation.
 *
 * @param G_real      Real part of linear TF at omega
 * @param G_imag      Imaginary part of linear TF at omega
 * @param df          Describing function at amplitude A
 * @param tolerance   Convergence tolerance for harmonic balance
 * @return            Limit cycle prediction
 *
 * Complexity: O(n_search * cost(eval))
 * Reference: Khalil (2002), section 7.2
 */
LimitCyclePrediction predict_limit_cycle_df(
    double (*G_real)(double omega, void *data),
    double (*G_imag)(double omega, void *data),
    DescribingFunction (*df_func)(double A, void *data),
    double A_min, double A_max, double omega_min, double omega_max,
    double tolerance, void *user_data);

/**
 * L2: DF-based stability analysis via extended Nyquist criterion.
 *
 * The closed-loop system with nonlinearity N(A) and linear part G(s)
 * is stable if the Nyquist plot of G(jw) does not encircle the critical
 * locus -1/N(A).
 *
 * @param G_real, G_imag  Linear plant frequency response
 * @param df              Describing function
 * @return                true if DF analysis predicts stability
 *
 * Reference: Slotine & Li, section 5.3
 */
bool df_stability_check(
    double (*G_real)(double omega, void *data),
    double (*G_imag)(double omega, void *data),
    double omega_min, double omega_max, int n_points,
    DescribingFunction (*df_func)(double A, void *data),
    double A, void *user_data);

/** Free a DescribingFunction. O(1). */
void free_describing_function(DescribingFunction *df);

 /* DESCRIBING_FUNCTION_H */

/** L6: Dual-input describing function (DIDF) for bias+sinusoid response. */
typedef struct {
    double complex gain_sinusoid;
    double bias_gain;
    double amplitude;
    double bias_level;
    double frequency;
} DualInputDF;

/** L5: DIDF for asymmetric nonlinearity with bias input. */
DualInputDF df_dual_input_relay(double M,double Delta,double A,double B);

/** L5: Random-input describing function (RIDF) for Gaussian input. */
DescribingFunction df_random_input(double (*nonlin)(double),double sigma);

/** L6: Limit cycle stability via Loeb criterion. */
bool loeb_stability_criterion(DescribingFunction (*df)(double,void*),
    double A,double omega,void *ud);

/** L6: Multi-loop describing function analysis. */
LimitCyclePrediction predict_multi_loop_df(
    double (*G_re)(double,void*),double (*G_im)(double,void*),
    DescribingFunction (*df1)(double,void*),DescribingFunction (*df2)(double,void*),
    double A_min,double A_max,double w_min,double w_max,double tol,void *ud);

/** Free DualInputDF. */
void free_dual_input_df(DualInputDF *df);

#endif /* DESCRIBING_FUNCTION_H */
