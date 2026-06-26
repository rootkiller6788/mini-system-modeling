/**
 * @file transfer_function.h
 * @brief Transfer Function - s-domain rational function representation of LTI systems
 *
 * Knowledge Coverage:
 *   L1: TransferFunction struct, proper/improper, type number, relative degree
 *   L2: Construction, destruction, cloning lifecycle management
 *   L3: Core representation using polynomial coefficients in ascending power order
 *   L4: Routh-Hurwitz stability via denominator polynomial analysis
 *   L5: System type classification, DC gain, pole/zero multiplicity counting
 *
 * Coefficient ordering (ascending powers):
 *   N(s) = num[0] + num[1]*s + num[2]*s^2 + ... + num[num_order]*s^num_order
 *   D(s) = den[0] + den[1]*s + den[2]*s^2 + ... + den[den_order]*s^den_order
 *
 * References:
 *   G.F. Franklin et al., Feedback Control of Dynamic Systems, 8th ed., Ch. 3
 *   K. Ogata, Modern Control Engineering, 5th ed., Ch. 2-3
 *   R.C. Dorf, R.H. Bishop, Modern Control Systems, 13th ed., Ch. 2
 *   MIT 6.302 / Stanford ENGR105 / Berkeley ME132 / ETH 151-0591
 */
#ifndef TRANSFER_FUNCTION_H
#define TRANSFER_FUNCTION_H

#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1: Core Data Types - Transfer Function Representations
 * ================================================================ */

/** Transfer function G(s)=N(s)/D(s) with optional time delay.
 *  Coefficients in ascending power order:
 *    num[0] + num[1]*s + ... + num[num_order]*s^num_order
 *    den[0] + den[1]*s + ... + den[den_order]*s^den_order
 */
typedef struct {
    double *num;         /**< Numerator coefficients [num_order+1] */
    double *den;         /**< Denominator coefficients [den_order+1] */
    int    num_order;    /**< Highest power of s in numerator */
    int    den_order;    /**< Highest power of s in denominator */
    double delay;        /**< Transport delay T_d: G(s)*exp(-T_d*s) */
    int    is_discrete;  /**< 1 if discrete-time (z-domain) */
    double sample_time;  /**< Sampling period Ts for discrete systems */
} TransferFunction;

/** System stability classification. */
typedef enum {
    TF_STABLE = 0,
    TF_MARGINALLY_STABLE = 1,
    TF_UNSTABLE = 2
} tf_stability_t;

/** System type: number of free integrators at s=0. */
typedef enum {
    TF_TYPE_0 = 0, TF_TYPE_1 = 1, TF_TYPE_2 = 2, TF_TYPE_3 = 3
} tf_system_type_t;

/** Properness classification. */
typedef enum {
    TF_STRICTLY_PROPER = 0, TF_PROPER = 1, TF_IMPROPER = 2
} tf_properness_t;

/** Phase classification. */
typedef enum {
    TF_MINIMUM_PHASE = 0, TF_NONMINIMUM_PHASE = 1
} tf_phase_type_t;

/** A single term in partial fraction expansion: r/(s - p)^m. */
typedef struct {
    double residue;
    double pole;
    int    multiplicity;
    int    is_complex;
} PFE_Term;

/** Partial fraction expansion result. */
typedef struct {
    PFE_Term *terms;
    int      num_terms;
    double   direct;
    double   delay;
} PartialFractionExpansion;

/** Zero-Pole-Gain representation. */
typedef struct {
    int     n_zeros, n_poles;
    double *zeros, *poles;
    double  gain;
    double  delay;
} ZeroPoleGain;

/** Single frequency-domain point (Bode plot datum). */
typedef struct {
    double frequency;
    double magnitude;
    double phase;
    double real;
    double imag;
} FreqPoint;

/** Time response point. */
typedef struct {
    double t;
    double y;
} TimePoint;

/* ================================================================
 * L1: Lifecycle
 * ================================================================ */

TransferFunction* tf_create(const double *num, int num_order,
                             const double *den, int den_order);
TransferFunction* tf_create_from_gain(double k);
TransferFunction* tf_create_integrator(void);
TransferFunction* tf_create_first_order_lag(double K, double tau);
TransferFunction* tf_create_second_order(double K, double wn, double zeta);
TransferFunction* tf_create_lead(double K, double T, double alpha);
TransferFunction* tf_create_lag(double K, double T, double beta);
TransferFunction* tf_create_pid(double Kp, double Ki, double Kd);
TransferFunction* tf_create_delay(double Td);
TransferFunction* tf_clone(const TransferFunction *src);
void tf_destroy(TransferFunction *tf);

/* ================================================================
 * L2: Core Properties
 * ================================================================ */

tf_properness_t tf_properness(const TransferFunction *tf);
int tf_relative_degree(const TransferFunction *tf);
tf_system_type_t tf_system_type(const TransferFunction *tf);
tf_phase_type_t tf_phase_type(const TransferFunction *tf);
int tf_has_rhp_zeros(const TransferFunction *tf);
double tf_dc_gain(const TransferFunction *tf);
double tf_hf_gain(const TransferFunction *tf);
int tf_num_poles(const TransferFunction *tf);
int tf_num_zeros(const TransferFunction *tf);

/* ================================================================
 * L3: I/O and Utility
 * ================================================================ */

void tf_print(const TransferFunction *tf);
void tf_print_matlab(const TransferFunction *tf);
void tf_print_zpk(const TransferFunction *tf);
int tf_equals(const TransferFunction *a, const TransferFunction *b, double tol);
void tf_normalize(TransferFunction *tf);
void tf_normalize_monic(TransferFunction *tf);
void tf_trim(TransferFunction *tf);

/* ================================================================
 * L5: Canonical Transformations
 * ================================================================ */

TransferFunction* tf_inverse(const TransferFunction *tf);
TransferFunction* tf_derivative(const TransferFunction *tf);
TransferFunction* tf_integral(const TransferFunction *tf);
TransferFunction* tf_frequency_shift(const TransferFunction *tf, double a);
TransferFunction* tf_time_scale(const TransferFunction *tf, double alpha);

#ifdef __cplusplus
}
#endif
#endif /* TRANSFER_FUNCTION_H */
