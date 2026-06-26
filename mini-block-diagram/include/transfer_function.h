/**
 * @file transfer_function.h
 * @brief Transfer Function �� s-domain rational functions for LTI systems
 *
 * Knowledge Coverage:
 *   L1: TransferFunction = N(s)/D(s), Polynomial, ZeroPoleGain, FreqPoint
 *   L2: Arithmetic (add, multiply, divide), interconnections (series, parallel, feedback)
 *   L3: Complex evaluation, frequency response, DC gain, Horner's method
 *   L4: Routh-Hurwitz stability criterion �� sign changes in first column
 *   L5: Bilinear (Tustin) transform s<->z, step response metrics
 *
 * The transfer function G(s) = Y(s)/U(s) is the ratio of the Laplace
 * transform of the output to the input, assuming zero initial conditions.
 * It fully characterizes the forced response of an LTI system.
 *
 * References:
 *   G.F. Franklin et al., "Feedback Control of Dynamic Systems", 8th ed., Ch. 3
 *   K. Ogata, "Modern Control Engineering", 5th ed., Ch. 2
 *   MIT 6.302 / Stanford ENGR105 / Berkeley ME132 / ETH 151-0591
 */
#ifndef TRANSFER_FUNCTION_H
#define TRANSFER_FUNCTION_H
#include <stddef.h>
#include "blockdiagram.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double *coeff; int order; } Polynomial;
typedef struct { int n_zeros,n_poles; double *zeros,*poles,gain,delay; } ZeroPoleGain;
typedef enum { TF_STABLE=0, TF_MARGINALLY_STABLE=1, TF_UNSTABLE=2 } tf_stability_t;
typedef struct { double frequency,magnitude,phase,real,imag; } FreqPoint;

Polynomial* poly_create(const double *c, int order);
void poly_destroy(Polynomial *p);
Polynomial* poly_clone(const Polynomial *s);
void poly_eval_complex(const Polynomial *p, double re, double im, double *ore, double *oim);
Polynomial* poly_add(const Polynomial *a, const Polynomial *b);
Polynomial* poly_subtract(const Polynomial *a, const Polynomial *b);
Polynomial* poly_multiply(const Polynomial *a, const Polynomial *b);
Polynomial* poly_scale(const Polynomial *p, double k);
int poly_degree(const Polynomial *p);
int poly_is_zero(const Polynomial *p);
void poly_print(const Polynomial *p, const char *var);

TransferFunction* tf_create(const double *num, int no, const double *den, int dno);
TransferFunction* tf_create_from_gain(double k);
void tf_destroy(TransferFunction *tf);
TransferFunction* tf_clone(const TransferFunction *s);
TransferFunction* tf_add(const TransferFunction *a, const TransferFunction *b);
TransferFunction* tf_subtract(const TransferFunction *a, const TransferFunction *b);
TransferFunction* tf_multiply(const TransferFunction *a, const TransferFunction *b);
TransferFunction* tf_divide(const TransferFunction *a, const TransferFunction *b);
TransferFunction* tf_inverse(const TransferFunction *tf);
TransferFunction* tf_series(const TransferFunction *g1, const TransferFunction *g2);
TransferFunction* tf_parallel(const TransferFunction *g1, const TransferFunction *g2);
TransferFunction* tf_feedback(const TransferFunction *G, const TransferFunction *H);
TransferFunction* tf_positive_feedback(const TransferFunction *G, const TransferFunction *H);
int tf_is_proper(const TransferFunction *tf);
TransferFunction* tf_make_proper(TransferFunction *tf);
double tf_dc_gain(const TransferFunction *tf);
double tf_magnitude_at(const TransferFunction *tf, double w);
double tf_phase_at(const TransferFunction *tf, double w);
void tf_eval_at(const TransferFunction *tf, double sr, double si, double *or, double *oi);
int tf_frequency_response(const TransferFunction *tf, double ws, double we, int n, FreqPoint *out);
int tf_find_roots_real(const double *c, int order, double *r, int maxr);
tf_stability_t tf_stability_routh(const TransferFunction *tf);
TransferFunction* tf_s_to_z(const TransferFunction *tf_s, double Ts);
TransferFunction* tf_z_to_s(const TransferFunction *tf_z, double Ts);
double tf_rise_time(const TransferFunction *tf);
double tf_settling_time(const TransferFunction *tf, double tol);
double tf_overshoot(const TransferFunction *tf);
double tf_steady_state_error(const TransferFunction *tf);
void tf_print(const TransferFunction *tf);

#ifdef __cplusplus
}
#endif
#endif
