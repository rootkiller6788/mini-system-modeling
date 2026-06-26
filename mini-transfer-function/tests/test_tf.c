/**
 * @file test_tf.c
 * @brief Unit tests for mini-transfer-function library
 *
 * Tests cover: L1-L6 knowledge levels
 * Uses standard assert() for all tests.
 */
#include "transfer_function.h"
#include "tf_polynomial.h"
#include "tf_analysis.h"
#include "tf_interconnections.h"
#include "tf_conversion.h"
#include "tf_advanced.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_run = 0, tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %s: ", #name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL - %s\n", msg); } while(0)
#define ASSERT_TRUE(cond, msg) do { if (!(cond)) { FAIL(msg); return; } } while(0)

/* ================================================================
 * L1: Transfer Function Construction
 * ================================================================ */

static void test_tf_create_gain(void) {
    TEST(tf_create_gain);
    TransferFunction *g = tf_create_from_gain(5.0);
    ASSERT_TRUE(g != NULL, "create failed");
    ASSERT_TRUE(fabs(tf_dc_gain(g) - 5.0) < 1e-12, "wrong DC gain");
    ASSERT_TRUE(tf_num_poles(g) == 0, "gain should have 0 poles");
    ASSERT_TRUE(tf_num_zeros(g) == 0, "gain should have 0 zeros");
    tf_destroy(g);
    PASS();
}

static void test_tf_create_integrator(void) {
    TEST(tf_create_integrator);
    TransferFunction *g = tf_create_integrator();
    ASSERT_TRUE(g != NULL, "create failed");
    ASSERT_TRUE(tf_system_type(g) == TF_TYPE_1, "integrator is type 1");
    ASSERT_TRUE(isinf(tf_dc_gain(g)), "integrator has infinite DC gain");
    tf_destroy(g);
    PASS();
}

static void test_tf_create_first_order_lag(void) {
    TEST(tf_create_first_order_lag);
    TransferFunction *g = tf_create_first_order_lag(2.0, 0.5);
    ASSERT_TRUE(g != NULL, "create failed");
    ASSERT_TRUE(fabs(tf_dc_gain(g) - 2.0) < 1e-12, "wrong DC gain");
    ASSERT_TRUE(tf_num_poles(g) == 1, "1st order has 1 pole");
    ASSERT_TRUE(tf_properness(g) == TF_STRICTLY_PROPER, "1st order lag is strictly proper");
    tf_destroy(g);
    PASS();
}

static void test_tf_create_second_order(void) {
    TEST(tf_create_second_order);
    TransferFunction *g = tf_create_second_order(1.0, 2.0, 0.5);
    ASSERT_TRUE(g != NULL, "create failed");
    ASSERT_TRUE(fabs(tf_dc_gain(g) - 1.0) < 1e-12, "DC gain should be 1");
    ASSERT_TRUE(tf_num_poles(g) == 2, "2nd order has 2 poles");
    tf_destroy(g);
    PASS();
}

static void test_tf_create_pid(void) {
    TEST(tf_create_pid);
    TransferFunction *g = tf_create_pid(1.0, 0.5, 0.1);
    ASSERT_TRUE(g != NULL, "create failed");
    ASSERT_TRUE(tf_system_type(g) == TF_TYPE_1, "PID has type 1 (integrator)");
    /* PID: (Kd*s^2 + Kp*s + Ki)/s - numerator order 2, den order 1 */
    ASSERT_TRUE(g->num_order == 2 && g->den_order == 1, "wrong PID orders");
    tf_destroy(g);
    PASS();
}

/* ================================================================
 * L2: Core Properties
 * ================================================================ */

static void test_tf_properness(void) {
    TEST(tf_properness);
    TransferFunction *strict = tf_create_first_order_lag(1.0, 1.0);
    ASSERT_TRUE(tf_properness(strict) == TF_STRICTLY_PROPER, "lag should be strictly proper");

    double num[2] = {1.0, 2.0}, den[2] = {1.0, 1.0};
    TransferFunction *proper_tf = tf_create(num, 1, den, 1);
    ASSERT_TRUE(tf_properness(proper_tf) == TF_PROPER, "equal order is proper");

    double num2[3] = {1.0, 2.0, 3.0}, den2[2] = {1.0, 1.0};
    TransferFunction *improper = tf_create(num2, 2, den2, 1);
    ASSERT_TRUE(tf_properness(improper) == TF_IMPROPER, "more zeros is improper");

    tf_destroy(strict); tf_destroy(proper_tf); tf_destroy(improper);
    PASS();
}

static void test_tf_system_type(void) {
    TEST(tf_system_type);
    TransferFunction *t0 = tf_create_from_gain(1.0);
    ASSERT_TRUE(tf_system_type(t0) == TF_TYPE_0, "gain is type 0");

    TransferFunction *t1 = tf_create_integrator();
    ASSERT_TRUE(tf_system_type(t1) == TF_TYPE_1, "integrator is type 1");

    double num[1] = {1.0}, den[3] = {0.0, 0.0, 1.0};
    TransferFunction *t2 = tf_create(num, 0, den, 2);
    ASSERT_TRUE(tf_system_type(t2) == TF_TYPE_2, "1/s^2 is type 2");

    tf_destroy(t0); tf_destroy(t1); tf_destroy(t2);
    PASS();
}

/* ================================================================
 * L3: Polynomial Operations
 * ================================================================ */

static void test_poly_create_eval(void) {
    TEST(poly_create_eval);
    double c[3] = {1.0, 2.0, 3.0}; /* 1 + 2s + 3s^2 */
    Polynomial *p = poly_create(c, 2);
    ASSERT_TRUE(p != NULL, "create failed");
    ASSERT_TRUE(fabs(poly_eval(p, 0.0) - 1.0) < 1e-12, "P(0)=1");
    ASSERT_TRUE(fabs(poly_eval(p, 1.0) - 6.0) < 1e-12, "P(1)=1+2+3=6");
    ASSERT_TRUE(poly_degree(p) == 2, "degree should be 2");
    poly_destroy(p);
    PASS();
}

static void test_poly_arithmetic(void) {
    TEST(poly_arithmetic);
    double c1[2] = {1.0, 1.0}; /* 1 + s */
    double c2[2] = {2.0, 1.0}; /* 2 + s */
    Polynomial *a = poly_create(c1, 1);
    Polynomial *b = poly_create(c2, 1);

    Polynomial *sum = poly_add(a, b);
    ASSERT_TRUE(fabs(poly_eval(sum, 0.0) - 3.0) < 1e-12, "add const");
    ASSERT_TRUE(fabs(poly_eval(sum, 1.0) - 5.0) < 1e-12, "add at s=1");

    Polynomial *prod = poly_multiply(a, b);
    ASSERT_TRUE(fabs(poly_eval(prod, 0.0) - 2.0) < 1e-12, "mul const");
    ASSERT_TRUE(fabs(poly_eval(prod, 1.0) - 6.0) < 1e-12, "mul at s=1");

    poly_destroy(a); poly_destroy(b); poly_destroy(sum); poly_destroy(prod);
    PASS();
}

static void test_poly_derivative(void) {
    TEST(poly_derivative);
    double c[3] = {3.0, 2.0, 1.0}; /* 3 + 2s + s^2 */
    Polynomial *p = poly_create(c, 2);
    Polynomial *dp = poly_derivative(p);

    ASSERT_TRUE(fabs(poly_eval(dp, 0.0) - 2.0) < 1e-12, "dP/ds(0)=2");
    ASSERT_TRUE(fabs(poly_eval(dp, 1.0) - 4.0) < 1e-12, "dP/ds(1)=4");

    poly_destroy(p); poly_destroy(dp);
    PASS();
}

/* ================================================================
 * L4: Stability Theorems
 * ================================================================ */

static void test_routh_stable(void) {
    TEST(routh_stable);
    /* Stable 2nd-order: s^2 + 2*s + 1 = (s+1)^2 */
    TransferFunction *g = tf_create_second_order(1.0, 1.0, 1.0);
    /* G(s) = 1/(s^2 + 2s + 1) - critically damped */
    double num[1] = {1.0}, den[3] = {1.0, 2.0, 1.0};
    TransferFunction *g2 = tf_create(num, 0, den, 2);
    ASSERT_TRUE(tf_stability_routh(g2) == TF_STABLE, "should be stable");
    tf_destroy(g); tf_destroy(g2);
    PASS();
}

static void test_routh_unstable(void) {
    TEST(routh_unstable);
    /* Unstable: s^2 - 1 = (s+1)(s-1) */
    double num[1] = {1.0}, den[3] = {-1.0, 0.0, 1.0};
    TransferFunction *g = tf_create(num, 0, den, 2);
    ASSERT_TRUE(tf_stability_routh(g) == TF_UNSTABLE, "should be unstable");
    tf_destroy(g);
    PASS();
}

static void test_final_value_theorem(void) {
    TEST(final_value_theorem);
    TransferFunction *g = tf_create_first_order_lag(3.0, 0.5);
    double result;
    int ok = tf_final_value_theorem(g, &result);
    ASSERT_TRUE(ok, "FVT should apply");
    ASSERT_TRUE(fabs(result - 3.0) < 1e-12, "final value should be 3");
    tf_destroy(g);
    PASS();
}

/* ================================================================
 * L4: Frequency Domain
 * ================================================================ */

static void test_freq_response(void) {
    TEST(freq_response);
    TransferFunction *g = tf_create_first_order_lag(1.0, 1.0);
    ASSERT_TRUE(g != NULL, "create failed");

    /* At s=0: G(0) = 1 */
    double re, im;
    tf_eval_at(g, 0.0, 0.0, &re, &im);
    ASSERT_TRUE(fabs(re - 1.0) < 1e-10, "G(0) real=1");
    ASSERT_TRUE(fabs(im) < 1e-10, "G(0) imag=0");

    /* At high frequency: G(j*inf) -> 0 */
    double mag_high = tf_magnitude_at(g, 1000.0);
    ASSERT_TRUE(mag_high < 0.01, "high freq gain should be small");

    tf_destroy(g);
    PASS();
}

static void test_step_response_chars(void) {
    TEST(step_response_chars);
    TransferFunction *g = tf_create_second_order(1.0, 1.0, 0.5);
    /* PO = 100*exp(-pi*0.5/sqrt(0.75)) ~ 16.3% */
    double po = tf_overshoot(g);
    ASSERT_TRUE(po > 0.0 && po < 100.0, "overshoot in range");
    double tr = tf_rise_time(g);
    ASSERT_TRUE(tr > 0.0, "rise time positive");
    double ts = tf_settling_time(g, 0.02);
    ASSERT_TRUE(ts > 0.0, "settling time positive");
    tf_destroy(g);
    PASS();
}

/* ================================================================
 * L5: Interconnections and Conversions
 * ================================================================ */

static void test_feedback_connection(void) {
    TEST(feedback_connection);
    TransferFunction *G = tf_create_first_order_lag(2.0, 1.0);
    TransferFunction *H = tf_create_from_gain(1.0);
    TransferFunction *cl = tf_feedback(G, H);

    ASSERT_TRUE(cl != NULL, "feedback failed");
    /* G_cl = G/(1+G) = (2/(s+1))/(1+2/(s+1)) = 2/(s+3) */
    ASSERT_TRUE(fabs(tf_dc_gain(cl) - 2.0/3.0) < 1e-10, "wrong closed-loop gain");

    tf_destroy(G); tf_destroy(H); tf_destroy(cl);
    PASS();
}

static void test_tf_to_ss(void) {
    TEST(tf_to_ss);
    TransferFunction *g = tf_create_first_order_lag(2.0, 1.0);
    StateSpace *ss = tf_to_ss_controllable(g);
    ASSERT_TRUE(ss != NULL, "conversion failed");
    ASSERT_TRUE(ss->n == 1, "1st order = 1 state");

    TransferFunction *g2 = ss_to_tf(ss);
    ASSERT_TRUE(g2 != NULL, "reverse conversion failed");
    ASSERT_TRUE(fabs(tf_dc_gain(g2) - 2.0) < 1e-8, "DC gain preserved");

    tf_destroy(g); ss_destroy(ss); tf_destroy(g2);
    PASS();
}

static void test_series_parallel(void) {
    TEST(series_parallel);
    TransferFunction *g1 = tf_create_first_order_lag(1.0, 1.0);
    TransferFunction *g2 = tf_create_first_order_lag(2.0, 1.0);

    TransferFunction *ser = tf_series(g1, g2);
    ASSERT_TRUE(ser != NULL, "series failed");
    ASSERT_TRUE(fabs(tf_dc_gain(ser) - 2.0) < 1e-10, "series DC gain");

    TransferFunction *par = tf_parallel(g1, g2);
    ASSERT_TRUE(par != NULL, "parallel failed");
    ASSERT_TRUE(fabs(tf_dc_gain(par) - 3.0) < 1e-10, "parallel DC gain");

    tf_destroy(g1); tf_destroy(g2); tf_destroy(ser); tf_destroy(par);
    PASS();
}

/* ================================================================
 * L6: Canonical Systems
 * ================================================================ */

static void test_dc_motor(void) {
    TEST(dc_motor);
    /* Realistic DC motor parameters */
    TransferFunction *motor = tf_dc_motor_speed(0.01, 0.1, 0.01, 0.01, 1.0, 0.5);
    ASSERT_TRUE(motor != NULL, "DC motor TF failed");
    ASSERT_TRUE(tf_num_poles(motor) == 2, "DC motor is 2nd order");
    ASSERT_TRUE(tf_stability_routh(motor) != TF_UNSTABLE, "motor should be stable");
    tf_destroy(motor);
    PASS();
}

static void test_mass_spring_damper(void) {
    TEST(mass_spring_damper);
    TransferFunction *msd = tf_mass_spring_damper(1.0, 0.5, 10.0);
    ASSERT_TRUE(msd != NULL, "MSD TF failed");
    ASSERT_TRUE(tf_num_poles(msd) == 2, "MSD is 2nd order");
    tf_destroy(msd);
    PASS();
}

static void test_rlc_lowpass(void) {
    TEST(rlc_lowpass);
    TransferFunction *rlc = tf_rlc_lowpass(100.0, 0.001, 0.000001);
    ASSERT_TRUE(rlc != NULL, "RLC TF failed");
    ASSERT_TRUE(tf_num_poles(rlc) == 2, "RLC is 2nd order");
    tf_destroy(rlc);
    PASS();
}

/* ================================================================
 * L7: Applications
 * ================================================================ */

static void test_suspension(void) {
    TEST(quarter_car_suspension);
    /* Tesla Model S-like parameters */
    TransferFunction *susp = tf_quarter_car_suspension(400.0, 50.0, 20000.0, 1500.0, 200000.0);
    ASSERT_TRUE(susp != NULL, "suspension TF failed");
    ASSERT_TRUE(tf_num_poles(susp) == 4, "quarter-car is 4th order");
    tf_destroy(susp);
    PASS();
}

static void test_inverted_pendulum(void) {
    TEST(inverted_pendulum);
    TransferFunction *invp = tf_inverted_pendulum(1.0, 0.1, 0.5);
    ASSERT_TRUE(invp != NULL, "inverted pendulum TF failed");
    ASSERT_TRUE(tf_stability_routh(invp) == TF_UNSTABLE, "inverted pendulum is unstable");
    tf_destroy(invp);
    PASS();
}

/* ================================================================
 * L8: Advanced Methods
 * ================================================================ */

static void test_pade_approximation(void) {
    TEST(pade_approximation);
    /* Pade [1/1] of exp(-s): (2-s)/(2+s) */
    TransferFunction *pade = tf_pade_approximation(1.0, 1, 1);
    ASSERT_TRUE(pade != NULL, "Pade failed");
    ASSERT_TRUE(tf_num_poles(pade) == 1, "Pade [1/1] has 1 pole");
    ASSERT_TRUE(tf_stability_routh(pade) != TF_UNSTABLE, "Pade should be stable");
    tf_destroy(pade);
    PASS();
}

static void test_fractional_integrator(void) {
    TEST(fractional_integrator);
    /* Half-integrator: 1/s^0.5 */
    TransferFunction *frac = tf_fractional_integrator(0.5, 3, 0.01, 100.0);
    ASSERT_TRUE(frac != NULL, "fractional integrator failed");
    ASSERT_TRUE(tf_num_poles(frac) >= 1, "should have poles from approximation");
    tf_destroy(frac);
    PASS();
}

/* ================================================================
 * Main
 * ================================================================ */

int main(void) {
    printf("=== mini-transfer-function Unit Tests ===\n\n");

    printf("-- L1: Definitions --\n");
    test_tf_create_gain();
    test_tf_create_integrator();
    test_tf_create_first_order_lag();
    test_tf_create_second_order();
    test_tf_create_pid();

    printf("\n-- L2: Core Properties --\n");
    test_tf_properness();
    test_tf_system_type();

    printf("\n-- L3: Polynomial Algebra --\n");
    test_poly_create_eval();
    test_poly_arithmetic();
    test_poly_derivative();

    printf("\n-- L4: Stability and Frequency --\n");
    test_routh_stable();
    test_routh_unstable();
    test_final_value_theorem();
    test_freq_response();
    test_step_response_chars();

    printf("\n-- L5: Interconnections and Conversions --\n");
    test_feedback_connection();
    test_tf_to_ss();
    test_series_parallel();

    printf("\n-- L6: Canonical Systems --\n");
    test_dc_motor();
    test_mass_spring_damper();
    test_rlc_lowpass();

    printf("\n-- L7: Applications --\n");
    test_suspension();
    test_inverted_pendulum();

    printf("\n-- L8: Advanced Methods --\n");
    test_pade_approximation();
    test_fractional_integrator();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
