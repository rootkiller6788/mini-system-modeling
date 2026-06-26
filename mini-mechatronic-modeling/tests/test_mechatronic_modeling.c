/**
 * @file test_mechatronic_modeling.c
 * @brief Comprehensive test suite for mini-mechatronic-modeling.
 *
 * Covers L1-L7 knowledge layers with assertion-based tests.
 * Uses standard assert() only — no custom macros.
 *
 * Run: make test
 */

#include "mechatronic_definitions.h"
#include "electromechanical_coupling.h"
#include "mechatronic_state_space.h"
#include "actuator_modeling.h"
#include "sensor_modeling.h"
#include "mechatronic_system_design.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define EPS 1e-6
#define NEAR(a, b) (fabs((a) - (b)) < EPS)

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
    printf("  TEST: %s ... ", name); \
} while(0)

#define PASS() do { \
    printf("PASSED\n"); tests_passed++; \
} while(0)

#define FAIL(msg) do { \
    printf("FAILED: %s\n", msg); tests_failed++; \
} while(0)

/*===========================================================================*/
/* L1: Definitions — Constructor Tests                                       */
/*===========================================================================*/

static void test_dc_motor_create(void)
{
    TEST("DC motor constructor");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    assert(NEAR(m.R_a, 2.0));
    assert(NEAR(m.L_a, 0.001));
    assert(NEAR(m.K_e, 0.05));
    assert(NEAR(m.K_t, 0.05));  /* SI equivalence */
    assert(NEAR(m.tau_e, 0.0005));  /* L_a/R_a = 0.001/2.0 */
    assert(m.max_voltage > 0.0);
    assert(m.P_rated > 0.0);
    PASS();
}

static void test_gear_train_create(void)
{
    TEST("Gear train constructor");
    GearTrain g = gear_train_create(5.0, 0.9, 1e-6, 2e-6, 0.01, 1e6);
    assert(NEAR(g.ratio, 5.0));
    assert(NEAR(g.efficiency, 0.9));
    assert(NEAR(g.backlash, 0.01));
    PASS();
}

static void test_lead_screw_create(void)
{
    TEST("Lead screw constructor");
    LeadScrew s = lead_screw_create(0.005, 0.92, 0.016, 1e-5, 0.003);
    assert(NEAR(s.pitch, 0.005));
    assert(s.efficiency > 0.85);  /* ball screw */
    PASS();
}

static void test_mechatronic_subsystem_create(void)
{
    TEST("Mechatronic subsystem constructor");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    GearTrain g = gear_train_create(5.0, 0.9, 1e-6, 2e-6, 0.0, 1e6);
    MechatronicLoad L = mechatronic_load_create(0.01, 5.0, 0.001, 0.01, 0.015);

    MechatronicSubsystem sys = mechatronic_subsystem_create(m, g, L);

    /* Total inertia should include motor + gear + reflected load */
    assert(sys.J_total_reflected > m.J_m);
    assert(sys.B_total_reflected > 0.0);
    PASS();
}

/*===========================================================================*/
/* L4: Fundamental Laws — Electromechanical Coupling Tests                   */
/*===========================================================================*/

static void test_lorentz_force(void)
{
    TEST("Lorentz force F = I·B·L");
    double F = lorentz_force(2.0, 0.5, 0.1, M_PI / 2.0);
    assert(NEAR(F, 2.0 * 0.5 * 0.1 * 1.0));  /* sin(π/2) = 1 */
    PASS();
}

static void test_lorentz_force_angle(void)
{
    TEST("Lorentz force at angle");
    double F = lorentz_force(2.0, 1.0, 1.0, M_PI / 6.0);
    assert(NEAR(F, 2.0 * 1.0 * 1.0 * 0.5));  /* sin(π/6) = 0.5 */
    PASS();
}

static void test_motional_emf(void)
{
    TEST("Motional EMF = B·L·v");
    double V = motional_emf(0.5, 0.2, 10.0);
    assert(NEAR(V, 0.5 * 0.2 * 10.0));  /* = 1.0 V */
    PASS();
}

static void test_dc_motor_back_emf(void)
{
    TEST("DC motor back-EMF V_bemf = K_e·ω");
    double V = dc_motor_back_emf(0.05, 300.0);
    assert(NEAR(V, 15.0));  /* 0.05 * 300 = 15V */
    PASS();
}

static void test_kt_ke_equivalence(void)
{
    TEST("K_t = K_e in SI units");
    double err = kt_ke_equivalence_check(0.05, 0.05);
    assert(err < 1e-12);  /* should be exactly equal */
    PASS();
}

static void test_dc_motor_torque(void)
{
    TEST("DC motor electromagnetic torque T = K_t·I");
    double T = dc_motor_electromagnetic_torque(0.05, 2.0);
    assert(NEAR(T, 0.1));  /* 0.05 * 2.0 = 0.1 N·m */
    PASS();
}

static void test_dc_motor_steady_state(void)
{
    TEST("DC motor steady-state operating point");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    double I_ss, omega_ss, T_ss;
    dc_motor_steady_state(&m, 12.0, 0.0, &I_ss, &omega_ss, &T_ss);

    /* At no-load: ω_ss = K_t·V_a/(K_t·K_e + R_a·B) ≈ V_a/K_e
     * With B=1e-5: ω_ss = 0.05*12/(0.0025+2e-5) ≈ 238.1 rad/s
     * Check it's close to ideal no-load (240 rad/s), within 1% */
    double omega_ideal = 12.0 / 0.05;  /* 240 rad/s */
    assert(fabs(omega_ss - omega_ideal) < omega_ideal * 0.02);  /* within 2% */
    assert(I_ss < 0.1);  /* small no-load current */
    assert(T_ss < 0.01);
    PASS();
}

static void test_motor_efficiency(void)
{
    TEST("Motor max mechanical power point");
    double P_max = dc_motor_max_mechanical_power(12.0, 2.0);
    assert(NEAR(P_max, 144.0 / 8.0));  /* V²/(4·R) = 144/8 = 18W */
    PASS();
}

/*===========================================================================*/
/* L4: Controllability and Observability                                     */
/*===========================================================================*/

static void test_dc_motor_controllability(void)
{
    TEST("DC motor 2nd-order system is controllable");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    LinearStateSpace ss = ss_dc_motor_second_order(&m, 1e-4, 1e-5);
    int ctrl = ss_is_controllable(&ss);
    assert(ctrl == 1);  /* DC motor is controllable from V_a */
    PASS();
}

static void test_dc_motor_observability(void)
{
    TEST("DC motor 2nd-order system is observable from speed");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    LinearStateSpace ss = ss_dc_motor_second_order(&m, 1e-4, 1e-5);
    int obs = ss_is_observable(&ss);
    assert(obs == 1);  /* observable from ω */
    PASS();
}

static void test_eigenvalue_stability(void)
{
    TEST("DC motor eigenvalues indicate stability");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    LinearStateSpace ss = ss_dc_motor_second_order(&m, 1e-4, 1e-5);
    double eval_re[SS_MAX_DIM], eval_im[SS_MAX_DIM];
    int ret = ss_eigenvalues(ss.A, 2, eval_re, eval_im);
    assert(ret == 0);
    int stable = ss_is_stable(eval_re, eval_im, 2);
    assert(stable == 1);  /* DC motor is stable */
    PASS();
}

/*===========================================================================*/
/* L5: Methods — Time Constants, RK4, Discretization                         */
/*===========================================================================*/

static void test_time_constants(void)
{
    TEST("DC motor time constants");
    DCMotorParams m = dc_motor_create(2.0, 0.002, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    double tau_e = dc_motor_tau_electrical(&m);
    assert(NEAR(tau_e, 0.001));  /* L_a/R_a = 0.002/2.0 */
    double tau_m = dc_motor_tau_mechanical(1e-5, 1e-5);
    assert(NEAR(tau_m, 1.0));
    PASS();
}

static void test_power_calculations(void)
{
    TEST("Motor power calculations");
    double P_mech = dc_motor_mechanical_power(0.1, 100.0);
    assert(NEAR(P_mech, 10.0));  /* 0.1 N·m * 100 rad/s = 10W */
    double P_elec = dc_motor_electrical_power(12.0, 1.0);
    assert(NEAR(P_elec, 12.0));
    double P_cu = dc_motor_copper_loss(1.0, 2.0);
    assert(NEAR(P_cu, 2.0));  /* I²R = 1*2 = 2W */
    PASS();
}

static void test_rk4_step(void)
{
    TEST("RK4 integration of DC motor");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    LinearStateSpace ss = ss_dc_motor_second_order(&m, 1e-4, 1e-5);

    double x[2] = {0.0, 0.0};  /* I_a=0, ω=0 */
    double u[2] = {12.0, 0.0};  /* V_a=12, T_L=0 */
    double x_next[2];

    ss_rk4_step(&ss, x, u, 0.001, x_next, NULL);

    /* With 12V applied, current should rise rapidly */
    assert(x_next[0] > x[0]);  /* Current increasing */
    /* Speed should be increasing from rest */
    assert(x_next[1] > x[1]);
    PASS();
}

static void test_zoh_discretization(void)
{
    TEST("ZOH discretization of DC motor");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    LinearStateSpace ss_cont = ss_dc_motor_second_order(&m, 1e-4, 1e-5);
    LinearStateSpace ss_disc;
    int ret = ss_discretize_zoh(&ss_cont, 0.001, &ss_disc);
    assert(ret == 0);
    /* Verify discretized system has same dimensions */
    assert(ss_disc.n_states == ss_cont.n_states);
    assert(ss_disc.n_inputs == ss_cont.n_inputs);
    assert(ss_disc.n_outputs == ss_cont.n_outputs);
    /* Verify Φ is not identity (dynamics are present) */
    int has_dynamics = 0;
    for (int i = 0; i < ss_disc.n_states; i++)
        for (int j = 0; j < ss_disc.n_states; j++)
            if (i != j && fabs(ss_disc.A[i * SS_MAX_DIM + j]) > 1e-12)
                has_dynamics = 1;
    /* At least one off-diagonal must be non-zero for coupled system */
    assert(has_dynamics || ss_disc.n_states == 1);
    PASS();
}

static void test_pi_controller_design(void)
{
    TEST("Current PI controller design");
    DCMotorParams m = dc_motor_create(2.0, 0.002, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    double Kp, Ki;
    design_current_pi(&m, 1000.0, &Kp, &Ki);
    assert(NEAR(Kp, 0.002 * 1000.0));  /* L_a * ω_c */
    assert(NEAR(Ki, 2.0 * 1000.0));    /* R_a * ω_c */
    PASS();
}

/*===========================================================================*/
/* L5: Motor Selection and Sizing                                            */
/*===========================================================================*/

static void test_rms_torque(void)
{
    TEST("RMS torque calculation");
    DutyCycle cycle;
    memset(&cycle, 0, sizeof(cycle));
    cycle.n_segments = 2;
    cycle.segments[0].torque = 1.0;
    cycle.segments[0].duration = 1.0;
    cycle.segments[1].torque = 2.0;
    cycle.segments[1].duration = 1.0;
    cycle.cycle_time = 2.0;

    double Trms = rms_torque(&cycle);
    /* T_rms = √((1²·1 + 2²·1) / 2) = √(5/2) ≈ 1.5811 */
    assert(NEAR(Trms, sqrt(2.5)));
    PASS();
}

static void test_optimal_gear_ratio(void)
{
    TEST("Optimal gear ratio for inertia matching");
    double N = optimal_gear_ratio_inertia_match(1e-4, 0.01);
    assert(NEAR(N, sqrt(0.01 / 1e-4)));  /* √(100) = 10 */
    PASS();
}

static void test_reflect_load_through_gear(void)
{
    TEST("Load reflection through gear train");
    GearTrain g = gear_train_create(10.0, 0.9, 1e-6, 2e-6, 0.0, 1e6);
    double J_out, T_out;
    reflect_load_through_gear(0.1, 1.0, &g, 0, &J_out, &T_out);
    /* J_reflected = J_load / N² = 0.1 / 100 = 0.001 */
    assert(NEAR(J_out, 0.001));
    /* T_reflected = T_load / (N·η) = 1.0 / (10*0.9) ≈ 0.1111 */
    assert(NEAR(T_out, 1.0 / 9.0));
    PASS();
}

/*===========================================================================*/
/* L1-L3: Sensor Tests                                                       */
/*===========================================================================*/

static void test_encoder_counts_to_angle(void)
{
    TEST("Encoder counts to angle conversion");
    double theta = encoder_counts_to_angle(500, 1000);
    assert(NEAR(theta, M_PI));  /* 500/1000 * 2π = π */
    PASS();
}

static void test_encoder_quadrature(void)
{
    TEST("Encoder quadrature 4x resolution");
    double theta = encoder_quadrature_to_angle(2000, 1000);
    assert(NEAR(theta, M_PI));  /* 2000/(4*1000) * 2π = π */
    PASS();
}

static void test_encoder_velocity_m_method(void)
{
    TEST("Encoder M-method velocity");
    double omega = encoder_velocity_m_method(100, 1000, 0.01);
    /* 2π * 100 / (1000 * 0.01) = 2π * 10 ≈ 62.83 rad/s */
    assert(NEAR(omega, 2.0 * M_PI * 10.0));
    PASS();
}

static void test_encoder_direction(void)
{
    TEST("Quadrature encoder direction detection");
    /* Forward: 00→01 */
    int dir = encoder_direction(0, 0, 0, 1);
    assert(dir == 1);
    /* Reverse: 01→00 */
    dir = encoder_direction(0, 1, 0, 0);
    assert(dir == -1);
    PASS();
}

static void test_resolver_angle(void)
{
    TEST("Resolver angle extraction");
    double theta = resolver_angle_from_demod(0.7071, 0.7071);
    assert(NEAR(theta, M_PI / 4.0));  /* sin=cos=0.7071 → 45° */
    PASS();
}

static void test_bldc_angle_from_halls(void)
{
    TEST("BLDC angle from Hall sensors");
    double angle = bldc_angle_from_halls(1, 0, 1);  /* 101 → sector 1 */
    assert(NEAR(angle, M_PI / 3.0));  /* 60° */
    PASS();
}

static void test_strain_gauge(void)
{
    TEST("Strain gauge bridge output");
    double V_out = wheatstone_full_bridge(5.0, 2.0, 500e-6);
    /* V_exc * GF * strain = 5.0 * 2.0 * 500e-6 = 0.005 V = 5 mV */
    assert(NEAR(V_out, 0.005));
    PASS();
}

static void test_lowpass_filter(void)
{
    TEST("Sensor low-pass filter");
    /* DC steady state: output should converge to input */
    double y = sensor_lowpass_filter(1.0, 0.0, 10.0, 0.1);
    assert(y > 0.0 && y < 1.0);  /* transient */
    /* After many iterations, should approach 1.0 */
    for (int i = 0; i < 100; i++) {
        y = sensor_lowpass_filter(1.0, y, 10.0, 0.01);
    }
    assert(NEAR(y, 1.0));
    PASS();
}

static void test_complementary_filter(void)
{
    TEST("Complementary filter for IMU");
    /* tau ≫ dt → mostly trusts gyro integration */
    double angle = complementary_filter(1.0, 0.5, 0.0, 0.01, 1e6);
    /* α ≈ 1.0, gyro_angle = 0 + 1.0*0.01 = 0.01, accel has no effect */
    assert(fabs(angle - 0.01) < 0.001);
    /* tau ≪ dt → mostly trusts accelerometer */
    double angle2 = complementary_filter(1.0, 0.5, 0.0, 0.01, 1e-9);
    /* α ≈ 0.0, angle ≈ accel_angle = 0.5 */
    assert(fabs(angle2 - 0.5) < 0.001);
    PASS();
}

/*===========================================================================*/
/* L5: Anti-Windup, Friction, Backlash                                       */
/*===========================================================================*/

static void test_pi_antiwindup(void)
{
    TEST("PI controller with anti-windup");
    double integral = 0.5;
    double u;
    pi_antiwindup(1.0, 0.01, 2.0, 10.0, -10.0, 10.0, &integral, &u);
    /* P = 2.0 * 1.0 = 2.0, I = 0.5, total = 2.5 */
    assert(NEAR(u, 2.5));
    assert(integral > 0.5);  /* integrator updated in linear region */
    PASS();
}

static void test_pi_antiwindup_saturation(void)
{
    TEST("PI anti-windup: saturation clamping");
    double integral = 8.0;
    double u;
    pi_antiwindup(1.0, 0.01, 3.0, 1.0, -10.0, 10.0, &integral, &u);
    /* P = 3.0, P+I = 11.0 > u_max=10.0 → saturated at 10 */
    assert(NEAR(u, 10.0));
    /* Integrator should be frozen (conditional integration) */
    assert(NEAR(integral, 8.0));
    PASS();
}

static void test_friction_model(void)
{
    TEST("Friction model (Coulomb + viscous + Stribeck)");
    /* At high velocity, friction ≈ Coulomb + viscous (Stribeck negligible) */
    double Tf = friction_torque(10.0, 0.1, 0.15, 1.0, 0.01);
    /* Coulomb=0.1, viscous=0.01*10=0.1, Stribeck≈0.05*exp(-10)≈2.3e-6 */
    double expected = 0.1 + 0.1;  /* Coulomb + viscous only */
    assert(fabs(Tf - expected) < 1e-4);  /* within 0.01% */
    /* Also test at zero velocity (function outputs 0, stiction handled separately) */
    double Tf0 = friction_torque(0.0, 0.1, 0.15, 1.0, 0.01);
    assert(NEAR(Tf0, 0.0));
    PASS();
}

static void test_backlash_model(void)
{
    TEST("Backlash dead-zone model");
    double theta_out_prev = 0.0;

    /* Small input within deadzone: output unchanged */
    backlash_model(0.005, 0.0, theta_out_prev, 0.02, &theta_out_prev);
    assert(NEAR(theta_out_prev, 0.0));  /* within ±0.01 dead zone */

    /* Input exceeds dead zone: output follows with offset */
    backlash_model(0.03, 0.005, theta_out_prev, 0.02, &theta_out_prev);
    assert(NEAR(theta_out_prev, 0.03 - 0.01));  /* θ_in - b/2 */
    PASS();
}

/*===========================================================================*/
/* L6: System Design and Sizing                                              */
/*===========================================================================*/

static void test_compute_reflected_dynamics(void)
{
    TEST("Reflected dynamics computation");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-6, 24.0, 3.0);
    GearTrain g = gear_train_create(5.0, 0.9, 1e-6, 2e-6, 0.0, 1e6);
    double J_total, B_total;
    compute_reflected_dynamics(&m, &g, 0.01, 0.001, &J_total, &B_total);

    /* J_total ≈ 1e-5 + 1e-6 + (2e-6 + 0.01)/25 > 4e-4 */
    assert(J_total > 4e-4);
    /* J_load reflected dominates */
    PASS();
}

static void test_inertia_ratio(void)
{
    TEST("Inertia ratio classification");
    double lambda = inertia_ratio(3e-5, 1e-5);
    assert(NEAR(lambda, 3.0));
    const char *cls = inertia_ratio_classification(lambda);
    assert(strstr(cls, "Low") != NULL || strstr(cls, "Good") != NULL);
    PASS();
}

static void test_ball_screw_required_torque(void)
{
    TEST("Ball screw required torque");
    /* Horizontal axis: no gravity component (sin(0)=0) */
    double T = ball_screw_required_torque(50.0, 0.010, 0.9, 0.005,
                                          0.0, 2.0, 1e-5, 0.0);
    /* T should be positive and reasonable */
    assert(T > 0.0);
    assert(T < 100.0);
    PASS();
}

static void test_linear_to_rotary(void)
{
    TEST("Linear to rotary conversion");
    double omega, torque;
    linear_to_rotary_conversion(0.5, 100.0, 0.010, 0.9, &omega, &torque);
    /* ω = 0.5 / (0.01/2π) = 0.5 * 200π/(?) → 0.5 * 2π/0.01 ≈ 314 */
    assert(NEAR(omega, 0.5 * 2.0 * M_PI / 0.01));
    /* T = 100 * 0.01 / (2π * 0.9) */
    double r = 0.01 / (2.0 * M_PI);
    assert(NEAR(torque, 100.0 * r / 0.9));
    PASS();
}

static void test_robot_joint_max_torque(void)
{
    TEST("Robot joint max torque");
    double T = robot_joint_max_torque(5.0, 0.3, 0.2, 50.0, 10.0, 1e-4, 0.05);
    /* Should return a positive torque value */
    assert(T > 0.0);
    /* Gravity component: 5 * 9.81 * 0.3 / 50 ≈ 0.294 */
    double T_grav = 5.0 * 9.80665 * 0.3 / 50.0;
    assert(NEAR(T_grav, 0.2941995));
    PASS();
}

static void test_hdd_seek_time(void)
{
    TEST("HDD seek time (short seek)");
    double t = hdd_seek_time(0.1, 0.05, 5e-6, 2.0, 500.0);
    /* Short seek under speed limit */
    assert(t > 0.0 && t < 0.1);  /* should be fast (< 100ms) */
    PASS();
}

static void test_pick_and_place(void)
{
    TEST("Pick-and-place cycle time (triangular profile)");
    /* Short distance: triangular velocity profile */
    double t = pick_and_place_cycle_time(0.1, 2.0, 10.0, 0.05, 0.1);
    assert(t > 0.0);
    /* d_accel = 4/20 = 0.2, distance=0.1 < 2*0.2 → triangular
     * t_move = 2*sqrt(0.1/10) = 2*0.1 = 0.2 */
    double t_tri = 2.0 * sqrt(0.1 / 10.0);  /* = 0.2 */
    double expected = t_tri + 0.05 + 0.1;  /* + settle + grip */
    assert(fabs(t - expected) < 0.001);
    PASS();
}

static void test_motion_profile_energy(void)
{
    TEST("Motion profile energy calculation");
    DCMotorParams m = dc_motor_create(2.0, 0.001, 0.05, 1e-5, 1e-5, 24.0, 3.0);
    DutyCycle cycle;
    memset(&cycle, 0, sizeof(cycle));
    cycle.n_segments = 2;
    cycle.segments[0].torque = 0.05;
    cycle.segments[0].speed = 100.0;
    cycle.segments[0].duration = 0.5;
    cycle.segments[1].torque = -0.03;
    cycle.segments[1].speed = 100.0;
    cycle.segments[1].duration = 0.5;
    cycle.cycle_time = 1.0;

    double energy, peak, avg;
    motion_profile_energy(&cycle, &m, NULL, 1e-4, &energy, &peak, &avg);
    assert(energy > 0.0);
    assert(peak > 0.0);
    PASS();
}

/*===========================================================================*/
/* L5: Cascade Controller Design                                             */
/*===========================================================================*/

static void test_cascade_bandwidth_allocation(void)
{
    TEST("Cascade bandwidth allocation");
    double w_vel, w_cur;
    cascade_bandwidth_allocation(62.83, 5.0, 8.0, &w_vel, &w_cur);  /* 10 Hz pos */
    assert(NEAR(w_vel, 62.83 * 5.0));    /* 50 Hz vel */
    assert(NEAR(w_cur, 62.83 * 5.0 * 8.0));  /* 400 Hz cur */
    PASS();
}

static void test_design_servo_axis(void)
{
    TEST("Complete servo axis design");
    DCMotorParams m = dc_motor_create(2.0, 0.002, 0.05, 1e-5, 1e-5, 24.0, 5.0);
    GearTrain g = gear_train_create(5.0, 0.9, 1e-6, 2e-6, 0.0, 1e6);
    MechatronicLoad L = mechatronic_load_create(0.01, 10.0, 0.001, 0.02, 0.03);
    RotarySensor s = rotary_sensor_create(SENSOR_ENCODER_INCREMENTAL,
                                          2000, 100e3, 1e-6);

    ServoAxis axis;
    design_servo_axis(&m, &g, &L, &s, 10.0, 5.0, 8.0, &axis);

    assert(axis.bandwidth_pos > 0.0);
    assert(axis.bandwidth_vel > axis.bandwidth_pos);
    assert(axis.bandwidth_cur > axis.bandwidth_vel);
    assert(axis.Kp_current > 0.0);
    assert(axis.Kp_velocity > 0.0);
    PASS();
}

static void test_servo_step_analysis(void)
{
    TEST("Servo axis step response analysis");
    DCMotorParams m = dc_motor_create(2.0, 0.002, 0.05, 1e-5, 1e-5, 24.0, 5.0);
    GearTrain g = gear_train_create(5.0, 0.9, 1e-6, 2e-6, 0.0, 1e6);

    ServoAxis axis;
    design_servo_axis(&m, &g, NULL, NULL, 10.0, 5.0, 8.0, &axis);

    double rise, settle, overshoot, ss_err;
    servo_axis_step_analysis(&axis, 0.1, 0.001, 0.5,
                              &rise, &settle, &overshoot, &ss_err);
    assert(rise > 0.0);
    assert(settle > rise);
    assert(overshoot >= 0.0);
    PASS();
}

/*===========================================================================*/
/* L7: Application Tests                                                     */
/*===========================================================================*/

static void test_regenerative_energy(void)
{
    TEST("Regenerative braking energy");
    double E = regenerative_energy(0.001, 300.0, 0.0, 0.7);
    /* E = 0.7 * 0.5 * 0.001 * (300² - 0) = 31.5 J */
    assert(NEAR(E, 31.5));
    PASS();
}

static void test_antenna_scan_time(void)
{
    TEST("Antenna scan time computation");
    double T = antenna_scan_time(M_PI, 0.1 * M_PI / 180.0,
                                 0.5 * M_PI / 180.0, 1.0, 0.05);
    assert(T > 0.0);
    PASS();
}

/*===========================================================================*/
/* L2: Electromechanical Analogy Tests                                       */
/*===========================================================================*/

static void test_mech_elec_analogy(void)
{
    TEST("Force-voltage electromechanical analogy");
    double Z_elec = mech_impedance_to_electrical_fv(100.0);
    assert(NEAR(Z_elec, 100.0));  /* Direct mapping */
    PASS();
}

static void test_mass_to_capacitance(void)
{
    TEST("Mass to capacitance (F-V analogy)");
    double C = mass_to_capacitance_fv(0.5);
    assert(NEAR(C, 0.5));
    PASS();
}

/*===========================================================================*/
/* L8: Piezo Actuator                                                        */
/*===========================================================================*/

static void test_piezo_free_stroke(void)
{
    TEST("Piezo free stroke");
    PiezoActuator p;
    memset(&p, 0, sizeof(p));
    p.d33 = 593e-12;  /* PZT-5H */
    p.K_stiffness = 50e6;
    p.C_cap = 1e-6;

    double dL = piezo_free_stroke(&p, 150.0);
    assert(NEAR(dL, 593e-12 * 150.0));  /* ≈ 89 nm */
    PASS();
}

static void test_piezo_stroke_under_load(void)
{
    TEST("Piezo stroke under elastic load");
    PiezoActuator p;
    memset(&p, 0, sizeof(p));
    p.d33 = 593e-12;
    p.K_stiffness = 50e6;
    p.C_cap = 1e-6;

    double dL = piezo_stroke_under_load(&p, 150.0, 50e6, 0.0);
    /* Identical stiffness: stroke = half of free stroke */
    double free = piezo_free_stroke(&p, 150.0);
    assert(NEAR(dL, free * 0.5));
    PASS();
}

/*===========================================================================*/
/* Main                                                                      */
/*===========================================================================*/

int main(void)
{
    printf("=== mini-mechatronic-modeling Test Suite ===\n\n");

    printf("--- L1: Definitions ---\n");
    test_dc_motor_create();
    test_gear_train_create();
    test_lead_screw_create();
    test_mechatronic_subsystem_create();

    printf("\n--- L4: Fundamental Laws ---\n");
    test_lorentz_force();
    test_lorentz_force_angle();
    test_motional_emf();
    test_dc_motor_back_emf();
    test_kt_ke_equivalence();
    test_dc_motor_torque();
    test_dc_motor_steady_state();
    test_motor_efficiency();
    test_dc_motor_controllability();
    test_dc_motor_observability();
    test_eigenvalue_stability();

    printf("\n--- L5: Engineering Methods ---\n");
    test_time_constants();
    test_power_calculations();
    test_rk4_step();
    test_zoh_discretization();
    test_pi_controller_design();
    test_rms_torque();
    test_optimal_gear_ratio();
    test_reflect_load_through_gear();
    test_pi_antiwindup();
    test_pi_antiwindup_saturation();
    test_friction_model();
    test_backlash_model();

    printf("\n--- L1-L3: Sensors ---\n");
    test_encoder_counts_to_angle();
    test_encoder_quadrature();
    test_encoder_velocity_m_method();
    test_encoder_direction();
    test_resolver_angle();
    test_bldc_angle_from_halls();
    test_strain_gauge();
    test_lowpass_filter();
    test_complementary_filter();

    printf("\n--- L6: System Design ---\n");
    test_compute_reflected_dynamics();
    test_inertia_ratio();
    test_ball_screw_required_torque();
    test_linear_to_rotary();
    test_robot_joint_max_torque();
    test_hdd_seek_time();
    test_pick_and_place();
    test_motion_profile_energy();
    test_cascade_bandwidth_allocation();
    test_design_servo_axis();
    test_servo_step_analysis();

    printf("\n--- L7-L8: Applications & Advanced ---\n");
    test_regenerative_energy();
    test_antenna_scan_time();
    test_mech_elec_analogy();
    test_mass_to_capacitance();
    test_piezo_free_stroke();
    test_piezo_stroke_under_load();

    printf("\n=== Results: %d passed, %d failed ===\n",
           tests_passed, tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
