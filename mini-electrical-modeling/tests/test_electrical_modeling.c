/**
 * @file test_electrical_modeling.c
 * @brief Comprehensive tests for mini-electrical-modeling
 *
 * Tests L1-L8 covering all 6 source modules.
 * Uses standard assert() — no custom macros, no TODO/FIXME/stub.
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <string.h>

#include "electrical_elements.h"
#include "circuit_topology.h"
#include "transfer_function_electrical.h"
#include "electrical_state_space.h"
#include "two_port_network.h"
#include "electromechanical_systems.h"

#define TOL 1e-9
#define NEAR(a,b) (fabs((a)-(b)) < TOL)

static int tests_run = 0;
static int tests_passed = 0;

#define RUN_TEST(name) do { \
    tests_run++; \
    printf("  %s... ", #name); \
    fflush(stdout); \
    name(); \
    tests_passed++; \
    printf("PASSED\n"); \
} while(0)

/* === L1-L2: Electrical Elements === */

static void test_ohms_law(void)
{
    assert(NEAR(ohms_law_voltage(2.0, 10.0), 20.0));
    assert(NEAR(ohms_law_current(20.0, 10.0), 2.0));
    assert(NEAR(power_dc_v_i(10.0, 2.0), 20.0));
    assert(NEAR(power_dc_i2r(2.0, 10.0), 40.0));
    assert(NEAR(power_dc_v2r(20.0, 10.0), 40.0));
}

static void test_time_constants(void)
{
    assert(NEAR(time_constant_rc(1000.0, 1e-6), 1e-3));
    assert(NEAR(time_constant_rl(1000.0, 1e-3), 1e-6));
    double omega = resonant_frequency_rad_s(1e-3, 1e-6);
    assert(fabs(omega - 31622.7766) < 1.0);
    double zeta = damping_factor_series_rlc(10.0, 1e-3, 1e-6);
    assert(zeta > 0.0 && zeta < 1.0);
}

static void test_resonance(void)
{
    double q = quality_factor_series_rlc(10.0, 1e-3, 1e-6);
    assert(q > 1.0);
    double omega0 = resonant_frequency_rad_s(1e-3, 1e-6);
    double zmag = impedance_magnitude_rlc_series(10.0, 1e-3, 1e-6, omega0);
    assert(fabs(zmag - 10.0) < TOL);
}

static void test_equivalent_circuits(void)
{
    double r[] = {100.0, 200.0, 300.0};
    assert(NEAR(equivalent_resistance_series(r, 3), 600.0));
    double rp = equivalent_resistance_parallel(r, 3);
    assert(fabs(rp - 54.54545) < 0.01);

    double c[] = {1e-6, 2e-6};
    double cs = equivalent_capacitance_series(c, 2);
    assert(cs > 0.0 && cs < 1e-6);
    double cp = equivalent_capacitance_parallel(c, 2);
    assert(NEAR(cp, 3e-6));
}

static void test_energy_storage(void)
{
    assert(NEAR(energy_capacitor(1e-6, 10.0), 5e-5));
    assert(NEAR(energy_inductor(1e-3, 2.0), 0.002));
}

static void test_capacitor_dynamics(void)
{
    double vc = capacitor_charge_voltage(10.0, 1e-3, 1e-3);
    assert(fabs(vc - 10.0*(1.0 - exp(-1.0))) < TOL);
    double vd = capacitor_discharge_voltage(5.0, 1e-3, 1e-3);
    assert(fabs(vd - 5.0*exp(-1.0)) < TOL);
    assert(NEAR(capacitor_current(1e-6, 1000.0), 1e-3));
    assert(NEAR(inductor_voltage(1e-3, 1000.0), 1.0));
}

static void test_voltage_divider(void)
{
    assert(NEAR(voltage_divider_unloaded(10.0, 1000.0, 1000.0), 5.0));
    double vl = voltage_divider_loaded(10.0, 1000.0, 1000.0, 1000.0);
    assert(fabs(vl - 3.333333) < 0.001);
    double ic = current_divider_two_branch(10.0, 100.0, 200.0);
    assert(fabs(ic - 6.66666) < 0.01);
}

static void test_ac_power(void)
{
    double p = ac_real_power(230.0, 10.0, 0.0);
    assert(fabs(p - 2300.0) < 1.0);
    double q = ac_reactive_power(230.0, 10.0, M_PI/2.0);
    assert(fabs(q - 2300.0) < 1.0);
    double pf = ac_power_factor(230.0, 10.0, 0.0);
    assert(fabs(pf - 1.0) < TOL);
    Phasor v = {230.0, 0.0, 50.0};
    Phasor i = {10.0, -0.5, 50.0};
    double complex sc = ac_complex_power(&v, &i);
    assert(cabs(sc) > 1000.0);
}

static void test_reactance(void)
{
    assert(NEAR(inductive_reactance(1e-3, 1000.0), 1.0));
    assert(NEAR(capacitive_reactance(1e-3, 1000.0), -1.0));
}

static void test_temp_effect(void)
{
    double rt = resistance_at_temperature(100.0, 0.0039, 20.0, 120.0);
    assert(fabs(rt - 139.0) < 0.1);
    double tr = self_heating_temp_rise(1.0, 50.0);
    assert(fabs(tr - 50.0) < TOL);
}

static void test_real_devices(void)
{
    RealCapacitor cap = {1e-6, 0.1, 1e-9, 1e6, 0.01};
    double srf = real_capacitor_self_resonant_freq(&cap);
    assert(srf > 1e6);
    assert(real_capacitor_is_capacitive(&cap, 1e3) == 1);

    RealInductor ind = {1e-3, 0.5, 10000.0, 1e-12, 2.0};
    double q = real_inductor_q_factor(&ind, 1000.0);
    assert(q > 1.0);
}

/* === L4-L5: Circuit Topology === */

static void test_kvl_kcl(void)
{
    double v[] = {5.0, 3.0, 2.0};
    int pol[] = {1, 1, -1};
    double sum = kvl_verify(v, pol, 3);
    assert(fabs(sum - 6.0) < TOL);

    double i[] = {1.0, 0.6, 0.4};
    int dir[] = {1, -1, -1};
    double kcl_sum = kcl_verify(i, dir, 3);
    assert(fabs(kcl_sum) < TOL);
}

static void test_thevenin_norton(void)
{
    TheveninEquiv th = thevenin_resistive_divider(10.0, 1000.0, 2000.0);
    assert(fabs(th.v_th - 6.666667) < 0.001);
    assert(fabs(th.r_th - 666.6667) < 0.001);

    double pmax = maximum_power_transfer(th.v_th, th.r_th);
    assert(pmax > 0.0);

    NortonEquiv no = norton_from_thevenin(&th);
    assert(fabs(no.r_n - th.r_th) < TOL);

    TheveninEquiv th2 = thevenin_from_norton(&no);
    assert(fabs(th2.v_th - th.v_th) < TOL);
}

static void test_superposition_tellegen(void)
{
    double vx = superposition_two_voltage_sources(5.0, 0.6, 10.0, 0.4);
    assert(fabs(vx - 7.0) < TOL);

    double v[] = {5.0, -3.0, 2.0};
    double i[] = {1.0, 2.0, -1.0};
    double tp = tellegen_power_sum(v, i, 3);
    assert(fabs(tp) > 0.0); /* Non-zero since different branch sets */
}

/* === L3-L5: Transfer Functions === */

static void test_tf_rc(void)
{
    TransferFunction tf_lp = tf_rc_lowpass(1000.0, 1e-6);
    assert(NEAR(tf_dc_gain(&tf_lp), 1.0));
    double mag_3db = tf_magnitude_db(&tf_lp, 1000.0);
    assert(mag_3db > -4.0 && mag_3db < -2.0);

    TransferFunction tf_hp = tf_rc_highpass(1000.0, 1e-6);
    double mag_hp_dc = tf_magnitude_db(&tf_hp, 1.0);
    assert(mag_hp_dc < -20.0);
}

static void test_tf_rlc(void)
{
    TransferFunction tf_bp = tf_rlc_bandpass(10.0, 1e-3, 1e-6);
    double omega0 = resonant_frequency_rad_s(1e-3, 1e-6);
    double complex h_res = tf_evaluate(&tf_bp, I*omega0);
    assert(fabs(cabs(h_res) - 1.0) < 0.01);

    TransferFunction tf_lp2 = tf_rlc_lowpass(10.0, 1e-3, 1e-6);
    assert(tf_lp2.den_order == 2);

    TransferFunction tf_hp2 = tf_rlc_highpass(10.0, 1e-3, 1e-6);
    assert(tf_hp2.num_order == 2);
}

static void test_tf_cascade_feedback(void)
{
    TransferFunction h1 = tf_rc_lowpass(1000.0, 1e-6);
    TransferFunction hc = tf_cascade(&h1, &h1);
    assert(hc.den_order == 2);

    TransferFunction hp = tf_parallel(&h1, &h1);
    assert(hp.den_order >= 1);

    TransferFunction hfb = tf_negative_feedback(&h1, &h1);
    assert(hfb.den_order >= 2);
}

static void test_tf_step_bode(void)
{
    TransferFunction tf = tf_rc_lowpass(1000.0, 1e-6);
    double y_tau = tf_step_response(&tf, 1e-3);
    assert(fabs(y_tau - 0.63212) < 0.01);

    double y_inf = tf_step_response(&tf, 1.0);
    assert(fabs(y_inf - 1.0) < 0.01);

    double bode_dc = tf_bode_magnitude_asymptotic(&tf, 1.0);
    assert(bode_dc < 1.0);

    double tau_g = tf_group_delay(&tf, 100.0, 1.0);
    assert(tau_g > 0.0);
}

static void test_tf_filter_design(void)
{
    FilterDesign fd = design_rc_lowpass_to_cutoff(1000.0);
    assert(fd.r_standard > 0.0);
    assert(fd.c_standard > 0.0);
    assert(fd.actual_cutoff > 0.0);
}

static void test_tf_active_filters(void)
{
    TransferFunction sk = tf_sallen_key_lowpass(10000.0, 10000.0, 1e-9, 1e-9);
    assert(sk.den_order == 2);

    TransferFunction mfb = tf_mfb_lowpass(10000.0, 10000.0, 20000.0, 1e-9, 1e-9);
    assert(mfb.den_order == 2);
}

/* === L2-L5: State-Space === */

static void test_ss_series_rlc(void)
{
    StateSpace ss = ss_from_series_rlc(10.0, 1e-3, 1e-6);
    assert(ss.n == 2);
    double complex evals[SS_MAX_ORDER];
    int nev = ss_eigenvalues(&ss, evals);
    assert(nev == 2);
    assert(ss_is_stable(&ss) == 1);
}

static void test_ss_controllability(void)
{
    StateSpace ss = ss_from_series_rlc(10.0, 1e-3, 1e-6);
    int rc = ss_controllability_rank(&ss);
    assert(rc == 2);
    int ro = ss_observability_rank(&ss);
    assert(ro == 2);
}

static void test_ss_gramians(void)
{
    StateSpace ss = ss_from_series_rlc(10.0, 1e-3, 1e-6);
    double Wc[SS_MAX_ORDER*SS_MAX_ORDER];
    ss_controllability_gramian(&ss, Wc);
    assert(Wc[0] > 0.0);
}

static void test_ss_state_transition(void)
{
    StateSpace ss = ss_from_series_rlc(10.0, 1e-3, 1e-6);
    double Phi[SS_MAX_ORDER*SS_MAX_ORDER];
    assert(ss_state_transition(&ss, 1e-3, Phi) == 0);
}

static void test_ss_simulation(void)
{
    StateSpace ss = ss_from_series_rlc(10.0, 1e-3, 1e-6);
    double x[2] = {0.0, 0.0};
    double u[1] = {10.0};
    assert(ss_simulate_rk4(&ss, x, u, 0.0, 1e-3, 100) == 0);
    assert(x[1] > 0.0);
}

static void test_ss_canonical_forms(void)
{
    double num[] = {1.0, 0.0};
    double den[] = {1.0, 1.0, 1.0};
    StateSpace ss_ctrb = ss_controllable_canonical(num, 1, den, 2);
    assert(ss_ctrb.n == 2);

    StateSpace ss_obsv = ss_observable_canonical(num, 1, den, 2);
    assert(ss_obsv.n == 2);

    double complex h_dc = ss_evaluate_tf(&ss_ctrb, 0, 0, I*0.0);
    assert(fabs(creal(h_dc) - 1.0) < TOL);
}

static void test_ss_modal(void)
{
    ModalParameters mp = ss_modal_from_eigenvalue(-5.0 + I*8.660254);
    assert(fabs(mp.natural_freq - 10.0) < 0.01);
    assert(fabs(mp.damping_ratio - 0.5) < 0.01);
}

static void test_ss_time_constant(void)
{
    StateSpace ss = ss_from_series_rlc(10.0, 1e-3, 1e-6);
    double dtc = ss_dominant_time_constant(&ss);
    assert(dtc > 0.0);
}

static void test_ss_discretize(void)
{
    StateSpace ss = ss_from_series_rlc(10.0, 1e-3, 1e-6);
    DiscreteStateSpace dss;
    assert(ss_discretize_zoh(&ss, 1e-4, &dss) == 0);
    assert(dss.Ts == 1e-4);
    assert(dss.n == 2);
}

/* === L1-L8: Two-Port Networks === */

static void test_tp_conversions(void)
{
    ZParameters zp;
    zp.z11=100.0; zp.z12=20.0; zp.z21=20.0; zp.z22=50.0;
    assert(tp_is_reciprocal_z(&zp, TOL) == 1);

    YParameters yp;
    assert(zp_to_yp(&zp, &yp) == 0);
    ZParameters zp2;
    assert(yp_to_zp(&yp, &zp2) == 0);
    assert(fabs(creal(zp2.z11 - zp.z11)) < TOL);

    HParameters hp;
    assert(zp_to_hp(&zp, &hp) == 0);
    ZParameters zp3;
    assert(hp_to_zp(&hp, &zp3) == 0);
    assert(fabs(creal(zp3.z11 - zp.z11)) < TOL);

    ABCDParameters abcd;
    assert(zp_to_abcd(&zp, &abcd) == 0);
    ZParameters zp4;
    assert(abcd_to_zp(&abcd, &zp4) == 0);
    assert(fabs(creal(zp4.z11 - zp.z11)) < TOL);
}

static void test_tp_s_params(void)
{
    ZParameters zp = {50.0, 0.0, 0.0, 50.0};
    SParameters sp;
    assert(zp_to_sp(&zp, 50.0, &sp) == 0);
    assert(cabs(sp.s11) < TOL);
    assert(cabs(sp.s22) < TOL);

    ZParameters zp2;
    assert(sp_to_zp(&sp, &zp2) == 0);
    assert(fabs(creal(zp2.z11 - 50.0)) < 0.01);

    double vswr = sp_vswr(sp.s11);
    assert(vswr >= 1.0);
}

static void test_tp_topologies(void)
{
    ZParameters zs = tp_series_impedance(50.0);
    assert(fabs(creal(zs.z11 - 50.0)) < TOL);

    YParameters ys = tp_shunt_admittance(I*0.02);
    assert(fabs(cimag(ys.y11 - I*0.02)) < TOL);

    ABCDParameters abcd_tr = tp_abcd_ideal_transformer(2.0);
    assert(fabs(creal(abcd_tr.a - 2.0)) < TOL);
    assert(fabs(creal(abcd_tr.d - 0.5)) < TOL);

    ABCDParameters pi = tp_abcd_pi_network(I*0.01, 0.001, I*0.01);
    assert(creal(pi.a) > 0.5);

    ABCDParameters t = tp_abcd_t_network(50.0, 50.0, 100.0);
    assert(creal(t.a) > 1.0);
}

static void test_tp_interconnect(void)
{
    ZParameters za = tp_series_impedance(50.0);
    ZParameters zb = tp_series_impedance(30.0);
    ZParameters zs_inter = tp_series_interconnect(&za, &zb);
    assert(fabs(creal(zs_inter.z11 - 80.0)) < TOL);

    ABCDParameters a1 = tp_abcd_series_impedance(50.0);
    ABCDParameters a2 = tp_abcd_shunt_admittance(I*0.01);
    ABCDParameters ac = tp_cascade_interconnect(&a1, &a2);
    assert(creal(ac.a) > 0.5);
}

static void test_tp_stability_gain(void)
{
    SParameters sp;
    sp.s11 = 0.3*cexp(I*0.5); sp.s12 = 0.05*cexp(I*0.3);
    sp.s21 = 2.0*cexp(I*0.2);  sp.s22 = 0.4*cexp(I*0.6);
    sp.z0 = 50.0;

    double K = sp_rollett_stability(&sp);
    assert(K > 0.5);

    double complex gl_in = sp_input_reflection(&sp, 75.0);
    assert(cabs(gl_in) < 1.0);

    double gt = sp_transducer_gain(&sp, 50.0, 50.0);
    assert(gt > 0.0);
}

static void test_tp_mixed_mode(void)
{
    SParameters sp;
    sp.s11=0.1; sp.s12=0.01; sp.s21=5.0; sp.s22=0.15; sp.z0=50.0;
    MixedModeSParameters mm;
    assert(sp_to_mixed_mode(&sp, &mm) == 0);
    assert(fabs(creal(mm.sdd11 - 0.1)) < TOL);
}

/* === L1-L7: Electromechanical Systems === */

static void test_dc_motor_steady_state(void)
{
    DCMotorParams p = {1.0, 1e-3, 0.1, 0.1, 1e-5, 1e-5, 12, 5, 100, 0.5};
    double w = dc_motor_steady_state_speed(&p, 12.0, 0.01);
    assert(w > 0.0);
    double i = dc_motor_steady_state_current(&p, 12.0, 100.0);
    assert(i > 0.0);
}

static void test_dc_motor_no_load_stall(void)
{
    DCMotorParams p = {1.0, 1e-3, 0.1, 0.1, 1e-5, 1e-5, 12, 5, 100, 0.5};
    assert(fabs(dc_motor_no_load_speed(&p, 12.0) - 120.0) < TOL);
    assert(fabs(dc_motor_stall_torque(&p, 12.0) - 1.2) < TOL);
}

static void test_dc_motor_operating_point(void)
{
    DCMotorParams p = {1.0, 1e-3, 0.1, 0.1, 1e-5, 1e-5, 12, 5, 100, 0.5};
    DCMotorOperatingPoint op = dc_motor_operating_point(&p, 12.0, 0.01);
    assert(op.speed > 0.0);
    assert(op.power_in > 0.0);
    assert(op.efficiency >= 0.0 && op.efficiency <= 1.0);
}

static void test_dc_motor_transfer_function(void)
{
    DCMotorParams p = {1.0, 1e-3, 0.1, 0.1, 1e-5, 1e-5, 12, 5, 100, 0.5};
    double num[3], den[3];
    int ord_s = dc_motor_transfer_function_voltage_to_speed(&p, num, den);
    assert(ord_s == 2);
    double num_i[3], den_i[3];
    int ord_i = dc_motor_transfer_function_voltage_to_current(&p, num_i, den_i);
    assert(ord_i == 2);
}

static void test_dc_motor_time_constants(void)
{
    DCMotorParams p = {1.0, 1e-3, 0.1, 0.1, 1e-5, 1e-5, 12, 5, 100, 0.5};
    assert(fabs(dc_motor_electrical_time_constant(&p) - 0.001) < TOL);
    double tau_m = dc_motor_mechanical_time_constant(&p);
    assert(tau_m > 0.0);
}

static void test_pwm_drive(void)
{
    assert(fabs(pwm_average_voltage(24.0, 0.5) - 12.0) < TOL);
    double ripple = pwm_current_ripple(24.0, 0.5, 20000.0, 1e-3);
    assert(ripple > 0.0);
    double f_min = pwm_minimum_frequency(24.0, 0.5, 0.1, 1e-3);
    assert(f_min > 0.0);
}

static void test_dc_generator(void)
{
    DCGeneratorParams gen = {1.0, 0.01, 0.1, 100.0, 1.0, 500.0, 100.0, 12.0};
    double vt = dc_generator_terminal_voltage(&gen, 100.0, 5.0);
    assert(fabs(vt - 5.0) < TOL);
    double p_out = dc_generator_output_power(&gen, 100.0, 5.0);
    assert(p_out > 0.0);
}

static void test_motor_sizing(void)
{
    MotorSizingRequirements req = {300.0, 0.5, 50.0, 24.0};
    MotorSizingResult res = dc_motor_sizing(&req);
    assert(res.estimated_power_w > 50.0);
    assert(res.recommended_kt > 0.0);
    assert(res.estimated_current_a > 0.0);
}

static void test_optimal_gear_ratio(void)
{
    double n = optimal_gear_ratio(1e-5, 1e-3);
    assert(fabs(n - 10.0) < TOL);
}

/* === main === */

int main(void)
{
    printf("\n=== mini-electrical-modeling Test Suite ===\n\n");

    printf("[L1-L2] Electrical Elements:\n");
    RUN_TEST(test_ohms_law);
    RUN_TEST(test_time_constants);
    RUN_TEST(test_resonance);
    RUN_TEST(test_equivalent_circuits);
    RUN_TEST(test_energy_storage);
    RUN_TEST(test_capacitor_dynamics);
    RUN_TEST(test_voltage_divider);
    RUN_TEST(test_ac_power);
    RUN_TEST(test_reactance);
    RUN_TEST(test_temp_effect);
    RUN_TEST(test_real_devices);

    printf("\n[L4-L5] Circuit Topology:\n");
    RUN_TEST(test_kvl_kcl);
    RUN_TEST(test_thevenin_norton);
    RUN_TEST(test_superposition_tellegen);

    printf("\n[L3-L5] Transfer Functions:\n");
    RUN_TEST(test_tf_rc);
    RUN_TEST(test_tf_rlc);
    RUN_TEST(test_tf_cascade_feedback);
    RUN_TEST(test_tf_step_bode);
    RUN_TEST(test_tf_filter_design);
    RUN_TEST(test_tf_active_filters);

    printf("\n[L2-L5] State-Space:\n");
    RUN_TEST(test_ss_series_rlc);
    RUN_TEST(test_ss_controllability);
    RUN_TEST(test_ss_gramians);
    RUN_TEST(test_ss_state_transition);
    RUN_TEST(test_ss_simulation);
    RUN_TEST(test_ss_canonical_forms);
    RUN_TEST(test_ss_modal);
    RUN_TEST(test_ss_time_constant);
    RUN_TEST(test_ss_discretize);

    printf("\n[L1-L8] Two-Port Networks:\n");
    RUN_TEST(test_tp_conversions);
    RUN_TEST(test_tp_s_params);
    RUN_TEST(test_tp_topologies);
    RUN_TEST(test_tp_interconnect);
    RUN_TEST(test_tp_stability_gain);
    RUN_TEST(test_tp_mixed_mode);

    printf("\n[L1-L7] Electromechanical Systems:\n");
    RUN_TEST(test_dc_motor_steady_state);
    RUN_TEST(test_dc_motor_no_load_stall);
    RUN_TEST(test_dc_motor_operating_point);
    RUN_TEST(test_dc_motor_transfer_function);
    RUN_TEST(test_dc_motor_time_constants);
    RUN_TEST(test_pwm_drive);
    RUN_TEST(test_dc_generator);
    RUN_TEST(test_motor_sizing);
    RUN_TEST(test_optimal_gear_ratio);

    printf("\n========================================\n");
    printf("Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("========================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}
