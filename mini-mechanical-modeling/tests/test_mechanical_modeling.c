/**
 * @file test_mechanical_modeling.c
 * @brief Comprehensive tests for mini-mechanical-modeling
 *
 * Tests L1-L8 covering all 6 source modules.
 * Uses standard assert() --- no custom macros, no TODO/FIXME/stub.
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <string.h>

#include "mechanical_elements.h"
#include "mechanical_transfer_function.h"
#include "mechanical_state_space.h"
#include "rotational_systems.h"
#include "multi_dof_systems.h"
#include "mechanical_applications.h"

#define TOL 1e-9
#define NEAR(a,b) (fabs((a)-(b)) < TOL)

static int tests_run = 0;
static int tests_passed = 0;

#define RUN_TEST(name) do {     tests_run++;     printf("  %s... ", #name);     fflush(stdout);     name();     tests_passed++;     printf("PASSED\n"); } while(0)

/* === L1-L2: Mechanical Elements === */

static void test_newton_second_law(void)
{
    assert(NEAR(newton_second_law_force(2.0, 10.0), 20.0));
    assert(NEAR(newton_second_law_acceleration(20.0, 2.0), 10.0));
    assert(newton_second_law_acceleration(10.0, 0.0) > 1e100);
    assert(newton_second_law_acceleration(-10.0, 0.0) < -1e100);
}

static void test_hookes_law(void)
{
    assert(NEAR(hookes_law_force(100.0, 0.1), -10.0));
    assert(NEAR(hookes_law_displacement(100.0, 10.0), 0.1));
    assert(NEAR(spring_potential_energy(100.0, 0.1), 0.5));
    double k[] = {100.0, 200.0};
    assert(fabs(spring_series_equivalent(k, 2) - 66.666667) < 1e-5); /* 1/(1/100+1/200) */
    assert(NEAR(spring_parallel_equivalent(k, 2), 300.0));
}

static void test_friction(void)
{
    assert(NEAR(coulomb_friction_magnitude(0.3, 10.0), 3.0));
    double Fc = coulomb_friction_force(0.3, 10.0, 5.0);
    assert(Fc < -0.01); /* opposes motion */
    double Fc_neg = coulomb_friction_force(0.3, 10.0, -5.0);
    assert(Fc_neg > 0.01);
    assert(NEAR(coulomb_friction_force(0.3, 10.0, 0.0), 0.0));
}

static void test_stribeck_friction(void)
{
    FrictionModel fm = {FRICTION_STRIBECK, 0.5, 0.3, 0.01, 0.1, 10.0};
    double f1 = stribeck_friction_force(&fm, 2.0);
    assert(f1 < 0.0); /* opposes motion */
    double f_zero = stribeck_friction_force(&fm, 0.0);
    assert(NEAR(f_zero, 0.0));
    assert(static_friction_holds(&fm, 2.0) == 1);  /* 2 < 0.5*10=5 */
    assert(static_friction_holds(&fm, 6.0) == 0);  /* 6 > 5 */
}

static void test_energy(void)
{
    assert(NEAR(kinetic_energy_translational(2.0, 3.0), 9.0)); /* 0.5*2*9=9 */
    assert(NEAR(gravitational_potential_energy(2.0, 5.0, 9.81), 98.1));
    double E = total_mechanical_energy(2.0, 3.0, 100.0, 0.1, 5.0, 9.81);
    assert(E > 0.0);
    double Wd = energy_dissipated_per_cycle(10.0, 20.0, 0.05);
    assert(Wd > 0.0);
}

static void test_natural_frequency(void)
{
    double wn = natural_frequency_undamped(1000.0, 10.0);
    assert(NEAR(wn, 10.0)); /* sqrt(1000/10)=10 */
    double fn = natural_frequency_hz(1000.0, 10.0);
    assert(NEAR(fn, 10.0/(2.0*M_PI)));
}

static void test_damping_params(void)
{
    /* m=1, k=100, omega_n=10, zeta = b/(2*10) = b/20 */
    double zeta = damping_ratio(5.0, 100.0, 1.0);
    assert(NEAR(zeta, 0.25));
    double b_cr = critical_damping_coefficient(100.0, 1.0);
    assert(NEAR(b_cr, 20.0));
    double wd = damped_natural_frequency(100.0, 1.0, 5.0);
    assert(wd > 0.0 && wd < 10.0);
}

static void test_log_decrement(void)
{
    /* x0=1.0, x10=0.535 (approx), n=10 => delta=0.0625 => zeta~=0.01 */
    double delta = logarithmic_decrement(1.0, 0.535, 10);
    assert(delta > 0.05 && delta < 0.07);
    double zeta = damping_ratio_from_logdec(delta);
    assert(zeta > 0.005 && zeta < 0.015);
}

static void test_q_factor(void)
{
    double Q = quality_factor_from_damping(0.05);
    assert(NEAR(Q, 10.0)); /* Q=1/(2*0.05)=10 */
    double z = damping_ratio_from_q_factor(10.0);
    assert(NEAR(z, 0.05));
    assert(NEAR(loss_factor_from_damping(0.05), 0.1));
}

static void test_static_deflection(void)
{
    double ds = static_deflection(10.0, 1000.0, 9.81);
    assert(NEAR(ds, 0.0981));
    double k = stiffness_from_static_deflection(10.0, 0.0981, 9.81);
    assert(NEAR(k, 1000.0));
    double wn = natural_freq_from_static_deflection(0.0981, 9.81);
    assert(NEAR(wn, 10.0));
}

static void test_mechanical_impedance(void)
{
    double complex Z = mechanical_impedance_sdof(1.0, 5.0, 100.0, 10.0);
    /* At resonance w=10: Z = b + j*(10*1 - 100/10) = 5 + j*(10-10) = 5 */
    assert(NEAR(creal(Z), 5.0));
    assert(fabs(cimag(Z)) < TOL);
    double mag = mechanical_impedance_magnitude(1.0, 5.0, 100.0, 10.0);
    assert(NEAR(mag, 5.0));
    double complex Y = mechanical_mobility(1.0, 5.0, 100.0, 10.0);
    assert(NEAR(creal(Y), 1.0/5.0));
}

static void test_transmissibility(void)
{
    /* r=2, zeta=0.1: TR = sqrt((1+0.16)/(9+0.16)) = sqrt(1.16/9.16) ~= 0.356 */
    double TR = force_transmissibility(2.0, 0.1);
    assert(TR < 1.0 && TR > 0.3);
    double TRd = displacement_transmissibility(2.0, 0.1);
    assert(NEAR(TR, TRd));
    double eff = isolation_efficiency(2.0, 0.1);
    assert(eff > 60.0);
}

static void test_moments_of_inertia(void)
{
    assert(NEAR(inertia_point_mass(2.0, 3.0), 18.0));  /* 2*9 */
    assert(NEAR(inertia_rod_end(6.0, 2.0), 8.0));       /* (1/3)*6*4 */
    assert(NEAR(inertia_rod_center(6.0, 2.0), 2.0));    /* (1/12)*6*4 */
    assert(NEAR(inertia_solid_cylinder(10.0, 0.5), 1.25)); /* 0.5*10*0.25 */
    assert(NEAR(inertia_solid_sphere(10.0, 0.5), 1.0));   /* 0.4*10*0.25 */
    /* Parallel axis: I_cm=2, m=6, d=2 => I=2+6*4=26 */
    assert(NEAR(parallel_axis_theorem(2.0, 6.0, 2.0), 26.0));
    assert(NEAR(radius_of_gyration(18.0, 2.0), 3.0));
}

static void test_time_domain_response(void)
{
    double zeta = 0.3, wn = 10.0;
    double tp = peak_time_underdamped(zeta, wn);
    assert(tp > 0.3 && tp < 0.35);
    double OS = percent_overshoot(zeta);
    assert(OS > 35.0 && OS < 40.0);
    double ts2 = settling_time_2percent(zeta, wn);
    assert(ts2 > 1.3 && ts2 < 1.35);
    double tr = rise_time_underdamped(zeta, wn);
    assert(tr > 0.1 && tr < 0.2);
}

static void test_step_response(void)
{
    /* At t large, step response ~= static deflection */
    double x = sdof_step_response_displacement(0.1, 10.0, 0.3, 5.0);
    assert(x > 0.09 && x < 0.12);
    double x0 = sdof_step_response_displacement(0.1, 10.0, 0.3, 0.0);
    assert(NEAR(x0, 0.0));
}

static void test_resonance(void)
{
    double Xmax = resonance_amplitude_force(1.0, 100.0, 0.05);
    assert(Xmax > 0.09 && Xmax < 0.11); /* ~ 1/(2*k*zeta) = 1/(2*100*0.05)=0.1 */
    double phi = phase_angle_sdof(1.0, 0.1);
    assert(fabs(phi - M_PI/2.0) < TOL);
    double bw = half_power_bandwidth(0.05, 100.0);
    assert(NEAR(bw, 10.0));
}

static void test_structural_stiffness(void)
{
    /* Circular beam: d=0.02, I=pi*d^4/64, E=200e9, L=1.0, k=3*E*I/L^3 */
    double I_circ = area_moment_circular(0.02);
    assert(I_circ > 7e-9 && I_circ < 8e-9);
    double J_p = polar_moment_circular(0.02);
    assert(J_p > 1.5e-8 && J_p < 1.6e-8);
    double k_cant = cantilever_beam_stiffness(200e9, I_circ,1.0);
    assert(k_cant > 4000.0 && k_cant < 5000.0);
    double k_ss = simply_supported_beam_stiffness(200e9, I_circ,1.0);
    assert(k_ss > 70000.0);
    double I_rect = area_moment_rectangular(0.05, 0.02);
    assert(I_rect > 3.2e-8 && I_rect < 3.5e-8);
}

/* === L3-L5: Transfer Functions === */

static void test_tf_sdof(void)
{
    MechanicalTF tf_x = tf_sdof_displacement(1.0, 5.0, 100.0);
    assert(tf_is_proper(&tf_x));
    double dc = tf_dc_gain(&tf_x);
    assert(NEAR(dc, 0.01)); /* 1/k = 1/100 */

    MechanicalTF tf_v = tf_sdof_velocity(1.0, 5.0, 100.0);
    assert(tf_v.is_strictly_proper);

    /* Poles for m=1,b=5,k=100: s = -2.5 +/- j*9.68 */
    double complex p1, p2;
    tf_sdof_poles(1.0, 5.0, 100.0, &p1, &p2);
    assert(fabs(creal(p1) + 2.5) < 0.1);
    assert(tf_sdof_is_underdamped(1.0, 5.0, 100.0) == 1);
    assert(tf_sdof_is_overdamped(1.0, 30.0, 100.0) == 1);
    assert(tf_sdof_is_critically_damped(1.0, 20.0, 100.0) == 1);
}

static void test_tf_base_excitation(void)
{
    MechanicalTF tf = tf_base_excitation_displacement(1.0, 5.0, 100.0);
    assert(tf_is_proper(&tf));
}

static void test_tf_operations(void)
{
    MechanicalTF h1 = tf_sdof_displacement(1.0, 5.0, 100.0);
    MechanicalTF h2 = tf_sdof_displacement(1.0, 5.0, 100.0);

    MechanicalTF hc = tf_cascade(&h1, &h2);
    assert(hc.den_order == 4);

    MechanicalTF hp = tf_parallel(&h1, &h2);
    assert(hp.den_order >= 2);

    MechanicalTF hfb = tf_neg_feedback(&h1, &h2);
    assert(hfb.den_order >= 2);
}

static void test_tf_evaluation(void)
{
    MechanicalTF tf = tf_sdof_displacement(1.0, 5.0, 100.0);
    double complex h0 = tf_eval(&tf, 0.0);
    assert(NEAR(creal(h0), 0.01));

    double mag_db = tf_mag_db(&tf, 10.0);
    assert(mag_db > -40.0); /* resonant frequency ~10 rad/s */

    double mag_dc = tf_mag(&tf, 0.01);
    assert(fabs(mag_dc - 0.01) < 1e-6); /* low-freq magnitude ~ DC gain */
}

static void test_tf_damping_estimate(void)
{
    MechanicalTF tf = tf_sdof_displacement(1.0, 5.0, 100.0);
    double zeta = tf_estimate_damping(&tf);
    assert(NEAR(zeta, 0.25));
    double wn = tf_estimate_natural_freq(&tf);
    assert(NEAR(wn, 10.0));
}

static void test_tf_steady_state(void)
{
    MechanicalTF tf = tf_sdof_displacement(1.0, 5.0, 100.0);
    double amp = tf_steady_state_amp(&tf, 1.0, 10.0);
    assert(amp > 0.0);
    double phase = tf_steady_state_phase(&tf, 10.0);
    assert(phase < 0.0); /* phase lag */
}

static void test_tf_freq_response(void)
{
    MechanicalTF tf = tf_sdof_displacement(1.0, 5.0, 100.0);
    FreqResp fr;
    tf_freq_response(&tf, 0.1, 100.0, 50, &fr);
    assert(fr.n_points == 50);
    assert(fr.mag_db != NULL);
    freq_resp_free(&fr);
}

/* === L2-L5: State-Space === */

static void test_ss_sdof(void)
{
    MechSS ss = ss_sdof(1.0, 5.0, 100.0);
    assert(ss.n == 2);
    assert(NEAR(ss.A[1], 1.0));
    assert(NEAR(ss.A[2], -100.0));
    assert(NEAR(ss.A[3], -5.0));
}

static void test_ss_stability(void)
{
    MechSS ss = ss_sdof(1.0, 5.0, 100.0);
    assert(ss_is_stable(&ss) == 1);
}

static void test_ss_controllability(void)
{
    MechSS ss = ss_sdof(1.0, 5.0, 100.0);
    assert(ss_is_controllable(&ss) == 1);
    assert(ss_is_observable(&ss) == 1);
}

static void test_ss_simulation(void)
{
    MechSS ss = ss_sdof(1.0, 5.0, 100.0);
    double x[2] = {0.0, 0.0};
    double u = 1.0;
    ss_rk4_step(&ss, x, &u, 0.001, x);
    assert(x[1] > 0.0); /* velocity increased */
}

static void test_ss_2dof(void)
{
    MechSS ss = ss_2dof(10.0, 5.0, 500.0, 300.0, 100.0, 5.0, 3.0, 1.0);
    assert(ss.n == 4);
    double complex evals[SS_MAX_ORDER];
    int nev = ss_eigenvalues(&ss, evals);
    assert(nev > 0);
}

static void test_ss_inverted_pendulum(void)
{
    MechSS ss = ss_inverted_pendulum(1.0, 0.1, 0.5, 9.81);
    assert(ss.n == 4);
    assert(ss_is_stable(&ss) == 0); /* unstable */
}

static void test_ss_discretize(void)
{
    MechSS ss = ss_sdof(1.0, 5.0, 100.0);
    MechDSS dss;
    assert(ss_discretize_zoh(&ss, 0.01, &dss) == 0);
    assert(dss.Ts == 0.01);
    assert(dss.n == 2);
}

static void test_ss_modal(void)
{
    MechSS ss = ss_sdof(1.0, 5.0, 100.0);
    double wn, zeta;
    ss_dominant_mode(&ss, &wn, &zeta);
    assert(wn > 5.0 && wn < 15.0);
    assert(zeta > 0.1 && zeta < 0.4);
}

static void test_ss_gramians(void)
{
    MechSS ss = ss_sdof(1.0, 5.0, 100.0);
    double Wc[SS_MAX_ORDER*SS_MAX_ORDER];
    assert(ss_gramian_ctrb(&ss, Wc) == 0);
    assert(Wc[0] > 0.0);
}

/* === L1-L7: Rotational Systems === */

static void test_rotational_newton(void)
{
    assert(NEAR(rotational_torque(0.1, 100.0), 10.0));
    assert(NEAR(rotational_accel(10.0, 0.1), 100.0));
    assert(NEAR(angular_momentum(0.1, 50.0), 5.0));
}

static void test_torsional_spring(void)
{
    assert(NEAR(torsional_spring_torque(1000.0, 0.1), -100.0));
    assert(NEAR(torsional_potential_energy(1000.0, 0.1), 5.0));
}

static void test_gear_ratio(void)
{
    double N = gear_ratio_calc(20, 60); /* N=3 */
    assert(NEAR(N, 3.0));
    assert(NEAR(gear_output_speed(300.0, 3.0), 100.0));
    assert(NEAR(gear_output_torque(10.0, 3.0), 30.0));
    assert(NEAR(reflected_inertia(0.09, 3.0), 0.01));
}

static void test_optimal_gear_ratio(void)
{
    double N_opt = optimal_gear_ratio(0.01, 0.09);
    assert(NEAR(N_opt, 3.0)); /* sqrt(0.09/0.01)=3 */
}

static void test_rack_pinion(void)
{
    double x = rack_pinion_position(0.1, M_PI/2.0);
    assert(fabs(x - 0.15708) < 0.001);
    double v = rack_pinion_velocity(0.1, 10.0);
    assert(NEAR(v, 1.0));
    double F = rack_pinion_force(0.1, 5.0);
    assert(NEAR(F, 50.0));
}

static void test_leadscrew(void)
{
    double x = leadscrew_position(0.01, 2.0*M_PI); /* one rev */
    assert(NEAR(x, 0.01));
    double F = leadscrew_force(0.01, 2.0, 0.9);
    assert(F > 1100.0); /* 2*pi*2*0.9/0.01 ~= 1131 */
}

static void test_rotating_unbalance(void)
{
    double Fc = centrifugal_force(0.01, 0.05, 100.0);
    assert(NEAR(Fc, 5.0)); /* 0.01*0.05*10000=5 */
    double n_cr = critical_speed_rpm(1e6, 10.0);
    assert(n_cr > 2900.0 && n_cr < 3100.0);
}

static void test_flywheel(void)
{
    double E = flywheel_energy(0.1, 100.0);
    assert(NEAR(E, 500.0)); /* 0.5*0.1*10000=500 */
}

/* === L5-L8: Multi-DOF Systems === */

static void test_mdof_2dof(void)
{
    MDOFSystem sys;
    mdof_build_2dof(&sys, 10.0, 5.0, 500.0, 300.0, 100.0, 5.0, 3.0, 1.0);
    assert(sys.n_dof == 2);

    double evals[MDOF_MAX_DOF], evecs[MDOF_MAX_DOF*MDOF_MAX_DOF];
    int n_modes = mdof_eigen_solve(&sys, evals, evecs);
    assert(n_modes == 2);

    double freqs[MDOF_MAX_DOF];
    int nf;
    mdof_natural_freqs(&sys, freqs, &nf);
    assert(nf == 2);
    assert(freqs[0] > 0.0);
    assert(freqs[1] > 0.0);
}

static void test_mdof_chain(void)
{
    double m[] = {10.0, 10.0, 10.0};
    double k[] = {1000.0, 1000.0, 1000.0, 1000.0}; /* 4 springs */
    MDOFSystem sys;
    mdof_build_chain(&sys, m, k, 3);
    assert(sys.n_dof == 3);
    assert(NEAR(sys.K[0], 2000.0));
}

static void test_mdof_modal_params(void)
{
    MDOFSystem sys;
    mdof_build_2dof(&sys, 10.0, 5.0, 500.0, 300.0, 100.0, 5.0, 3.0, 1.0);
    double evals[MDOF_MAX_DOF], evecs[MDOF_MAX_DOF*MDOF_MAX_DOF];
    mdof_eigen_solve(&sys, evals, evecs);

    ModalParams mp;
    mdof_modal_params(&sys, evecs, 0, &mp);
    assert(mp.natural_freq > 0.0);

    mdof_mass_normalize(&sys, evecs);
    assert(mdof_verify_ortho(&sys, evecs, 2, 0.01) == 1); /* M-orthogonal within tolerance */
}

static void test_mdof_orthogonality(void)
{
    MDOFSystem sys;
    mdof_build_2dof(&sys, 10.0, 5.0, 500.0, 300.0, 100.0, 5.0, 3.0, 1.0);
    double evals[MDOF_MAX_DOF], evecs[MDOF_MAX_DOF*MDOF_MAX_DOF];
    mdof_eigen_solve(&sys, evals, evecs);

    double ortho = mdof_mass_ortho(&sys, evecs, &evecs[2]);
    assert(fabs(ortho) < 0.01);
}

static void test_mdof_mac(void)
{
    double phi_a[] = {1.0, 2.0, 3.0};
    double phi_b[] = {1.0, 2.0, 3.0};
    double mac_val = mdof_mac(phi_a, phi_b, 3);
    /* Perfect correlation: MAC = 1 */
    assert(fabs(mac_val - 1.0) < TOL);
}

static void test_rayleigh_damping(void)
{
    double alpha, beta;
    rayleigh_damping_coeffs(10.0, 100.0, 0.02, 0.05, &alpha, &beta);
    /* Verify: zeta at w1 should be ~0.02 */
    double z1 = rayleigh_modal_damping(alpha, beta, 10.0);
    assert(fabs(z1 - 0.02) < 0.01);
}

static void test_mdof_energy(void)
{
    MDOFSystem sys;
    mdof_build_2dof(&sys, 10.0, 5.0, 500.0, 300.0, 100.0, 5.0, 3.0, 1.0);
    double vel[] = {1.0, 0.5};
    double disp[] = {0.1, 0.05};
    double T = mdof_kinetic_energy(&sys, vel);
    assert(T > 0.0);
    double U = mdof_potential_energy(&sys, disp);
    assert(U > 0.0);
    double E = mdof_total_energy(&sys, disp, vel);
    assert(NEAR(E, T + U));
}

/* === L6-L7: Applications === */

static void test_quarter_car(void)
{
    QuarterCarModel qcm = {350.0, 40.0, 25000.0, 2000.0, 200000.0, 50.0};
    double body_hz, wheel_hz;
    quarter_car_natural_freqs(&qcm, &body_hz, &wheel_hz);
    assert(body_hz > 1.0 && body_hz < 2.0);
    assert(wheel_hz > 10.0 && wheel_hz < 15.0);

    double defl = quarter_car_static_deflection(&qcm, 9.81);
    assert(defl > 0.1 && defl < 0.2);

    double zeta = quarter_car_damping_ratio(&qcm);
    assert(zeta > 0.2 && zeta < 0.5);
}

static void test_vibration_isolator(void)
{
    IsolatorReq req = {100.0, 100.0, 1.0, 0.1, 250.0, 4};
    IsolatorDesign design = design_vibration_isolator(&req);
    assert(design.natural_freq_hz < 20.0);
    assert(design.freq_ratio > 2.0);
    assert(design.achieved_transmissibility < 0.3);
}

static void test_robot_joint(void)
{
    RobotJoint j = {0.5, 10.0, 1.0, 100.0, 0.001, 50000.0, 10.0, 200.0};
    double J_eff = robot_joint_effective_inertia(&j, 1.0);
    assert(J_eff > 0.001);

    double wn = robot_joint_natural_freq(&j, 1.0);
    assert(wn > 0.0);

    double Kp, Kd;
    robot_joint_pd_gains(&j, 50.0, 0.7, &Kp, &Kd);
    assert(Kp > 0.0);
    assert(Kd > 0.0);
}

static void test_precision_machine(void)
{
    double k_frame = machine_frame_stiffness(10000.0, 0.00001);
    assert(fabs(k_frame - 1e9) < 1e-3); /* floating point: 10000/0.00001 */

    double dL = thermal_expansion(1.0, 12e-6, 50.0);
    assert(NEAR(dL, 0.0006));
}

static void test_aerospace(void)
{
    double bw = actuator_bandwidth_for_control(20.0, 6.0, 45.0);
    assert(bw >= 100.0);

    double fd = solar_panel_deploy_freq(5.0, 2.0, 100.0);
    assert(fd > 0.0);
}

/* === main === */

int main(void)
{
    printf("\n=== mini-mechanical-modeling Test Suite ===\n\n");

    printf("[L1-L2] Mechanical Elements:\n");
    RUN_TEST(test_newton_second_law);
    RUN_TEST(test_hookes_law);
    RUN_TEST(test_friction);
    RUN_TEST(test_stribeck_friction);
    RUN_TEST(test_energy);
    RUN_TEST(test_natural_frequency);
    RUN_TEST(test_damping_params);
    RUN_TEST(test_log_decrement);
    RUN_TEST(test_q_factor);
    RUN_TEST(test_static_deflection);
    RUN_TEST(test_mechanical_impedance);
    RUN_TEST(test_transmissibility);
    RUN_TEST(test_moments_of_inertia);
    RUN_TEST(test_time_domain_response);
    RUN_TEST(test_step_response);
    RUN_TEST(test_resonance);
    RUN_TEST(test_structural_stiffness);

    printf("\n[L3-L5] Transfer Functions:\n");
    RUN_TEST(test_tf_sdof);
    RUN_TEST(test_tf_base_excitation);
    RUN_TEST(test_tf_operations);
    RUN_TEST(test_tf_evaluation);
    RUN_TEST(test_tf_damping_estimate);
    RUN_TEST(test_tf_steady_state);
    RUN_TEST(test_tf_freq_response);

    printf("\n[L2-L5] State-Space:\n");
    RUN_TEST(test_ss_sdof);
    RUN_TEST(test_ss_stability);
    RUN_TEST(test_ss_controllability);
    RUN_TEST(test_ss_simulation);
    RUN_TEST(test_ss_2dof);
    RUN_TEST(test_ss_inverted_pendulum);
    RUN_TEST(test_ss_discretize);
    RUN_TEST(test_ss_modal);
    RUN_TEST(test_ss_gramians);

    printf("\n[L1-L7] Rotational Systems:\n");
    RUN_TEST(test_rotational_newton);
    RUN_TEST(test_torsional_spring);
    RUN_TEST(test_gear_ratio);
    RUN_TEST(test_optimal_gear_ratio);
    RUN_TEST(test_rack_pinion);
    RUN_TEST(test_leadscrew);
    RUN_TEST(test_rotating_unbalance);
    RUN_TEST(test_flywheel);

    printf("\n[L5-L8] Multi-DOF Systems:\n");
    RUN_TEST(test_mdof_2dof);
    RUN_TEST(test_mdof_chain);
    RUN_TEST(test_mdof_modal_params);
    RUN_TEST(test_mdof_orthogonality);
    RUN_TEST(test_mdof_mac);
    RUN_TEST(test_rayleigh_damping);
    RUN_TEST(test_mdof_energy);

    printf("\n[L6-L7] Applications:\n");
    RUN_TEST(test_quarter_car);
    RUN_TEST(test_vibration_isolator);
    RUN_TEST(test_robot_joint);
    RUN_TEST(test_precision_machine);
    RUN_TEST(test_aerospace);

    printf("\n========================================\n");
    printf("Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("========================================\n");

    return (tests_passed == tests_run) ? 0 : 1;
}
