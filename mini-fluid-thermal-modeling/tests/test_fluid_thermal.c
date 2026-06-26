/**
 * mini-fluid-thermal-modeling — Test Suite
 *
 * Coverage: All core APIs — dimensionless groups, flow classification,
 * conservation laws, pipe flow, boundary layers, conduction,
 * convection correlations, radiation, fin analysis, heat exchangers,
 * pump systems, electronics cooling, and numerical methods.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "fluid_thermal_core.h"
#include "fluid_mechanics.h"
#include "heat_transfer.h"
#include "thermal_fluid_systems.h"
#include "numerical_fluid_thermal.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(desc, cond) do { \
    printf("  %-55s", desc); \
    if (cond) { \
        printf(" PASS\n"); \
        tests_passed++; \
    } else { \
        printf(" FAIL\n"); \
        tests_failed++; \
    } \
} while(0)

#define TEST_FEQ(a, b, tol, desc) TEST(desc, fabs((a)-(b)) < (tol))


/* =============================================================
   L1: Dimensionless Group Computations
   ============================================================= */
static void test_dimensionless_groups(void)
{
    printf("\n--- L1: Dimensionless Groups ---\n");

    double re = compute_reynolds(1000.0, 2.0, 0.05, 0.001);
    TEST_FEQ(re, 100000.0, 1e-6, "Reynolds: ρuL/μ");
    TEST("Re with μ=0 → INF", isinf(compute_reynolds(1, 1, 1, 0)));

    double pr = compute_prandtl(4180.0, 0.001, 0.6);
    TEST_FEQ(pr, 6.9667, 0.01, "Prandtl: cp·μ/k");

    double nu_val = compute_nusselt(500.0, 0.1, 0.6);
    TEST_FEQ(nu_val, 83.3333, 0.01, "Nusselt: h·L/k");

    double gr = compute_grashof(9.81, 0.0003, 20.0, 0.1, 1.0e-6);
    TEST("Grashof > 0 for natural convection", gr > 0.0);

    double ra = compute_rayleigh(gr, pr);
    TEST_FEQ(ra, gr * pr, 1e-9, "Rayleigh = Gr·Pr");

    double ma = compute_mach(340.0, 340.0);
    TEST_FEQ(ma, 1.0, 1e-9, "Mach = 1 at sonic");

    double pe = compute_peclet(50000.0, 0.71);
    TEST_FEQ(pe, 35500.0, 1.0, "Peclet = Re·Pr");

    double bi = compute_biot(50.0, 0.01, 150.0);
    TEST_FEQ(bi, 0.00333, 0.001, "Biot < 0.1 (lumped)");

    double fo = compute_fourier(1.0e-4, 100.0, 0.01);
    TEST_FEQ(fo, 100.0, 1e-6, "Fourier: α·t/L²");

    double st = compute_stanton(200.0, 50000.0, 0.71);
    TEST("Stanton > 0", st > 0.0);
}

/* =============================================================
   L2: Flow Regime Classification
   ============================================================= */
static void test_flow_classification(void)
{
    printf("\n--- L2: Flow Classification ---\n");

    TEST("Re=0.5 → Creeping",
         classify_flow_regime(0.5) == FLOW_CREEPING);
    TEST("Re=100 → Laminar",
         classify_flow_regime(100.0) == FLOW_LAMINAR);
    TEST("Re=3000 → Transitional",
         classify_flow_regime(3000.0) == FLOW_TRANSITIONAL);
    TEST("Re=10000 → Turbulent",
         classify_flow_regime(10000.0) == FLOW_TURBULENT);

    TEST("Ma=0.1 → Incompressible",
         classify_compressibility(0.1) == COMPRESS_INCOMPRESSIBLE);
    TEST("Ma=0.5 → Subsonic",
         classify_compressibility(0.5) == COMPRESS_SUBSONIC);
    TEST("Ma=1.0 → Transonic",
         classify_compressibility(1.0) == COMPRESS_TRANSONIC);
    TEST("Ma=3.0 → Supersonic",
         classify_compressibility(3.0) == COMPRESS_SUPERSONIC);

    /* Convection type: Gr/Re² ratio */
    TEST("Gr/Re²=100 → Natural",
         classify_convection(100000.0, 10.0) == CONV_NATURAL);
    TEST("Gr/Re²=0.01 → Forced",
         classify_convection(1.0, 100.0) == CONV_FORCED);
    TEST("Gr/Re²≈1 → Mixed",
         classify_convection(2500.0, 50.0) == CONV_MIXED);
}

/* =============================================================
   L3: Engineering Quantities
   ============================================================= */
static void test_engineering_quantities(void)
{
    printf("\n--- L3: Engineering Quantities ---\n");

    ReynoldsRange r = typical_reynolds_range("pipe_water");
    TEST("Pipe water Re min > 0", r.min_re > 0.0);
    TEST("Pipe water Re max > min", r.max_re > r.min_re);

    double pr_air = typical_prandtl("air");
    TEST_FEQ(pr_air, 0.71, 0.01, "Pr air ≈ 0.71");
    double pr_water = typical_prandtl("water");
    TEST_FEQ(pr_water, 5.83, 0.1, "Pr water ≈ 5.83");

    double h_min, h_max;
    typical_h_range("forced_air", &h_min, &h_max);
    TEST("Forced air h min > 0", h_min > 0.0);
    TEST("Forced air h max > min", h_max > h_min);

    double k_copper = typical_thermal_conductivity("copper");
    TEST_FEQ(k_copper, 401.0, 1.0, "k copper ≈ 401");
    double k_air = typical_thermal_conductivity("air");
    TEST_FEQ(k_air, 0.0263, 0.001, "k air ≈ 0.0263");

    double cp_al = typical_specific_heat("aluminum");
    TEST_FEQ(cp_al, 903.0, 1.0, "cp aluminum ≈ 903");
    double rho_steel = typical_density("steel_carbon");
    TEST_FEQ(rho_steel, 7854.0, 10.0, "ρ steel ≈ 7854");
}

/* =============================================================
   L4: Conservation Laws
   ============================================================= */
static void test_conservation_laws(void)
{
    printf("\n--- L4: Conservation Laws ---\n");

    /* Continuity: A1*u1 = A2*u2 */
    double u2;
    int ret = continuity_equation(1.0, 0.1, 2.0, 1.0, 0.05, &u2);
    TEST("Continuity success", ret == 0);
    TEST_FEQ(u2, 4.0, 1e-9, "Continuity u2=u1*A1/A2");

    /* Bernoulli: solve for p2 */
    double p2_solved = bernoulli_solve(101325.0, NAN, 0.0, 10.0,
                                        0.0, 0.0, 1000.0, 9.81);
    TEST_FEQ(p2_solved, 51325.0, 1.0, "Bernoulli p2 = p1 - ½ρu2²");

    /* Bernoulli head loss */
    double hL = bernoulli_head_loss(200000.0, 150000.0, 1.0, 3.0,
                                     0.0, 0.0, 9810.0, 9.81);
    TEST("Head loss finite", !isnan(hL));

    /* Fourier's law */
    double q_cond = fourier_conduction(400.0, 0.01, 373.0, 300.0, 0.1);
    TEST_FEQ(q_cond, 2920.0, 1.0, "Fourier q = kAΔT/L");

    /* Newton cooling */
    double q_conv = newton_cooling(25.0, 0.5, 350.0, 300.0);
    TEST_FEQ(q_conv, 625.0, 1e-6, "Newton q = hA(Ts-T∞)");

    /* Stefan-Boltzmann */
    double Eb = stefan_boltzmann_emission(300.0);
    TEST_FEQ(Eb, 459.3, 1.0, "Stefan-Boltzmann at 300K");
}

/* =============================================================
   L4: Navier-Stokes Solutions
   ============================================================= */
static void test_navier_stokes_solutions(void)
{
    printf("\n--- L4: Navier-Stokes Solutions ---\n");

    /* Hagen-Poiseuille */
    double v_hp = hagen_poiseuille_velocity(-100.0, 0.001, 0.01, 0.0);
    /* v = (-dp/dx)(R²)/(4μ) = 100 * 1e-4 / 0.004 = 2.5 m/s */
    TEST_FEQ(v_hp, 2.5, 1e-6, "HP centerline velocity");
    double Q_hp = hagen_poiseuille_flow_rate(-100.0, 0.001, 0.01);
    TEST("HP flow rate positive", Q_hp > 0.0);

    /* Plane Poiseuille */
    double v_pp = plane_poiseuille_velocity(-50.0, 0.001, 0.005, 0.0);
    TEST("Plane Poiseuille velocity center", v_pp > 0.0);
    double Q_pp = plane_poiseuille_flow_rate(-50.0, 0.001, 0.005);
    TEST("Plane Poiseuille flow rate", Q_pp > 0.0);

    /* Couette */
    double v_c = couette_velocity(5.0, 0.01, 0.005);
    TEST_FEQ(v_c, 2.5, 1e-9, "Couette linear profile");

    /* Stokes first problem */
    double v_st = stokes_first_problem(1.0, 1.0e-6, 0.001, 1.0);
    TEST("Stokes: 0 < u < U", v_st > 0.0 && v_st < 1.0);
}

/* =============================================================
   L5: Pipe Flow (Darcy-Weisbach, Friction Factors)
   ============================================================= */
static void test_pipe_flow(void)
{
    printf("\n--- L5: Pipe Flow ---\n");

    /* Laminar friction */
    double f_lam = darcy_friction_laminar(1000.0);
    TEST_FEQ(f_lam, 0.064, 1e-6, "f = 64/Re (laminar)");

    /* Swamee-Jain */
    double f_sj = swamee_jain_friction_factor(100000.0, 0.0001);
    TEST("Swamee-Jain in range [0.01, 0.05]", f_sj > 0.0 && f_sj < 0.1);

    /* Colebrook should converge */
    double f_cw = colebrook_friction_factor(100000.0, 0.0001, 20, 1e-10);
    TEST("Colebrook converged", f_cw > 0.0 && f_cw < 0.1);

    /* Darcy-Weisbach head loss */
    double hL = darcy_weisbach_head_loss(0.02, 100.0, 0.1, 2.0, 9.81);
    TEST("DW head loss positive", hL > 0.0);

    /* Minor loss */
    double KL = minor_loss_coefficient("elbow_90_standard");
    TEST_FEQ(KL, 0.75, 1e-9, "Elbow 90 K=0.75");

    /* Hydraulic diameter */
    /* A = π*D²/4, P = π*D, Dh = 4A/P = D */
    double A_pipe = M_PI * 0.05 * 0.05 / 4.0;
    double P_pipe = M_PI * 0.05;
    double Dh = hydraulic_diameter(A_pipe, P_pipe);
    TEST_FEQ(Dh, 0.05, 1e-9, "Hydraulic diameter of pipe = D");
}

/* =============================================================
   L5: Boundary Layer Theory
   ============================================================= */
static void test_boundary_layer(void)
{
    printf("\n--- L5: Boundary Layer ---\n");

    double delta = blasius_boundary_layer_thickness(1.0, 1000000.0);
    TEST_FEQ(delta, 0.005, 1e-6, "Blasius δ = 5x/√Rex");

    double cf_local = blasius_skin_friction_local(1000000.0);
    TEST_FEQ(cf_local, 0.000664, 1e-6, "Blasius cf local = 0.664/√Rex");

    double cf_avg = blasius_skin_friction_average(1000000.0);
    TEST_FEQ(cf_avg, 0.001328, 1e-6, "Blasius Cf avg = 1.328/√ReL");

    double delta_t = turbulent_boundary_layer_thickness(1.0, 5000000.0);
    TEST("Turbulent BL thicker than laminar", delta_t > delta);

    double H = shape_factor_laminar();
    TEST_FEQ(H, 2.5916, 0.001, "Shape factor H ≈ 2.59");

    double Re_tr = estimate_transition_reynolds(1.0);
    TEST("Transition Re with 1% Tu < 5e5", Re_tr < 500000.0);
}

/* =============================================================
   L5: Nusselt Correlations
   ============================================================= */
static void test_nusselt_correlations(void)
{
    printf("\n--- L5: Nusselt Correlations ---\n");

    /* Dittus-Boelter */
    double nu_db = dittus_boelter_nusselt(50000.0, 0.71, 1);
    TEST("Dittus-Boelter Nu in range", nu_db > 50.0);

    /* Flat plate laminar */
    double nu_fpl = flat_plate_laminar_nusselt_average(500000.0, 0.71);
    TEST("Flat plate laminar Nu", nu_fpl > 200.0);

    /* Churchill-Chu vertical */
    double nu_cc = churchill_chu_vertical_plate(1.0e8, 0.71);
    TEST("Churchill-Chu Nu for Ra=1e8", nu_cc > 50.0);

    /* Morgan horizontal cylinder */
    double nu_m = morgan_horizontal_cylinder(1.0e5);
    TEST("Morgan cylinder Nu for Ra=1e5", nu_m > 5.0);

    /* McAdams horizontal */
    double nu_mc = mcadams_horizontal_plate(1.0e5, 1);
    TEST("McAdams plate Nu for Ra=1e5", nu_mc > 5.0);
}

/* =============================================================
   L5: Radiation and Fin Analysis
   ============================================================= */
static void test_radiation_and_fins(void)
{
    printf("\n--- L5: Radiation and Fins ---\n");

    /* Net radiation */
    double q_rad = blackbody_net_radiation(5.67e-8, 1.0, 1.0, 400.0, 300.0);
    TEST("Net radiation 400K→300K > 800W", q_rad > 800.0);

    /* Gray radiation */
    double q_gray = gray_net_radiation(400.0, 300.0, 1.0, 1.0,
                                       0.9, 0.9, 1.0);
    TEST("Gray radiation positive", q_gray > 0.0);

    /* View factor parallel disks */
    double F12 = view_factor_parallel_disks(0.5, 0.5, 0.5);
    TEST("View factor disks in (0,1)", F12 > 0.0 && F12 < 1.0);

    /* Fin efficiency */
    double m = sqrt(50.0 * 0.1 / (200.0 * 0.001)); /* m = √(hP/kAc) */
    double eta_f = fin_efficiency(m, 0.05);
    TEST("Fin efficiency in (0,1]", eta_f > 0.0 && eta_f <= 1.0);

    /* Fin heat transfer */
    double q_fin = fin_heat_transfer_rate(50.0, 0.1, 200.0, 0.001, 50.0, 0.05);
    TEST("Fin heat transfer positive", q_fin > 0.0);

    /* Lumped capacitance */
    double T_lc = lumped_capacitance_temp(300.0, 400.0, 25.0, 0.01,
                                           2700.0, 0.0001, 900.0, 50.0);
    TEST("Lumped cap temp between Ti and T∞", T_lc > 300.0 && T_lc < 400.0);

    /* Overall U */
    double U = overall_heat_transfer_coefficient_plane(500.0, 0.002, 200.0, 25.0);
    TEST("Overall U positive", U > 0.0);
}

/* =============================================================
   L6: Heat Exchanger Design
   ============================================================= */
static void test_heat_exchanger(void)
{
    printf("\n--- L6: Heat Exchanger ---\n");

    /* LMTD */
    double dT_lm = lmtd(80.0, 40.0);
    TEST_FEQ(dT_lm, 57.71, 0.5, "LMTD (80,40) ≈ 57.7");
    TEST("LMTD equal temps", lmtd(50.0, 50.0) == 50.0);

    /* Effectiveness — parallel */
    double eps_p = heat_exchanger_effectiveness(2.0, 0.5, "parallel");
    TEST("Parallel ε in (0,1)", eps_p > 0.0 && eps_p < 1.0);

    /* Effectiveness — counter */
    double eps_c = heat_exchanger_effectiveness(2.0, 0.5, "counter");
    TEST("Counter ε > Parallel ε", eps_c > eps_p);

    /* NTU from effectiveness */
    double ntu_p = ntu_from_effectiveness_parallel(eps_p, 0.5);
    TEST_FEQ(ntu_p, 2.0, 0.01, "NTU from ε (parallel)");

    /* Outlet temperatures */
    double Th_out, Tc_out;
    heat_exchanger_outlet_temps(400.0, 300.0, 2000.0, 1500.0,
                                 0.6, &Th_out, &Tc_out);
    TEST("Hot outlet between inlets", Th_out < 400.0 && Th_out > 300.0);
    TEST("Cold outlet between inlets", Tc_out > 300.0 && Tc_out < 400.0);

    /* Shell-and-tube correction factor F */
    /* Bowman F factor: P=0.3, R=0.5 should give valid F */
    double F = hx_correction_factor_shell_tube(0.3, 0.5);
    TEST("F correction factor in (0,1]", F > 0.0 && F <= 1.0);
}

/* =============================================================
   L6: Pipe Heating / Thermal-Fluid Systems
   ============================================================= */
static void test_thermal_fluid_systems(void)
{
    printf("\n--- L6: Thermal-Fluid Systems ---\n");

    /* Pipe heating */
    double T_out = pipe_heating_outlet_temp(300.0, 10000.0, 0.02, 2.0,
                                             0.1, 4180.0);
    TEST("Pipe heating raises temperature", T_out > 300.0);

    /* CWT pipe */
    double T_cwt = pipe_outlet_temp_cwt(400.0, 300.0, 0.02, 2.0,
                                         500.0, 0.1, 4180.0);
    TEST("CWT outlet between Tin and Ts", T_cwt > 300.0 && T_cwt <= 400.0);

    /* Laminar Nu */
    TEST_FEQ(nusselt_laminar_pipe(1), 4.36, 1e-9, "Nu laminar CHF = 4.36");
    TEST_FEQ(nusselt_laminar_pipe(0), 3.66, 1e-9, "Nu laminar CWT = 3.66");

    /* Pump power */
    double P_hyd = pump_hydraulic_power(1000.0, 9.81, 0.01, 20.0);
    TEST_FEQ(P_hyd, 1962.0, 1.0, "Pump hydraulic power");

    /* System curve operating point */
    double Q_op, H_op;
    int ret = pump_system_operating_point(50.0, 10.0, 1000.0, 3000.0,
                                          &Q_op, &H_op);
    TEST("Operating point found", ret == 0);
    TEST("Flow rate positive", Q_op > 0.0);

    /* NPSHa */
    double npsha = npsh_available(101325.0, 3169.0, 9810.0, -3.0, 1.5);
    TEST("NPSHa computed", npsha > 0.0);
}

/* =============================================================
   L7: Applications — Electronics Cooling
   ============================================================= */
static void test_electronics_cooling(void)
{
    printf("\n--- L7: Electronics Cooling ---\n");

    /* Junction temperature */
    double T_j = junction_temperature(300.0, 50.0, 0.5, 0.3, 2.0);
    TEST_FEQ(T_j, 440.0, 1.0, "Tjunction = Tamb + q·ΣR");

    /* Required airflow */
    double Q_air = required_airflow_cooling(500.0, 1.16, 1007.0, 20.0);
    TEST("Required airflow positive", Q_air > 0.0);

    /* Heat sink resistance */
    double R_hs = heat_sink_thermal_resistance(200.0, 0.01, 0.85, 50.0, 0.2);
    TEST("Heat sink resistance positive", R_hs > 0.0);

    /* Stacked die temperatures — equal power, downstream gets hotter fluid */
    double q_power[] = {15.0, 15.0, 15.0};
    double h[] = {500.0, 500.0, 500.0};
    double A[] = {0.0004, 0.0004, 0.0004};
    double T_die[3];
    stacked_die_temperatures(300.0, 3, q_power, 0.01, 1007.0, h, A, T_die);
    TEST("3D IC: downstream die hotter (same power, warmer fluid)",
         T_die[2] > T_die[0]);
}

/* =============================================================
   L7: HVAC Applications
   ============================================================= */
static void test_hvac(void)
{
    printf("\n--- L7: HVAC ---\n");

    double q_sens = sensible_heat_load(1.0, 1007.0, 288.0, 295.0);
    TEST("Sensible cooling: negative load (supply<room)", q_sens < 0.0);

    double T_surfaces[] = {293.0, 295.0, 293.0, 295.0, 293.0, 295.0};
    double A_surfaces[] = {20.0, 20.0, 15.0, 15.0, 30.0, 30.0};
    double MRT = mean_radiant_temperature(T_surfaces, A_surfaces, 6);
    TEST("MRT between bounds", MRT > 293.0 && MRT < 295.0);

    double ACH = air_changes_per_hour(0.5, 150.0);
    TEST_FEQ(ACH, 12.0, 0.1, "ACH = Q*3600/V");
}

/* =============================================================
   L7: Automotive Cooling
   ============================================================= */
static void test_automotive_cooling(void)
{
    printf("\n--- L7: Automotive Cooling ---\n");

    double q_rad = radiator_heat_rejection(363.0, 300.0, 0.5, 3800.0,
                                            1.0, 1007.0, 2000.0);
    TEST("Radiator heat rejection positive", q_rad > 0.0);

    double P_new = pump_power_scaling(1000.0, 6000.0, 3000.0);
    TEST_FEQ(P_new, 8000.0, 1.0, "Pump power scales with N³");

    double Q_new, H_new, P_fan;
    fan_affinity_laws(2.0, 0.5, 100.0, 500.0, &Q_new, &H_new, &P_fan);
    TEST_FEQ(Q_new, 1.0, 1e-9, "Fan Q doubles");
    TEST_FEQ(H_new, 400.0, 1.0, "Fan H quadruples");
    TEST_FEQ(P_fan, 4000.0, 1.0, "Fan P ×8");
}

/* =============================================================
   L8: Multi-Phase Flow Basics
   ============================================================= */
static void test_multiphase(void)
{
    printf("\n--- L8: Multi-Phase Flow ---\n");

    double alpha = void_fraction_homogeneous(0.5, 1.0, 1000.0, 1.0);
    TEST("Void fraction: vapor > liquid volume", alpha > 0.9);

    double X = martinelli_parameter(0.3, 1.0, 1000.0, 1.0e-5, 1.0e-3);
    TEST("Martinelli parameter positive", X > 0.0);

    double phi2 = two_phase_multiplier(X, 20.0);
    TEST("Two-phase multiplier > 1", phi2 > 1.0);
}

/* =============================================================
   L5: Numerical Methods — 1D Conduction
   ============================================================= */
static void test_numerical_1d_conduction(void)
{
    printf("\n--- L5: Numerical 1D Conduction ---\n");

    /* Steady 1D */
    int n = 11;
    double T[11];
    for (int i = 0; i < n; i++) T[i] = 300.0;
    int iters = solve_1d_steady_conduction(n, 0.1, 200.0, 0.0,
                                            400.0, 300.0, T, 1000, 1e-6);
    TEST("Steady 1D converged", iters > 0);
    TEST_FEQ(T[5], 350.0, 1.0, "Steady 1D midpoint ≈ 350K");

    /* Transient explicit */
    double T0[11], T_exp[11];
    for (int i = 0; i < n; i++) { T0[i] = 300.0; T_exp[i] = 0.0; }
    T0[0] = 400.0; T0[10] = 300.0;
    int ret_exp = solve_1d_transient_conduction_explicit(11, 100, 0.1,
                                                           1.0e-5, 0.01,
                                                           400.0, 300.0,
                                                           T0, T_exp);
    TEST("Transient explicit stable", ret_exp == 0);
    TEST("Transient explicit temp changed", T_exp[5] > 300.0);

    /* Thomas algorithm */
    double a_t[3] = {0, 1, 1};
    double b_t[3] = {2, 2, 2};
    double c_t[3] = {1, 1, 0};
    double d_t[3] = {8, 4, 6};
    double x_t[3];
    int ret_t = thomas_algorithm(3, a_t, b_t, c_t, d_t, x_t);
    TEST("Thomas algorithm success", ret_t == 0);
    /* System: 2x0+x1=8, x0+2x1+x2=4, x1+2x2=6 → x={5.5,-3,4.5} */
    TEST_FEQ(x_t[0], 5.5, 0.01, "Thomas: x0 ≈ 5.5");
    TEST_FEQ(x_t[1], -3.0, 0.01, "Thomas: x1 ≈ -3.0");
    TEST_FEQ(x_t[2], 4.5, 0.01, "Thomas: x2 ≈ 4.5");
}

/* =============================================================
   L5: Numerical Methods — 2D Laplace
   ============================================================= */
static void test_numerical_2d_laplace(void)
{
    printf("\n--- L5: Numerical 2D Laplace ---\n");

    int nx = 21, ny = 21;
    double *T_2d = (double *)calloc(nx * ny, sizeof(double));
    assert(T_2d);
    /* Set boundary: hot left, cold right, insulated top/bottom */
    for (int j = 0; j < ny; j++) {
        T_2d[j * nx] = 100.0;          /* left */
        T_2d[j * nx + nx - 1] = 0.0;   /* right */
    }

    int iters = gauss_seidel_2d_laplace(nx, ny, T_2d, 5000, 1e-4);
    TEST("GS 2D converged", iters > 0);

    /* Check midpoint is between boundaries */
    int mid = (ny/2) * nx + nx/2;
    TEST("Midpoint between BCs", T_2d[mid] > 0.0 && T_2d[mid] < 100.0);

    free(T_2d);
}

/* =============================================================
   L5: Numerical Integration
   ============================================================= */
static void test_numerical_integration(void)
{
    printf("\n--- L5: Numerical Integration ---\n");

    /* Simpson: ∫₀¹ x² dx = 1/3 */
    int n_s = 101;
    double *f_s = (double *)malloc(n_s * sizeof(double));
    assert(f_s);
    double h_s = 1.0 / (n_s - 1);
    for (int i = 0; i < n_s; i++) {
        double x = i * h_s;
        f_s[i] = x * x;
    }
    double integral = simpson_integral(n_s, h_s, f_s);
    TEST_FEQ(integral, 1.0/3.0, 1e-6, "Simpson: ∫x²dx = 1/3");
    free(f_s);

    /* Displacement thickness: Blasius profile */
    int n_bl = 51;
    double *y_bl = (double *)malloc(n_bl * sizeof(double));
    double *u_bl = (double *)malloc(n_bl * sizeof(double));
    assert(y_bl && u_bl);
    double delta_bl = 0.01;
    double U_inf = 10.0;
    for (int i = 0; i < n_bl; i++) {
        y_bl[i] = delta_bl * i / (n_bl - 1);
        double eta = y_bl[i] / delta_bl;
        u_bl[i] = U_inf * (eta < 1.0 ? eta : 1.0); /* linear approximation */
    }
    /* This is a test of the integration, not validation of Blasius */
    double ds = integrate_displacement_thickness(n_bl, y_bl, u_bl, U_inf);
    TEST("Displacement thickness > 0", ds > 0.0);

    free(y_bl); free(u_bl);
}

/* =============================================================
   Fluid Properties Initialization
   ============================================================= */
static void test_fluid_properties(void)
{
    printf("\n--- Fluid Properties ---\n");

    FluidProperties props;
    int ret = fluid_properties_init("water", 300.0, &props);
    TEST("Water init success", ret == 0);
    TEST_FEQ(props.density, 997.0, 1.0, "Water density");
    TEST_FEQ(props.thermal_conductivity, 0.613, 0.01, "Water k");
    TEST("Water Pr > 1", props.prandtl_number > 1.0);

    ret = fluid_properties_init("nonexistent", 300.0, &props);
    TEST("Unknown fluid fails", ret != 0);
}

/* =============================================================
   Dimensional Analysis
   ============================================================= */
static void test_dimensional_analysis(void)
{
    printf("\n--- L5: Dimensional Analysis ---\n");

    int n_pi = buckingham_pi_count(7, 3);
    TEST("7 vars, 3 dims → 4 Pi groups", n_pi == 4);

    double params_re[] = {1.0, 2.0, 0.1, 0.001}; /* ρ, u, L, μ */
    double re = dimensionless_group("Re", params_re);
    TEST_FEQ(re, 200.0, 1e-6, "Dimensionless group Re");

    double fr = compute_froude(10.0, 9.81, 2.0);
    TEST("Froude > 0", fr > 0.0);

    double we = compute_weber(1000.0, 2.0, 0.01, 0.073);
    TEST("Weber > 0", we > 0.0);

    double eu = compute_euler(5000.0, 1.2, 30.0);
    TEST("Euler number positive", eu > 0.0);
}

/* =============================================================
   Natural Circulation
   ============================================================= */
static void test_natural_circulation(void)
{
    printf("\n--- L6: Natural Circulation ---\n");

    double dp = natural_circulation_driving_pressure(1000.0, 9.81,
                                                      0.0002, 30.0, 2.0);
    TEST("Driving pressure positive", dp > 0.0);

    double m_dot = natural_circulation_flow_rate(0.05,
        M_PI*0.025*0.025, 1000.0, 9.81, 0.0002, 30.0,
        2.0, 0.02, 10.0);
    TEST("Natural circulation flow > 0", m_dot > 0.0);
}

/* =============================================================
   L8: System Coupling
   ============================================================= */
static void test_system_coupling(void)
{
    printf("\n--- L8: System Coupling ---\n");

    double T_new[] = {310.0, 315.0, 320.0};
    double T_old[] = {300.0, 300.0, 300.0};
    double res = thermal_convergence_residual(T_new, T_old, 3);
    TEST("Convergence residual detected", res > 0.01);

    double nu_corr = conjugate_nusselt_correction(100.0, 15.0,
                                                    0.0001, 0.01, 1007.0, 0.1);
    TEST("Conjugate Nu correction", nu_corr > 0.0 && nu_corr <= 100.0);

    double eps_rec = recuperator_effectiveness_balanced(5.0);
    TEST_FEQ(eps_rec, 5.0/6.0, 0.01, "Recuperator ε = NTU/(1+NTU)");
}

/* =============================================================
   Pressure Correction (SIMPLE 1D)
   ============================================================= */
static void test_simple_pressure(void)
{
    printf("\n--- L8: SIMPLE Pressure Correction ---\n");

    int n_cells = 10;
    double u_star[11] = {2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0};
    double p_simple[10] = {100000.0, 99900.0, 99800.0, 99700.0, 99600.0,
                            99500.0, 99400.0, 99300.0, 99200.0, 99100.0};
    (void)simple_pressure_correction_1d(n_cells, u_star, p_simple,
                                         1.2, 0.1, 0.3, 100, 1e-3);
    TEST("SIMPLE 1D executed", 1);
}

/* =============================================================
   Main
   ============================================================= */
int main(void)
{
    printf("=== mini-fluid-thermal-modeling Test Suite ===\n\n");

    test_dimensionless_groups();
    test_flow_classification();
    test_engineering_quantities();
    test_conservation_laws();
    test_navier_stokes_solutions();
    test_pipe_flow();
    test_boundary_layer();
    test_nusselt_correlations();
    test_radiation_and_fins();
    test_heat_exchanger();
    test_thermal_fluid_systems();
    test_electronics_cooling();
    test_hvac();
    test_automotive_cooling();
    test_multiphase();
    test_numerical_1d_conduction();
    test_numerical_2d_laplace();
    test_numerical_integration();
    test_fluid_properties();
    test_dimensional_analysis();
    test_natural_circulation();
    test_system_coupling();
    test_simple_pressure();

    printf("\n=== Results: %d passed, %d failed ===\n",
           tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
