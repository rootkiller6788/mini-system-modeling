/**
 * mini-fluid-thermal-modeling — Coupled Thermal-Fluid Systems Implementation
 *
 * Pipe flow with heat transfer, pump/fan systems, electronics cooling,
 * HVAC analysis, automotive thermal management, natural circulation,
 * multi-phase flow, and system-level coupling.
 *
 * Knowledge Levels: L6-L8
 */

#include "thermal_fluid_systems.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =====================================================================
 * L6: Pipe Flow with Heat Transfer
 * ===================================================================== */

double pipe_heating_outlet_temp(double T_in, double q_w,
                                double D, double L,
                                double m_dot, double cp)
{
    if (m_dot <= 0.0 || cp <= 0.0) return T_in;
    double q_total = q_w * M_PI * D * L;
    return T_in + q_total / (m_dot * cp);
}

double pipe_wall_temperature_chf(double T_in, double q_w, double D,
                                 double L, double m_dot, double cp,
                                 double x, double h)
{
    if (m_dot <= 0.0 || cp <= 0.0) return T_in + q_w / ((h > 0) ? h : 1e10);
    if (h <= 0.0) return INFINITY;
    double T_fluid_at_x = T_in + (q_w * M_PI * D * L / (m_dot * cp)) * (x / L);
    return T_fluid_at_x + q_w / h;
}

double pipe_outlet_temp_cwt(double Ts, double T_in, double D, double L,
                            double h, double m_dot, double cp)
{
    if (m_dot <= 0.0 || cp <= 0.0 || h <= 0.0) return T_in;
    double exponent = -M_PI * D * L * h / (m_dot * cp);
    return Ts - (Ts - T_in) * exp(exponent);
}

double pipe_length_for_temp_change(double Ts, double T_in, double T_out,
                                   double D, double h,
                                   double m_dot, double cp)
{
    if (h <= 0.0 || m_dot <= 0.0 || cp <= 0.0 || D <= 0.0) return INFINITY;
    if ((T_out - Ts) == 0.0 || (T_in - Ts) == 0.0) return 0.0;
    double ratio = (T_out - Ts) / (T_in - Ts);
    if (ratio <= 0.0) return INFINITY; /* cannot reach T_out */
    return -(m_dot * cp / (M_PI * D * h)) * log(ratio);
}

double nusselt_laminar_pipe(int constant_heat_flux)
{
    return constant_heat_flux ? 4.36 : 3.66;
}

/* =====================================================================
 * L6: Pump and Fan System Curves
 * ===================================================================== */

double pump_hydraulic_power(double rho, double g, double Q, double H)
{
    return rho * g * Q * H;
}

double pump_shaft_power(double hydraulic_power, double efficiency)
{
    if (efficiency <= 0.0 || efficiency > 1.0) return INFINITY;
    return hydraulic_power / efficiency;
}

double system_curve_head(double H_static, double K_sys, double Q)
{
    return H_static + K_sys * Q * Q;
}

double pump_curve_head(double H_shutoff, double K_pump, double Q)
{
    return H_shutoff - K_pump * Q * Q;
}

int pump_system_operating_point(double H_shutoff, double H_static,
                                double K_pump, double K_sys,
                                double *Q_op, double *H_op)
{
    if (!Q_op || !H_op) return -1;
    double K_total = K_pump + K_sys;
    if (K_total <= 0.0 || H_shutoff < H_static) return -1;
    *Q_op = sqrt((H_shutoff - H_static) / K_total);
    *H_op = H_shutoff - K_pump * (*Q_op) * (*Q_op);
    return 0;
}

double npsh_available(double p_atm, double p_vapor, double gamma,
                      double z_suction, double hL_suction)
{
    if (gamma <= 0.0) return INFINITY;
    return (p_atm - p_vapor) / gamma + z_suction - hL_suction;
}

/* =====================================================================
 * L6: Natural Circulation Loops
 * ===================================================================== */

double natural_circulation_driving_pressure(double rho0, double g,
                                            double beta, double deltaT,
                                            double H_loop)
{
    return rho0 * g * beta * fabs(deltaT) * H_loop;
}

double natural_circulation_flow_rate(double D, double A, double rho,
                                     double g, double beta, double deltaT,
                                     double H, double f, double L_eff)
{
    if (f <= 0.0 || L_eff <= 0.0 || D <= 0.0) return 0.0;
    double numerator = 2.0 * D * A * A * rho * rho
                       * g * beta * fabs(deltaT) * H;
    double denominator = f * L_eff;
    if (denominator <= 0.0) return 0.0;
    if (numerator < 0.0) return 0.0;
    return sqrt(numerator / denominator);
}

/* =====================================================================
 * L7: Electronics Cooling Applications
 * ===================================================================== */

double heat_sink_thermal_resistance(double k_base, double A_base,
                                    double eta_fin, double h,
                                    double A_fin_total)
{
    double R_spreading = 0.0;
    double R_convection = 0.0;

    if (k_base > 0.0 && A_base > 0.0) {
        /* Simplified spreading resistance for square heat source
         * on square base. Semi-infinite approximation. */
        R_spreading = 1.0 / (2.0 * k_base * sqrt(A_base));
    }
    if (eta_fin > 0.0 && h > 0.0 && A_fin_total > 0.0) {
        R_convection = 1.0 / (eta_fin * h * A_fin_total);
    }
    return R_spreading + R_convection;
}

double junction_temperature(double T_ambient, double q,
                            double R_jc, double R_cs, double R_sa)
{
    return T_ambient + q * (R_jc + R_cs + R_sa);
}

double required_airflow_cooling(double q_total, double rho_air,
                                double cp_air, double deltaT)
{
    if (rho_air <= 0.0 || cp_air <= 0.0 || deltaT <= 0.0) return INFINITY;
    return q_total / (rho_air * cp_air * deltaT);
}

void stacked_die_temperatures(double T_in, int n_dies,
                              const double q_power[],
                              double m_dot, double cp,
                              const double h[], const double A[],
                              double T_die[])
{
    if (!q_power || !h || !A || !T_die || n_dies <= 0) return;
    double cumulative_q_dot_mcp = 0.0;

    for (int i = 0; i < n_dies; i++) {
        cumulative_q_dot_mcp += q_power[i] / (m_dot * cp);
        double T_fluid = T_in + cumulative_q_dot_mcp - q_power[i]/(m_dot*cp);
        T_die[i] = T_fluid + q_power[i] / (h[i] * A[i]);
    }
}

/* =====================================================================
 * L7: HVAC and Building Thermal Analysis
 * ===================================================================== */

double sensible_heat_load(double m_dot, double cp,
                          double T_supply, double T_room)
{
    return fabs(m_dot) * cp * (T_supply - T_room);
}

double latent_heat_load(double m_dot, double h_fg,
                        double w_room, double w_supply)
{
    return fabs(m_dot) * h_fg * (w_room - w_supply);
}

double air_changes_per_hour(double Q_air, double V_room)
{
    if (V_room <= 0.0) return INFINITY;
    return Q_air * 3600.0 / V_room;
}

double mean_radiant_temperature(const double T_surfaces[],
                                const double A_surfaces[],
                                int n_surfaces)
{
    if (!T_surfaces || !A_surfaces || n_surfaces <= 0) return NAN;
    double sum_TA = 0.0;
    double sum_A = 0.0;
    for (int i = 0; i < n_surfaces; i++) {
        sum_TA += T_surfaces[i] * A_surfaces[i];
        sum_A += A_surfaces[i];
    }
    if (sum_A == 0.0) return NAN;
    return sum_TA / sum_A;
}

/* =====================================================================
 * L7: Automotive / Industrial Cooling
 * ===================================================================== */

double radiator_heat_rejection(double T_coolant_in, double T_air_in,
                               double m_dot_coolant, double cp_coolant,
                               double m_dot_air, double cp_air,
                               double UA)
{
    double C_coolant = m_dot_coolant * cp_coolant;
    double C_air = m_dot_air * cp_air;
    double C_min = (C_coolant < C_air) ? C_coolant : C_air;
    double C_max = (C_coolant > C_air) ? C_coolant : C_air;
    if (C_min <= 0.0) return 0.0;
    double Cr = C_min / C_max;
    double NTU = UA / C_min;
    /* Cross-flow both unmixed effectiveness */
    double epsilon;
    if (Cr == 0.0) {
        epsilon = 1.0 - exp(-NTU);
    } else {
        double inner = exp(-Cr * pow(NTU, 0.78)) - 1.0;
        epsilon = 1.0 - exp((1.0 / Cr) * pow(NTU, 0.22) * inner);
    }
    return epsilon * C_min * (T_coolant_in - T_air_in);
}

double pump_power_scaling(double P_ref, double N_new, double N_ref)
{
    if (N_ref == 0.0) return INFINITY;
    double ratio = N_new / N_ref;
    return P_ref * ratio * ratio * ratio;
}

void fan_affinity_laws(double N_ratio,
                       double Q_ref, double H_ref, double P_ref,
                       double *Q_new, double *H_new, double *P_new)
{
    if (Q_new) *Q_new = Q_ref * N_ratio;
    if (H_new) *H_new = H_ref * N_ratio * N_ratio;
    if (P_new) *P_new = P_ref * N_ratio * N_ratio * N_ratio;
}

/* =====================================================================
 * L8: System-Level Coupling
 * ===================================================================== */

double thermal_convergence_residual(const double T_new[],
                                    const double T_old[], int n)
{
    if (!T_new || !T_old || n <= 0) return INFINITY;
    double max_res = 0.0;
    for (int i = 0; i < n; i++) {
        if (T_old[i] == 0.0) continue;
        double res = fabs(T_new[i] - T_old[i]) / fabs(T_old[i]);
        if (res > max_res) max_res = res;
    }
    return max_res;
}

double conjugate_nusselt_correction(double Nu, double k_wall,
                                    double A_wall, double m_dot,
                                    double cp, double L)
{
    if (m_dot <= 0.0 || cp <= 0.0 || L <= 0.0) return Nu;
    double correction = 1.0 + (k_wall * A_wall / (m_dot * cp * L)) * Nu;
    if (correction == 0.0) return Nu;
    return Nu / correction;
}

double recuperator_effectiveness_balanced(double ntu)
{
    return ntu / (1.0 + ntu);
}

double tube_bank_pressure_drop(int N_L, double chi, double rho,
                               double u_max, double f)
{
    return N_L * chi * (0.5 * rho * u_max * u_max) * f;
}

/* =====================================================================
 * L8: Multi-Phase Flow Basics
 * ===================================================================== */

double void_fraction_homogeneous(double quality, double rho_v,
                                 double rho_l, double slip_ratio)
{
    if (rho_v <= 0.0 || rho_l <= 0.0) return NAN;
    if (quality <= 0.0) return 0.0;
    if (quality >= 1.0) return 1.0;
    double denom = quality + (1.0 - quality) * (rho_v / rho_l) * slip_ratio;
    if (denom == 0.0) return NAN;
    return quality / denom;
}

double two_phase_multiplier(double X, double C)
{
    if (X <= 0.0) return INFINITY;
    return 1.0 + C / X + 1.0 / (X * X);
}

double martinelli_parameter(double quality, double rho_v, double rho_l,
                            double mu_v, double mu_l)
{
    if (quality <= 0.0) return INFINITY;  /* all liquid */
    if (quality >= 1.0) return 0.0;       /* all vapor */
    if (rho_v <= 0.0 || rho_l <= 0.0 || mu_v <= 0.0 || mu_l <= 0.0)
        return NAN;
    double term1 = pow((1.0 - quality) / quality, 0.9);
    double term2 = sqrt(rho_v / rho_l);
    double term3 = pow(mu_l / mu_v, 0.1);
    return term1 * term2 * term3;
}
