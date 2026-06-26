/**
 * Example 3: Electronics Cooling — CPU Heat Sink Design
 *
 * Designs an air-cooled aluminum heat sink for a 95W CPU.
 * Computes required fin area, airflow, junction temperature,
 * and evaluates thermal resistance budget.
 *
 * Application: Desktop CPU cooling with forced air convection.
 * Demonstrates L7 (electronics cooling) coupled analysis.
 */

#include <stdio.h>
#include <math.h>
#include "../include/fluid_thermal_core.h"
#include "../include/heat_transfer.h"
#include "../include/thermal_fluid_systems.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Example 3: CPU Heat Sink Design ===\n\n");

    /* CPU specifications */
    double q_cpu = 95.0;        /* W — CPU TDP */
    double T_amb = 308.15;      /* K — ambient (35°C inside case) */
    double T_j_max = 368.15;    /* K — max junction temp (95°C) */
    double T_j_target = 358.15; /* K — target (85°C with margin) */

    printf("CPU TDP:             %.0f W\n", q_cpu);
    printf("Ambient temperature: %.1f°C\n", T_amb - 273.15);
    printf("Max junction temp:   %.1f°C\n", T_j_max - 273.15);
    printf("Target junction:     %.1f°C\n\n", T_j_target - 273.15);

    /* Thermal budget */
    double total_budget = T_j_target - T_amb;
    printf("Total ΔT budget:    %.1f K\n\n", total_budget);

    /* Package resistances (typical values) */
    double R_jc = 0.33;  /* K/W — junction-to-case (IHS) */
    double R_cs = 0.05;  /* K/W — case-to-sink (TIM) */
    double R_sa_max = total_budget / q_cpu - R_jc - R_cs;

    printf("--- Thermal Resistance Budget ---\n");
    printf("R_jc (die→IHS):     %.2f K/W\n", R_jc);
    printf("R_cs (TIM):         %.2f K/W\n", R_cs);
    printf("R_sa max allowed:   %.3f K/W\n\n", R_sa_max);

    /* Heat sink geometry */
    double k_al = typical_thermal_conductivity("aluminum"); /* W/(m·K) */
    double A_base = 0.0016;    /* m² — 40×40 mm base */
    double L_fin = 0.04;       /* m — fin length */
    double t_fin = 0.001;      /* m — fin thickness */
    int n_fins = 20;
    double W_hs = 0.04;        /* m — heat sink width */

    printf("--- Heat Sink Geometry ---\n");
    printf("Base area:          %.0f mm²\n", A_base * 1e6);
    printf("Material:           Aluminum (k=%.0f W/m·K)\n", k_al);
    printf("Number of fins:     %d\n", n_fins);
    printf("Fin length:         %.0f mm\n", L_fin * 1000.0);
    printf("Fin thickness:      %.1f mm\n\n", t_fin * 1000.0);

    /* Total fin surface area */
    double A_single_fin = 2.0 * L_fin * W_hs; /* both sides */
    double A_fins_total = n_fins * A_single_fin;
    double A_base_exposed = A_base - n_fins * t_fin * W_hs;
    double A_total = A_fins_total + A_base_exposed;

    printf("Total fin area:     %.4f m²\n", A_fins_total);
    printf("Exposed base area:  %.4f m²\n", A_base_exposed);
    printf("Total area:         %.4f m²\n\n", A_total);

    /* Forced convection coefficient */
    /* Using fan: assume u_air ≈ 2 m/s over fins */
    double u_air = 2.0;
    double rho_air = 1.16;
    double mu_air = 1.85e-5;
    double k_air = 0.026;
    double cp_air = 1007.0;
    double D_h = 4.0 * (W_hs * 0.002) / (2.0 * (W_hs + 0.002));

    double Re_air = compute_reynolds(rho_air, u_air, D_h, mu_air);
    double Pr_air = compute_prandtl(cp_air, mu_air, k_air);

    printf("--- Airflow Analysis ---\n");
    printf("Air velocity:       %.1f m/s\n", u_air);
    printf("Hydraulic diameter: %.4f m\n", D_h);
    printf("Re_air:             %.0f\n", Re_air);
    printf("Pr_air:             %.3f\n", Pr_air);

    /* Nusselt — Dittus-Boelter for channel flow */
    double Nu_air;
    if (Re_air >= 10000.0) {
        Nu_air = dittus_boelter_nusselt(Re_air, Pr_air, 1);
    } else {
        double Re_L = Re_air; /* use flat plate laminar for simplicity */
        Nu_air = flat_plate_laminar_nusselt_average(Re_L, Pr_air);
    }
    double h_forced = compute_nusselt(Nu_air, D_h, k_air)
                      * k_air / D_h;
    printf("Nusselt number:     %.1f\n", Nu_air);
    printf("h (convection):     %.1f W/(m²·K)\n\n", h_forced);

    /* Fin efficiency */
    double m_fin = sqrt(h_forced * 2.0 * (W_hs + t_fin)
                        / (k_al * W_hs * t_fin));
    double eta_single = fin_efficiency(m_fin, L_fin);
    double eta_overall = 1.0 - (A_fins_total / A_total) * (1.0 - eta_single);

    printf("--- Fin Performance ---\n");
    printf("Fin parameter m:    %.1f 1/m\n", m_fin);
    printf("Single fin η:       %.3f\n", eta_single);
    printf("Overall surface η:  %.3f\n\n", eta_overall);

    /* Heat sink thermal resistance */
    double R_hs = heat_sink_thermal_resistance(k_al, A_base,
                                                eta_overall, h_forced,
                                                A_total);
    printf("--- Total Performance ---\n");
    printf("R_heat_sink:        %.4f K/W\n", R_hs);

    /* Junction temperature */
    double T_j = junction_temperature(T_amb, q_cpu, R_jc, R_cs, R_hs);
    printf("T_junction:         %.1f°C\n\n", T_j - 273.15);

    /* Required fan airflow */
    double deltaT_air = 15.0; /* allowable air temperature rise */
    double Q_air = required_airflow_cooling(q_cpu, rho_air, cp_air,
                                             deltaT_air);
    printf("Required airflow:   %.4f m³/s (%.1f CFM)\n",
           Q_air, Q_air * 2118.88);

    /* Fan affinity: if current fan at 2000 RPM, what at 3000? */
    double Q_new, H_new, P_new;
    fan_affinity_laws(1.5, Q_air, 50.0, 2.0, &Q_new, &H_new, &P_new);
    printf("Fan at 1.5× speed:  Q=%.4f m³/s, P=%.1f W\n\n", Q_new, P_new);

    /* Check if design meets requirements */
    if (T_j <= T_j_target) {
        printf("✓ Design PASS: Tj=%.1f°C ≤ %.1f°C target\n",
               T_j - 273.15, T_j_target - 273.15);
    } else {
        printf("✗ Design FAIL: Tj=%.1f°C > %.1f°C target\n",
               T_j - 273.15, T_j_target - 273.15);
        printf("  Recommendation: increase fin count or airflow.\n");
    }

    return 0;
}
