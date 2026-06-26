/**
 * Example 2: Shell-and-Tube Heat Exchanger Design
 *
 * Designs a counter-flow shell-and-tube heat exchanger using
 * the ε-NTU method. Compares counter-flow, parallel-flow, and
 * shell-and-tube (1 shell pass) configurations.
 *
 * Problem: Cool 0.5 kg/s of engine oil from 120°C to 80°C
 * using 0.3 kg/s of water entering at 25°C. Determine
 * required heat transfer area and compare arrangements.
 */

#include <stdio.h>
#include <math.h>
#include "../include/fluid_thermal_core.h"
#include "../include/heat_transfer.h"

int main(void)
{
    printf("=== Example 2: Heat Exchanger Design (ε-NTU Method) ===\n\n");

    /* Given */
    double Th_in = 393.15;     /* K — oil inlet (120°C) */
    double Th_out_req = 353.15; /* K — desired oil outlet (80°C) */
    double Tc_in = 298.15;     /* K — water inlet (25°C) */

    double m_dot_h = 0.5;      /* kg/s — oil mass flow */
    double cp_h = 1900.0;      /* J/(kg·K) — oil specific heat */
    double m_dot_c = 0.3;      /* kg/s — water mass flow */
    double cp_c = 4179.0;      /* J/(kg·K) — water specific heat */

    /* Capacity rates */
    double Ch = heat_capacity_rate(m_dot_h, cp_h);
    double Cc = heat_capacity_rate(m_dot_c, cp_c);
    double C_min = (Ch < Cc) ? Ch : Cc;
    double C_max = (Ch > Cc) ? Ch : Cc;
    double Cr = C_min / C_max;

    printf("Hot fluid (oil):  ṁh=%.2f kg/s, cp=%.0f J/(kg·K)\n",
           m_dot_h, cp_h);
    printf("Cold fluid (H2O): ṁc=%.2f kg/s, cp=%.0f J/(kg·K)\n",
           m_dot_c, cp_c);
    printf("Ch = %.0f W/K, Cc = %.0f W/K\n", Ch, Cc);
    printf("Cmin = %.0f W/K, Cr = %.3f\n\n", C_min, Cr);

    /* Required heat duty */
    double q_req = Ch * (Th_in - Th_out_req);
    printf("Required heat duty: q = %.2f kW\n\n", q_req / 1000.0);

    /* Maximum possible heat transfer */
    double q_max = C_min * (Th_in - Tc_in);
    double epsilon_req = q_req / q_max;
    printf("q_max = %.2f kW\n", q_max / 1000.0);
    printf("Required effectiveness: ε = %.3f\n\n", epsilon_req);

    /* Compare arrangements */
    const char *arrangements[] = {"counter", "parallel", "shell_tube_1"};
    const char *arr_names[] = {"Counter-flow", "Parallel-flow",
                                "1 Shell / 2 Tube"};

    double U_est = 500.0; /* W/(m²·K) — estimated overall U */

    printf("--- Arrangement Comparison (U=%.0f W/m²K) ---\n\n", U_est);

    for (int i = 0; i < 3; i++) {
        /* Compute NTU needed */
        double NTU;
        if (i == 0) {
            NTU = ntu_from_effectiveness_counter(epsilon_req, Cr);
        } else if (i == 1) {
            NTU = ntu_from_effectiveness_parallel(epsilon_req, Cr);
        } else {
            /* For shell-and-tube, iteratively solve ε-NTU relation.
             * Using the counter relation as approximation for shell-tube
             * (exact for n shells in series). */
            NTU = ntu_from_effectiveness_counter(epsilon_req, Cr);
        }

        double A_req = NTU * C_min / U_est;

        /* Verify effectiveness */
        double eps_actual;
        if (i < 2) {
            eps_actual = heat_exchanger_effectiveness(
                NTU, Cr, arrangements[i]);
        } else {
            eps_actual = heat_exchanger_effectiveness(
                NTU, Cr, "shell_tube_1");
        }

        /* Outlet temperatures */
        double Th_out, Tc_out;
        heat_exchanger_outlet_temps(Th_in, Tc_in, Ch, Cc,
                                     eps_actual, &Th_out, &Tc_out);

        printf("  %s:\n", arr_names[i]);
        printf("    NTU required:  %.2f\n", NTU);
        printf("    Area required: %.2f m²\n", A_req);
        printf("    Effectiveness: %.4f\n", eps_actual);
        printf("    Th_out: %.1f°C, Tc_out: %.1f°C\n\n",
               Th_out - 273.15, Tc_out - 273.15);
    }

    /* LMTD method for the counter-flow case */
    double dT1_c = Th_in - Tc_in;  /* (Th,i - Tc,o) for counter — unknown yet */
    /* We know: Th,o = Th_out_req, need Tc,o from energy balance */
    double Tc_out_c = Tc_in + q_req / Cc;
    dT1_c = Th_in - Tc_out_c;
    double dT2_c = Th_out_req - Tc_in;

    double dT_lm_c = lmtd(dT1_c, dT2_c);
    double A_lmtd = q_req / (U_est * dT_lm_c);

    printf("--- LMTD Method (counter-flow) ---\n\n");
    printf("  ΔT1 = %.1f K, ΔT2 = %.1f K\n", dT1_c, dT2_c);
    printf("  LMTD = %.1f K\n", dT_lm_c);
    printf("  A_LMTD = %.2f m²\n\n", A_lmtd);

    printf("Conclusion: Counter-flow requires the least area (%.2f m²).\n",
           ntu_from_effectiveness_counter(epsilon_req, Cr) * C_min / U_est);

    return 0;
}
