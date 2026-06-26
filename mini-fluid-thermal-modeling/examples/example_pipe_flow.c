/**
 * Example 1: Pipe Flow Analysis — Water Distribution System
 *
 * Designs a water supply pipe for a 100 m horizontal run,
 * computing Reynolds number, friction factor, head loss,
 * and pressure drop. Demonstrates Darcy-Weisbach with
 * Colebrook friction factor for turbulent flow.
 *
 * Problem: Supply 0.001 m³/s of water through a 100 m pipe
 * with 50 mm diameter. Determine head loss, pressure drop,
 * and required pump power.
 */

#include <stdio.h>
#include <math.h>
#include "../include/fluid_thermal_core.h"
#include "../include/fluid_mechanics.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
    printf("=== Example 1: Pipe Flow Analysis ===\n\n");

    /* Input parameters */
    double Q = 0.001;        /* m³/s — domestic water flow */
    double D = 0.05;         /* m — pipe diameter */
    double L = 100.0;        /* m — pipe length */
    double eps = 0.00015;    /* m — roughness (commercial steel) */
    double rho = 997.0;      /* kg/m³ — water density at 25°C */
    double mu = 8.55e-4;     /* Pa·s — water viscosity */
    double g = 9.81;         /* m/s² */

    /* Compute flow parameters */
    double A = M_PI * D * D / 4.0;
    double u = Q / A;
    double Re = compute_reynolds(rho, u, D, mu);
    FlowRegime regime = classify_flow_regime(Re);
    double eps_D = eps / D;

    printf("Pipe diameter:       %.3f m\n", D);
    printf("Flow area:           %.6f m²\n", A);
    printf("Flow velocity:       %.2f m/s\n", u);
    printf("Reynolds number:     %.0f\n", Re);
    printf("Relative roughness:  %.5f\n", eps_D);

    const char *regime_str[] = {"Creeping", "Laminar", "Transitional",
                                 "Turbulent", "Hypersonic"};
    printf("Flow regime:         %s\n", regime_str[regime]);

    /* Compute friction factor */
    double f;
    if (Re < 2300.0) {
        f = darcy_friction_laminar(Re);
        printf("Friction method:     Darcy (laminar, f=64/Re)\n");
    } else {
        f = colebrook_friction_factor(Re, eps_D, 30, 1e-10);
        printf("Friction method:     Colebrook-White (iterative)\n");
    }
    printf("Darcy friction f:    %.5f\n", f);

    /* Head loss */
    double hL = darcy_weisbach_head_loss(f, L, D, u, g);
    printf("Head loss:           %.2f m\n", hL);

    /* Pressure drop */
    double dP = darcy_weisbach_pressure_drop(f, L, D, rho, u);
    printf("Pressure drop:       %.2f Pa (%.3f bar)\n", dP, dP/1e5);

    /* Pump power required */
    double P_hyd = rho * g * Q * hL;
    printf("Hydraulic power:     %.2f W\n", P_hyd);

    /* Include minor losses */
    double K_entrance = minor_loss_coefficient("pipe_entrance");
    double K_exit = minor_loss_coefficient("pipe_exit");
    double K_elbow = minor_loss_coefficient("elbow_90_standard");
    double hL_minor = (K_entrance + K_exit + 2.0 * K_elbow)
                      * (u * u / (2.0 * g));
    printf("\nWith minor losses (entrance + exit + 2 elbows):\n");
    printf("Minor head loss:     %.2f m\n", hL_minor);
    printf("Total head loss:     %.2f m\n", hL + hL_minor);

    double P_total = rho * g * Q * (hL + hL_minor);
    printf("Total pump power:    %.2f W\n", P_total);

    printf("\nConclusion: %.2f W pump needed for this water supply.\n",
           P_total);

    return 0;
}
