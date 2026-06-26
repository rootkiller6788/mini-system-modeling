/**
 * Example 3: DC Motor Speed Control Block Diagram
 *
 * L6: Canonical Problem — DC motor modeling and control
 * L7: Application — DC motor speed control with PI controller
 *
 * Models a DC motor with:
 *   - Electrical dynamics: L*di/dt + R*i = V - Kb*w
 *   - Mechanical dynamics: J*dw/dt + B*w = Kt*i - TL
 *
 * Block diagram structure:
 *   V_ref ->(+) -> [PI] -> [1/(Ls+R)] -> [Kt] ->(+) -> [1/(Js+B)] -> w
 *            ^-                                  ^- TL              |
 *            |--- [Kb] <-------------------------------------------|
 *
 * This is a classic control problem appearing in:
 *   MIT 6.302, Stanford ENGR105, Berkeley ME132, ETH 151-0591
 */
#include "../include/blockdiagram.h"
#include "../include/transfer_function.h"
#include "../include/block_algebra.h"
#include <stdio.h>

int main(void) {
    printf("=== Example 3: DC Motor Speed Control ===\n\n");

    /* Motor parameters (typical small DC motor, e.g., maxon) */
    double R = 1.0;    /* Armature resistance [Ohm] */
    double L = 0.001;  /* Armature inductance [H] */
    double Kt = 0.05;  /* Torque constant [Nm/A] */
    double Kb = 0.05;  /* Back-EMF constant [Vs/rad] */
    double J = 1e-5;   /* Rotor inertia [kg*m^2] */
    double B = 1e-6;   /* Viscous friction [Nm*s/rad] */

    /* Electrical TF: I(s)/Va(s) = 1/(Ls + R) */
    double n_elec[] = {1.0};       double d_elec[] = {R, L};
    TransferFunction *G_elec = tf_create(n_elec, 0, d_elec, 1);

    /* Mechanical TF: w(s)/Tm(s) = 1/(Js + B) */
    double n_mech[] = {1.0};       double d_mech[] = {B, J};
    TransferFunction *G_mech = tf_create(n_mech, 0, d_mech, 1);

    /* PI Controller: C(s) = Kp + Ki/s = (Kp*s + Ki)/s */
    double Kp = 0.1, Ki = 10.0;
    double n_pi[] = {Ki, Kp};      double d_pi[] = {0.0, 1.0};
    TransferFunction *C_pi = tf_create(n_pi, 1, d_pi, 1);

    printf("DC Motor Parameters:\n");
    printf("  R=%.4f, L=%.6f, Kt=%.4f, Kb=%.4f, J=%.2e, B=%.2e\n",
           R, L, Kt, Kb, J, B);
    printf("  PI: Kp=%.2f, Ki=%.2f\n\n", Kp, Ki);

    printf("Electrical TF:  G_elec(s) = "); tf_print(G_elec);
    printf("Mechanical TF:  G_mech(s) = "); tf_print(G_mech);

    /* Build block diagram:
     * V_ref -> summer1 -> [C_pi] -> [G_elec*Kt] -> summer2 -> [G_mech] -> w
     *            ^-  [Kb] <---------------------------------|
     */
    BlockDiagram *bd = bd_create("DC_Motor", 12, 16);
    int vr   = bd_add_node(bd, BD_NODE_INPUT, "V_ref");
    int s1   = bd_add_node(bd, BD_NODE_SUMMER, "Err");
    int ctl  = bd_add_node(bd, BD_NODE_BLOCK, "C(s)");
    int el   = bd_add_node(bd, BD_NODE_BLOCK, "1/(Ls+R)");
    int kt   = bd_add_node(bd, BD_NODE_GAIN, "Kt");
    int s2   = bd_add_node(bd, BD_NODE_SUMMER, "Tm+TL");
    int tl   = bd_add_node(bd, BD_NODE_INPUT, "TL");
    int mech = bd_add_node(bd, BD_NODE_BLOCK, "1/(Js+B)");
    int w    = bd_add_node(bd, BD_NODE_OUTPUT, "w");
    int tko  = bd_add_node(bd, BD_NODE_TAKEOFF, "T");
    int kb   = bd_add_node(bd, BD_NODE_GAIN, "Kb");

    bd_set_transfer_function(bd, ctl, C_pi);
    bd_set_transfer_function(bd, el, G_elec);
    bd_set_gain(bd, kt, Kt);
    bd_set_transfer_function(bd, mech, G_mech);
    bd_set_gain(bd, kb, Kb);

    double s1s[] = {1.0, -1.0}; bd_set_summer_signs(bd, s1, s1s, 2);
    double s2s[] = {1.0, -1.0}; bd_set_summer_signs(bd, s2, s2s, 2);

    bd_add_edge(bd, vr,  s1,  1.0);
    bd_add_edge(bd, s1,  ctl, 1.0);
    bd_add_edge(bd, ctl, el,  1.0);
    bd_add_edge(bd, el,  kt,  1.0);
    bd_add_edge(bd, kt,  s2,  1.0);
    bd_add_edge(bd, tl,  s2,  1.0);
    bd_add_edge(bd, s2,  mech,1.0);
    bd_add_edge(bd, mech,w,   1.0);
    bd_add_edge(bd, mech,tko, 1.0);
    bd_add_edge(bd, tko, kb,  1.0);
    bd_add_edge(bd, kb,  s1,  1.0);

    printf("\nBlock Diagram:\n");
    bd_print(bd, 0);

    /* DC gain analysis */
    double dc_elec = tf_dc_gain(G_elec);
    double dc_mech = tf_dc_gain(G_mech);
    double dc_open = dc_elec * Kt * dc_mech;
    printf("\nDC Analysis:\n");
    printf("  G_elec(0) = %.4f, G_mech(0) = %.2f\n", dc_elec, dc_mech);
    printf("  Open-loop DC gain = %.6f\n", dc_open);
    printf("  Closed-loop with unity FB: T(0) = %.6f\n",
           dc_open / (1.0 + dc_open * Kb));

    /* Cleanup */
    tf_destroy(G_elec); tf_destroy(G_mech); tf_destroy(C_pi);
    bd_destroy(bd);

    printf("\n=== Done ===\n");
    return 0;
}
