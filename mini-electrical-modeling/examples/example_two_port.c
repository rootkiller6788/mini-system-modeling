/**
 * @file example_two_port.c
 * @brief End-to-end example: Two-port network analysis of a cascaded amplifier stage
 *
 * Demonstrates: Z/Y/S/ABCD parameter conversions, cascade interconnection,
 * S-parameter stability analysis, maximum gain, VSWR, matching.
 *
 * Problem: Analyze a two-stage RF amplifier cascade. Each stage is characterized
 * by its S-parameters at 2.4 GHz (WiFi band). Determine overall gain,
 * stability, and input/output match.
 *
 * Reference: Pozar (2012) Chapter 12 — Microwave Amplifier Design
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "two_port_network.h"

int main(void)
{
    printf("=== Two-Port Network Analysis: RF Amplifier Cascade ===\n\n");

    double freq_ghz = 2.4;
    double z0 = 50.0;
    printf("Frequency: %.1f GHz, Reference Z0: %.0f ohm\n\n", freq_ghz, z0);

    /* Stage 1 S-parameters (typical LNA at 2.4 GHz) */
    SParameters stage1;
    stage1.s11 = 0.3 * cexp(I * (-1.2));    /* input reflection */
    stage1.s12 = 0.02 * cexp(I * 0.8);       /* reverse isolation */
    stage1.s21 = 3.16 * cexp(I * (-2.5));    /* forward gain (~10 dB) */
    stage1.s22 = 0.25 * cexp(I * (-0.9));    /* output reflection */
    stage1.z0 = z0;

    printf("=== Stage 1 (LNA) S-Parameters ===\n");
    printf("  S11: %.3f angle(%.1f deg)\n", cabs(stage1.s11), carg(stage1.s11)*180/M_PI);
    printf("  S12: %.3f angle(%.1f deg)\n", cabs(stage1.s12), carg(stage1.s12)*180/M_PI);
    printf("  S21: %.3f angle(%.1f deg)  -> Gain: %.1f dB\n",
           cabs(stage1.s21), carg(stage1.s21)*180/M_PI,
           20*log10(cabs(stage1.s21)));
    printf("  S22: %.3f angle(%.1f deg)\n", cabs(stage1.s22), carg(stage1.s22)*180/M_PI);

    /* Stage 2 S-parameters (typical driver amp) */
    SParameters stage2;
    stage2.s11 = 0.35 * cexp(I * (-1.5));
    stage2.s12 = 0.03 * cexp(I * 1.2);
    stage2.s21 = 2.0 * cexp(I * (-1.8));     /* ~6 dB gain */
    stage2.s22 = 0.28 * cexp(I * (-1.1));
    stage2.z0 = z0;

    printf("\n=== Stage 2 (Driver) S-Parameters ===\n");
    printf("  S11: %.3f, S12: %.3f, S21: %.3f (%.1f dB), S22: %.3f\n",
           cabs(stage2.s11), cabs(stage2.s12),
           cabs(stage2.s21), 20*log10(cabs(stage2.s21)),
           cabs(stage2.s22));

    /* Stability analysis of each stage */
    printf("\n=== Stability Analysis ===\n");
    double K1 = sp_rollett_stability(&stage1);
    double K2 = sp_rollett_stability(&stage2);
    printf("  Stage 1 Rollett K = %.3f %s\n", K1,
           K1 > 1.0 ? "(Unconditionally Stable)" : "(Potentially Unstable)");
    printf("  Stage 2 Rollett K = %.3f %s\n", K2,
           K2 > 1.0 ? "(Unconditionally Stable)" : "(Potentially Unstable)");

    /* Convert S to Z for series interconnection computation */
    ZParameters zp1, zp2;
    sp_to_zp(&stage1, &zp1);
    sp_to_zp(&stage2, &zp2);

    /* Cascade via ABCD parameters */
    ABCDParameters abcd1, abcd2;
    zp_to_abcd(&zp1, &abcd1);
    zp_to_abcd(&zp2, &abcd2);
    ABCDParameters abcd_cascade = tp_cascade_interconnect(&abcd1, &abcd2);

    /* Convert back to S-parameters for the cascade */
    ZParameters zp_cascade;
    abcd_to_zp(&abcd_cascade, &zp_cascade);
    SParameters sp_cascade;
    zp_to_sp(&zp_cascade, z0, &sp_cascade);

    printf("\n=== Cascaded Amplifier S-Parameters ===\n");
    printf("  S11: %.3f  -> Return Loss: %.1f dB, VSWR: %.2f\n",
           cabs(sp_cascade.s11),
           -20*log10(cabs(sp_cascade.s11)),
           sp_vswr(sp_cascade.s11));
    printf("  S12: %.3f  -> Reverse Isolation: %.1f dB\n",
           cabs(sp_cascade.s12),
           -20*log10(cabs(sp_cascade.s12)));
    printf("  S21: %.3f  -> Total Gain: %.1f dB\n",
           cabs(sp_cascade.s21),
           20*log10(cabs(sp_cascade.s21)));
    printf("  S22: %.3f  -> Output Return Loss: %.1f dB, VSWR: %.2f\n",
           cabs(sp_cascade.s22),
           -20*log10(cabs(sp_cascade.s22)),
           sp_vswr(sp_cascade.s22));

    /* Cascade stability */
    double K_cascade = sp_rollett_stability(&sp_cascade);
    printf("\n  Cascade Rollett K = %.3f %s\n", K_cascade,
           K_cascade > 1.0 ? "(Unconditionally Stable)" : "(Potentially Unstable)");

    /* Maximum available gain */
    double mag = sp_max_available_gain(&sp_cascade);
    printf("  Max Available Gain: %.1f dB\n",
           (mag > 0) ? 10*log10(mag) : 0.0);

    /* Transducer gain with 50 ohm terminations */
    double gt = sp_transducer_gain(&sp_cascade, 50.0, 50.0);
    printf("  Transducer Gain (50 ohm): %.1f dB\n",
           (gt > 0) ? 10*log10(gt) : -999.0);

    /* Check reciprocity */
    int is_recip = tp_is_reciprocal_z(&zp_cascade, 1e-6);
    printf("\n  Reciprocity: %s\n", is_recip ? "YES" : "NO (active device)");

    /* Passivity check */
    int is_passive = tp_is_passive_z(&zp_cascade, 1e-6);
    printf("  Passivity: %s\n", is_passive ? "YES" : "NO (amplifier has gain)");

    /* Mixed-mode S-parameters for differential analysis */
    MixedModeSParameters mm;
    sp_to_mixed_mode(&sp_cascade, &mm);
    printf("\n=== Mixed-Mode S-Parameters (Differential) ===\n");
    printf("  Sdd11: %.3f, Sdd21: %.3f\n",
           cabs(mm.sdd11), cabs(mm.sdd21));
    printf("  Common-mode rejection: Sdc21 = %.6f (ideally 0)\n",
           cabs(mm.sdc21));

    /* Pi-network equivalent */
    printf("\n=== Pi-Network Equivalent ===\n");
    double complex Y1 = I * 2.0 * M_PI * 2.4e9 * 0.5e-12; /* 0.5 pF */
    double complex Y2 = 1.0/1000.0 + I*0.0;               /* 1 kohm */
    double complex Y3 = I * 2.0 * M_PI * 2.4e9 * 0.3e-12; /* 0.3 pF */
    ABCDParameters abcd_pi = tp_abcd_pi_network(Y1, Y2, Y3);

    ZParameters zp_pi;
    abcd_to_zp(&abcd_pi, &zp_pi);
    SParameters sp_pi;
    zp_to_sp(&zp_pi, z0, &sp_pi);

    printf("  Pi-network S21 gain at 2.4 GHz: %.1f dB\n",
           20*log10(cabs(sp_pi.s21)));

    printf("\n=== Two-Port Analysis Complete ===\n");
    return 0;
}
