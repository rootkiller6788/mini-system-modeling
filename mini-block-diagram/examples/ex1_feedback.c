/**
 * Example 1: Feedback System Block Diagram Reduction
 *
 * Demonstrates standard negative feedback configuration:
 *   R(s) -> (+) -> [G(s)] -> Y(s)
 *            ^- [H(s)] <-|
 *
 * L6: Canonical Problem — Feedback loop reduction
 * L7: Application — DC motor speed control
 *
 * This example solves for the closed-loop transfer function
 * T(s) = G(s) / (1 + G(s)*H(s)) using both direct formula
 * and block diagram reduction.
 */
#include "../include/blockdiagram.h"
#include "../include/transfer_function.h"
#include "../include/block_algebra.h"
#include "../include/block_reduction.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 1: Feedback System Reduction ===\n\n");

    /* Define plant G(s) = 10/(s+2) and sensor H(s) = 1/(s+5) */
    double nG[] = {10.0};           double dG[] = {2.0, 1.0};
    double nH[] = {1.0};            double dH[] = {5.0, 1.0};
    TransferFunction *G = tf_create(nG, 0, dG, 1);
    TransferFunction *H = tf_create(nH, 0, dH, 1);

    printf("Plant G(s) = "); tf_print(G);
    printf("Sensor H(s) = "); tf_print(H);

    /* Method 1: Direct formula */
    TransferFunction *T_direct = tf_feedback(G, H);
    printf("\nMethod 1 (direct formula): T(s) = G/(1+GH)\n");
    printf("  T(s) = "); tf_print(T_direct);

    /* Method 2: Build block diagram and reduce */
    BlockDiagram *bd = ba_build_standard_feedback(G, H);
    printf("\nMethod 2 (block diagram reduction):\n");
    printf("  Original diagram: "); bd_print(bd, 0);
    ReduceOptions opts = br_options_default();
    br_reduce(bd, &opts, NULL);
    TransferFunction *T_reduced = br_extract_transfer(bd);
    printf("  Reduced T(s) = "); tf_print(T_reduced);

    /* Verify equivalence */
    double dc1 = tf_dc_gain(T_direct);
    double dc2 = tf_dc_gain(T_reduced);
    printf("\nDC gain (direct):  %.4f\n", dc1);
    printf("DC gain (reduced): %.4f\n", dc2);
    printf("Match: %s\n", fabs(dc1 - dc2) < 1e-9 ? "YES" : "NO");

    /* Cleanup */
    tf_destroy(G); tf_destroy(H);
    tf_destroy(T_direct); tf_destroy(T_reduced);
    bd_destroy(bd);

    printf("\n=== Done ===\n");
    return 0;
}
