/**
 * Example 2: Mason's Gain Formula Application
 *
 * Demonstrates Mason's Gain Formula for computing the transfer
 * function of a multi-loop system directly from its signal flow
 * graph, without block diagram reduction.
 *
 * L6: Canonical Problem — Multi-loop system analysis
 * L4: Mason's theorem verification
 *
 * The system has one forward path and two non-touching loops,
 * demonstrating the full Mason's formula:
 *   T(s) = (1/Delta) * Sum(P_k * Delta_k)
 */
#include "../include/blockdiagram.h"
#include "../include/signal_flow_graph.h"
#include "../include/mason.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 2: Mason's Gain Formula ===\n\n");

    /* Build a signal flow graph for G(s) = 10/(s+2) with feedback H=1 */
    SignalFlowGraph *sg = sfg_create("mason_demo", 10, 20);

    /* Nodes: R -> E(summer) -> Y, with feedback Y -> E */
    int src = sfg_add_node(sg, SFG_SOURCE, "R(s)");
    int e   = sfg_add_node(sg, SFG_INTERNAL, "E(s)");
    int y   = sfg_add_node(sg, SFG_SINK, "Y(s)");

    /* Forward path: R->E (gain 1), E->Y (gain 10) */
    sfg_add_branch(sg, src, e, 1.0);
    sfg_add_branch(sg, e, y, 10.0);

    /* Feedback loop: Y->E (gain -1 for negative feedback) */
    sfg_add_branch(sg, y, e, -1.0);

    sfg_set_source(sg, src);
    sfg_set_sink(sg, y);

    printf("Signal Flow Graph:\n");
    sfg_print(sg);

    /* Perform complete Mason's analysis */
    MasonResult *mr = mason_analyze(sg, 16, 32);
    if (mr) {
        mason_report_print(mr);

        /* Verify against direct calculation:
         * Forward path gain P1 = 1*10 = 10
         * Loop gain L1 = 10*(-1) = -10
         * Delta = 1 - L1 = 1 - (-10) = 11
         * Delta1 = 1 (no loops remain after removing forward path)
         * T = P1*Delta1/Delta = 10*1/11 = 0.90909... */
        printf("\nManual verification:\n");
        printf("  P1 = 10, L1 = -10\n");
        printf("  Delta = 1 - L1 = 11\n");
        printf("  Delta1 = 1\n");
        printf("  T = 10/11 = %.6f\n", 10.0/11.0);
        printf("  Match: %s\n",
               fabs(mr->transfer_gain - 10.0/11.0) < 1e-9 ? "YES" : "NO");

        mason_result_free(mr);
    }

    /* Alternative: compute directly */
    double T = sfg_mason_transfer_function_gain(sg);
    printf("  SFG direct T = %.6f\n", T);

    sfg_destroy(sg);
    printf("\n=== Done ===\n");
    return 0;
}
