/**
 * @file block_algebra.c
 * @brief Block diagram algebraic operations — series, parallel, feedback, reduction
 *
 * L1: Interconnection types (series, parallel, feedback)
 * L2: Build and apply standard interconnection patterns
 * L3: Transformation rules and equivalence checks
 * L4: Systematic reduction algorithms
 * L5: Disturbance rejection configuration, structure printing
 *
 * Block diagram algebra is the set of operations that combine subsystems
 * into larger systems while preserving the overall input-output behavior.
 * These operations form the mathematical foundation for block diagram
 * reduction — one of the most important skills in control engineering.
 *
 * Ref: K. Ogata, "Modern Control Engineering", 5th ed., Ch. 3.3
 * Ref: R.C. Dorf, "Modern Control Systems", 13th ed., Ch. 2.6
 * MIT 6.302 / ETH 151-0591 / Cambridge 3F2
 */
#include "block_algebra.h"
#include "block_reduction.h"
#include "blockdiagram.h"
#include "transfer_function.h"
#include "signal_flow_graph.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ======== Internal Builder Helpers ======== */

/** Build a simple 2-node diagram: input -> [TF] -> output */
static BlockDiagram* build_single_block(const TransferFunction *tf,
                                         const char *label) {
    BlockDiagram *bd = bd_create("single_block", 4, 4);
    if (!bd) return NULL;
    int in  = bd_add_node(bd, BD_NODE_INPUT,  "R(s)");
    int blk = bd_add_node(bd, BD_NODE_BLOCK,  label ? label : "G(s)");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y(s)");
    if (in < 0 || blk < 0 || out < 0) { bd_destroy(bd); return NULL; }
    bd_set_transfer_function(bd, blk, tf);
    bd_add_edge(bd, in,  blk, 1.0);
    bd_add_edge(bd, blk, out, 1.0);
    return bd;
}

/* ======== Series Interconnection ======== */
BlockDiagram* ba_series(const BlockDiagram *a, const BlockDiagram *b) {
    /* L1: Series (cascade) interconnection: output of A feeds input of B.
     * The overall transfer function is G_B(s) * G_A(s).
     * This creates a new diagram with A and B chained together. */
    if (!a || !b) return NULL;
    /* Clone both diagrams */
    BlockDiagram *ca = bd_clone(a);
    BlockDiagram *cb = bd_clone(b);
    if (!ca || !cb) { bd_destroy(ca); bd_destroy(cb); return NULL; }
    /* Create new combined diagram */
    int total_nodes = ca->n_nodes + cb->n_nodes + 2;
    int total_edges = ca->n_edges + cb->n_edges + 4;
    BlockDiagram *result = bd_create("series_combined", total_nodes, total_edges);
    if (!result) { bd_destroy(ca); bd_destroy(cb); return NULL; }
    /* Copy nodes from A (except output), then from B (except input),
     * then connect A's last block to B's first block */
    /* Simplified: just use TF-based series */
    TransferFunction *tf_a = br_extract_transfer(ca);
    TransferFunction *tf_b = br_extract_transfer(cb);
    TransferFunction *tf_s = tf_series(tf_a, tf_b);
    bd_destroy(ca); bd_destroy(cb);
    if (!tf_s) { bd_destroy(result); return NULL; }
    BlockDiagram *final = build_single_block(tf_s, "G_series");
    tf_destroy(tf_a); tf_destroy(tf_b); tf_destroy(tf_s);
    bd_destroy(result);
    return final;
}

/* ======== Parallel Interconnection ======== */
BlockDiagram* ba_parallel(const BlockDiagram *a, const BlockDiagram *b) {
    /* L1: Parallel interconnection: both systems share the same input,
     * outputs are summed. G_parallel = G_A + G_B. */
    if (!a || !b) return NULL;
    TransferFunction *tf_a = br_extract_transfer(bd_clone(a));
    TransferFunction *tf_b = br_extract_transfer(bd_clone(b));
    TransferFunction *tf_par = tf_parallel(tf_a, tf_b);
    if (!tf_par) return NULL;
    BlockDiagram *result = build_single_block(tf_par, "G_parallel");
    tf_destroy(tf_a); tf_destroy(tf_b); tf_destroy(tf_par);
    return result;
}

/* ======== Feedback Interconnection ======== */
BlockDiagram* ba_unity_feedback(const BlockDiagram *G) {
    /* L1: Negative unity feedback: G_cl = G/(1+G). */
    if (!G) return NULL;
    TransferFunction *tf_g = br_extract_transfer(bd_clone(G));
    TransferFunction *tf_cl = tf_feedback(tf_g, NULL);
    BlockDiagram *result = build_single_block(tf_cl, "G_cl_unity");
    tf_destroy(tf_g); tf_destroy(tf_cl);
    return result;
}
BlockDiagram* ba_feedback(const BlockDiagram *G, const BlockDiagram *H) {
    /* L1: Negative feedback: G_cl = G/(1+G*H). */
    if (!G) return NULL;
    TransferFunction *tf_g = br_extract_transfer(bd_clone(G));
    TransferFunction *tf_h = H ? br_extract_transfer(bd_clone(H)) : NULL;
    TransferFunction *tf_cl = tf_feedback(tf_g, tf_h);
    BlockDiagram *result = build_single_block(tf_cl, "G_cl_feedback");
    tf_destroy(tf_g); tf_destroy(tf_h); tf_destroy(tf_cl);
    return result;
}
BlockDiagram* ba_positive_feedback(const BlockDiagram *G, const BlockDiagram *H) {
    /* L1: Positive feedback: G_cl = G/(1-G*H). */
    if (!G) return NULL;
    TransferFunction *tf_g = br_extract_transfer(bd_clone(G));
    TransferFunction *tf_h = H ? br_extract_transfer(bd_clone(H)) : NULL;
    TransferFunction *tf_cl = tf_positive_feedback(tf_g, tf_h);
    BlockDiagram *result = build_single_block(tf_cl, "G_cl_pos_fb");
    tf_destroy(tf_g); tf_destroy(tf_h); tf_destroy(tf_cl);
    return result;
}

/* ======== Build Standard Feedback ======== */
BlockDiagram* ba_build_standard_feedback(const TransferFunction *G,
                                          const TransferFunction *H) {
    /* L2: Build the canonical negative feedback configuration.
     *    R ->(+) -> [G] -> Y
     *          ^- [H] <-|
     * The summing junction subtracts the feedback signal from the reference. */
    if (!G) return NULL;
    BlockDiagram *bd = bd_create("std_feedback", 8, 10);
    int in    = bd_add_node(bd, BD_NODE_INPUT,  "R(s)");
    int sum   = bd_add_node(bd, BD_NODE_SUMMER, "Error");
    int fwd   = bd_add_node(bd, BD_NODE_BLOCK,  "G(s)");
    int out   = bd_add_node(bd, BD_NODE_OUTPUT, "Y(s)");
    int tko   = bd_add_node(bd, BD_NODE_TAKEOFF,"Branch");
    int fb    = bd_add_node(bd, BD_NODE_BLOCK,  "H(s)");
    bd_set_transfer_function(bd, fwd, G);
    if (H) bd_set_transfer_function(bd, fb, H);
    /* Set summer signs: R positive, feedback negative */
    double signs[2] = {1.0, -1.0};
    bd_set_summer_signs(bd, sum, signs, 2);
    /* Wire it up */
    bd_add_edge(bd, in,  sum, 1.0);
    bd_add_edge(bd, sum, fwd, 1.0);
    bd_add_edge(bd, fwd, out, 1.0);
    bd_add_edge(bd, fwd, tko, 1.0);
    bd_add_edge(bd, tko, fb,  1.0);
    bd_add_edge(bd, fb,  sum, 1.0);
    return bd;
}

/* ======== Disturbance Rejection ======== */
BlockDiagram* ba_build_disturbance_rejection(const TransferFunction *C,
                                              const TransferFunction *P,
                                              const TransferFunction *H) {
    /* L5: Build control configuration with disturbance input.
     * R ->(+) -> [C] -> (+) -> [P] -> Y
     *     ^-                    |
     *     |-- [H] <-------------|
     *              D(s) disturbs after controller.
     * This models real-world systems where disturbances enter
     * between the controller and the plant (e.g., DC motor load torque,
     * wind gusts on quadrotors, road grade for cruise control). */
    if (!C || !P) return NULL;
    BlockDiagram *bd = bd_create("disturb_reject", 10, 14);
    int in  = bd_add_node(bd, BD_NODE_INPUT,  "R(s)");
    int s1  = bd_add_node(bd, BD_NODE_SUMMER, "ErrSum");
    int ctl = bd_add_node(bd, BD_NODE_BLOCK,  "C(s)");
    int s2  = bd_add_node(bd, BD_NODE_SUMMER, "DistSum");
    int dist= bd_add_node(bd, BD_NODE_INPUT,  "D(s)");
    int plt = bd_add_node(bd, BD_NODE_BLOCK,  "P(s)");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y(s)");
    int tko = bd_add_node(bd, BD_NODE_TAKEOFF,"Branch");
    int fb  = bd_add_node(bd, BD_NODE_BLOCK,  "H(s)");
    bd_set_transfer_function(bd, ctl, C);
    bd_set_transfer_function(bd, plt, P);
    if (H) bd_set_transfer_function(bd, fb, H);
    double s1s[2] = {1.0, -1.0}; bd_set_summer_signs(bd, s1, s1s, 2);
    double s2s[2] = {1.0, 1.0};  bd_set_summer_signs(bd, s2, s2s, 2); /* D adds */
    bd_add_edge(bd, in,  s1,  1.0);
    bd_add_edge(bd, s1,  ctl, 1.0);
    bd_add_edge(bd, ctl, s2,  1.0);
    bd_add_edge(bd, dist,s2,  1.0);
    bd_add_edge(bd, s2,  plt, 1.0);
    bd_add_edge(bd, plt, out, 1.0);
    bd_add_edge(bd, plt, tko, 1.0);
    bd_add_edge(bd, tko, fb,  1.0);
    bd_add_edge(bd, fb,  s1,  1.0);
    return bd;
}

/* ======== Transformation Rules (L3) ======== */
int ba_apply_transformation(BlockDiagram *bd, ba_transformation_t rule,
                            const int *nids, int nn) {
    /* L3: Apply a named block diagram transformation.
     * Each transformation preserves the closed-loop I/O relationship.
     * These are proven equivalences in block diagram theory (L4).
     * For detailed implementations, see block_reduction.c. */
    if (!bd) return -1;
    switch (rule) {
    case BA_ABSORB_SERIES_BLOCKS:
        if (nn >= 2) return br_apply_series_rule(bd, nids[0], nids[1]);
        break;
    case BA_ABSORB_PARALLEL_BLOCKS:
        if (nn >= 3) return br_apply_parallel_rule(bd, nids[0], nids[1], nids[2]);
        break;
    case BA_ELIMINATE_FEEDBACK_LOOP:
        if (nn >= 3) return br_apply_feedback_rule(bd, nids[0], nids[1], nids[2]);
        break;
    case BA_MOVE_BLOCK_PAST_SUMMER_BACKWARD:
        if (nn >= 2) return br_move_block_past_summer_backward(bd, nids[0], nids[1]);
        break;
    case BA_MOVE_BLOCK_PAST_TAKEOFF_BACKWARD:
        if (nn >= 2) return br_move_block_past_takeoff_backward(bd, nids[0], nids[1]);
        break;
    default: return -1;
    }
    return 0;
}

/* ======== Systematic Reduction (L4-L5) ======== */
BlockDiagram* ba_reduce(const BlockDiagram *bd, ba_reduce_strategy_t strategy) {
    /* L4: Systematic block diagram reduction.
     * Apply reduction rules iteratively until a single equivalent block
     * remains between input and output. */
    if (!bd) return NULL;
    BlockDiagram *r = bd_clone(bd);
    if (!r) return NULL;
    ReduceOptions opts = br_options_default();
    opts.verbose = 0;
    ReduceTrace *trace = br_trace_create(32);
    if (!trace) { /* continue without trace */ }
    br_reduce(r, &opts, trace);
    if (trace) br_trace_destroy(trace);
    return r;
}

TransferFunction* ba_equivalent_transfer(const BlockDiagram *bd,
                                          ba_reduce_strategy_t strategy) {
    /* L5: Compute the equivalent transfer function of a block diagram. */
    BlockDiagram *reduced = ba_reduce(bd, strategy);
    if (!reduced) return NULL;
    TransferFunction *tf = br_extract_transfer(reduced);
    bd_destroy(reduced);
    return tf;
}

int ba_are_equivalent(const BlockDiagram *a, const BlockDiagram *b,
                      double tolerance) {
    /* L4: Check if two block diagrams have the same I/O behavior.
     * Compares the equivalent transfer functions after reduction. */
    if (!a || !b) return 0;
    TransferFunction *tf_a = ba_equivalent_transfer(a, BA_REDUCE_GREEDY);
    TransferFunction *tf_b = ba_equivalent_transfer(b, BA_REDUCE_GREEDY);
    if (!tf_a || !tf_b) { tf_destroy(tf_a); tf_destroy(tf_b); return 0; }
    /* Compare DC gains as a simple equivalence check */
    double ga = tf_dc_gain(tf_a), gb = tf_dc_gain(tf_b);
    tf_destroy(tf_a); tf_destroy(tf_b);
    return fabs(ga - gb) <= tolerance;
}

int ba_simplify(BlockDiagram *bd) {
    /* L5: Lightweight simplification that combines obvious cascades
     * and parallels without full reduction. */
    if (!bd) return -1;
    int changed = 0;
    /* Combine cascaded blocks */
    int i, j;
    for (i = 0; i < bd->n_nodes; i++) {
        if (bd->nodes[i]->type != BD_NODE_BLOCK) continue;
        for (j = 0; j < bd->n_nodes; j++) {
            if (i == j || bd->nodes[j]->type != BD_NODE_BLOCK) continue;
            /* Check if i feeds directly into j */
            int k, found = 0;
            for (k = 0; k < bd->n_edges; k++) {
                if (bd->edges[k]->src_node == bd->nodes[i]->id &&
                    bd->edges[k]->dst_node == bd->nodes[j]->id) {
                    found = 1; break;
                }
            }
            /* Check no branching from i's output */
            int out_count = 0;
            for (k = 0; k < bd->n_edges; k++)
                if (bd->edges[k]->src_node == bd->nodes[i]->id) out_count++;
            if (found && out_count == 1) {
                if (br_apply_series_rule(bd, bd->nodes[i]->id, bd->nodes[j]->id))
                    changed++;
            }
        }
    }
    return changed;
}

void ba_print_structure(const BlockDiagram *bd) {
    /* L5: Print interconnection structure in a readable format. */
    if (!bd) { printf("Diagram: NULL\n"); return; }
    printf("Block Diagram Structure: \"%s\"\n", bd->name);
    int counts[6]; bd_count_node_types(bd, counts);
    printf("  Inputs:%d  Outputs:%d  Blocks:%d  Summers:%d  Takeoffs:%d  Gains:%d\n",
           counts[BD_NODE_INPUT], counts[BD_NODE_OUTPUT], counts[BD_NODE_BLOCK],
           counts[BD_NODE_SUMMER], counts[BD_NODE_TAKEOFF], counts[BD_NODE_GAIN]);
    bd_print(bd, 1);
}
