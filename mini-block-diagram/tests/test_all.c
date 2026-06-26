#include "../include/blockdiagram.h"
#include "../include/transfer_function.h"
#include "../include/signal_flow_graph.h"
#include "../include/block_algebra.h"
#include "../include/mason.h"
#include "../include/block_reduction.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EPS 1e-9
static int tests_run = 0, tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  TEST %s... ", name);
#define OK() tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)

void test_bd_create_destroy(void) {
    TEST("bd_create/destroy")
    BlockDiagram *bd = bd_create("test", 10, 10);
    assert(bd != NULL);
    assert(bd->n_nodes == 0);
    assert(bd->n_edges == 0);
    bd_destroy(bd);
    OK();
}
void test_bd_add_nodes(void) {
    TEST("bd_add_node")
    BlockDiagram *bd = bd_create("test", 10, 10);
    int in = bd_add_node(bd, BD_NODE_INPUT, "R");
    int blk = bd_add_node(bd, BD_NODE_BLOCK, "G");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y");
    assert(in > 0 && blk > 0 && out > 0);
    assert(bd->n_nodes == 3);
    assert(bd->input_node_id == in);
    assert(bd->output_node_id == out);
    bd_destroy(bd);
    OK();
}
void test_bd_add_edges(void) {
    TEST("bd_add_edge")
    BlockDiagram *bd = bd_create("test", 10, 10);
    int in = bd_add_node(bd, BD_NODE_INPUT, "R");
    int blk = bd_add_node(bd, BD_NODE_BLOCK, "G");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y");
    int e1 = bd_add_edge(bd, in, blk, 1.0);
    int e2 = bd_add_edge(bd, blk, out, 2.5);
    assert(e1 > 0 && e2 > 0);
    assert(bd->n_edges == 2);
    bd_destroy(bd);
    OK();
}
void test_bd_find_forward_paths(void) {
    TEST("bd_find_forward_paths")
    BlockDiagram *bd = bd_create("test", 10, 10);
    int in = bd_add_node(bd, BD_NODE_INPUT, "R");
    int blk = bd_add_node(bd, BD_NODE_BLOCK, "G");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y");
    bd_add_edge(bd, in, blk, 2.0);
    bd_add_edge(bd, blk, out, 3.0);
    bd_path_t paths[10];
    int np = bd_find_forward_paths(bd, 10, paths);
    assert(np == 1);
    assert(fabs(paths[0].gain - 6.0) < EPS);
    bd_path_free(&paths[0]);
    bd_destroy(bd);
    OK();
}
void test_bd_validate(void) {
    TEST("bd_validate")
    BlockDiagram *bd = bd_create("test", 10, 10);
    int in = bd_add_node(bd, BD_NODE_INPUT, "R");
    int blk = bd_add_node(bd, BD_NODE_BLOCK, "G");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y");
    bd_add_edge(bd, in, blk, 1.0);
    bd_add_edge(bd, blk, out, 1.0);
    assert(bd_validate(bd) == 0);
    bd_destroy(bd);
    OK();
}
void test_tf_create(void) {
    TEST("tf_create")
    double n[] = {1.0}; double d[] = {1.0, 2.0, 1.0};
    TransferFunction *tf = tf_create(n, 0, d, 2);
    assert(tf != NULL);
    assert(tf->num_order == 0);
    assert(tf->den_order == 2);
    assert(tf_is_proper(tf));
    tf_destroy(tf);
    OK();
}
void test_tf_dc_gain(void) {
    TEST("tf_dc_gain")
    double n[] = {5.0}; double d[] = {1.0, 3.0};
    TransferFunction *tf = tf_create(n, 0, d, 1);
    assert(fabs(tf_dc_gain(tf) - 5.0) < EPS);
    tf_destroy(tf);
    OK();
}
void test_tf_series(void) {
    TEST("tf_series")
    double n[] = {2.0}; double d[] = {1.0};
    TransferFunction *g1 = tf_create(n, 0, d, 0);
    TransferFunction *g2 = tf_create(n, 0, d, 0);
    TransferFunction *gs = tf_series(g1, g2);
    assert(fabs(tf_dc_gain(gs) - 4.0) < EPS);
    tf_destroy(g1); tf_destroy(g2); tf_destroy(gs);
    OK();
}
void test_tf_feedback(void) {
    TEST("tf_feedback")
    double n[] = {10.0}; double d[] = {1.0};
    TransferFunction *G = tf_create(n, 0, d, 0);
    TransferFunction *cl = tf_feedback(G, NULL);
    assert(fabs(tf_dc_gain(cl) - 10.0/11.0) < EPS);
    tf_destroy(G); tf_destroy(cl);
    OK();
}
void test_tf_routh_stable(void) {
    TEST("tf_stability_routh stable")
    double d[] = {1.0, 3.0, 2.0}; double n[] = {1.0};
    TransferFunction *tf = tf_create(n, 0, d, 2);
    assert(tf_stability_routh(tf) == TF_STABLE);
    tf_destroy(tf);
    OK();
}
void test_tf_routh_unstable(void) {
    TEST("tf_stability_routh unstable")
    double d[] = {1.0, -1.0, 1.0}; double n[] = {1.0};
    TransferFunction *tf = tf_create(n, 0, d, 2);
    assert(tf_stability_routh(tf) == TF_UNSTABLE);
    tf_destroy(tf);
    OK();
}
void test_sfg_create(void) {
    TEST("sfg_create")
    SignalFlowGraph *sg = sfg_create("test", 10, 10);
    assert(sg != NULL);
    sfg_destroy(sg);
    OK();
}
void test_sfg_mason_simple(void) {
    TEST("sfg_mason_simple")
    SignalFlowGraph *sg = sfg_create("test", 10, 20);
    int src = sfg_add_node(sg, SFG_SOURCE, "X");
    int n1  = sfg_add_node(sg, SFG_INTERNAL, "N1");
    int snk = sfg_add_node(sg, SFG_SINK, "Y");
    sfg_add_branch(sg, src, n1, 2.0);
    sfg_add_branch(sg, n1, snk, 3.0);
    sfg_set_source(sg, src); sfg_set_sink(sg, snk);
    double T = sfg_mason_transfer_function_gain(sg);
    printf("(T=%.6f) ", T);
    assert(T > 0.0); /* Mason's gain should be positive for this graph */
    sfg_destroy(sg);
    OK();
}
void test_sfg_mason_feedback(void) {
    TEST("sfg_mason_feedback")
    SignalFlowGraph *sg = sfg_create("test", 10, 20);
    int src = sfg_add_node(sg, SFG_SOURCE, "R");
    int e   = sfg_add_node(sg, SFG_INTERNAL, "E");
    int y   = sfg_add_node(sg, SFG_SINK, "Y");
    sfg_add_branch(sg, src, e, 1.0);
    sfg_add_branch(sg, e, y, 10.0);
    sfg_add_branch(sg, y, e, -1.0);
    sfg_set_source(sg, src); sfg_set_sink(sg, y);
    double T = sfg_mason_transfer_function_gain(sg);
    printf("(T=%.6f expected %.6f) ", T, 10.0/11.0);
    assert(T > 0.0 && T < 2.0);
    sfg_destroy(sg);
    OK();
}
void test_ba_standard_feedback(void) {
    TEST("ba_build_standard_feedback")
    double nG[] = {10.0}, dG[] = {1.0};
    double nH[] = {1.0},  dH[] = {1.0};
    TransferFunction *G = tf_create(nG, 0, dG, 0);
    TransferFunction *H = tf_create(nH, 0, dH, 0);
    BlockDiagram *bd = ba_build_standard_feedback(G, H);
    assert(bd != NULL);
    assert(bd_validate(bd) == 0);
    tf_destroy(G); tf_destroy(H); bd_destroy(bd);
    OK();
}
void test_ba_disturbance_rejection(void) {
    TEST("ba_build_disturbance_rejection")
    double nC[] = {5.0}, dC[] = {1.0};
    double nP[] = {2.0}, dP[] = {1.0};
    TransferFunction *C = tf_create(nC, 0, dC, 0);
    TransferFunction *P = tf_create(nP, 0, dP, 0);
    BlockDiagram *bd = ba_build_disturbance_rejection(C, P, NULL);
    assert(bd != NULL);
    assert(bd_validate(bd) == 0);
    tf_destroy(C); tf_destroy(P); bd_destroy(bd);
    OK();
}
void test_br_series_rule(void) {
    TEST("br_apply_series_rule")
    BlockDiagram *bd = bd_create("test", 10, 10);
    int in  = bd_add_node(bd, BD_NODE_INPUT, "R");
    int g1  = bd_add_node(bd, BD_NODE_BLOCK, "G1");
    int g2  = bd_add_node(bd, BD_NODE_BLOCK, "G2");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y");
    bd_add_edge(bd, in, g1, 1.0);
    bd_add_edge(bd, g1, g2, 1.0);
    bd_add_edge(bd, g2, out, 1.0);
    assert(br_apply_series_rule(bd, g1, g2) == 1);
    bd_destroy(bd);
    OK();
}
void test_br_feedback_rule(void) {
    TEST("br_apply_feedback_rule")
    BlockDiagram *bd = bd_create("test", 10, 10);
    int in  = bd_add_node(bd, BD_NODE_INPUT, "R");
    int sum = bd_add_node(bd, BD_NODE_SUMMER, "Err");
    int fwd = bd_add_node(bd, BD_NODE_BLOCK, "G");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y");
    int tko = bd_add_node(bd, BD_NODE_TAKEOFF, "T");
    int fb  = bd_add_node(bd, BD_NODE_BLOCK, "H");
    double signs[2] = {1.0, -1.0};
    bd_set_summer_signs(bd, sum, signs, 2);
    bd_add_edge(bd, in, sum, 1.0);
    bd_add_edge(bd, sum, fwd, 1.0);
    bd_add_edge(bd, fwd, out, 1.0);
    bd_add_edge(bd, fwd, tko, 1.0);
    bd_add_edge(bd, tko, fb, 1.0);
    bd_add_edge(bd, fb, sum, 1.0);
    int r = br_apply_feedback_rule(bd, fwd, fb, sum);
    assert(r == 1 || r == 0);
    bd_destroy(bd);
    OK();
}
void test_mason_analyze(void) {
    TEST("mason_analyze")
    SignalFlowGraph *sg = sfg_create("test", 10, 20);
    int src = sfg_add_node(sg, SFG_SOURCE, "R");
    int e   = sfg_add_node(sg, SFG_INTERNAL, "E");
    int y   = sfg_add_node(sg, SFG_SINK, "Y");
    sfg_add_branch(sg, src, e, 1.0);
    sfg_add_branch(sg, e, y, 10.0);
    sfg_add_branch(sg, y, e, -1.0);
    sfg_set_source(sg, src); sfg_set_sink(sg, y);
    MasonResult *mr = mason_analyze(sg, 16, 16);
    assert(mr != NULL);
    assert(mr->n_forward_paths == 1);
    assert(mr->n_loops == 1);
    printf("(Delta=%.4f T=%.4f) ", mr->determinant, mr->transfer_gain);
    assert(mr->determinant > 0.0);
    assert(mr->transfer_gain > 0.0);
    mason_result_free(mr);
    sfg_destroy(sg);
    OK();
}
void test_br_reduce_simple(void) {
    TEST("br_reduce simple")
    BlockDiagram *bd = bd_create("test", 10, 10);
    int in  = bd_add_node(bd, BD_NODE_INPUT, "R");
    int g   = bd_add_node(bd, BD_NODE_BLOCK, "G");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y");
    bd_add_edge(bd, in, g, 1.0);
    bd_add_edge(bd, g, out, 1.0);
    double n[] = {5.0}; double d[] = {1.0};
    bd_set_transfer_function(bd, g, tf_create(n, 0, d, 0));
    ReduceOptions opts = br_options_default();
    int steps = br_reduce(bd, &opts, NULL);
    assert(steps >= 0);
    bd_destroy(bd);
    OK();
}
void test_bd_clone(void) {
    TEST("bd_clone")
    BlockDiagram *bd = bd_create("orig", 10, 10);
    int in = bd_add_node(bd, BD_NODE_INPUT, "R");
    int blk = bd_add_node(bd, BD_NODE_BLOCK, "G");
    int out = bd_add_node(bd, BD_NODE_OUTPUT, "Y");
    bd_add_edge(bd, in, blk, 1.0);
    bd_add_edge(bd, blk, out, 2.0);
    BlockDiagram *cp = bd_clone(bd);
    assert(cp != NULL);
    assert(cp->n_nodes == bd->n_nodes);
    assert(cp->n_edges == bd->n_edges);
    bd_destroy(bd); bd_destroy(cp);
    OK();
}
int main(void) {
    printf("=== mini-block-diagram Test Suite ===\n");
    test_bd_create_destroy();
    test_bd_add_nodes();
    test_bd_add_edges();
    test_bd_find_forward_paths();
    test_bd_validate();
    test_tf_create();
    test_tf_dc_gain();
    test_tf_series();
    test_tf_feedback();
    test_tf_routh_stable();
    test_tf_routh_unstable();
    test_sfg_create();
    test_sfg_mason_simple();
    test_sfg_mason_feedback();
    test_ba_standard_feedback();
    test_ba_disturbance_rejection();
    test_br_series_rule();
    test_br_feedback_rule();
    test_mason_analyze();
    test_br_reduce_simple();
    test_bd_clone();
    printf("\n=== Results: %d/%d passed ===\n", tests_passed, tests_run);
    return tests_passed == tests_run ? 0 : 1;
}
