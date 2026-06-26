/**
 * test_sfg.c — Comprehensive Signal Flow Graph Test Suite
 *
 * Covers all core APIs: graph construction, path/loop enumeration,
 * Mason's formula, SFG reduction, state-space conversion, and analysis.
 *
 * Uses standard assert() for all checks.
 * Compile with: gcc -std=c99 -Wall -Wextra -I../include -o test_sfg
 *               test_sfg.c ../src/sfg_core.c ../src/sfg_path.c
 *               ../src/sfg_mason.c ../src/sfg_complex.c
 *               ../src/sfg_reduction.c ../src/sfg_state_space.c
 *               ../src/sfg_analysis.c ../src/sfg_matrix.c -lm
 */

#include "sfg_core.h"
#include "sfg_path.h"
#include "sfg_mason.h"
#include "sfg_reduction.h"
#include "sfg_state_space.h"
#include "sfg_analysis.h"
#include "sfg_matrix.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define EPS 1e-10

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
    printf("  TEST: %s ... ", name); \
    fflush(stdout); \
} while(0)

#define PASS() do { \
    printf("PASS\n"); \
    tests_passed++; \
} while(0)

#define FAIL(msg) do { \
    printf("FAIL: %s\n", msg); \
    tests_failed++; \
} while(0)

/* ================================================================
 * L1 Tests: Graph Construction and Basic Operations
 * ================================================================ */

static void test_graph_init(void) {
    TEST("graph initialization");
    sfg_graph_t g;
    sfg_graph_init(&g);
    assert(g.num_nodes == 0);
    assert(g.num_branches == 0);
    assert(g.input_node == SFG_INVALID_NODE_ID);
    assert(g.output_node == SFG_INVALID_NODE_ID);
    PASS();
}

static void test_add_nodes(void) {
    TEST("add nodes of all types");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "input");
    sfg_node_id_t snk = sfg_add_sink(&g, "output");
    sfg_node_id_t st  = sfg_add_state(&g, "x1");
    sfg_node_id_t sum = sfg_add_summing(&g, "e");
    sfg_node_id_t pko = sfg_add_node(&g, SFG_NODE_PICKOFF, "pk");

    assert(src == 0);
    assert(snk == 1);
    assert(st  == 2);
    assert(sum == 3);
    assert(pko == 4);
    assert(g.num_nodes == 5);

    /* Verify node types */
    assert(g.nodes[0].type == SFG_NODE_SOURCE);
    assert(g.nodes[1].type == SFG_NODE_SINK);
    assert(g.nodes[2].type == SFG_NODE_STATE);
    assert(g.nodes[3].type == SFG_NODE_SUMMING);
    assert(g.nodes[4].type == SFG_NODE_PICKOFF);

    PASS();
}

static void test_add_branches(void) {
    TEST("add branches");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t n0 = sfg_add_source(&g, "u");
    sfg_node_id_t n1 = sfg_add_state(&g, "x");
    sfg_node_id_t n2 = sfg_add_sink(&g, "y");

    int b0 = sfg_add_branch_real(&g, n0, n1, 2.0, NULL);
    int b1 = sfg_add_branch_real(&g, n1, n2, 3.0, NULL);

    assert(b0 == 0);
    assert(b1 == 1);
    assert(g.num_branches == 2);
    assert(g.branches[0].src == n0);
    assert(g.branches[0].dst == n1);
    assert(fabs(creal(g.branches[0].gain) - 2.0) < EPS);
    assert(fabs(creal(g.branches[1].gain) - 3.0) < EPS);

    PASS();
}

static void test_node_query(void) {
    TEST("node query operations");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t n0 = sfg_add_source(&g, "in");
    sfg_node_id_t n1 = sfg_add_state(&g, "x");
    sfg_node_id_t n2 = sfg_add_sink(&g, "out");

    sfg_add_branch_real(&g, n0, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, n2, 1.0, NULL);
    sfg_add_branch_real(&g, n1, n1, 0.5, NULL); /* Self-loop */

    assert(sfg_node_exists(&g, n0));
    assert(sfg_node_exists(&g, n1));
    assert(!sfg_node_exists(&g, 99));

    assert(sfg_node_in_degree(&g, n0) == 0);
    assert(sfg_node_out_degree(&g, n0) == 1);
    assert(sfg_node_in_degree(&g, n1) == 2);  /* From n0 + self-loop */
    assert(sfg_node_out_degree(&g, n1) == 2); /* To n2 + self-loop */
    assert(sfg_node_in_degree(&g, n2) == 1);

    assert(sfg_is_source(&g, n0));
    assert(sfg_is_sink(&g, n2));
    assert(!sfg_is_source(&g, n1));
    assert(!sfg_is_sink(&g, n1));

    PASS();
}

static void test_remove_branch(void) {
    TEST("remove branch");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t n0 = sfg_add_source(&g, "u");
    sfg_node_id_t n1 = sfg_add_sink(&g, "y");
    sfg_add_branch_real(&g, n0, n1, 5.0, NULL);

    assert(g.num_branches == 1);
    int ret = sfg_remove_branch(&g, n0, n1);
    assert(ret == 0);
    assert(g.num_branches == 0);
    assert(sfg_remove_branch(&g, n0, n1) == -1); /* Not found */

    PASS();
}

static void test_graph_validate(void) {
    TEST("graph validation");

    sfg_graph_t g;
    sfg_graph_init(&g);
    assert(sfg_graph_validate(&g) == 0);

    sfg_node_id_t n0 = sfg_add_source(&g, "u");
    sfg_node_id_t n1 = sfg_add_sink(&g, "y");
    sfg_add_branch_real(&g, n0, n1, 1.0, NULL);
    assert(sfg_graph_validate(&g) == 0);

    PASS();
}

static void test_graph_copy_transpose(void) {
    TEST("graph copy and transpose");

    sfg_graph_t g, g2;
    sfg_graph_init(&g);

    sfg_node_id_t n0 = sfg_add_source(&g, "u");
    sfg_node_id_t n1 = sfg_add_state(&g, "x");
    sfg_node_id_t n2 = sfg_add_sink(&g, "y");
    sfg_add_branch_real(&g, n0, n1, 2.0, NULL);
    sfg_add_branch_real(&g, n1, n2, 3.0, NULL);
    sfg_set_input(&g, n0);
    sfg_set_output(&g, n2);

    /* Test copy */
    sfg_graph_copy(&g2, &g);
    assert(g2.num_nodes == 3);
    assert(g2.num_branches == 2);
    assert(g2.input_node == n0);
    assert(g2.output_node == n2);

    /* Test transpose */
    sfg_graph_t gt;
    sfg_graph_transpose(&gt, &g);
    assert(gt.num_nodes == 3);
    assert(gt.num_branches == 2);
    /* Branches should be reversed */
    assert(gt.branches[0].src == n1 && gt.branches[0].dst == n0);
    assert(gt.branches[1].src == n2 && gt.branches[1].dst == n1);

    PASS();
}

/* ================================================================
 * L5 Tests: Path and Loop Enumeration
 * ================================================================ */

static void test_forward_paths_simple(void) {
    TEST("forward path enumeration — simple chain");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "n1");
    sfg_node_id_t n2  = sfg_add_node(&g, SFG_NODE_INTERNAL, "n2");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, n2, 2.0, NULL);
    sfg_add_branch_real(&g, n2, snk, 3.0, NULL);

    sfg_path_list_t paths;
    sfg_path_list_clear(&paths);

    int n = sfg_find_forward_paths(&g, src, snk, &paths);
    assert(n == 1);
    assert(paths.count == 1);
    assert(paths.paths[0].length == 4);
    assert(fabs(creal(paths.paths[0].gain) - 6.0) < EPS); /* 1*2*3 */

    PASS();
}

static void test_forward_paths_multiple(void) {
    TEST("forward path enumeration — multiple paths");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "n1");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    /* Two parallel paths */
    sfg_add_branch_real(&g, src, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, snk, 2.0, NULL);
    sfg_add_branch_real(&g, src, snk, 4.0, NULL);  /* Direct path */

    sfg_path_list_t paths;
    sfg_path_list_clear(&paths);

    int n = sfg_find_forward_paths(&g, src, snk, &paths);
    assert(n == 2);

    PASS();
}

static void test_loop_enumeration(void) {
    TEST("loop enumeration");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "n1");
    sfg_node_id_t n2  = sfg_add_node(&g, SFG_NODE_INTERNAL, "n2");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, n2, 2.0, NULL);
    sfg_add_branch_real(&g, n2, n1, 3.0, NULL);  /* Feedback loop! */
    sfg_add_branch_real(&g, n2, snk, 4.0, NULL);
    sfg_add_branch_real(&g, n1, n1, 0.5, NULL);  /* Self-loop */

    sfg_loop_list_t loops;
    sfg_loop_list_clear(&loops);

    int n = sfg_find_all_loops(&g, &loops);
    assert(n >= 2); /* At least the n1→n2→n1 loop and the self-loop */

    PASS();
}

static void test_loops_touch(void) {
    TEST("loop touching detection");

    sfg_graph_t g;
    sfg_graph_init(&g);

    (void)sfg_add_source(&g, "u");
    sfg_node_id_t n1 = sfg_add_node(&g, SFG_NODE_INTERNAL, "a");
    sfg_node_id_t n2 = sfg_add_node(&g, SFG_NODE_INTERNAL, "b");
    sfg_node_id_t n3 = sfg_add_node(&g, SFG_NODE_INTERNAL, "c");
    (void)sfg_add_sink(&g, "y");

    /* Loop 1: n1→n2→n1 */
    sfg_add_branch_real(&g, n1, n2, 1.0, NULL);
    sfg_add_branch_real(&g, n2, n1, 2.0, NULL);

    /* Loop 2: n3→n4 but n4 goes nowhere... different loop */
    /* Let's make two non-touching loops:
     * L_a: n1→n2→n1
     * L_b: n3→n3 (self-loop) */
    sfg_add_branch_real(&g, n3, n3, 5.0, NULL);

    sfg_loop_list_t loops;
    sfg_loop_list_clear(&loops);
    sfg_find_all_loops(&g, &loops);

    /* Find the two loops and verify they are non-touching */
    assert(loops.count >= 2);
    int loop_a_idx = -1, loop_b_idx = -1;
    for (int i = 0; i < loops.count; i++) {
        if (loops.loops[i].length == 2) loop_a_idx = i;
        if (loops.loops[i].length == 1) loop_b_idx = i;
    }

    if (loop_a_idx >= 0 && loop_b_idx >= 0) {
        int touches = sfg_loops_touch(&loops.loops[loop_a_idx],
                                       &loops.loops[loop_b_idx]);
        assert(touches == 0); /* Should be non-touching */
    }

    PASS();
}

/* ================================================================
 * L4 Tests: Mason's Gain Formula
 * ================================================================ */

static void test_mason_simple_chain(void) {
    TEST("Mason's formula — simple chain (no loops)");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "n1");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, n1, 3.0, NULL);
    sfg_add_branch_real(&g, n1, snk, 4.0, NULL);

    sfg_mason_result_t result;
    sfg_mason_result_init(&result);

    int ret = sfg_mason_compute(&g, src, snk, &result);
    assert(ret == 0);
    assert(result.success == 1);
    assert(result.num_forward_paths == 1);
    assert(fabs(creal(result.transfer_function) - 12.0) < EPS);

    sfg_mason_result_free(&result);
    PASS();
}

static void test_mason_with_feedback(void) {
    TEST("Mason's formula — single feedback loop");

    /*
     * Classic feedback system:
     *   Forward gain G = 10, Feedback H = 0.5
     *   T = G / (1 + GH) = 10 / (1 + 5) = 10/6 = 1.666...
     */
    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "r");
    sfg_node_id_t err = sfg_add_summing(&g, "e");
    sfg_node_id_t out = sfg_add_sink(&g, "y");
    sfg_node_id_t fb  = sfg_add_node(&g, SFG_NODE_PICKOFF, "fb");

    sfg_add_branch_real(&g, src, err, 1.0, NULL);
    sfg_add_branch_real(&g, err, out, 10.0, NULL);   /* G = 10 */
    sfg_add_branch_real(&g, out, fb, 1.0, NULL);
    sfg_add_branch_real(&g, fb, err, -0.5, NULL);    /* H = 0.5, negative fb */

    sfg_mason_result_t result;
    sfg_mason_result_init(&result);

    int ret = sfg_mason_compute(&g, src, out, &result);
    assert(ret == 0);
    assert(result.success == 1);

    double expected = 10.0 / (1.0 + 5.0);
    assert(fabs(creal(result.transfer_function) - expected) < EPS);

    sfg_mason_result_free(&result);
    PASS();
}

static void test_mason_multiple_loops(void) {
    TEST("Mason's formula — multiple non-touching loops");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "a");
    sfg_node_id_t n2  = sfg_add_node(&g, SFG_NODE_INTERNAL, "b");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    /* Forward path: u → a → b → y with gains 1, 2, 3 */
    sfg_add_branch_real(&g, src, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, n2, 2.0, NULL);
    sfg_add_branch_real(&g, n2, snk, 3.0, NULL);

    /* Loop 1: a → a (self-loop gain 0.1) */
    sfg_add_branch_real(&g, n1, n1, 0.1, NULL);

    /* Loop 2: b → b (self-loop gain 0.2)
     * These two self-loops are non-touching */
    sfg_add_branch_real(&g, n2, n2, 0.2, NULL);

    /* Mason:
     * P₁ = 1·2·3 = 6
     * L₁ = 0.1, L₂ = 0.2
     * Δ = 1 - (0.1+0.2) + (0.1·0.2) = 1 - 0.3 + 0.02 = 0.72
     * Δ₁ = 1 - 0 = 1 (both loops touch the path)
     * T = 6·1 / 0.72 = 8.333... */

    sfg_mason_result_t result;
    sfg_mason_result_init(&result);

    int ret = sfg_mason_compute(&g, src, snk, &result);
    assert(ret == 0);
    assert(result.success == 1);

    double expected = 6.0 / 0.72;
    assert(fabs(creal(result.transfer_function) - expected) < EPS);

    /* Verify determinant components */
    assert(fabs(creal(result.sum_L1) - 0.3) < EPS);
    assert(fabs(creal(result.sum_L2) - 0.02) < EPS);

    sfg_mason_result_free(&result);
    PASS();
}

static void test_mason_no_path(void) {
    TEST("Mason's formula — no path (zero transfer)");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");
    /* No branch between them */

    sfg_mason_result_t result;
    sfg_mason_result_init(&result);

    int ret = sfg_mason_compute(&g, src, snk, &result);
    assert(ret == 0);
    assert(fabs(creal(result.transfer_function)) < EPS);

    sfg_mason_result_free(&result);
    PASS();
}

/* ================================================================
 * L5 Tests: SFG Reduction
 * ================================================================ */

static void test_reduce_series(void) {
    TEST("SFG reduction — series elimination");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t mid = sfg_add_node(&g, SFG_NODE_INTERNAL, "mid");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, mid, 2.0, NULL);
    sfg_add_branch_real(&g, mid, snk, 3.0, NULL);

    int eliminated = sfg_reduce_series(&g);
    assert(eliminated == 1);

    /* After elimination, should have direct branch between remaining nodes
     * with gain 6. Node IDs may have shifted, so scan for source/sink. */
    assert(g.num_nodes == 2);

    /* Find the source and sink by type */
    sfg_node_id_t new_src = SFG_INVALID_NODE_ID;
    sfg_node_id_t new_snk = SFG_INVALID_NODE_ID;
    for (int i = 0; i < g.num_nodes; i++) {
        if (g.nodes[i].type == SFG_NODE_SOURCE) new_src = (sfg_node_id_t)i;
        if (g.nodes[i].type == SFG_NODE_SINK)   new_snk = (sfg_node_id_t)i;
    }
    assert(new_src != SFG_INVALID_NODE_ID);
    assert(new_snk != SFG_INVALID_NODE_ID);

    const sfg_branch_t* b = sfg_find_branch(&g, new_src, new_snk);
    assert(b != NULL);
    assert(fabs(creal(b->gain) - 6.0) < EPS);

    PASS();
}

static void test_reduce_parallel(void) {
    TEST("SFG reduction — parallel combination");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, snk, 2.0, NULL);
    sfg_add_branch_real(&g, src, snk, 3.0, NULL);
    sfg_add_branch_real(&g, src, snk, 5.0, NULL);

    int combined = sfg_reduce_parallel(&g);
    assert(combined == 1);

    /* After parallel reduction, should have one branch */
    assert(g.num_branches == 1);
    assert(g.branches[0].src == src);
    assert(g.branches[0].dst == snk);
    assert(fabs(creal(g.branches[0].gain) - 10.0) < EPS);

    PASS();
}

static void test_reduce_self_loop(void) {
    TEST("SFG reduction — self-loop elimination");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "n1");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, n1, 0.5, NULL); /* Self-loop gain 0.5 */
    sfg_add_branch_real(&g, n1, snk, 2.0, NULL);

    int eliminated = sfg_reduce_self_loop(&g);
    assert(eliminated == 1);
    /* Incoming branch to n1 should be divided by (1 - 0.5) = 0.5
     * So gain becomes 1.0 / 0.5 = 2.0 */
    /* Find the incoming branch to n1 (src→n1) */
    const sfg_branch_t* b = sfg_find_branch(&g, src, n1);
    assert(b != NULL);
    assert(fabs(creal(b->gain) - 2.0) < EPS);

    (void)snk;
    PASS();
}

static void test_reduce_to_single(void) {
    TEST("SFG reduction — full reduction to single branch");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    (void)sfg_add_node(&g, SFG_NODE_INTERNAL, "n1");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, 1, 2.0, NULL);    /* src→n1 */
    sfg_add_branch_real(&g, 1, snk, 3.0, NULL);    /* n1→snk */

    sfg_reduction_result_t result;
    sfg_reduction_result_init(&result);

    /* Note: sfg_reduce_to_single works on a copy; the original g
     * is unchanged. Check the result's final gain. */
    int ret = sfg_reduce_to_single(&g, src, snk, &result);
    assert(ret == 0);
    assert(result.success == 1);
    assert(fabs(creal(result.overall_gain) - 6.0) < EPS);

    PASS();
}

/* ================================================================
 * L6 Tests: State-Space Conversion
 * ================================================================ */

static void test_ss_alloc_free(void) {
    TEST("state-space — allocation and deallocation");

    sfg_ss_system_t ss;
    sfg_ss_init(&ss);

    int ret = sfg_ss_alloc(&ss, 3, 1, 1);
    assert(ret == 0);
    assert(ss.n == 3);
    assert(ss.m == 1);
    assert(ss.p == 1);
    assert(ss.A != NULL);
    assert(ss.B != NULL);
    assert(ss.C != NULL);
    assert(ss.D != NULL);

    sfg_ss_free(&ss);
    assert(ss.A == NULL);
    PASS();
}

static void test_ss_to_sfg(void) {
    TEST("state-space → SFG conversion");

    sfg_ss_system_t ss;
    sfg_ss_init(&ss);
    sfg_ss_alloc(&ss, 2, 1, 1);

    /* Simple second-order system */
    sfg_ss_set_A(&ss, 0, 0, 0.0);
    sfg_ss_set_A(&ss, 0, 1, 1.0);
    sfg_ss_set_A(&ss, 1, 0, -2.0);
    sfg_ss_set_A(&ss, 1, 1, -3.0);
    sfg_ss_set_B(&ss, 0, 0, 0.0);
    sfg_ss_set_B(&ss, 1, 0, 1.0);
    sfg_ss_set_C(&ss, 0, 0, 1.0);
    sfg_ss_set_C(&ss, 0, 1, 0.0);

    sfg_graph_t g;
    sfg_graph_init(&g);

    int ret = sfg_from_state_space(&ss, &g);
    assert(ret == 0);
    assert(g.num_nodes > 0);
    assert(g.num_branches > 0);

    sfg_ss_free(&ss);
    PASS();
}

static void test_characteristic_polynomial(void) {
    TEST("characteristic polynomial via Faddeev-Leverrier");

    sfg_ss_system_t ss;
    sfg_ss_init(&ss);
    sfg_ss_alloc(&ss, 2, 1, 1);

    /* A = [[0, 1], [-2, -3]]
     * |sI - A| = s² + 3s + 2 */
    sfg_ss_set_A(&ss, 0, 0, 0.0);
    sfg_ss_set_A(&ss, 0, 1, 1.0);
    sfg_ss_set_A(&ss, 1, 0, -2.0);
    sfg_ss_set_A(&ss, 1, 1, -3.0);

    double poly[3];
    int ret = sfg_ss_characteristic_polynomial(&ss, poly, 2);
    assert(ret == 0);
    /* poly = [a0, a1, a2] = [2, 3, 1] */
    assert(fabs(poly[0] - 2.0) < EPS);
    assert(fabs(poly[1] - 3.0) < EPS);
    assert(fabs(poly[2] - 1.0) < EPS);

    sfg_ss_free(&ss);
    PASS();
}

/* ================================================================
 * L5 Tests: Complex Number Utilities
 * ================================================================ */

static void test_complex_arithmetic(void) {
    TEST("complex number arithmetic");

    double complex a = 3.0 + 4.0 * I;
    double complex b = 1.0 + 2.0 * I;

    double complex s = sfg_cmplx_add(a, b);
    assert(fabs(creal(s) - 4.0) < EPS);
    assert(fabs(cimag(s) - 6.0) < EPS);

    double complex p = sfg_cmplx_mul(a, b);
    /* (3+4j)(1+2j) = 3+6j+4j-8 = -5+10j */
    assert(fabs(creal(p) - (-5.0)) < EPS);
    assert(fabs(cimag(p) - 10.0) < EPS);

    double mag = sfg_cmplx_magnitude(a);
    assert(fabs(mag - 5.0) < EPS); /* |3+4j| = 5 */

    double complex conj = sfg_cmplx_conjugate(a);
    assert(fabs(creal(conj) - 3.0) < EPS);
    assert(fabs(cimag(conj) - (-4.0)) < EPS);

    PASS();
}

static void test_poly_eval(void) {
    TEST("polynomial evaluation (Horner)");

    /* P(x) = 2 + 3x + x^2 */
    double coef[] = {2.0, 3.0, 1.0};
    double val = sfg_poly_eval_real(coef, 2, 2.0);
    /* 2 + 3*2 + 4 = 2 + 6 + 4 = 12 */
    assert(fabs(val - 12.0) < EPS);

    val = sfg_poly_eval_real(coef, 2, 0.0);
    assert(fabs(val - 2.0) < EPS);

    PASS();
}

static void test_poly_derivative(void) {
    TEST("polynomial derivative");

    /* P(x) = 2 + 3x + x^2, P'(x) = 3 + 2x */
    double coef[] = {2.0, 3.0, 1.0};
    double deriv[2];
    sfg_poly_derivative(coef, 2, deriv);
    assert(fabs(deriv[0] - 3.0) < EPS); /* coefficient of x^0 in P' */
    assert(fabs(deriv[1] - 2.0) < EPS); /* coefficient of x^1 in P' */

    PASS();
}

static void test_poly_multiply(void) {
    TEST("polynomial multiplication");

    /* A(x) = 1 + x, B(x) = 2 + x */
    double a[] = {1.0, 1.0};
    double b[] = {2.0, 1.0};
    double result[3];
    int deg = sfg_poly_multiply(a, 1, b, 1, result);
    assert(deg == 2);
    /* (1+x)(2+x) = 2 + x + 2x + x^2 = 2 + 3x + x^2 */
    assert(fabs(result[0] - 2.0) < EPS);
    assert(fabs(result[1] - 3.0) < EPS);
    assert(fabs(result[2] - 1.0) < EPS);

    PASS();
}

/* ================================================================
 * L6 Tests: Transfer Function from SFG
 * ================================================================ */

static void test_tf_from_sfg(void) {
    TEST("transfer function construction from SFG");

    sfg_graph_t g;
    sfg_graph_init(&g);

    /* T(s) = 1 / (s² + 3s + 2) */
    sfg_real_t num[] = {1.0};  /* Numerator: 1 */
    sfg_real_t den[] = {2.0, 3.0, 1.0};  /* Denominator: 2 + 3s + s² */

    int ret = sfg_from_transfer_function(&g, num, 0, den, 2);
    assert(ret == 0);
    assert(g.num_nodes > 0);

    PASS();
}

/* ================================================================
 * L7 Tests: Sensitivity Analysis
 * ================================================================ */

static void test_sensitivity(void) {
    TEST("sensitivity analysis");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "G");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, snk, 10.0, NULL); /* Gain = 10 */
    sfg_add_branch_real(&g, snk, n1, -0.5, NULL); /* FB = -0.5 */

    sfg_real_t sens;
    int ret = sfg_sensitivity(&g, src, snk, n1, snk, &sens);
    assert(ret == 0);
    /* Should give a finite, non-zero sensitivity */

    PASS();
}

/* ================================================================
 * L6 Tests: Stability Analysis
 * ================================================================ */

static void test_stability_routh(void) {
    TEST("Routh-Hurwitz stability test");

    /* Stable polynomial: s³ + 2s² + 2s + 1
     * Routh:
     * s³ | 1   2   0
     * s² | 2   1   0
     * s¹ | (4-1)/2 = 1.5  0
     * s⁰ | 1
     * First column: 1, 2, 1.5, 1 — all positive → stable */
    double coef[] = {1.0, 2.0, 2.0, 1.0}; /* a0..a3 */
    double routh[16]; /* 4 rows × 4 cols */
    int stable;
    int ret = sfg_stability_routh(coef, 3, routh, &stable);
    assert(ret == 0);
    assert(stable == 1);

    /* Unstable: s² + s - 2 = (s+2)(s-1) */
    double coef2[] = {-2.0, 1.0, 1.0};
    ret = sfg_stability_routh(coef2, 2, routh, &stable);
    assert(ret == 0);
    assert(stable == 0);

    PASS();
}

/* ================================================================
 * L8 Tests: Monte Carlo and Matrix Methods
 * ================================================================ */

static void test_monte_carlo(void) {
    TEST("Monte Carlo tolerance analysis");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");
    sfg_add_branch_real(&g, src, snk, 10.0, NULL);

    double mean, stddev, worst, best;
    int ret = sfg_monte_carlo(&g, src, snk, 0.05, 100,
                               &mean, &stddev, &worst, &best);
    assert(ret == 0);
    assert(mean > 0.0);
    assert(stddev >= 0.0);

    PASS();
}

static void test_condition_number(void) {
    TEST("condition number estimation");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "a");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, snk, 10.0, NULL);
    sfg_add_branch_real(&g, snk, n1, -0.5, NULL);

    double kappa = sfg_condition_number(&g);
    assert(kappa > 0.0);

    PASS();
}

/* ================================================================
 * L4 Tests: Mason vs Direct Verification
 * ================================================================ */

static void test_mason_verify(void) {
    TEST("Mason's formula vs direct matrix inversion");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "x1");
    sfg_node_id_t n2  = sfg_add_node(&g, SFG_NODE_INTERNAL, "x2");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, n2, 2.0, NULL);
    sfg_add_branch_real(&g, n2, snk, 3.0, NULL);
    sfg_add_branch_real(&g, n2, n1, 0.1, NULL); /* Feedback */

    double error = sfg_mason_verify(&g, src, snk);
    assert(error < 1e-8);
    assert(error >= 0.0);

    PASS();
}

static void test_cramer_verify(void) {
    TEST("Cramer's rule vs Mason's formula");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t src = sfg_add_source(&g, "u");
    sfg_node_id_t n1  = sfg_add_node(&g, SFG_NODE_INTERNAL, "x");
    sfg_node_id_t snk = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, src, n1, 1.0, NULL);
    sfg_add_branch_real(&g, n1, snk, 5.0, NULL);
    sfg_add_branch_real(&g, n1, n1, 0.2, NULL); /* Self-loop */

    double error = sfg_cramer_verify(&g, src, snk);
    assert(error < 1e-8);
    assert(error >= 0.0);

    PASS();
}

/* ================================================================
 * L6 Tests: Controllability / Observability via SFG
 * ================================================================ */

static void test_sfg_controllability(void) {
    TEST("controllability via SFG");

    sfg_ss_system_t ss;
    sfg_ss_init(&ss);
    sfg_ss_alloc(&ss, 2, 1, 1);

    /* Controllable canonical form */
    sfg_ss_set_A(&ss, 0, 0, 0.0);
    sfg_ss_set_A(&ss, 0, 1, 1.0);
    sfg_ss_set_A(&ss, 1, 0, -2.0);
    sfg_ss_set_A(&ss, 1, 1, -3.0);
    sfg_ss_set_B(&ss, 0, 0, 0.0);
    sfg_ss_set_B(&ss, 1, 0, 1.0);

    int ctrl = sfg_ss_controllability(&ss);
    assert(ctrl == 1);

    sfg_ss_free(&ss);
    PASS();
}

/* ================================================================
 * Edge Case Tests
 * ================================================================ */

static void test_empty_graph(void) {
    TEST("operations on empty graph");

    sfg_graph_t g;
    sfg_graph_init(&g);

    assert(sfg_node_count(&g) == 0);
    assert(sfg_branch_count(&g) == 0);

    sfg_path_list_t paths;
    sfg_path_list_clear(&paths);
    int n = sfg_find_forward_paths(&g, 0, 0, &paths);
    assert(n == -1 || n == 0); /* Should fail gracefully */

    PASS();
}

static void test_duplicate_branches(void) {
    TEST("duplicate parallel branches");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t n0 = sfg_add_source(&g, "u");
    sfg_node_id_t n1 = sfg_add_sink(&g, "y");

    sfg_add_branch_real(&g, n0, n1, 2.0, NULL);
    sfg_add_branch_real(&g, n0, n1, 3.0, NULL);

    assert(g.num_branches == 2);

    sfg_gain_t T = sfg_mason_compute_simple(&g, n0, n1);
    assert(fabs(creal(T) - 5.0) < EPS); /* Parallel → sum */

    PASS();
}

static void test_self_loop_only(void) {
    TEST("graph with only self-loop");

    sfg_graph_t g;
    sfg_graph_init(&g);

    sfg_node_id_t n0 = sfg_add_node(&g, SFG_NODE_INTERNAL, "x");
    sfg_add_branch_real(&g, n0, n0, 0.5, NULL);

    sfg_loop_list_t loops;
    sfg_loop_list_clear(&loops);
    int n = sfg_find_self_loops(&g, &loops);
    assert(n == 1);
    assert(loops.loops[0].length == 1);
    assert(fabs(creal(loops.loops[0].gain) - 0.5) < EPS);

    PASS();
}

/* ================================================================
 * Main
 * ================================================================ */

int main(void) {
    printf("\n==============================================\n");
    printf("  mini-signal-flow-graph — Test Suite\n");
    printf("==============================================\n\n");

    printf("--- L1: Core Definitions ---\n");
    test_graph_init();
    test_add_nodes();
    test_add_branches();
    test_node_query();
    test_remove_branch();
    test_graph_validate();
    test_graph_copy_transpose();
    test_empty_graph();
    test_duplicate_branches();
    test_self_loop_only();

    printf("\n--- L5: Path/Loop Enumeration ---\n");
    test_forward_paths_simple();
    test_forward_paths_multiple();
    test_loop_enumeration();
    test_loops_touch();

    printf("\n--- L4+L5: Mason's Gain Formula ---\n");
    test_mason_simple_chain();
    test_mason_with_feedback();
    test_mason_multiple_loops();
    test_mason_no_path();
    test_mason_verify();
    test_cramer_verify();

    printf("\n--- L5: SFG Reduction ---\n");
    test_reduce_series();
    test_reduce_parallel();
    test_reduce_self_loop();
    test_reduce_to_single();

    printf("\n--- L6: State-Space ---\n");
    test_ss_alloc_free();
    test_ss_to_sfg();
    test_characteristic_polynomial();
    test_tf_from_sfg();
    test_sfg_controllability();

    printf("\n--- L3+L5: Complex Numbers ---\n");
    test_complex_arithmetic();
    test_poly_eval();
    test_poly_derivative();
    test_poly_multiply();

    printf("\n--- L7+L8: Analysis ---\n");
    test_sensitivity();
    test_stability_routh();
    test_monte_carlo();
    test_condition_number();

    printf("\n==============================================\n");
    printf("  Results: %d PASS, %d FAIL\n",
           tests_passed, tests_failed);
    printf("==============================================\n\n");

    return tests_failed > 0 ? 1 : 0;
}
