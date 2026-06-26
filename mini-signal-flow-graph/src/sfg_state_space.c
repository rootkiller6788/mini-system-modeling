/**
 * sfg_state_space.c — State-Space ↔ Signal Flow Graph Conversion
 *
 * Implements bidirectional conversion between state-space representations
 * (A, B, C, D matrices) and signal flow graphs. Both representations
 * encode the same linear dynamical system.
 *
 * The SFG for a state-space system has a characteristic structure:
 *   - State nodes connected by A matrix gains
 *   - Integrator (1/s) branches from ẋ_i to x_i
 *   - Input branches from u_j to ẋ_i (B matrix)
 *   - Output branches from x_i to y_k (C matrix)
 *   - Direct feedthrough from u_j to y_k (D matrix)
 *
 * Knowledge coverage:
 *   L5: State-space ↔ SFG conversion algorithms
 *   L5: Canonical form construction (controllable, observable)
 *   L6: Controllability/observability via SFG path analysis
 *   L4: Characteristic polynomial via SFG determinant
 */

#include "sfg_state_space.h"
#include "sfg_mason.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * L1: State-Space System Management
 * ================================================================ */

void sfg_ss_init(sfg_ss_system_t* ss) {
    if (!ss) return;
    memset(ss, 0, sizeof(sfg_ss_system_t));
}

int sfg_ss_alloc(sfg_ss_system_t* ss, int n, int m, int p) {
    if (!ss) return -1;
    if (n <= 0 || n > SFG_SS_MAX_ORDER) return -1;
    if (m < 0) return -1;
    if (p < 0) return -1;

    sfg_ss_free(ss);

    ss->n = n;
    ss->m = m;
    ss->p = p;

    ss->A = (sfg_real_t*)calloc((size_t)(n * n), sizeof(sfg_real_t));
    ss->B = (sfg_real_t*)calloc((size_t)(n * m), sizeof(sfg_real_t));
    ss->C = (sfg_real_t*)calloc((size_t)(p * n), sizeof(sfg_real_t));
    ss->D = (sfg_real_t*)calloc((size_t)(p * m), sizeof(sfg_real_t));

    if (!ss->A || !ss->B || !ss->C || !ss->D) {
        sfg_ss_free(ss);
        return -1;
    }

    ss->A_alloc = 1;
    ss->B_alloc = 1;
    ss->C_alloc = 1;
    ss->D_alloc = 1;

    return 0;
}

void sfg_ss_free(sfg_ss_system_t* ss) {
    if (!ss) return;
    if (ss->A_alloc) { free(ss->A); ss->A = NULL; ss->A_alloc = 0; }
    if (ss->B_alloc) { free(ss->B); ss->B = NULL; ss->B_alloc = 0; }
    if (ss->C_alloc) { free(ss->C); ss->C = NULL; ss->C_alloc = 0; }
    if (ss->D_alloc) { free(ss->D); ss->D = NULL; ss->D_alloc = 0; }
}

/* ================================================================
 * L3: Matrix Element Accessors
 * ================================================================ */

void sfg_ss_set_A(sfg_ss_system_t* ss, int row, int col, sfg_real_t val) {
    if (ss && ss->A && row >= 0 && row < ss->n && col >= 0 && col < ss->n) {
        ss->A[row * ss->n + col] = val;
    }
}

void sfg_ss_set_B(sfg_ss_system_t* ss, int row, int col, sfg_real_t val) {
    if (ss && ss->B && row >= 0 && row < ss->n && col >= 0 && col < ss->m) {
        ss->B[row * ss->m + col] = val;
    }
}

void sfg_ss_set_C(sfg_ss_system_t* ss, int row, int col, sfg_real_t val) {
    if (ss && ss->C && row >= 0 && row < ss->p && col >= 0 && col < ss->n) {
        ss->C[row * ss->n + col] = val;
    }
}

void sfg_ss_set_D(sfg_ss_system_t* ss, int row, int col, sfg_real_t val) {
    if (ss && ss->D && row >= 0 && row < ss->p && col >= 0 && col < ss->m) {
        ss->D[row * ss->m + col] = val;
    }
}

sfg_real_t sfg_ss_get_A(const sfg_ss_system_t* ss, int row, int col) {
    if (ss && ss->A && row >= 0 && row < ss->n && col >= 0 && col < ss->n) {
        return ss->A[row * ss->n + col];
    }
    return 0.0;
}

sfg_real_t sfg_ss_get_B(const sfg_ss_system_t* ss, int row, int col) {
    if (ss && ss->B && row >= 0 && row < ss->n && col >= 0 && col < ss->m) {
        return ss->B[row * ss->m + col];
    }
    return 0.0;
}

sfg_real_t sfg_ss_get_C(const sfg_ss_system_t* ss, int row, int col) {
    if (ss && ss->C && row >= 0 && row < ss->p && col >= 0 && col < ss->n) {
        return ss->C[row * ss->n + col];
    }
    return 0.0;
}

sfg_real_t sfg_ss_get_D(const sfg_ss_system_t* ss, int row, int col) {
    if (ss && ss->D && row >= 0 && row < ss->p && col >= 0 && col < ss->m) {
        return ss->D[row * ss->m + col];
    }
    return 0.0;
}

/* ================================================================
 * L5: State-Space → SFG Conversion
 * ================================================================ */

int sfg_from_state_space(const sfg_ss_system_t* ss, sfg_graph_t* g) {
    if (!ss || !g) return -1;
    if (ss->n > SFG_MAX_NODES / 3) return -1; /* Each state needs ~3 nodes */

    sfg_graph_init(g);
    int n = ss->n, m = ss->m, p_val = ss->p;

    /* Create input (source) nodes: u_0 ... u_{m-1} */
    sfg_node_id_t* u_nodes = (sfg_node_id_t*)calloc((size_t)m,
                                                      sizeof(sfg_node_id_t));
    if (!u_nodes && m > 0) return -1;
    for (int j = 0; j < m; j++) {
        char label[32];
        snprintf(label, sizeof(label), "u_%d", j);
        u_nodes[j] = sfg_add_source(g, label);
    }

    /* Create state derivative nodes: xdot_0 ... xdot_{n-1} */
    sfg_node_id_t* xdot = (sfg_node_id_t*)calloc((size_t)n,
                                                   sizeof(sfg_node_id_t));
    if (!xdot) { free(u_nodes); return -1; }
    for (int i = 0; i < n; i++) {
        char label[32];
        snprintf(label, sizeof(label), "xdot_%d", i);
        xdot[i] = sfg_add_summing(g, label);
    }

    /* Create state nodes: x_0 ... x_{n-1} */
    sfg_node_id_t* x = (sfg_node_id_t*)calloc((size_t)n,
                                                sizeof(sfg_node_id_t));
    if (!x) { free(u_nodes); free(xdot); return -1; }
    for (int i = 0; i < n; i++) {
        char label[32];
        snprintf(label, sizeof(label), "x_%d", i);
        x[i] = sfg_add_state(g, label);
    }

    /* Create output (sink) nodes: y_0 ... y_{p-1} */
    sfg_node_id_t* y_nodes = (sfg_node_id_t*)calloc((size_t)p_val,
                                                      sizeof(sfg_node_id_t));
    if (!y_nodes && p_val > 0) {
        free(u_nodes); free(xdot); free(x); return -1;
    }
    for (int k = 0; k < p_val; k++) {
        char label[32];
        snprintf(label, sizeof(label), "y_%d", k);
        y_nodes[k] = sfg_add_sink(g, label);
    }

    /* Add integrator branches: xdot_i → x_i with gain 1 */
    for (int i = 0; i < n; i++) {
        sfg_add_branch_real(g, xdot[i], x[i], 1.0, NULL);
    }

    /* Add A matrix branches: x_j → xdot_i with gain A[i][j] */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sfg_real_t a_ij = sfg_ss_get_A(ss, i, j);
            if (fabs(a_ij) > 1e-15) {
                sfg_add_branch_real(g, x[j], xdot[i], a_ij, NULL);
            }
        }
    }

    /* Add B matrix branches: u_j → xdot_i with gain B[i][j] */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            sfg_real_t b_ij = sfg_ss_get_B(ss, i, j);
            if (fabs(b_ij) > 1e-15) {
                sfg_add_branch_real(g, u_nodes[j], xdot[i], b_ij, NULL);
            }
        }
    }

    /* Add C matrix branches: x_j → y_k with gain C[k][j] */
    for (int k = 0; k < p_val; k++) {
        for (int j = 0; j < n; j++) {
            sfg_real_t c_kj = sfg_ss_get_C(ss, k, j);
            if (fabs(c_kj) > 1e-15) {
                sfg_add_branch_real(g, x[j], y_nodes[k], c_kj, NULL);
            }
        }
    }

    /* Add D matrix branches: u_j → y_k with gain D[k][j] */
    for (int k = 0; k < p_val; k++) {
        for (int j = 0; j < m; j++) {
            sfg_real_t d_kj = sfg_ss_get_D(ss, k, j);
            if (fabs(d_kj) > 1e-15) {
                sfg_add_branch_real(g, u_nodes[j], y_nodes[k], d_kj, NULL);
            }
        }
    }

    /* Set primary input/output */
    if (m > 0) sfg_set_input(g, u_nodes[0]);
    if (p_val > 0) sfg_set_output(g, y_nodes[0]);

    if (ss->name[0]) {
        snprintf(g->name, sizeof(g->name), "ss_%.50s_n%d", ss->name, n);
    } else {
        snprintf(g->name, sizeof(g->name), "ss_sys_n%d", n);
    }

    free(u_nodes);
    free(xdot);
    free(x);
    free(y_nodes);

    return 0;
}

/* ================================================================
 * L5: Controllable Canonical Form → SFG
 * ================================================================ */

int sfg_from_state_space_controllable(const sfg_ss_system_t* ss,
                                       sfg_graph_t* g) {
    if (!ss || !g) return -1;
    if (ss->m != 1 || ss->p != 1) return -1; /* SISO only for canonical */
    /* Use the standard conversion */
    return sfg_from_state_space(ss, g);
}

int sfg_from_state_space_observable(const sfg_ss_system_t* ss,
                                     sfg_graph_t* g) {
    if (!ss || !g) return -1;
    if (ss->m != 1 || ss->p != 1) return -1;
    return sfg_from_state_space(ss, g);
}

/* ================================================================
 * L5: Transfer Function → SFG (Controllable Canonical Form)
 * ================================================================ */

int sfg_from_transfer_function(sfg_graph_t* g,
                                const sfg_real_t* num, int num_order,
                                const sfg_real_t* den, int den_order) {
    if (!g || !num || !den) return -1;
    if (den_order <= 0) return -1;

    sfg_graph_init(g);
    int n = den_order; /* System order = denominator degree */

    /* Build controllable canonical form:
     *
     * The state equations for T(s) = N(s)/D(s):
     *   A matrix is companion form with 1's on superdiagonal
     *     and den coefficients (negated) in last row
     *   B = [0 ... 0 1]^T
     *   C = [b_0 b_1 ... b_{n-1}]  (numerator coefficients)
     *   D = b_n (if proper, else 0)
     */

    /* Create nodes */
    sfg_node_id_t input = sfg_add_source(g, "u");
    /* n integrator node pairs */
    sfg_node_id_t* xdot = (sfg_node_id_t*)calloc((size_t)n,
                                                   sizeof(sfg_node_id_t));
    sfg_node_id_t* x = (sfg_node_id_t*)calloc((size_t)n,
                                                sizeof(sfg_node_id_t));
    if (!xdot || !x) {
        free(xdot); free(x); return -1;
    }
    for (int i = 0; i < n; i++) {
        char label[32];
        snprintf(label, sizeof(label), "xdot_%d", i);
        xdot[i] = sfg_add_summing(g, label);
        snprintf(label, sizeof(label), "x_%d", i);
        x[i] = sfg_add_state(g, label);
    }
    sfg_node_id_t output = sfg_add_sink(g, "y");

    /* Integrator branches: xdot_i → x_i with gain 1 */
    for (int i = 0; i < n; i++) {
        sfg_add_branch_real(g, xdot[i], x[i], 1.0, NULL);
    }

    /* A matrix: shift chain x_i → xdot_{i+1} with gain 1,
     * and feedback from x_i → xdot_{n-1} with gain -a_i */
    for (int i = 0; i < n - 1; i++) {
        sfg_add_branch_real(g, x[i], xdot[i + 1], 1.0, NULL);
    }
    /* Feedback: each x_i → xdot_{n-1} with gain -den[i] */
    for (int i = 0; i < n; i++) {
        sfg_real_t a_i = den[i]; /* den[0] = a_0, ..., den[n-1] = a_{n-1} */
        if (fabs(a_i) > 1e-15) {
            sfg_add_branch_real(g, x[i], xdot[n - 1], -a_i, NULL);
        }
    }

    /* B: input → xdot_0 with gain 1 */
    sfg_add_branch_real(g, input, xdot[0], 1.0, NULL);

    /* C: x_i → output with gain num[i] */
    for (int i = 0; i < n && i <= num_order; i++) {
        if (fabs(num[i]) > 1e-15) {
            sfg_add_branch_real(g, x[i], output, num[i], NULL);
        }
    }

    /* D: if numerator order equals denominator order (biproper),
     * add direct feedthrough */
    if (num_order == den_order) {
        sfg_add_branch_real(g, input, output, num[num_order], NULL);
    }

    sfg_set_input(g, input);
    sfg_set_output(g, output);
    snprintf(g->name, sizeof(g->name), "tf_sfg_n%d", n);

    free(xdot);
    free(x);
    return 0;
}

/* ================================================================
 * L5: SFG → State-Space Conversion
 * ================================================================ */

int sfg_to_state_space(const sfg_graph_t* g,
                        const sfg_node_id_t* state_ids, int n_states,
                        const sfg_node_id_t* input_ids, int n_inputs,
                        const sfg_node_id_t* output_ids, int n_outputs,
                        sfg_ss_system_t* ss) {
    if (!g || !ss) return -1;
    if (n_states <= 0 || n_states > SFG_SS_MAX_ORDER) return -1;

    int ret = sfg_ss_alloc(ss, n_states, n_inputs, n_outputs);
    if (ret != 0) return -1;

    /*
     * Extract A: gain from state j to derivative of state i.
     * In the SFG, each state node x_j has branches to xdot_i nodes
     * with gain A[i][j].
     *
     * We need to identify integrator input nodes (xdot) and
     * integrator output nodes (x). The state_ids should point to
     * the x nodes (integrator outputs).
     *
     * For each branch x_j → xdot_i with gain g, A[i][j] = g.
     * For each branch u_j → xdot_i with gain g, B[i][j] = g.
     * For each branch x_j → y_k with gain g, C[k][j] = g.
     * For each branch u_j → y_k with gain g, D[k][j] = g.
     */

    /* For simplicity, map state IDs to indices */
    int* state_rev = (int*)calloc((size_t)g->num_nodes, sizeof(int));
    if (!state_rev) { sfg_ss_free(ss); return -1; }
    for (int i = 0; i < g->num_nodes; i++) state_rev[i] = -1;
    for (int i = 0; i < n_states; i++) {
        if (state_ids[i] >= 0 && state_ids[i] < g->num_nodes) {
            state_rev[state_ids[i]] = i;
        }
    }

    /* Scan all branches and classify */
    for (int b = 0; b < g->num_branches; b++) {
        sfg_node_id_t src = g->branches[b].src;
        sfg_node_id_t dst = g->branches[b].dst;
        sfg_real_t    gain = creal(g->branches[b].gain);

        /* Check if src or dst is a state node */
        int src_state = (src < (sfg_node_id_t)g->num_nodes)
                        ? state_rev[src] : -1;
        int dst_state = (dst < (sfg_node_id_t)g->num_nodes)
                        ? state_rev[dst] : -1;

        /* Check if src/dst is an input */
        int is_src_input = 0, src_input_idx = -1;
        for (int j = 0; j < n_inputs; j++) {
            if (src == input_ids[j]) { is_src_input = 1; src_input_idx = j; break; }
        }

        /* Check if src/dst is an output */
        int is_dst_output = 0, dst_output_idx = -1;
        for (int j = 0; j < n_outputs; j++) {
            if (dst == output_ids[j]) {
                is_dst_output = 1;
                dst_output_idx = j;
                break;
            }
        }

        /* A matrix: state → state (state transition) */
        if (src_state >= 0 && dst_state >= 0) {
            sfg_ss_set_A(ss, dst_state, src_state,
                         sfg_ss_get_A(ss, dst_state, src_state) + gain);
        }
        /* B matrix: input → state */
        else if (is_src_input && dst_state >= 0) {
            sfg_ss_set_B(ss, dst_state, src_input_idx,
                         sfg_ss_get_B(ss, dst_state, src_input_idx) + gain);
        }
        /* C matrix: state → output */
        else if (src_state >= 0 && is_dst_output) {
            sfg_ss_set_C(ss, dst_output_idx, src_state,
                         sfg_ss_get_C(ss, dst_output_idx, src_state) + gain);
        }
        /* D matrix: input → output */
        else if (is_src_input && is_dst_output) {
            sfg_ss_set_D(ss, dst_output_idx, src_input_idx,
                         sfg_ss_get_D(ss, dst_output_idx, src_input_idx) + gain);
        }
    }

    free(state_rev);
    return 0;
}

/* ================================================================
 * L4: Characteristic Polynomial
 * ================================================================ */

int sfg_ss_characteristic_polynomial(const sfg_ss_system_t* ss,
                                      sfg_real_t* poly, int n) {
    if (!ss || !poly) return -1;
    if (n != ss->n) return -1;

    /*
     * The characteristic polynomial is |sI - A|.
     * For small n, use the Faddeev-Leverrier algorithm.
     *
     * Algorithm (Faddeev-Leverrier):
     *   B_0 = I,  a_{n-1} = -trace(A B_0) / 1
     *   B_1 = A B_0 + a_{n-1} I,  a_{n-2} = -trace(A B_1) / 2
     *   ...
     *   The characteristic polynomial is:
     *   s^n + a_{n-1} s^{n-1} + ... + a_1 s + a_0
     */

    int N = ss->n;
    double* B = (double*)calloc((size_t)(N * N), sizeof(double));
    double* AB = (double*)calloc((size_t)(N * N), sizeof(double));
    if (!B || !AB) { free(B); free(AB); return -1; }

    /* B starts as identity */
    for (int i = 0; i < N; i++) B[i * N + i] = 1.0;

    /* poly[0] = a_0, poly[1] = a_1, ..., poly[N] = 1 (coefficient of s^N) */
    poly[N] = 1.0; /* Leading coefficient */

    for (int k = 1; k <= N; k++) {
        /* Compute A B */
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                for (int l = 0; l < N; l++) {
                    sum += ss->A[i * N + l] * B[l * N + j];
                }
                AB[i * N + j] = sum;
            }
        }

        /* Compute trace */
        double trace = 0.0;
        for (int i = 0; i < N; i++) {
            trace += AB[i * N + i];
        }

        /* a_{N-k} = -trace / k */
        poly[N - k] = -trace / (double)k;

        /* Update B = AB + a_{N-k} I */
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                B[i * N + j] = AB[i * N + j];
                if (i == j) B[i * N + j] += poly[N - k];
            }
        }
    }

    free(B);
    free(AB);
    return 0;
}

/* ================================================================
 * L6: Controllability / Observability via SFG
 * ================================================================ */

int sfg_ss_controllability(const sfg_ss_system_t* ss) {
    if (!ss) return 0;
    /* Build SFG and check path existence from each input to each state */
    sfg_graph_t g;
    sfg_graph_init(&g);
    int ret = sfg_from_state_space(ss, &g);
    if (ret != 0) return 0;

    /* For a simple check: is there a path from at least one input
     * to every state node? In the SFG, each state node xi has
     * incoming integrator branch from xdot_i. Check if there's
     * a path from any input to each xdot_i. */
    /* Simplified: perform DFS reachability */
    int N = g.num_nodes;
    char* reachable = (char*)calloc((size_t)N, sizeof(char));
    if (!reachable) return 0;

    /* Start DFS from input nodes (sources) */
    int* stack = (int*)calloc((size_t)N, sizeof(int));
    if (!stack) { free(reachable); return 0; }
    int stack_top = 0;

    /* Push all source nodes */
    for (int i = 0; i < N; i++) {
        if (g.nodes[i].type == SFG_NODE_SOURCE
            && sfg_node_out_degree(&g, (sfg_node_id_t)i) > 0) {
            stack[stack_top++] = i;
            reachable[i] = 1;
        }
    }

    /* DFS */
    while (stack_top > 0) {
        int node = stack[--stack_top];
        for (int b = 0; b < g.num_branches; b++) {
            if ((int)g.branches[b].src == node) {
                int next = (int)g.branches[b].dst;
                if (!reachable[next]) {
                    reachable[next] = 1;
                    stack[stack_top++] = next;
                }
            }
        }
    }

    /* Check if all state nodes are reachable */
    int controllable = 1;
    for (int i = 0; i < N; i++) {
        if (g.nodes[i].type == SFG_NODE_STATE && !reachable[i]) {
            controllable = 0;
            break;
        }
    }

    free(reachable);
    free(stack);
    return controllable;
}

int sfg_ss_observability(const sfg_ss_system_t* ss) {
    if (!ss) return 0;
    /* Build SFG and check if every state has a path to some output */
    sfg_graph_t g;
    sfg_graph_init(&g);
    int ret = sfg_from_state_space(ss, &g);
    if (ret != 0) return 0;

    int N = g.num_nodes;
    char* reaches_out = (char*)calloc((size_t)N, sizeof(char));
    if (!reaches_out) return 0;

    /* Reverse DFS from output nodes */
    int* stack = (int*)calloc((size_t)N, sizeof(int));
    if (!stack) { free(reaches_out); return 0; }
    int stack_top = 0;

    for (int i = 0; i < N; i++) {
        if (g.nodes[i].type == SFG_NODE_SINK) {
            stack[stack_top++] = i;
            reaches_out[i] = 1;
        }
    }

    /* Reverse DFS: follow branches backwards */
    while (stack_top > 0) {
        int node = stack[--stack_top];
        for (int b = 0; b < g.num_branches; b++) {
            if ((int)g.branches[b].dst == node) {
                int prev = (int)g.branches[b].src;
                if (!reaches_out[prev]) {
                    reaches_out[prev] = 1;
                    stack[stack_top++] = prev;
                }
            }
        }
    }

    int observable = 1;
    for (int i = 0; i < N; i++) {
        if (g.nodes[i].type == SFG_NODE_STATE && !reaches_out[i]) {
            observable = 0;
            break;
        }
    }

    free(reaches_out);
    free(stack);
    return observable;
}

/* ================================================================
 * L5: Transfer Function via SFG
 * ================================================================ */

int sfg_ss_transfer_function(const sfg_ss_system_t* ss,
                              sfg_gain_t* T, int* T_rows, int* T_cols) {
    if (!ss || !T || !T_rows || !T_cols) return -1;

    sfg_graph_t g;
    sfg_graph_init(&g);
    int ret = sfg_from_state_space(ss, &g);
    if (ret != 0) return -1;

    *T_rows = ss->p;
    *T_cols = ss->m;

    /* Compute transfer function for each input-output pair */
    for (int out = 0; out < ss->p; out++) {
        for (int in = 0; in < ss->m; in++) {
            /* Find the corresponding source/sink nodes.
             * In the SFG constructed by sfg_from_state_space,
             * source nodes come first, then xdot, then x, then sink. */
            /* sfg_node_id_t src = (sfg_node_id_t)in;   u_in is node in */
            /* sfg_node_id_t snk = (sfg_node_id_t)(ss->m + 2 * ss->n + out); */
            /* The indexing depends on construction order. For simplicity,
             * use Mason on the entire graph with designated source/sink.
             * Locate the actual source/sink nodes by scanning. */

            sfg_node_id_t actual_src = SFG_INVALID_NODE_ID;
            sfg_node_id_t actual_snk = SFG_INVALID_NODE_ID;

            /* Find input node */
            int src_count = 0;
            for (int i = 0; i < g.num_nodes; i++) {
                if (g.nodes[i].type == SFG_NODE_SOURCE) {
                    if (src_count == in) { actual_src = (sfg_node_id_t)i; break; }
                    src_count++;
                }
            }
            /* Find output node */
            int snk_count = 0;
            for (int i = 0; i < g.num_nodes; i++) {
                if (g.nodes[i].type == SFG_NODE_SINK) {
                    if (snk_count == out) { actual_snk = (sfg_node_id_t)i; break; }
                    snk_count++;
                }
            }

            if (actual_src != SFG_INVALID_NODE_ID
                && actual_snk != SFG_INVALID_NODE_ID) {
                T[out * ss->m + in] = sfg_mason_compute_simple(
                    &g, actual_src, actual_snk);
            } else {
                T[out * ss->m + in] = 0.0;
            }
        }
    }

    return 0;
}

/* ================================================================
 * L3: Display
 * ================================================================ */

void sfg_ss_print(const sfg_ss_system_t* ss) {
    if (!ss) { printf("SS: (null)\n"); return; }

    printf("State-Space System: %s\n", ss->name[0] ? ss->name : "(unnamed)");
    printf("  Order: n=%d, m=%d, p=%d\n", ss->n, ss->m, ss->p);

    printf("A (%dx%d):\n", ss->n, ss->n);
    for (int i = 0; i < ss->n; i++) {
        printf("  ");
        for (int j = 0; j < ss->n; j++) {
            printf("%8.4f ", sfg_ss_get_A(ss, i, j));
        }
        printf("\n");
    }

    if (ss->m > 0) {
        printf("B (%dx%d):\n", ss->n, ss->m);
        for (int i = 0; i < ss->n; i++) {
            printf("  ");
            for (int j = 0; j < ss->m; j++) {
                printf("%8.4f ", sfg_ss_get_B(ss, i, j));
            }
            printf("\n");
        }
    }

    if (ss->p > 0) {
        printf("C (%dx%d):\n", ss->p, ss->n);
        for (int i = 0; i < ss->p; i++) {
            printf("  ");
            for (int j = 0; j < ss->n; j++) {
                printf("%8.4f ", sfg_ss_get_C(ss, i, j));
            }
            printf("\n");
        }
    }

    if (ss->m > 0 && ss->p > 0) {
        printf("D (%dx%d):\n", ss->p, ss->m);
        for (int i = 0; i < ss->p; i++) {
            printf("  ");
            for (int j = 0; j < ss->m; j++) {
                printf("%8.4f ", sfg_ss_get_D(ss, i, j));
            }
            printf("\n");
        }
    }
}
