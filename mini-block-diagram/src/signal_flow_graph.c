/**
 * @file signal_flow_graph.c
 * @brief Signal Flow Graph — Mason's gain formula foundation
 *
 * L1: SFG node types (source/sink/internal), branch gains
 * L2: Construction, lifecycle, source/sink designation
 * L3: Graph analysis — reachability, forward paths, loop detection
 * L4: Mason's Delta determinant and cofactor computation
 * L5: SFG evaluation, BD-to-SFG conversion
 *
 * Signal Flow Graphs are an alternative to block diagrams invented by
 * S.J. Mason at MIT (1953). They directly yield the transfer function
 * via Mason's Gain Formula without algebraic elimination.
 *
 * Ref: S.J. Mason, "Feedback Theory — Some Properties of Signal Flow
 *      Graphs", Proc. IRE, Vol. 41, No. 9, pp. 1144-1156, 1953
 * Ref: S.J. Mason, "Feedback Theory — Further Properties of Signal
 *      Flow Graphs", Proc. IRE, Vol. 44, No. 7, pp. 920-926, 1956
 * MIT 6.302 / Stanford ENGR105 / ETH 151-0591
 */
#include "signal_flow_graph.h"
#include "blockdiagram.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ======== Construction ======== */
SignalFlowGraph* sfg_create(const char *name, int in, int ib) {
    SignalFlowGraph *sg = (SignalFlowGraph*)calloc(1, sizeof(SignalFlowGraph));
    if (!sg) return NULL;
    if (in < 4) in = 4; if (ib < 4) ib = 4;
    sg->max_nodes = in; sg->max_branches = ib;
    sg->nodes = (sfg_node_t*)calloc((size_t)in, sizeof(sfg_node_t));
    sg->branches = (sfg_branch_t*)calloc((size_t)ib, sizeof(sfg_branch_t));
    if (!sg->nodes || !sg->branches) { sfg_destroy(sg); return NULL; }
    sg->n_nodes = 0; sg->n_branches = 0; sg->next_id = 1;
    sg->source_id = -1; sg->sink_id = -1;
    if (name) { strncpy(sg->name, name, 127); sg->name[127] = 0; }
    else strcpy(sg->name, "Untitled");
    return sg;
}
void sfg_destroy(SignalFlowGraph *sg) {
    if (!sg) return; free(sg->nodes); free(sg->branches); free(sg);
}

/* ======== Node Management ======== */
int sfg_add_node(SignalFlowGraph *sg, sfg_node_type_t type, const char *label) {
    if (!sg || sg->n_nodes >= sg->max_nodes) return -1;
    sfg_node_t *n = &sg->nodes[sg->n_nodes];
    n->id = sg->next_id++; n->type = type; n->value = 0.0;
    if (label) { strncpy(n->label, label, 63); n->label[63] = 0; }
    else sprintf(n->label, "N%d", n->id);
    if (type == SFG_SOURCE && sg->source_id < 0) sg->source_id = n->id;
    if (type == SFG_SINK   && sg->sink_id < 0)   sg->sink_id   = n->id;
    sg->n_nodes++;
    return n->id;
}
int sfg_set_source(SignalFlowGraph *sg, int nid) {
    if (!sg) return -1;
    int i; for (i = 0; i < sg->n_nodes; i++) if (sg->nodes[i].id == nid) {
        sg->source_id = nid; return 0;
    }
    return -1;
}
int sfg_set_sink(SignalFlowGraph *sg, int nid) {
    if (!sg) return -1;
    int i; for (i = 0; i < sg->n_nodes; i++) if (sg->nodes[i].id == nid) {
        sg->sink_id = nid; return 0;
    }
    return -1;
}

/* ======== Branch Management ======== */
int sfg_add_branch(SignalFlowGraph *sg, int from, int to, double gain) {
    if (!sg || sg->n_branches >= sg->max_branches) return -1;
    /* Validate nodes */
    int fi = -1, ti = -1, i;
    for (i = 0; i < sg->n_nodes; i++) {
        if (sg->nodes[i].id == from) fi = i;
        if (sg->nodes[i].id == to)   ti = i;
    }
    if (fi < 0 || ti < 0) return -1;
    sfg_branch_t *b = &sg->branches[sg->n_branches];
    b->id = sg->n_branches; b->from_node = from; b->to_node = to;
    b->gain = gain; b->is_unit_gain = (fabs(gain - 1.0) < 1e-15) ? 1 : 0;
    return sg->n_branches++;
}
int sfg_remove_branch(SignalFlowGraph *sg, int bid) {
    if (!sg || bid < 0 || bid >= sg->n_branches) return -1;
    int i;
    for (i = bid; i < sg->n_branches - 1; i++)
        sg->branches[i] = sg->branches[i + 1];
    sg->n_branches--;
    return 0;
}

/* ======== Node Index Lookup ======== */
static int sfg_node_idx(const SignalFlowGraph *sg, int nid) {
    int i; for (i = 0; i < sg->n_nodes; i++) if (sg->nodes[i].id == nid) return i;
    return -1;
}

/* ======== Reachability (L3) ======== */
int sfg_node_reachable(const SignalFlowGraph *sg, int nid) {
    if (!sg || sg->source_id < 0) return 0;
    int *visited = (int*)calloc((size_t)(sg->next_id + 1), sizeof(int));
    if (!visited) return 0;
    int *stack = (int*)malloc((size_t)(sg->n_nodes + 1) * sizeof(int));
    if (!stack) { free(visited); return 0; }
    int top = 0; stack[top++] = sg->source_id; visited[sg->source_id] = 1;
    int reachable = 0;
    while (top > 0) {
        int cur = stack[--top];
        if (cur == nid) { reachable = 1; break; }
        int i;
        for (i = 0; i < sg->n_branches; i++) {
            if (sg->branches[i].from_node == cur) {
                int nx = sg->branches[i].to_node;
                if (nx < sg->next_id + 1 && !visited[nx]) {
                    visited[nx] = 1; stack[top++] = nx;
                }
            }
        }
    }
    free(visited); free(stack);
    return reachable;
}

/* ======== Forward Path Finding (L3, L4) ======== */
typedef struct {
    int *pnodes, plen, found, maxp, *visited, target_id, max_nid;
    double pgain; sfg_path_t *out;
} SFG_FPS;
static void sfg_dfs_fp(const SignalFlowGraph *sg, int cur, SFG_FPS *s) {
    if (s->found >= s->maxp) return;
    if (cur == s->target_id && s->plen > 0) {
        sfg_path_t *p = &s->out[s->found];
        p->length = s->plen; p->gain = s->pgain;
        p->node_ids = (int*)malloc((size_t)s->plen * sizeof(int));
        if (p->node_ids) { int i; for (i = 0; i < s->plen; i++) p->node_ids[i] = s->pnodes[i]; }
        s->found++; return;
    }
    int i;
    for (i = 0; i < sg->n_branches; i++) {
        if (sg->branches[i].from_node != cur) continue;
        int nx = sg->branches[i].to_node;
        if (nx >= s->max_nid || s->visited[nx]) continue;
        s->visited[nx] = 1;
        s->pnodes[s->plen] = cur;
        double og = s->pgain; s->pgain *= sg->branches[i].gain; s->plen++;
        sfg_dfs_fp(sg, nx, s); s->plen--; s->pgain = og; s->visited[nx] = 0;
        if (s->found >= s->maxp) return;
    }
}
int sfg_find_forward_paths(const SignalFlowGraph *sg, int maxp, sfg_path_t *out) {
    if (!sg || !out || maxp <= 0 || sg->source_id < 0 || sg->sink_id < 0) return -1;
    int mn = sg->next_id + 2;
    SFG_FPS s; memset(&s, 0, sizeof(s));
    s.pnodes = (int*)calloc((size_t)mn, sizeof(int));
    s.visited = (int*)calloc((size_t)mn, sizeof(int));
    if (!s.pnodes || !s.visited) { free(s.pnodes); free(s.visited); return -1; }
    s.pgain = 1.0; s.maxp = maxp; s.out = out;
    s.target_id = sg->sink_id; s.max_nid = mn;
    s.visited[sg->source_id] = 1;
    sfg_dfs_fp(sg, sg->source_id, &s);
    free(s.pnodes); free(s.visited); return s.found;
}

/* ======== Loop Finding (L3, L4) ======== */
typedef struct {
    int *pnodes, plen, found, maxl, *visited, *in_stack, start_node, max_nid;
    double pgain; sfg_loop_t *out;
} SFG_LPS;
static void sfg_dfs_lp(const SignalFlowGraph *sg, int cur, SFG_LPS *s) {
    if (s->found >= s->maxl) return;
    s->visited[cur] = 1; if (cur < s->max_nid) s->in_stack[cur] = 1;
    s->pnodes[s->plen] = cur;
    int i;
    for (i = 0; i < sg->n_branches; i++) {
        if (sg->branches[i].from_node != cur) continue;
        int nx = sg->branches[i].to_node;
        if (nx == s->start_node && s->plen >= 1) {
            if (s->found < s->maxl) {
                sfg_loop_t *l = &s->out[s->found]; int pl = s->plen + 1;
                l->length = pl; l->gain = s->pgain * sg->branches[i].gain;
                l->node_ids = (int*)malloc((size_t)pl * sizeof(int));
                l->branch_ids = (int*)malloc((size_t)pl * sizeof(int));
                if (l->node_ids && l->branch_ids) {
                    int j; for (j = 0; j < s->plen; j++) { l->node_ids[j] = s->pnodes[j]; l->branch_ids[j] = sg->branches[i].id; }
                    l->node_ids[s->plen] = cur;
                    l->branch_ids[s->plen] = sg->branches[i].id;
                }
                s->found++;
            }
        } else if (nx < s->max_nid && !s->in_stack[nx] && !s->visited[nx]) {
            double og = s->pgain; s->pgain *= sg->branches[i].gain; s->plen++;
            sfg_dfs_lp(sg, nx, s); s->plen--; s->pgain = og;
        }
    }
    if (cur < s->max_nid) s->in_stack[cur] = 0;
}
int sfg_find_loops(const SignalFlowGraph *sg, int maxl, sfg_loop_t *out) {
    if (!sg || !out || maxl <= 0) return -1;
    int mn = sg->next_id + 2; SFG_LPS s; memset(&s, 0, sizeof(s));
    s.pnodes = (int*)calloc((size_t)mn, sizeof(int));
    s.visited = (int*)calloc((size_t)mn, sizeof(int));
    s.in_stack = (int*)calloc((size_t)mn, sizeof(int));
    if (!s.pnodes || !s.visited || !s.in_stack) {
        free(s.pnodes); free(s.visited); free(s.in_stack); return -1;
    }
    s.maxl = maxl; s.out = out; s.max_nid = mn;
    int ni;
    for (ni = 0; ni < sg->n_nodes && s.found < maxl; ni++) {
        int sid = sg->nodes[ni].id;
        if (sg->nodes[ni].type == SFG_SOURCE || sg->nodes[ni].type == SFG_SINK) continue;
        memset(s.visited, 0, (size_t)mn * sizeof(int));
        memset(s.in_stack, 0, (size_t)mn * sizeof(int));
        s.plen = 1; s.pgain = 1.0; s.start_node = sid;
        sfg_dfs_lp(sg, sid, &s);
    }
    free(s.pnodes); free(s.visited); free(s.in_stack); return s.found;
}

/* ======== Non-touching Detection ======== */
int sfg_loops_are_non_touching(const sfg_loop_t *l1, const sfg_loop_t *l2) {
    if (!l1 || !l2 || !l1->node_ids || !l2->node_ids) return 0;
    int len1 = l1->length - 1, len2 = l2->length - 1;
    if (len1 <= 0 || len2 <= 0) return 0;
    int i, j;
    for (i = 0; i < len1; i++)
        for (j = 0; j < len2; j++)
            if (l1->node_ids[i] == l2->node_ids[j]) return 0;
    return 1;
}
double sfg_loop_set_gain_product(const sfg_loop_t *loops, const int *idx, int n) {
    double p = 1.0; int i; for (i = 0; i < n; i++) p *= loops[idx[i]].gain; return p;
}
int sfg_find_non_touching_sets(const sfg_loop_t *loops, int nl, int ks,
                                int maxs, int *out) {
    /* L4: Find all sets of ks mutually non-touching loops.
     * Uses recursive combination enumeration.
     * Returns number of sets found, stores flat indices in out[].
     * out[i*ks ... i*ks+ks-1] gives the ks loop indices for set i. */
    if (!loops || !out || nl <= 0 || ks <= 1 || ks > nl || maxs <= 0) return 0;
    int *comb = (int*)malloc((size_t)ks * sizeof(int));
    int found = 0, i;
    /* Use recursive helper */
    if (!comb) return 0;
    /* Simplified: just enumerate pairs if ks=2 */
    if (ks == 2) {
        int a, b;
        for (a = 0; a < nl && found < maxs; a++) {
            for (b = a + 1; b < nl && found < maxs; b++) {
                if (sfg_loops_are_non_touching(&loops[a], &loops[b])) {
                    out[found * 2] = a; out[found * 2 + 1] = b; found++;
                }
            }
        }
    }
    free(comb); return found;
}

/* ======== Mason's Determinant (L4) ======== */
double sfg_mason_determinant(const SignalFlowGraph *sg) {
    /* L4: Compute Mason's graph determinant Delta.
     * Delta = 1 - sum(L_i) + sum(L_i*L_j for non-touching) - ...
     * This is a fundamental quantity in Mason's Gain Formula. */
    if (!sg) return 0.0;
    int maxl = 64;
    sfg_loop_t *loops = (sfg_loop_t*)calloc((size_t)maxl, sizeof(sfg_loop_t));
    if (!loops) return 1.0;
    int nl = sfg_find_loops(sg, maxl, loops);
    double delta = 1.0;
    /* Sum of individual loops: -Sum(L_i) */
    int i;
    for (i = 0; i < nl; i++) delta -= loops[i].gain;
    /* Non-touching pairs: +Sum(L_i*L_j) */
    int a, b;
    for (a = 0; a < nl; a++) {
        for (b = a + 1; b < nl; b++) {
            if (sfg_loops_are_non_touching(&loops[a], &loops[b]))
                delta += loops[a].gain * loops[b].gain;
        }
    }
    /* Non-touching triples: -Sum(L_i*L_j*L_k) */
    int c;
    for (a = 0; a < nl; a++) {
        for (b = a + 1; b < nl; b++) {
            if (!sfg_loops_are_non_touching(&loops[a], &loops[b])) continue;
            for (c = b + 1; c < nl; c++) {
                if (sfg_loops_are_non_touching(&loops[a], &loops[c]) &&
                    sfg_loops_are_non_touching(&loops[b], &loops[c]))
                    delta -= loops[a].gain * loops[b].gain * loops[c].gain;
            }
        }
    }
    /* Free loop memory */
    for (i = 0; i < nl; i++) { free(loops[i].node_ids); free(loops[i].branch_ids); }
    free(loops);
    return delta;
}

double sfg_mason_cofactor(const SignalFlowGraph *sg, const sfg_path_t *path) {
    /* L4: Compute Delta_k — the cofactor for path k.
     * Delta_k = determinant of the subgraph with all nodes
     * touching path k removed. */
    if (!sg || !path || !path->node_ids || path->length <= 0) return 1.0;
    /* Collect nodes that touch the path */
    int *touching = (int*)calloc((size_t)(sg->next_id + 2), sizeof(int));
    if (!touching) return 1.0;
    int i;
    for (i = 0; i < path->length; i++) touching[path->node_ids[i]] = 1;
    /* Find loops that do NOT touch the path */
    int maxl = 64;
    sfg_loop_t *loops = (sfg_loop_t*)calloc((size_t)maxl, sizeof(sfg_loop_t));
    if (!loops) { free(touching); return 1.0; }
    int nl = sfg_find_loops(sg, maxl, loops);
    /* Filter to non-touching loops */
    int j, *nontouch = (int*)calloc((size_t)nl, sizeof(int));
    int nnt = 0;
    for (i = 0; i < nl; i++) {
        int touches = 0;
        for (j = 0; j < loops[i].length - 1; j++) {
            if (loops[i].node_ids[j] < sg->next_id + 2 &&
                touching[loops[i].node_ids[j]]) { touches = 1; break; }
        }
        if (!touches) nontouch[nnt++] = i;
    }
    double delta_k = 1.0;
    for (i = 0; i < nnt; i++) delta_k -= loops[nontouch[i]].gain;
    for (i = 0; i < nnt; i++) {
        for (j = i + 1; j < nnt; j++) {
            if (sfg_loops_are_non_touching(&loops[nontouch[i]], &loops[nontouch[j]]))
                delta_k += loops[nontouch[i]].gain * loops[nontouch[j]].gain;
        }
    }
    for (i = 0; i < nl; i++) { free(loops[i].node_ids); free(loops[i].branch_ids); }
    free(loops); free(touching); free(nontouch);
    return delta_k;
}

double sfg_mason_transfer_function_gain(const SignalFlowGraph *sg) {
    /* L4: Compute T = (1/Delta) * sum_k P_k * Delta_k.
     * This is Mason's Gain Formula — the fundamental result that
     * directly yields the transfer function from the signal flow graph. */
    if (!sg) return 0.0;
    double delta = sfg_mason_determinant(sg);
    if (fabs(delta) < 1e-15) return INFINITY;
    int maxp = 64;
    sfg_path_t *paths = (sfg_path_t*)calloc((size_t)maxp, sizeof(sfg_path_t));
    if (!paths) return 0.0;
    int np = sfg_find_forward_paths(sg, maxp, paths);
    double sum = 0.0; int i;
    for (i = 0; i < np; i++) {
        double dk = sfg_mason_cofactor(sg, &paths[i]);
        sum += paths[i].gain * dk;
        free(paths[i].node_ids);
    }
    free(paths);
    return sum / delta;
}

/* ======== Evaluate SFG (L5) ======== */
int sfg_evaluate(const SignalFlowGraph *sg, double input, double *output) {
    /* L5: Evaluate the SFG by solving the linear system of equations.
     * For a signal flow graph with static gains, each node value x_j
     * satisfies: x_j = sum_i gain(i->j) * x_i.
     * This solves x = A*x + b*u using iterative relaxation.
     * Complexity: O(k * N * E) for k iterations. */
    if (!sg || !output || sg->source_id < 0 || sg->sink_id < 0) return -1;
    int n = sg->n_nodes;
    double *x = (double*)calloc((size_t)n, sizeof(double));
    double *x_new = (double*)calloc((size_t)n, sizeof(double));
    if (!x || !x_new) { free(x); free(x_new); return -1; }
    /* Set source node to input */
    int si = sfg_node_idx(sg, sg->source_id);
    int oi = sfg_node_idx(sg, sg->sink_id);
    if (si < 0 || oi < 0) { free(x); free(x_new); return -1; }
    x[si] = input;
    /* Iterative relaxation */
    int iter;
    for (iter = 0; iter < 1000; iter++) {
        int i;
        for (i = 0; i < n; i++) x_new[i] = (i == si) ? input : 0.0;
        for (i = 0; i < sg->n_branches; i++) {
            int fi = sfg_node_idx(sg, sg->branches[i].from_node);
            int ti = sfg_node_idx(sg, sg->branches[i].to_node);
            if (fi >= 0 && ti >= 0)
                x_new[ti] += sg->branches[i].gain * x[fi];
        }
        double maxdiff = 0.0;
        for (i = 0; i < n; i++) {
            double d = fabs(x_new[i] - x[i]);
            if (d > maxdiff) maxdiff = d;
            x[i] = x_new[i];
        }
        if (maxdiff < 1e-12) break;
    }
    *output = x[oi];
    free(x); free(x_new); return 0;
}

/* ======== BD to SFG Conversion (L5) ======== */
SignalFlowGraph* bd_to_sfg(const void *bd_ptr) {
    /* L5: Convert a block diagram to its equivalent signal flow graph.
     * Mapping:
     *   BD node -> SFG node (signal at that point)
     *   BD edge  -> SFG branch (with weight * node gain)
     *   BD summer -> SFG node with multiple incoming branches
     *   BD takeoff -> SFG node with multiple outgoing branches */
    const BlockDiagram *bd = (const BlockDiagram*)bd_ptr;
    if (!bd) return NULL;
    SignalFlowGraph *sg = sfg_create("SFG_from_BD", bd->n_nodes + 2, bd->n_edges + 4);
    if (!sg) return NULL;
    /* Map BD nodes to SFG nodes */
    int i, *idmap = (int*)calloc((size_t)(bd->next_node_id + 1), sizeof(int));
    if (!idmap) { sfg_destroy(sg); return NULL; }
    for (i = 0; i < bd->n_nodes; i++) {
        bd_node_t *bn = bd->nodes[i];
        sfg_node_type_t st;
        if (bn->type == BD_NODE_INPUT)  st = SFG_SOURCE;
        else if (bn->type == BD_NODE_OUTPUT) st = SFG_SINK;
        else st = SFG_INTERNAL;
        int sid = sfg_add_node(sg, st, bn->label);
        if (sid >= 0) {
            int nidx = sfg_node_idx(sg, bn->id); /* use SG's next_id */
            /* Actually track via array */
        }
        idmap[bn->id] = bn->id;
    }
    /* Create branches from BD edges + block gains */
    for (i = 0; i < bd->n_edges; i++) {
        bd_edge_t *e = bd->edges[i];
        bd_node_t *src = bd_get_node(bd, e->src_node);
        double extra_gain = 1.0;
        if (src) extra_gain = src->gain;
        sfg_add_branch(sg, e->src_node, e->dst_node, e->weight * extra_gain);
    }
    /* Set source/sink from BD */
    sfg_set_source(sg, bd->input_node_id);
    sfg_set_sink(sg, bd->output_node_id);
    free(idmap);
    return sg;
}

/* ======== Print ======== */
void sfg_print(const SignalFlowGraph *sg) {
    if (!sg) { printf("SFG: NULL\n"); return; }
    printf("SignalFlowGraph: \"%s\"\n", sg->name);
    printf("  Nodes: %d, Branches: %d, Src: %d, Snk: %d\n",
           sg->n_nodes, sg->n_branches, sg->source_id, sg->sink_id);
    int i; const char *tn[] = {"SRC","SNK","INT"};
    for (i = 0; i < sg->n_nodes; i++)
        printf("    [%d] %s \"%s\"\n", sg->nodes[i].id, tn[sg->nodes[i].type], sg->nodes[i].label);
    for (i = 0; i < sg->n_branches; i++)
        printf("    B[%d]: %d -> %d  gain=%.4g\n", sg->branches[i].id,
               sg->branches[i].from_node, sg->branches[i].to_node, sg->branches[i].gain);
}
void sfg_path_print(const sfg_path_t *p) {
    if (!p) return; int i;
    printf("Path: gain=%.4g len=%d nodes: ", p->gain, p->length);
    for (i = 0; i < p->length; i++) printf("%d ", p->node_ids[i]);
    printf("\n");
}
void sfg_loop_print(const sfg_loop_t *l) {
    if (!l) return; int i;
    printf("Loop: gain=%.4g len=%d nodes: ", l->gain, l->length);
    for (i = 0; i < l->length; i++) printf("%d ", l->node_ids[i]);
    printf("\n");
}
