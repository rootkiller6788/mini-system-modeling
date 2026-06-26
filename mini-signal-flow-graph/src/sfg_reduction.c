/**
 * sfg_reduction.c — Signal Flow Graph Topological Reduction
 *
 * Implements the five fundamental SFG reduction rules (R1-R5) and
 * an automated reduction procedure that simplifies an SFG to a
 * single equivalent transmittance between input and output.
 *
 * Reduction Rules:
 *   R1. Series (cascade) branch elimination
 *   R2. Parallel branch combination
 *   R3. Self-loop elimination
 *   R4. Node absorption (intermediate node elimination)
 *   R5. Node splitting (untangling paths)
 *
 * Knowledge coverage:
 *   L5: Five reduction rules as engineering methods
 *   L4: Algebraic equivalence of reduction to Cramer's rule
 *   L6: Verification against Mason's formula
 */

#include "sfg_reduction.h"
#include "sfg_path.h"
#include "sfg_mason.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * L1: Result Structure Management
 * ================================================================ */

void sfg_reduction_result_init(sfg_reduction_result_t* result) {
    if (!result) return;
    memset(result, 0, sizeof(sfg_reduction_result_t));
    result->overall_gain = 0.0;
    result->success = 0;
}

/* ================================================================
 * L5: Rule R1 — Series (Cascade) Branch Elimination
 * ================================================================ */

int sfg_reduce_series(sfg_graph_t* g) {
    if (!g) return 0;
    int eliminated = 0;

    /*
     * R1: For each node B that has exactly one incoming branch A→B
     * and exactly one outgoing branch B→C, eliminate node B and
     * create a direct branch A→C with gain = g_AB * g_BC.
     *
     * Exception: Do not eliminate source or sink nodes.
     */
    for (int b_idx = 0; b_idx < g->num_nodes; b_idx++) {
        sfg_node_id_t mid = (sfg_node_id_t)b_idx;

        /* Skip source and sink nodes */
        if (mid == g->input_node || mid == g->output_node) continue;

        int in_deg  = sfg_node_in_degree(g, mid);
        int out_deg = sfg_node_out_degree(g, mid);

        if (in_deg == 1 && out_deg == 1) {
            /* Find the incoming and outgoing branches */
            sfg_node_id_t pred = SFG_INVALID_NODE_ID;
            sfg_node_id_t succ = SFG_INVALID_NODE_ID;
            sfg_gain_t    g_in = 0.0, g_out = 0.0;

            for (int k = 0; k < g->num_branches; k++) {
                if (g->branches[k].dst == mid) {
                    pred = g->branches[k].src;
                    g_in = g->branches[k].gain;
                }
                if (g->branches[k].src == mid) {
                    succ = g->branches[k].dst;
                    g_out = g->branches[k].gain;
                }
            }

            if (pred != SFG_INVALID_NODE_ID
                && succ != SFG_INVALID_NODE_ID) {
                /* Create or update branch pred → succ */
                const sfg_branch_t* existing =
                    sfg_find_branch(g, pred, succ);
                sfg_gain_t new_gain = g_in * g_out;

                if (existing) {
                    /* Parallel combination with existing */
                    sfg_add_branch(g, pred, succ,
                                   existing->gain + new_gain, NULL);
                } else {
                    sfg_add_branch(g, pred, succ, new_gain, NULL);
                }

                /* Remove the series node */
                sfg_remove_node(g, mid);
                eliminated++;
                /* Restart scan since indices shifted */
                b_idx = -1;
            }
        }
    }

    return eliminated;
}

/* ================================================================
 * L5: Rule R2 — Parallel Branch Combination
 * ================================================================ */

int sfg_reduce_parallel(sfg_graph_t* g) {
    if (!g) return 0;
    int combined = 0;

    /*
     * R2: For every pair of nodes (i, j), combine all parallel
     * branches i→j into a single branch with gain = Σ gains.
     */
    for (int i = 0; i < g->num_nodes; i++) {
        for (int j = 0; j < g->num_nodes; j++) {
            /* Collect all branches from i to j */
            sfg_gain_t sum_gain = 0.0;
            int       count     = 0;
            int       first_idx = -1;

            for (int k = 0; k < g->num_branches; k++) {
                if (g->branches[k].src == i && g->branches[k].dst == j) {
                    sum_gain += g->branches[k].gain;
                    if (count == 0) first_idx = k;
                    count++;
                }
            }

            if (count > 1 && first_idx >= 0) {
                /* Update the first branch with the sum and remove others */
                g->branches[first_idx].gain = sum_gain;

                /* Remove other parallel branches (scan backwards) */
                int k = g->num_branches - 1;
                while (k > first_idx) {
                    if (g->branches[k].src == i && g->branches[k].dst == j) {
                        if (k < g->num_branches - 1) {
                            g->branches[k] = g->branches[g->num_branches - 1];
                        }
                        g->num_branches--;
                    }
                    k--;
                }
                combined++;
            }
        }
    }
    return combined;
}

/* ================================================================
 * L5: Rule R3 — Self-Loop Elimination
 * ================================================================ */

int sfg_reduce_self_loop(sfg_graph_t* g) {
    if (!g) return 0;
    int eliminated = 0;

    for (int n = 0; n < g->num_nodes; n++) {
        /* Find self-loop on this node */
        sfg_gain_t self_gain = 0.0;
        int        has_self  = 0;

        /* Collect and remove all self-loops on node n */
        int k = 0;
        while (k < g->num_branches) {
            if (g->branches[k].src == n && g->branches[k].dst == n) {
                self_gain += g->branches[k].gain;
                has_self = 1;
                /* Remove this self-loop */
                if (k < g->num_branches - 1) {
                    g->branches[k] = g->branches[g->num_branches - 1];
                }
                g->num_branches--;
            } else {
                k++;
            }
        }

        if (has_self) {
            /*
             * Divide all incoming branches (except self-loops, already
             * removed) by (1 - self_gain).
             *
             * Rationale: x_n = Σ(g_in × x_src) + self_gain × x_n
             *           → x_n (1 - self_gain) = Σ(g_in × x_src)
             *           → x_n = Σ(g_in / (1 - self_gain)) × x_src
             */
            sfg_gain_t denom = 1.0 + 0.0 * I - self_gain;

            /* Check for division by zero (self_gain == 1 means infinite gain) */
            if (cabs(denom) < 1e-15) {
                /* Restore self-loop and skip */
                sfg_add_branch(g, (sfg_node_id_t)n, (sfg_node_id_t)n,
                               self_gain, NULL);
                continue;
            }

            for (int b = 0; b < g->num_branches; b++) {
                if (g->branches[b].dst == n) {
                    g->branches[b].gain /= denom;
                }
            }
            eliminated++;
        }
    }
    return eliminated;
}

/* ================================================================
 * L5: Rule R4 — Node Absorption
 * ================================================================ */

int sfg_reduce_node_absorb(sfg_graph_t* g, sfg_node_id_t node_id) {
    if (!g) return -1;
    if (node_id < 0 || node_id >= g->num_nodes) return -1;

    /* Cannot absorb source or sink nodes */
    if (node_id == g->input_node || node_id == g->output_node) return -1;

    /* Eliminate self-loops on this node first */
    sfg_graph_t temp;
    sfg_graph_copy(&temp, g);

    /* Find all incoming and outgoing branches for this node */
    /* We need to iterate carefully to handle dynamic changes */

    /* Collect incoming branches */
    struct { sfg_node_id_t src; sfg_gain_t gain; } incoming[SFG_MAX_BRANCHES];
    int n_in = 0;
    /* Collect outgoing branches */
    struct { sfg_node_id_t dst; sfg_gain_t gain; } outgoing[SFG_MAX_BRANCHES];
    int n_out = 0;

    for (int b = 0; b < g->num_branches; b++) {
        if (g->branches[b].dst == node_id
            && g->branches[b].src != node_id) {
            if (n_in < SFG_MAX_BRANCHES) {
                incoming[n_in].src  = g->branches[b].src;
                incoming[n_in].gain = g->branches[b].gain;
                n_in++;
            }
        }
        if (g->branches[b].src == node_id
            && g->branches[b].dst != node_id) {
            if (n_out < SFG_MAX_BRANCHES) {
                outgoing[n_out].dst  = g->branches[b].dst;
                outgoing[n_out].gain = g->branches[b].gain;
                n_out++;
            }
        }
    }

    /* For each (incoming, outgoing) pair, create a new branch */
    for (int i = 0; i < n_in; i++) {
        for (int j = 0; j < n_out; j++) {
            sfg_gain_t new_gain = incoming[i].gain * outgoing[j].gain;
            const sfg_branch_t* existing =
                sfg_find_branch(g, incoming[i].src, outgoing[j].dst);
            if (existing) {
                /* Combine with existing parallel branch */
                sfg_add_branch(g, incoming[i].src, outgoing[j].dst,
                               existing->gain + new_gain, NULL);
            } else {
                sfg_add_branch(g, incoming[i].src, outgoing[j].dst,
                               new_gain, NULL);
            }
        }
    }

    /* Remove the absorbed node (this also removes its incident branches) */
    sfg_remove_node(g, node_id);

    return 0;
}

/* ================================================================
 * L5: Rule R5 — Node Splitting
 * ================================================================ */

int sfg_reduce_node_split(sfg_graph_t* g, sfg_node_id_t node_id) {
    if (!g) return -1;
    if (node_id < 0 || node_id >= g->num_nodes) return -1;

    /*
     * R5: Duplicate a node to isolate non-touching paths.
     * Create a copy of the node. All incoming branches go to the
     * original. Select outgoing branches go to the duplicate.
     *
     * This is used when paths that are touching through this node
     * need to be made non-touching.
     */

    if (g->num_nodes >= SFG_MAX_NODES) return -1;

    /* Create duplicate node */
    sfg_node_id_t dup = sfg_add_node(g, g->nodes[node_id].type, NULL);

    /* Redirect half the outgoing branches to the duplicate.
     * A simple heuristic: redirect the first half. */
    int n_out = sfg_node_out_degree(g, node_id);
    int redirect_count = n_out / 2;
    int redirected = 0;

    for (int b = 0; b < g->num_branches && redirected < redirect_count; b++) {
        if (g->branches[b].src == node_id) {
            g->branches[b].src = dup;
            redirected++;
        }
    }

    return 0;
}

/* ================================================================
 * L5: Automated Reduction
 * ================================================================ */

int sfg_is_reducible(const sfg_graph_t* g, sfg_node_id_t source,
                      sfg_node_id_t sink) {
    if (!g) return 0;
    /* A graph is reducible if source and sink are connected */
    sfg_path_list_t paths;
    sfg_path_list_clear(&paths);
    int n = sfg_find_forward_paths(g, source, sink, &paths);
    return (n > 0);
}

int sfg_reduce_interactive(sfg_graph_t* g, sfg_reduction_step_t* step) {
    if (!g || !step) return -1;
    memset(step, 0, sizeof(sfg_reduction_step_t));

    /* Priority order: parallel → self-loop → series → node absorb */
    int count;

    count = sfg_reduce_parallel(g);
    if (count > 0) {
        step->operation = SFG_REDUCE_PARALLEL;
        snprintf(step->description, sizeof(step->description),
                 "Combined %d parallel branch group(s)", count);
        return 1;
    }

    count = sfg_reduce_self_loop(g);
    if (count > 0) {
        step->operation = SFG_REDUCE_SELF_LOOP;
        snprintf(step->description, sizeof(step->description),
                 "Eliminated %d self-loop(s)", count);
        return 1;
    }

    count = sfg_reduce_series(g);
    if (count > 0) {
        step->operation = SFG_REDUCE_SERIES;
        snprintf(step->description, sizeof(step->description),
                 "Eliminated %d series node(s)", count);
        return 1;
    }

    /* Find an internal node to absorb (not source/sink) */
    for (int n = 0; n < g->num_nodes; n++) {
        if ((sfg_node_id_t)n != g->input_node
            && (sfg_node_id_t)n != g->output_node) {
            int ret = sfg_reduce_node_absorb(g, (sfg_node_id_t)n);
            if (ret == 0) {
                step->operation = SFG_REDUCE_NODE_ABSORB;
                snprintf(step->description, sizeof(step->description),
                         "Absorbed internal node %d", n);
                return 1;
            }
        }
    }

    return 0; /* No reduction possible */
}

int sfg_reduce_to_single(sfg_graph_t* g, sfg_node_id_t source,
                          sfg_node_id_t sink,
                          sfg_reduction_result_t* result) {
    if (!g || !result) return -1;
    if (source < 0 || source >= g->num_nodes) return -1;
    if (sink < 0 || sink >= g->num_nodes) return -1;
    sfg_reduction_result_init(result);

    /* Copy the graph for reduction */
    sfg_graph_copy(&result->reduced_graph, g);
    sfg_graph_t* rg = &result->reduced_graph;

    /* Track source and sink as their types (IDs may shift) */
    sfg_node_id_t cur_source = source;
    sfg_node_id_t cur_sink   = sink;
    rg->input_node  = cur_source;
    rg->output_node = cur_sink;

    int max_iterations = 1000;
    for (int iter = 0; iter < max_iterations; iter++) {
        /* Refresh source/sink by scanning for their types
         * (in case IDs shifted due to node removal) */
        if (!sfg_node_exists(rg, cur_source)
            || rg->nodes[(int)cur_source].type != SFG_NODE_SOURCE) {
            cur_source = SFG_INVALID_NODE_ID;
            for (int i = 0; i < rg->num_nodes; i++) {
                if (rg->nodes[i].type == SFG_NODE_SOURCE) {
                    cur_source = (sfg_node_id_t)i;
                    break;
                }
            }
        }
        if (!sfg_node_exists(rg, cur_sink)
            || rg->nodes[(int)cur_sink].type != SFG_NODE_SINK) {
            cur_sink = SFG_INVALID_NODE_ID;
            for (int i = 0; i < rg->num_nodes; i++) {
                if (rg->nodes[i].type == SFG_NODE_SINK) {
                    cur_sink = (sfg_node_id_t)i;
                    break;
                }
            }
        }
        rg->input_node  = cur_source;
        rg->output_node = cur_sink;

        if (cur_source == SFG_INVALID_NODE_ID
            || cur_sink == SFG_INVALID_NODE_ID) break;

        /* Check if we're done: only source and sink remaining,
         * with a single branch between them */
        if (rg->num_nodes == 2
            && sfg_find_branch(rg, cur_source, cur_sink) != NULL) {
            const sfg_branch_t* final_branch =
                sfg_find_branch(rg, cur_source, cur_sink);
            result->overall_gain = final_branch->gain;
            result->success = 1;

            if (result->num_steps < SFG_MAX_REDUCTION_STEPS) {
                sfg_reduction_step_t* st =
                    &result->steps[result->num_steps++];
                st->operation = SFG_REDUCE_COMPLETE;
                snprintf(st->description, sizeof(st->description),
                         "Reduction complete: T = ");
                sfg_cmplx_format(st->description + 22,
                                 sizeof(st->description) - 22,
                                 result->overall_gain);
            }
            return 0;
        }

        if (result->num_steps >= SFG_MAX_REDUCTION_STEPS) break;

        sfg_reduction_step_t step;
        int progress = sfg_reduce_interactive(rg, &step);
        if (progress == 1) {
            result->steps[result->num_steps++] = step;
        } else {
            /* Stuck — try absorbing any internal node */
            int absorbed = 0;
            for (int n = 0; n < rg->num_nodes && !absorbed; n++) {
                sfg_node_id_t nid = (sfg_node_id_t)n;
                if (nid != cur_source && nid != cur_sink) {
                    int ret2 = sfg_reduce_node_absorb(rg, nid);
                    if (ret2 == 0) {
                        absorbed = 1;
                        if (result->num_steps < SFG_MAX_REDUCTION_STEPS) {
                            sfg_reduction_step_t* st =
                                &result->steps[result->num_steps++];
                            st->operation = SFG_REDUCE_NODE_ABSORB;
                            snprintf(st->description, sizeof(st->description),
                                     "Forced absorption of node %d", n);
                        }
                    }
                }
            }
            if (!absorbed) break; /* Truly stuck */
        }
    }

    result->success = 0;
    return -1; /* Could not fully reduce */
}

double sfg_reduce_mason_compare(const sfg_graph_t* g,
                                 sfg_node_id_t source,
                                 sfg_node_id_t sink) {
    if (!g) return -1.0;

    /* Compute via Mason */
    sfg_gain_t T_mason = sfg_mason_compute_simple(g, source, sink);

    /* Compute via reduction */
    sfg_reduction_result_t reduce_result;
    sfg_reduction_result_init(&reduce_result);

    sfg_graph_t g_copy;
    sfg_graph_init(&g_copy);
    sfg_graph_copy(&g_copy, g);

    int ret = sfg_reduce_to_single(&g_copy, source, sink, &reduce_result);
    sfg_gain_t T_reduce = ret == 0 ? reduce_result.overall_gain : 0.0;

    double error = cabs(T_mason - T_reduce) / (cabs(T_mason) + 1e-15);
    return error;
}

/* ================================================================
 * L5: Reduction Display
 * ================================================================ */

void sfg_reduction_result_print(const sfg_reduction_result_t* result) {
    if (!result) return;

    printf("\n==========================================\n");
    printf("  SFG Reduction — Step-by-Step\n");
    printf("==========================================\n\n");

    for (int i = 0; i < result->num_steps; i++) {
        const char* op_names[] = {
            "SERIES", "PARALLEL", "SELF-LOOP", "NODE-ABSORB",
            "NODE-SPLIT", "SOURCE-SHIFT", "COMPLETE"
        };
        printf("  Step %3d: %-12s — %s\n",
               i + 1,
               op_names[result->steps[i].operation % 7],
               result->steps[i].description);
    }

    if (result->success) {
        printf("\n  ✓ Final Transfer Function: ");
        sfg_cmplx_print(result->overall_gain, NULL);
    } else {
        printf("\n  ✗ Reduction incomplete\n");
    }
    printf("==========================================\n\n");
}
