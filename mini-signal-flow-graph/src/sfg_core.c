/**
 * sfg_core.c — Signal Flow Graph Core Implementation
 *
 * Implements the fundamental graph data structure and basic operations
 * for signal flow graphs. The SFG is a directed graph where nodes
 * represent system variables and branches represent signal transmission
 * with associated gains (transmittances).
 *
 * Knowledge coverage:
 *   L1: Graph construction, node/branch types, basic queries
 *   L2: Graph manipulation (copy, transpose, merge), connectivity
 *   L3: Connection matrix, gain matrix, graph display
 */

#include "sfg_core.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* ================================================================
 * L1: Graph Initialization and Cleanup
 * ================================================================ */

sfg_graph_t* sfg_graph_init(sfg_graph_t* g) {
    if (!g) return NULL;
    memset(g, 0, sizeof(sfg_graph_t));
    g->input_node  = SFG_INVALID_NODE_ID;
    g->output_node = SFG_INVALID_NODE_ID;
    g->name[0] = '\0';
    return g;
}

void sfg_graph_clear(sfg_graph_t* g) {
    if (!g) return;
    memset(g, 0, sizeof(sfg_graph_t));
    g->input_node  = SFG_INVALID_NODE_ID;
    g->output_node = SFG_INVALID_NODE_ID;
}

/* ================================================================
 * L1: Node Management
 * ================================================================ */

sfg_node_id_t sfg_add_node(sfg_graph_t* g, sfg_node_type_t type,
                            const char* label) {
    if (!g) return SFG_INVALID_NODE_ID;
    if (g->num_nodes >= SFG_MAX_NODES) return SFG_INVALID_NODE_ID;

    sfg_node_id_t id = (sfg_node_id_t)g->num_nodes;
    g->nodes[id].id   = id;
    g->nodes[id].type = type;
    g->nodes[id].is_input_tapped = 0;

    if (label && label[0]) {
        /* Safe strncpy with null termination guarantee */
        strncpy(g->nodes[id].label, label, sizeof(g->nodes[id].label) - 1);
        g->nodes[id].label[sizeof(g->nodes[id].label) - 1] = '\0';
    } else {
        /* Auto-generate label based on type */
        const char* type_names[] = {
            "src", "snk", "int", "sum", "pko", "st", "dmy"
        };
        snprintf(g->nodes[id].label, sizeof(g->nodes[id].label),
                 "%s_%d", type_names[type % 7], (int)id);
    }

    g->num_nodes++;
    return id;
}

sfg_node_id_t sfg_add_source(sfg_graph_t* g, const char* label) {
    return sfg_add_node(g, SFG_NODE_SOURCE, label);
}

sfg_node_id_t sfg_add_sink(sfg_graph_t* g, const char* label) {
    return sfg_add_node(g, SFG_NODE_SINK, label);
}

sfg_node_id_t sfg_add_state(sfg_graph_t* g, const char* label) {
    return sfg_add_node(g, SFG_NODE_STATE, label);
}

sfg_node_id_t sfg_add_summing(sfg_graph_t* g, const char* label) {
    return sfg_add_node(g, SFG_NODE_SUMMING, label);
}

int sfg_remove_node(sfg_graph_t* g, sfg_node_id_t node_id) {
    if (!g) return -1;
    if (node_id < 0 || node_id >= g->num_nodes) return -1;

    /* Remove all branches incident to this node */
    int i = 0;
    while (i < g->num_branches) {
        if (g->branches[i].src == node_id || g->branches[i].dst == node_id) {
            /* Compact array by moving last element here */
            if (i < g->num_branches - 1) {
                g->branches[i] = g->branches[g->num_branches - 1];
            }
            g->num_branches--;
        } else {
            i++;
        }
    }

    /* Remove node by compacting node array */
    if (node_id < g->num_nodes - 1) {
        /* Shift subsequent nodes down, update their IDs */
        for (int j = node_id; j < g->num_nodes - 1; j++) {
            g->nodes[j] = g->nodes[j + 1];
            g->nodes[j].id = (sfg_node_id_t)j;
        }
        /* Update branch src/dst references for shifted nodes */
        for (int k = 0; k < g->num_branches; k++) {
            if (g->branches[k].src > node_id) g->branches[k].src--;
            if (g->branches[k].dst > node_id) g->branches[k].dst--;
        }
        /* Update input/output node references */
        if (g->input_node > node_id) g->input_node--;
        if (g->output_node > node_id) g->output_node--;
    }
    g->num_nodes--;

    /* Invalidate input/output if they referenced the removed node */
    if (g->input_node == node_id)  g->input_node  = SFG_INVALID_NODE_ID;
    if (g->output_node == node_id) g->output_node = SFG_INVALID_NODE_ID;

    return 0;
}

/* ================================================================
 * L1: Branch Management
 * ================================================================ */

int sfg_add_branch(sfg_graph_t* g, sfg_node_id_t src, sfg_node_id_t dst,
                   sfg_gain_t gain, const char* label) {
    if (!g) return -1;
    if (src < 0 || src >= g->num_nodes) return -1;
    if (dst < 0 || dst >= g->num_nodes) return -1;
    if (g->num_branches >= SFG_MAX_BRANCHES) return -1;

    int idx = g->num_branches;
    g->branches[idx].src  = src;
    g->branches[idx].dst  = dst;
    g->branches[idx].gain = gain;

    if (label && label[0]) {
        strncpy(g->branches[idx].label, label,
                sizeof(g->branches[idx].label) - 1);
        g->branches[idx].label[sizeof(g->branches[idx].label) - 1] = '\0';
    } else {
        g->branches[idx].label[0] = '\0';
    }

    g->num_branches++;
    return idx;
}

int sfg_add_branch_real(sfg_graph_t* g, sfg_node_id_t src, sfg_node_id_t dst,
                        sfg_real_t gain, const char* label) {
    return sfg_add_branch(g, src, dst, gain + 0.0 * I, label);
}

int sfg_remove_branch(sfg_graph_t* g, sfg_node_id_t src, sfg_node_id_t dst) {
    if (!g) return -1;
    for (int i = 0; i < g->num_branches; i++) {
        if (g->branches[i].src == src && g->branches[i].dst == dst) {
            /* Compact by moving last branch here */
            if (i < g->num_branches - 1) {
                g->branches[i] = g->branches[g->num_branches - 1];
            }
            g->num_branches--;
            return 0;
        }
    }
    return -1;  /* Not found */
}

const sfg_branch_t* sfg_find_branch(const sfg_graph_t* g,
                                     sfg_node_id_t src, sfg_node_id_t dst) {
    if (!g) return NULL;
    for (int i = 0; i < g->num_branches; i++) {
        if (g->branches[i].src == src && g->branches[i].dst == dst) {
            return &g->branches[i];
        }
    }
    return NULL;
}

/* ================================================================
 * L2: Graph Configuration
 * ================================================================ */

void sfg_set_input(sfg_graph_t* g, sfg_node_id_t node_id) {
    if (g && node_id >= 0 && node_id < g->num_nodes) {
        g->input_node = node_id;
    }
}

void sfg_set_output(sfg_graph_t* g, sfg_node_id_t node_id) {
    if (g && node_id >= 0 && node_id < g->num_nodes) {
        g->output_node = node_id;
    }
}

/* ================================================================
 * L2: Graph Query Operations
 * ================================================================ */

int sfg_node_count(const sfg_graph_t* g) {
    return g ? g->num_nodes : 0;
}

int sfg_branch_count(const sfg_graph_t* g) {
    return g ? g->num_branches : 0;
}

const sfg_node_t* sfg_get_node(const sfg_graph_t* g, sfg_node_id_t id) {
    if (!g || id < 0 || id >= g->num_nodes) return NULL;
    return &g->nodes[id];
}

const sfg_branch_t* sfg_get_branch(const sfg_graph_t* g, int idx) {
    if (!g || idx < 0 || idx >= g->num_branches) return NULL;
    return &g->branches[idx];
}

int sfg_node_in_degree(const sfg_graph_t* g, sfg_node_id_t node_id) {
    if (!g || node_id < 0 || node_id >= g->num_nodes) return 0;
    int count = 0;
    for (int i = 0; i < g->num_branches; i++) {
        if (g->branches[i].dst == node_id) count++;
    }
    return count;
}

int sfg_node_out_degree(const sfg_graph_t* g, sfg_node_id_t node_id) {
    if (!g || node_id < 0 || node_id >= g->num_nodes) return 0;
    int count = 0;
    for (int i = 0; i < g->num_branches; i++) {
        if (g->branches[i].src == node_id) count++;
    }
    return count;
}

int sfg_is_source(const sfg_graph_t* g, sfg_node_id_t node_id) {
    if (!g || node_id < 0 || node_id >= g->num_nodes) return 0;
    int in_deg  = sfg_node_in_degree(g, node_id);
    int out_deg = sfg_node_out_degree(g, node_id);
    /* A source has only outgoing branches (or is isolated as pure source) */
    return (in_deg == 0 && out_deg > 0);
}

int sfg_is_sink(const sfg_graph_t* g, sfg_node_id_t node_id) {
    if (!g || node_id < 0 || node_id >= g->num_nodes) return 0;
    int in_deg  = sfg_node_in_degree(g, node_id);
    int out_deg = sfg_node_out_degree(g, node_id);
    return (in_deg > 0 && out_deg == 0);
}

int sfg_node_exists(const sfg_graph_t* g, sfg_node_id_t id) {
    return (g && id >= 0 && id < g->num_nodes);
}

int sfg_branch_exists(const sfg_graph_t* g,
                      sfg_node_id_t src, sfg_node_id_t dst) {
    return (sfg_find_branch(g, src, dst) != NULL);
}

/* ================================================================
 * L2: Graph Manipulation
 * ================================================================ */

sfg_graph_t* sfg_graph_copy(sfg_graph_t* dst, const sfg_graph_t* src) {
    if (!dst || !src) return NULL;
    memcpy(dst, src, sizeof(sfg_graph_t));
    return dst;
}

sfg_graph_t* sfg_graph_transpose(sfg_graph_t* dst, const sfg_graph_t* src) {
    if (!dst || !src) return NULL;

    /* Copy node structure */
    memcpy(dst, src, sizeof(sfg_graph_t));

    /* Swap source and destination on all branches */
    for (int i = 0; i < dst->num_branches; i++) {
        sfg_node_id_t tmp = dst->branches[i].src;
        dst->branches[i].src = dst->branches[i].dst;
        dst->branches[i].dst = tmp;
    }

    /* Swap input and output designations */
    sfg_node_id_t tmp_io = dst->input_node;
    dst->input_node  = dst->output_node;
    dst->output_node = tmp_io;

    return dst;
}

void sfg_graph_reverse(sfg_graph_t* g) {
    if (!g) return;
    for (int i = 0; i < g->num_branches; i++) {
        sfg_node_id_t tmp = g->branches[i].src;
        g->branches[i].src = g->branches[i].dst;
        g->branches[i].dst = tmp;
    }
    sfg_node_id_t tmp_io = g->input_node;
    g->input_node  = g->output_node;
    g->output_node = tmp_io;
}

int sfg_merge_nodes(sfg_graph_t* g, sfg_node_id_t keep,
                    sfg_node_id_t absorb) {
    if (!g) return -1;
    if (keep < 0 || keep >= g->num_nodes) return -1;
    if (absorb < 0 || absorb >= g->num_nodes) return -1;
    if (keep == absorb) return 0;

    /* Redirect all branches FROM absorb to FROM keep */
    for (int i = 0; i < g->num_branches; i++) {
        if (g->branches[i].src == absorb) {
            g->branches[i].src = keep;
        }
    }

    /* Redirect all branches TO absorb to TO keep */
    for (int i = 0; i < g->num_branches; i++) {
        if (g->branches[i].dst == absorb) {
            g->branches[i].dst = keep;
        }
    }

    /* Remove the absorbed node */
    return sfg_remove_node(g, absorb);
}

/* ================================================================
 * L3: Adjacency / Connection Matrix
 * ================================================================ */

int sfg_connection_matrix(const sfg_graph_t* g, sfg_gain_t* C, int N) {
    if (!g || !C) return -1;
    if (N != g->num_nodes) return -1;

    /* Initialize all entries to 0 */
    for (int i = 0; i < N * N; i++) {
        C[i] = 0.0;
    }

    /* Fill connection matrix: C[src*N + dst] = gain */
    for (int b = 0; b < g->num_branches; b++) {
        int src = (int)g->branches[b].src;
        int dst = (int)g->branches[b].dst;
        C[src * N + dst] = g->branches[b].gain;
    }

    return 0;
}

int sfg_gain_matrix(const sfg_graph_t* g, sfg_gain_t* G, int N) {
    if (!g || !G) return -1;
    if (N != g->num_nodes) return -1;

    /* G[i][j] is gain from node j to node i (column j → row i) */
    /* This is the transpose of the connection matrix convention */
    for (int i = 0; i < N * N; i++) {
        G[i] = 0.0;
    }
    for (int b = 0; b < g->num_branches; b++) {
        int src = (int)g->branches[b].src;
        int dst = (int)g->branches[b].dst;
        G[dst * N + src] = g->branches[b].gain;
    }

    return 0;
}

void sfg_node_index_map(const sfg_graph_t* g, int* map) {
    if (!g || !map) return;
    for (int i = 0; i < g->num_nodes; i++) {
        map[i] = (int)g->nodes[i].id;
    }
}

/* ================================================================
 * L3: Graph Display and Debugging
 * ================================================================ */

void sfg_graph_print(const sfg_graph_t* g) {
    if (!g) {
        printf("SFG: (null)\n");
        return;
    }

    printf("=== Signal Flow Graph: %s ===\n",
           g->name[0] ? g->name : "(unnamed)");
    printf("Nodes: %d, Branches: %d\n", g->num_nodes, g->num_branches);
    printf("Input Node: %d, Output Node: %d\n\n",
           (int)g->input_node, (int)g->output_node);

    const char* type_names[] = {
        "SOURCE", "SINK", "INTERNAL", "SUMMING",
        "PICKOFF", "STATE", "DUMMY"
    };

    printf("--- Nodes ---\n");
    for (int i = 0; i < g->num_nodes; i++) {
        printf("  [%d] %-10s  type=%-10s  in=%d  out=%d\n",
               g->nodes[i].id,
               g->nodes[i].label,
               type_names[g->nodes[i].type % 7],
               sfg_node_in_degree(g, g->nodes[i].id),
               sfg_node_out_degree(g, g->nodes[i].id));
    }

    printf("\n--- Branches ---\n");
    for (int i = 0; i < g->num_branches; i++) {
        printf("  [%d] (%d) --> (%d)  gain = ",
               i, (int)g->branches[i].src, (int)g->branches[i].dst);
        sfg_cmplx_print(g->branches[i].gain, NULL);
    }
}

void sfg_graph_to_dot(const sfg_graph_t* g, const char* filename) {
    if (!g) return;

    FILE* fp = filename ? fopen(filename, "w") : stdout;
    if (!fp) return;

    fprintf(fp, "digraph SFG {\n");
    fprintf(fp, "  rankdir=LR;\n");
    fprintf(fp, "  node [shape=circle, style=filled, fillcolor=lightblue];\n");

    /* Emit nodes */
    for (int i = 0; i < g->num_nodes; i++) {
        const char* color = "lightblue";
        switch (g->nodes[i].type) {
            case SFG_NODE_SOURCE:  color = "lightgreen"; break;
            case SFG_NODE_SINK:    color = "lightcoral"; break;
            case SFG_NODE_SUMMING: color = "lightyellow"; break;
            case SFG_NODE_STATE:   color = "lightcyan"; break;
            default: break;
        }
        fprintf(fp, "  n%d [label=\"%s\\n[%d]\", fillcolor=%s];\n",
                g->nodes[i].id, g->nodes[i].label, g->nodes[i].id, color);
    }

    /* Emit branches */
    for (int i = 0; i < g->num_branches; i++) {
        double mag = cabs(g->branches[i].gain);
        double ang = carg(g->branches[i].gain) * 180.0 / 3.14159265358979323846;
        fprintf(fp, "  n%d -> n%d [label=\"%.3f∠%.1f°\"];\n",
                (int)g->branches[i].src, (int)g->branches[i].dst,
                mag, ang);
    }

    fprintf(fp, "}\n");

    if (filename) fclose(fp);
}

int sfg_graph_validate(const sfg_graph_t* g) {
    if (!g) return -1;

    /* Check node count */
    if (g->num_nodes < 0 || g->num_nodes > SFG_MAX_NODES) return -2;
    if (g->num_branches < 0 || g->num_branches > SFG_MAX_BRANCHES) return -3;

    /* Verify all branch endpoints reference valid nodes */
    for (int i = 0; i < g->num_branches; i++) {
        if (g->branches[i].src < 0 ||
            g->branches[i].src >= g->num_nodes) return -4;
        if (g->branches[i].dst < 0 ||
            g->branches[i].dst >= g->num_nodes) return -5;
    }

    /* Verify node IDs are sequential 0..N-1 */
    for (int i = 0; i < g->num_nodes; i++) {
        if (g->nodes[i].id != i) return -6;
    }

    return 0;
}
