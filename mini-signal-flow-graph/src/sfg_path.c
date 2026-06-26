/**
 * sfg_path.c — Signal Flow Graph Path and Loop Enumeration
 *
 * Implements DFS-based algorithms for finding forward paths, loops
 * (cycles), and non-touching loop groups in signal flow graphs.
 * These are the computational primitives upon which Mason's Gain
 * Formula operates.
 *
 * Knowledge coverage:
 *   L1: Path and loop data structures
 *   L5: DFS forward path enumeration (backtracking)
 *   L5: Johnson's algorithm for cycle enumeration
 *   L5: Non-touching loop group enumeration
 *   L4: Connection to Mason's determinant structure
 */

#include "sfg_path.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* ================================================================
 * L1: Data Structure Management
 * ================================================================ */

void sfg_path_list_clear(sfg_path_list_t* pl) {
    if (!pl) return;
    memset(pl, 0, sizeof(sfg_path_list_t));
}

void sfg_loop_list_clear(sfg_loop_list_t* ll) {
    if (!ll) return;
    memset(ll, 0, sizeof(sfg_loop_list_t));
}

void sfg_nt_groups_clear(sfg_nt_groups_t* nt) {
    if (!nt) return;
    memset(nt, 0, sizeof(sfg_nt_groups_t));
}

/* ================================================================
 * L5: Forward Path Enumeration (DFS with backtracking)
 * ================================================================ */

/**
 * dfs_forward_paths — Recursive DFS helper for forward path enumeration
 *
 * Algorithm:
 *   Starting from the source node, recursively explore all outgoing
 *   branches, avoiding nodes already on the current path. When the
 *   sink node is reached, record the path.
 *
 * The visited array is passed as a bit-set for efficiency.
 */
static void dfs_forward_paths(const sfg_graph_t* g,
                               sfg_node_id_t current,
                               sfg_node_id_t sink,
                               char* visited,
                               sfg_node_id_t* current_path,
                               int path_len,
                               sfg_gain_t accum_gain,
                               sfg_path_list_t* result) {
    /* Mark current node as visited */
    visited[(int)current] = 1;
    current_path[path_len] = current;
    path_len++;

    /* Check if we reached the sink */
    if (current == sink) {
        if (result->count < SFG_MAX_PATHS) {
            sfg_path_t* p = &result->paths[result->count];
            p->length = path_len;
            p->is_loop = 0;

            /* Copy node sequence */
            for (int i = 0; i < path_len; i++) {
                p->nodes[i] = current_path[i];
            }
            p->gain = accum_gain;
            result->count++;
        }
        /* Backtrack: unmark sink so that other paths can reach it */
        visited[(int)current] = 0;
        return;
    }

    /* Explore outgoing branches */
    for (int b = 0; b < g->num_branches; b++) {
        if (g->branches[b].src == current) {
            sfg_node_id_t next = g->branches[b].dst;
            if (!visited[(int)next]) {
                sfg_gain_t new_gain = accum_gain * g->branches[b].gain;
                dfs_forward_paths(g, next, sink, visited,
                                  current_path, path_len, new_gain, result);
            }
        }
    }

    /* Backtrack: unmark current node */
    visited[(int)current] = 0;
}

int sfg_find_forward_paths(const sfg_graph_t* g,
                           sfg_node_id_t source, sfg_node_id_t sink,
                           sfg_path_list_t* paths) {
    if (!g || !paths) return -1;
    if (!sfg_node_exists(g, source)) return -1;
    if (!sfg_node_exists(g, sink)) return -1;

    sfg_path_list_clear(paths);

    char* visited = (char*)calloc((size_t)g->num_nodes, sizeof(char));
    if (!visited) return -1;

    sfg_node_id_t* current_path =
        (sfg_node_id_t*)calloc((size_t)SFG_MAX_NODES, sizeof(sfg_node_id_t));
    if (!current_path) {
        free(visited);
        return -1;
    }

    dfs_forward_paths(g, source, sink, visited, current_path, 0,
                      1.0 + 0.0 * I, paths);

    free(visited);
    free(current_path);
    return paths->count;
}

/* ================================================================
 * L5: Loop Enumeration (Johnson's Algorithm adaptation)
 * ================================================================ */

/**
 * johnson_cycles_from — Find all elementary cycles starting from a
 * subset of nodes, following Johnson's algorithm structure.
 *
 * This is a simplified version of Johnson's algorithm for directed
 * graphs. For each node, we perform DFS to find all simple cycles
 * that have this node as their minimum-index vertex.
 */
static void find_cycles_from(const sfg_graph_t* g,
                              sfg_node_id_t start,
                              sfg_node_id_t current,
                              char* blocked,
                              char* visited,
                              sfg_node_id_t* stack,
                              int stack_len,
                              sfg_loop_list_t* result) {
    visited[(int)current] = 1;
    stack[stack_len] = current;
    stack_len++;

    /* Explore outgoing branches */
    for (int b = 0; b < g->num_branches; b++) {
        if (g->branches[b].src == current) {
            sfg_node_id_t next = g->branches[b].dst;

            if (next == start && stack_len >= 2) {
                /* Found a cycle back to start */
                if (result->count < SFG_MAX_LOOPS) {
                    sfg_path_t* loop = &result->loops[result->count];
                    loop->length = stack_len;
                    loop->is_loop = 1;

                    for (int i = 0; i < stack_len; i++) {
                        loop->nodes[i] = stack[i];
                    }

                    /* Compute loop gain */
                    loop->gain = 1.0 + 0.0 * I;
                    for (int i = 0; i < stack_len; i++) {
                        sfg_node_id_t from = stack[i];
                        sfg_node_id_t to   = stack[(i + 1) % stack_len];
                        /* For the last step back to start, we need
                           the branch from current to start */
                        const sfg_branch_t* branch;
                        if (i == stack_len - 1) {
                            branch = sfg_find_branch(g, current, start);
                        } else {
                            branch = sfg_find_branch(g, from, to);
                        }
                        if (branch) {
                            loop->gain *= branch->gain;
                        }
                    }
                    result->count++;
                }
            } else if (!visited[(int)next] && !blocked[(int)next]
                       && next >= start) {
                /* Continue DFS, but only to nodes not less than start
                   (to avoid duplicate cycle enumeration) */
                find_cycles_from(g, start, next, blocked, visited,
                                 stack, stack_len, result);
            }
        }
    }

    /* Backtrack */
    visited[(int)current] = 0;
    stack_len--;

    /* Block current node if no cycle through it was found */
    if (current != start) {
        /* Johnson's rule: block the node */
        int has_unblocked_path = 0;
        for (int b = 0; b < g->num_branches && !has_unblocked_path; b++) {
            if (g->branches[b].src == current
                && g->branches[b].dst >= start
                && !blocked[(int)g->branches[b].dst]) {
                has_unblocked_path = 1;
            }
        }
        if (!has_unblocked_path) {
            blocked[(int)current] = 1;
        }
    }
}

int sfg_find_all_loops(const sfg_graph_t* g, sfg_loop_list_t* loops) {
    if (!g || !loops) return -1;

    sfg_loop_list_clear(loops);

    int N = g->num_nodes;
    char* blocked = (char*)calloc((size_t)N, sizeof(char));
    char* visited = (char*)calloc((size_t)N, sizeof(char));
    sfg_node_id_t* stack =
        (sfg_node_id_t*)calloc((size_t)SFG_MAX_NODES, sizeof(sfg_node_id_t));

    if (!blocked || !visited || !stack) {
        free(blocked);
        free(visited);
        free(stack);
        return -1;
    }

    /* For each node as potential cycle start */
    for (sfg_node_id_t s = 0; s < (sfg_node_id_t)N; s++) {
        /* Reset blocked and visited for this start node */
        memset(blocked, 0, (size_t)N * sizeof(char));
        memset(visited, 0, (size_t)N * sizeof(char));

        find_cycles_from(g, s, s, blocked, visited, stack, 0, loops);

        if (loops->count >= SFG_MAX_LOOPS) break;
    }

    free(blocked);
    free(visited);
    free(stack);

    /* Also check for self-loops explicitly */
    for (int b = 0; b < g->num_branches; b++) {
        if (g->branches[b].src == g->branches[b].dst) {
            /* Self-loop — already found by Johnson? Check and add if not */
            int found = 0;
            for (int l = 0; l < loops->count; l++) {
                if (loops->loops[l].length == 1
                    && loops->loops[l].nodes[0] == g->branches[b].src) {
                    found = 1;
                    break;
                }
            }
            if (!found && loops->count < SFG_MAX_LOOPS) {
                sfg_path_t* loop = &loops->loops[loops->count];
                loop->length = 1;
                loop->is_loop = 1;
                loop->nodes[0] = g->branches[b].src;
                loop->gain = g->branches[b].gain;
                loops->count++;
            }
        }
    }

    return loops->count;
}

int sfg_find_self_loops(const sfg_graph_t* g, sfg_loop_list_t* loops) {
    if (!g || !loops) return -1;

    sfg_loop_list_clear(loops);

    for (int b = 0; b < g->num_branches; b++) {
        if (g->branches[b].src == g->branches[b].dst) {
            if (loops->count < SFG_MAX_LOOPS) {
                sfg_path_t* loop = &loops->loops[loops->count];
                loop->length = 1;
                loop->is_loop = 1;
                loop->nodes[0] = g->branches[b].src;
                loop->gain = g->branches[b].gain;
                loops->count++;
            }
        }
    }
    return loops->count;
}

/* ================================================================
 * L5: Non-Touching Analysis
 * ================================================================ */

int sfg_loops_touch(const sfg_path_t* loop1, const sfg_path_t* loop2) {
    if (!loop1 || !loop2) return 0;
    for (int i = 0; i < loop1->length; i++) {
        for (int j = 0; j < loop2->length; j++) {
            if (loop1->nodes[i] == loop2->nodes[j]) {
                return 1; /* Touch — share a node */
            }
        }
    }
    return 0; /* Non-touching */
}

int sfg_path_loops_touch(const sfg_path_t* path, const sfg_path_t* loop) {
    if (!path || !loop) return 0;
    for (int i = 0; i < path->length; i++) {
        for (int j = 0; j < loop->length; j++) {
            if (path->nodes[i] == loop->nodes[j]) {
                return 1;
            }
        }
    }
    return 0;
}

int sfg_paths_touch(const sfg_path_t* p1, const sfg_path_t* p2) {
    if (!p1 || !p2) return 0;
    for (int i = 0; i < p1->length; i++) {
        for (int j = 0; j < p2->length; j++) {
            if (p1->nodes[i] == p2->nodes[j]) {
                return 1;
            }
        }
    }
    return 0;
}

int sfg_find_non_touching_groups(const sfg_graph_t* g,
                                  const sfg_loop_list_t* all_loops,
                                  sfg_nt_groups_t* groups) {
    if (!g || !all_loops || !groups) return -1;

    sfg_nt_groups_clear(groups);
    int L = all_loops->count;

    if (L == 0) return 0;

    /* Precompute touch matrix: touch[i][j] = 1 if loop i touches loop j */
    char (*touch)[256] = NULL;
    if (L <= 256) {
        touch = (char(*)[256])calloc((size_t)L, 256 * sizeof(char));
    }
    /* Fall back to dynamic allocation for larger L */
    char* touch_flat = NULL;
    if (!touch) {
        touch_flat = (char*)calloc((size_t)L * (size_t)L, sizeof(char));
        if (!touch_flat) return -1;
    }

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            char val = (char)sfg_loops_touch(&all_loops->loops[i],
                                              &all_loops->loops[j]);
            if (touch) {
                touch[i][j] = val;
            } else {
                touch_flat[i * L + j] = val;
            }
        }
    }

    /* Generate non-touching groups by combinatorial exploration.
     * We use a recursive approach: starting from each loop,
     * add non-touching loops to build groups of increasing size. */

    int group_idx = 0;

    /* For each possible group size k = 2 up to L */
    for (int k = 2; k <= L && group_idx < SFG_MAX_NT_GROUPS; k++) {
        /* Use a recursive combination builder */
        int* current_group = (int*)calloc((size_t)k, sizeof(int));
        if (!current_group) break;

        /* Recursive closure to build combinations */
        #define MAX_LOOP_IDX 256
        /* Simplified: enumerate all k-combinations and test non-touching */
        int* indices = (int*)calloc((size_t)L, sizeof(int));
        if (!indices) { free(current_group); break; }

        /* Use stack-based k-combination enumeration */
        for (int i = 0; i < k; i++) indices[i] = i;

        while (1) {
            /* Test if all k loops are mutually non-touching */
            int all_nt = 1;
            for (int i = 0; i < k && all_nt; i++) {
                for (int j = i + 1; j < k && all_nt; j++) {
                    int ti = indices[i], tj = indices[j];
                    int touches;
                    if (touch) {
                        touches = touch[ti][tj];
                    } else {
                        touches = touch_flat[ti * L + tj];
                    }
                    if (touches) all_nt = 0;
                }
            }

            if (all_nt && group_idx < SFG_MAX_NT_GROUPS) {
                sfg_loop_group_t* grp = &groups->groups[group_idx];
                grp->num_loops = k;
                grp->group_gain = 1.0 + 0.0 * I;
                for (int i = 0; i < k; i++) {
                    grp->loop_indices[i] = indices[i];
                    grp->group_gain *= all_loops->loops[indices[i]].gain;
                }
                group_idx++;
            }

            /* Generate next combination */
            int pos = k - 1;
            while (pos >= 0 && indices[pos] == L - k + pos) pos--;
            if (pos < 0) break; /* Done all combinations */
            indices[pos]++;
            for (int j = pos + 1; j < k; j++) {
                indices[j] = indices[j - 1] + 1;
            }
        }

        free(indices);
        free(current_group);

        if (group_idx >= SFG_MAX_NT_GROUPS) break;
    }

    groups->count = group_idx;

    /* Fill sizes array */
    memset(groups->sizes, 0, sizeof(groups->sizes));
    for (int i = 0; i < group_idx; i++) {
        int sz = groups->groups[i].num_loops;
        if (sz <= SFG_MAX_LOOPS) {
            groups->sizes[sz]++;
        }
    }

    if (touch) free(touch);
    if (touch_flat) free(touch_flat);
    return group_idx;
}

/* ================================================================
 * L5: Path/Loop Utility Functions
 * ================================================================ */

sfg_gain_t sfg_path_gain_compute(const sfg_graph_t* g,
                                  const sfg_path_t* path) {
    if (!g || !path || path->length < 2) return 0.0;

    sfg_gain_t gain = 1.0 + 0.0 * I;
    for (int i = 0; i < path->length - 1; i++) {
        const sfg_branch_t* branch =
            sfg_find_branch(g, path->nodes[i], path->nodes[i + 1]);
        if (branch) {
            gain *= branch->gain;
        } else {
            return 0.0;
        }
    }
    return gain;
}

int sfg_path_contains_node(const sfg_path_t* path, sfg_node_id_t node) {
    if (!path) return 0;
    for (int i = 0; i < path->length; i++) {
        if (path->nodes[i] == node) return 1;
    }
    return 0;
}

void sfg_path_print(const sfg_path_t* path) {
    if (!path) { printf("(null path)\n"); return; }
    printf("%s: ", path->is_loop ? "Loop" : "Path");
    for (int i = 0; i < path->length; i++) {
        printf("%d", (int)path->nodes[i]);
        if (i < path->length - 1) printf(" → ");
    }
    if (path->is_loop) printf(" → %d", (int)path->nodes[0]);
    printf("  [gain = ");
    sfg_cmplx_print(path->gain, NULL);
    printf("]\n");
}

void sfg_loop_print(const sfg_path_t* loop) {
    sfg_path_print(loop);
}

/* ================================================================
 * L5: Sorting Utilities
 * ================================================================ */

static int cmp_gain_asc(const void* a, const void* b) {
    const sfg_path_t* pa = (const sfg_path_t*)a;
    const sfg_path_t* pb = (const sfg_path_t*)b;
    double ma = cabs(pa->gain);
    double mb = cabs(pb->gain);
    if (ma < mb) return -1;
    if (ma > mb) return 1;
    return 0;
}

void sfg_path_list_sort_by_gain(sfg_path_list_t* paths) {
    if (!paths || paths->count <= 1) return;
    qsort(paths->paths, (size_t)paths->count, sizeof(sfg_path_t), cmp_gain_asc);
}

void sfg_loop_list_sort_by_gain(sfg_loop_list_t* loops) {
    if (!loops || loops->count <= 1) return;
    qsort(loops->loops, (size_t)loops->count, sizeof(sfg_path_t), cmp_gain_asc);
}
