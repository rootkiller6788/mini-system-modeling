/**
 * @file bd_graph_utils.c
 * @brief Additional block diagram graph utilities
 *
 * L5: Adjacency matrix construction, connected component analysis,
 *     node distance computation, diagram statistics.
 */
#include "blockdiagram.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/** Build adjacency matrix of the block diagram (1 if edge exists).
 *  Returns matrix as flat array [n_nodes * n_nodes], caller must free. */
int* bd_adjacency_matrix(const BlockDiagram *bd, int *n) {
    if (!bd || !n) return NULL;
    *n = bd->n_nodes;
    int *adj = (int*)calloc((size_t)(bd->n_nodes * bd->n_nodes), sizeof(int));
    if (!adj) return NULL;
    int i;
    for (i = 0; i < bd->n_edges; i++) {
        int si = -1, di = -1, j;
        for (j = 0; j < bd->n_nodes; j++) {
            if (bd->nodes[j]->id == bd->edges[i]->src_node) si = j;
            if (bd->nodes[j]->id == bd->edges[i]->dst_node) di = j;
        }
        if (si >= 0 && di >= 0) adj[si * bd->n_nodes + di] = 1;
    }
    return adj;
}

/** Count the total number of feedback loops detected in the diagram. */
int bd_count_loops(const BlockDiagram *bd) {
    if (!bd) return -1;
    bd_path_t loops[128];
    return bd_find_loops(bd, 128, loops);
}

/** Count the total number of forward paths from input to output. */
int bd_count_forward_paths(const BlockDiagram *bd) {
    if (!bd) return -1;
    bd_path_t paths[64];
    return bd_find_forward_paths(bd, 64, paths);
}

/** Compute the total number of nodes reachable from the input node via BFS. */
int bd_reachable_node_count(const BlockDiagram *bd) {
    if (!bd || bd->input_node_id < 0) return -1;
    int max_nid = bd->next_node_id + 2;
    int *visited = (int*)calloc((size_t)max_nid, sizeof(int));
    if (!visited) return -1;
    int *queue = (int*)malloc((size_t)max_nid * sizeof(int));
    if (!queue) { free(visited); return -1; }
    int front = 0, back = 0, cnt = 0;
    queue[back++] = bd->input_node_id;
    visited[bd->input_node_id] = 1;
    while (front < back) {
        int cur = queue[front++]; cnt++;
        int i;
        for (i = 0; i < bd->n_edges; i++) {
            if (bd->edges[i]->src_node == cur) {
                int nx = bd->edges[i]->dst_node;
                if (nx < max_nid && !visited[nx]) {
                    visited[nx] = 1;
                    queue[back++] = nx;
                }
            }
        }
    }
    free(visited); free(queue);
    return cnt;
}

/** Compute the longest path length (in edges) from input to any node. */
int bd_longest_path_from_input(const BlockDiagram *bd) {
    if (!bd || bd->input_node_id < 0) return -1;
    int max_nid = bd->next_node_id + 2;
    int *dist = (int*)malloc((size_t)max_nid * sizeof(int));
    if (!dist) return -1;
    int i;
    for (i = 0; i < max_nid; i++) dist[i] = -1;
    int *queue = (int*)malloc((size_t)max_nid * sizeof(int));
    if (!queue) { free(dist); return -1; }
    int front = 0, back = 0;
    queue[back++] = bd->input_node_id;
    dist[bd->input_node_id] = 0;
    int max_dist = 0;
    while (front < back) {
        int cur = queue[front++];
        if (dist[cur] > max_dist) max_dist = dist[cur];
        for (i = 0; i < bd->n_edges; i++) {
            if (bd->edges[i]->src_node == cur) {
                int nx = bd->edges[i]->dst_node;
                if (nx < max_nid && dist[nx] < 0) {
                    dist[nx] = dist[cur] + 1;
                    queue[back++] = nx;
                }
            }
        }
    }
    free(dist); free(queue);
    return max_dist;
}
