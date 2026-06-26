/**
 * @file signal_flow_graph.h
 * @brief Signal Flow Graph �� Mason's formula foundation
 *
 * Knowledge Coverage:
 *   L1: SFG node types (source/sink/internal), branch gains
 *   L2: Construction, source/sink designation
 *   L3: Reachability (BFS), forward path finding (DFS), loop detection
 *   L4: Mason's Delta determinant, cofactor Delta_k, transfer function
 *   L5: SFG evaluation (iterative relaxation), BD-to-SFG conversion
 *
 * Signal Flow Graphs (S.J. Mason, MIT, 1953) are weighted directed graphs
 * where nodes represent signals and branches represent gains. Mason's Gain
 * Formula directly yields T(s) without algebraic elimination.
 *
 * References:
 *   S.J. Mason, "Feedback Theory �� Some Properties of Signal Flow Graphs",
 *   Proc. IRE, Vol. 41, pp. 1144-1156, 1953
 *   S.J. Mason, "Feedback Theory �� Further Properties...", Proc. IRE, 1956
 *   MIT 6.302 / Stanford ENGR105 / Cambridge 3F2
 */
#ifndef SIGNAL_FLOW_GRAPH_H
#define SIGNAL_FLOW_GRAPH_H
#include <stddef.h>
#include "blockdiagram.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef enum { SFG_SOURCE=0, SFG_SINK=1, SFG_INTERNAL=2 } sfg_node_type_t;

typedef struct { int id; sfg_node_type_t type; char label[64]; double value; } sfg_node_t;
typedef struct { int id,from_node,to_node,is_unit_gain; double gain; } sfg_branch_t;
typedef struct { int *node_ids,length; double gain; } sfg_path_t;
typedef struct { int *node_ids,*branch_ids,length; double gain; } sfg_loop_t;

typedef struct { int n_nodes,max_nodes; sfg_node_t *nodes; int n_branches,max_branches;
    sfg_branch_t *branches; int source_id,sink_id,next_id; char name[128]; } SignalFlowGraph;

SignalFlowGraph* sfg_create(const char *name, int in, int ib);
void sfg_destroy(SignalFlowGraph *sg);
int sfg_add_node(SignalFlowGraph *sg, sfg_node_type_t t, const char *label);
int sfg_add_branch(SignalFlowGraph *sg, int from, int to, double gain);
int sfg_remove_branch(SignalFlowGraph *sg, int bid);
int sfg_set_source(SignalFlowGraph *sg, int nid);
int sfg_set_sink(SignalFlowGraph *sg, int nid);
int sfg_node_reachable(const SignalFlowGraph *sg, int nid);
int sfg_find_forward_paths(const SignalFlowGraph *sg, int maxp, sfg_path_t *out);
int sfg_find_loops(const SignalFlowGraph *sg, int maxl, sfg_loop_t *out);
int sfg_loops_are_non_touching(const sfg_loop_t *l1, const sfg_loop_t *l2);
int sfg_find_non_touching_sets(const sfg_loop_t *loops, int nl, int ks, int maxs, int *out);
double sfg_loop_set_gain_product(const sfg_loop_t *loops, const int *idx, int n);
double sfg_mason_determinant(const SignalFlowGraph *sg);
double sfg_mason_cofactor(const SignalFlowGraph *sg, const sfg_path_t *p);
double sfg_mason_transfer_function_gain(const SignalFlowGraph *sg);
int sfg_evaluate(const SignalFlowGraph *sg, double in, double *out);
void sfg_print(const SignalFlowGraph *sg);
SignalFlowGraph* bd_to_sfg(const void *bd);
void sfg_path_print(const sfg_path_t *p);
void sfg_loop_print(const sfg_loop_t *l);

#ifdef __cplusplus
}
#endif
#endif
