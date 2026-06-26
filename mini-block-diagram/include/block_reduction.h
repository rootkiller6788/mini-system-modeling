/**
 * @file block_reduction.h
 * @brief Systematic Block Diagram Reduction �� rules and algorithms
 *
 * Knowledge Coverage:
 *   L1: Reduction rule IDs, ReduceStep (trace element), ReduceOptions
 *   L2: Individual rules �� series, parallel, feedback, block movement
 *   L3: Reduction trace for step-by-step audit trail
 *   L4: Iterative reduction algorithm with convergence guarantee
 *   L5: Transfer function extraction from reduced diagram
 *
 * Block diagram reduction simplifies a complex multi-loop diagram into
 * a single equivalent block by applying proven algebraic equivalences.
 * The process is systematic (guaranteed convergence for LTI systems)
 * but not unique in the order of rule application.
 *
 * Key rules: Series (multiply), Parallel (add), Feedback (G/(1+GH)),
 * Movement (past summers and take-off points).
 *
 * References:
 *   R.C. Dorf, "Modern Control Systems", 13th ed., Ch. 2.7
 *   N.S. Nise, "Control Systems Engineering", 8th ed., Ch. 5
 *   ETH 151-0555 / Cambridge 3F2 / Berkeley ME132
 */
#ifndef BLOCK_REDUCTION_H
#define BLOCK_REDUCTION_H
#include <stddef.h>
#include "blockdiagram.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef enum { REDUCE_RULE_SERIES=0, REDUCE_RULE_PARALLEL=1, REDUCE_RULE_FEEDBACK=2,
    REDUCE_RULE_MOVE_SUMMER=3, REDUCE_RULE_MOVE_PICKOFF=4, REDUCE_RULE_ELIM_LOOP=5,
    REDUCE_RULE_COMBINE_GAIN=6 } reduce_rule_t;

typedef struct { reduce_rule_t rule_applied; int node_ids[4],n_nodes; char desc[256]; } ReduceStep;
typedef struct { int n_steps,max_steps; ReduceStep *steps; } ReduceTrace;
typedef struct { int max_iterations,verbose,preserve_structure,eliminate_unit_gains;
    double zero_tolerance; } ReduceOptions;

int br_reduce(BlockDiagram *bd, const ReduceOptions *opts, ReduceTrace *trace);
TransferFunction* br_extract_transfer(const BlockDiagram *bd);
TransferFunction* br_reduce_and_extract(BlockDiagram *bd, const ReduceOptions *opts);
int br_apply_series_rule(BlockDiagram *bd, int b1, int b2);
int br_apply_parallel_rule(BlockDiagram *bd, int b1, int b2, int sid);
int br_apply_feedback_rule(BlockDiagram *bd, int fwd, int fb, int sid);
int br_move_block_past_summer_backward(BlockDiagram *bd, int bid, int sid);
int br_move_block_past_summer_forward(BlockDiagram *bd, int bid, int sid);
int br_move_block_past_takeoff_backward(BlockDiagram *bd, int bid, int tid);
int br_move_block_past_takeoff_forward(BlockDiagram *bd, int bid, int tid);
int br_eliminate_gain_loop(BlockDiagram *bd, int gid, int sid);
ReduceTrace* br_trace_create(int cap);
void br_trace_destroy(ReduceTrace *t);
void br_trace_add_step(ReduceTrace *t, reduce_rule_t r, const int *nids, int n, const char *d);
void br_trace_print(const ReduceTrace *t);
ReduceOptions br_options_default(void);
int br_is_fully_reduced(const BlockDiagram *bd);

#ifdef __cplusplus
}
#endif
#endif
