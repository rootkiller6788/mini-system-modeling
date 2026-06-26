/**
 * @file block_algebra.h
 * @brief Block Diagram Algebra �� interconnections and reduction
 *
 * Knowledge Coverage:
 *   L1: Interconnection types (series, parallel, feedback)
 *   L2: Standard configurations (canonical feedback, disturbance rejection)
 *   L3: Transformation rules (move block past summer/takeoff)
 *   L4: Systematic reduction �� iterative application of equivalence rules
 *   L5: Equivalence checking, simplification, structure printing
 *
 * Block diagram algebra governs how subsystems combine. Each rule is a
 * proven equivalence that preserves the closed-loop I/O relationship.
 *
 * References:
 *   K. Ogata, "Modern Control Engineering", 5th ed., Ch. 3.3
 *   R.C. Dorf & R.H. Bishop, "Modern Control Systems", 13th ed., Ch. 2.6
 *   MIT 6.302 / ETH 151-0591 / Cambridge 3F2
 */
#ifndef BLOCK_ALGEBRA_H
#define BLOCK_ALGEBRA_H
#include <stddef.h>
#include "blockdiagram.h"
#ifdef __cplusplus
extern "C" {
#endif

BlockDiagram* ba_series(const BlockDiagram *a, const BlockDiagram *b);
BlockDiagram* ba_parallel(const BlockDiagram *a, const BlockDiagram *b);
BlockDiagram* ba_unity_feedback(const BlockDiagram *G);
BlockDiagram* ba_feedback(const BlockDiagram *G, const BlockDiagram *H);
BlockDiagram* ba_positive_feedback(const BlockDiagram *G, const BlockDiagram *H);

typedef enum { BA_MOVE_BLOCK_PAST_SUMMER_BACKWARD=0, BA_MOVE_BLOCK_PAST_SUMMER_FORWARD=1,
    BA_MOVE_BLOCK_PAST_TAKEOFF_BACKWARD=2, BA_MOVE_BLOCK_PAST_TAKEOFF_FORWARD=3,
    BA_MOVE_SUMMER_PAST_SUMMER=4, BA_MOVE_TAKEOFF_PAST_TAKEOFF=5,
    BA_ELIMINATE_FEEDBACK_LOOP=6, BA_ABSORB_PARALLEL_BLOCKS=7,
    BA_ABSORB_SERIES_BLOCKS=8 } ba_transformation_t;

typedef enum { BA_REDUCE_GREEDY=0, BA_REDUCE_FEEDBACK=1, BA_REDUCE_CASCADE=2 } ba_reduce_strategy_t;

int ba_apply_transformation(BlockDiagram *bd, ba_transformation_t rule,
                            const int *nids, int nn);
BlockDiagram* ba_reduce(const BlockDiagram *bd, ba_reduce_strategy_t s);
TransferFunction* ba_equivalent_transfer(const BlockDiagram *bd, ba_reduce_strategy_t s);
int ba_are_equivalent(const BlockDiagram *a, const BlockDiagram *b, double tol);
int ba_simplify(BlockDiagram *bd);
BlockDiagram* ba_build_standard_feedback(const TransferFunction *G, const TransferFunction *H);
BlockDiagram* ba_build_disturbance_rejection(const TransferFunction *C,
                                              const TransferFunction *P, const TransferFunction *H);
void ba_print_structure(const BlockDiagram *bd);

#ifdef __cplusplus
}
#endif
#endif
