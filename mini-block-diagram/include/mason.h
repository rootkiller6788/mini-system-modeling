/**
 * @file mason.h
 * @brief Mason's Gain Formula �� T(s) = (1/Delta)*Sum(P_k*Delta_k)
 *
 * Knowledge Coverage:
 *   L1: MasonResult (forward paths, loops, Delta, cofactors)
 *   L2: Determinant computation, cofactor for each forward path
 *   L3: Non-touching loop set enumeration
 *   L4: Mason's theorem verification against block reduction
 *   L5: Sensitivity analysis dT/dg, report printing
 *
 * Mason's Gain Formula (1953) provides the exact closed-form transfer
 * function of any linear system from its signal flow graph without
 * algebraic manipulation. It is one of the most elegant results in
 * control theory.
 *
 * Formula: T = (1/Delta) * Sum_k [ P_k * Delta_k ]
 * where Delta = 1 - Sum(L_i) + Sum(L_i*L_j) - ...
 *
 * References:
 *   S.J. Mason, Proc. IRE, Vol. 41, 1953 / Vol. 44, 1956
 *   MIT 6.302 / Berkeley ME232 / Cambridge 3F2
 */
#ifndef MASON_H
#define MASON_H
#include <stddef.h>
#include "signal_flow_graph.h"
#include "transfer_function.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct { int n_forward_paths,n_loops,n_non_touching_pairs;
    double determinant,*path_gains,*cofactors,transfer_gain; } MasonResult;

double mason_determinant(const SignalFlowGraph *sg, int maxl);
double mason_cofactor(const SignalFlowGraph *sg, const int *pnids, int plen, int maxl);
MasonResult* mason_analyze(const SignalFlowGraph *sg, int maxp, int maxl);
void mason_result_free(MasonResult *r);
int mason_verify_against_reduction(const SignalFlowGraph *sg,
                                    const TransferFunction *tf, double tol);
TransferFunction* mason_transfer_function(const SignalFlowGraph *sg, int maxp, int maxl);
void mason_report_print(const MasonResult *r);
int* mason_non_touching_k_sets(const SignalFlowGraph *sg, int k, int maxs, int *oc);
double mason_sensitivity(const SignalFlowGraph *sg, int bid);

#ifdef __cplusplus
}
#endif
#endif
