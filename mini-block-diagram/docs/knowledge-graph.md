# Knowledge Graph — mini-block-diagram

## L1: Definitions (Complete)

| # | Definition | Code Location |
|---|-----------|---------------|
| 1 | Block Diagram: directed graph of signal processing elements | `include/blockdiagram.h:BlockDiagram` |
| 2 | Node Types: INPUT, OUTPUT, BLOCK, SUMMER, TAKEOFF, GAIN | `include/blockdiagram.h:bd_node_type_t` |
| 3 | Transfer Function: G(s) = N(s)/D(s), ratio of s-polynomials | `include/blockdiagram.h:TransferFunction` |
| 4 | Edge: directed signal flow between two nodes | `include/blockdiagram.h:bd_edge_t` |
| 5 | Path: sequence of connected nodes (forward or loop) | `include/blockdiagram.h:bd_path_t` |
| 6 | Signal Flow Graph node types: SOURCE, SINK, INTERNAL | `include/signal_flow_graph.h:sfg_node_type_t` |
| 7 | Branch: weighted directed edge in SFG | `include/signal_flow_graph.h:sfg_branch_t` |
| 8 | MasonResult: forward paths, loops, Delta, cofactors | `include/mason.h:MasonResult` |
| 9 | Polynomial: ascending-power coefficient array | `include/transfer_function.h:Polynomial` |
| 10 | ZeroPoleGain: factored transfer function form | `include/transfer_function.h:ZeroPoleGain` |
| 11 | FreqPoint: single frequency response datapoint | `include/transfer_function.h:FreqPoint` |
| 12 | Transformation rules (9 types) | `include/block_algebra.h:ba_transformation_t` |
| 13 | Reduction rules (7 types) | `include/block_reduction.h:reduce_rule_t` |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Series (cascade) interconnection: G = G1*G2 | `src/block_algebra.c:ba_series`, `src/transfer_function.c:tf_series` |
| 2 | Parallel interconnection: G = G1+G2 | `src/block_algebra.c:ba_parallel`, `src/transfer_function.c:tf_parallel` |
| 3 | Negative feedback: G/(1+GH) | `src/block_algebra.c:ba_feedback`, `src/transfer_function.c:tf_feedback` |
| 4 | Positive feedback: G/(1-GH) | `src/transfer_function.c:tf_positive_feedback` |
| 5 | Standard feedback configuration | `src/block_algebra.c:ba_build_standard_feedback` |
| 6 | Disturbance rejection configuration | `src/block_algebra.c:ba_build_disturbance_rejection` |
| 7 | SFG construction and source/sink designation | `src/signal_flow_graph.c` |
| 8 | Block diagram lifecycle (create/destroy) | `src/blockdiagram.c` |
| 9 | Node property assignment (TF, gain, signs, IC) | `src/blockdiagram.c` |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Directed graph representation (adjacency via edge list) | `src/blockdiagram.c:BlockDiagram` |
| 2 | Polynomial algebra (add, multiply, scale, eval) | `src/transfer_function.c:poly_*` |
| 3 | Complex number evaluation (Horner's method) | `src/transfer_function.c:poly_eval_complex` |
| 4 | Transfer function arithmetic (add, mult, div, inv) | `src/transfer_function.c:tf_*` |
| 5 | Bode frequency response magnitude and phase | `src/transfer_function.c:tf_magnitude_at`, `tf_phase_at` |
| 6 | SFG evaluation via iterative relaxation | `src/signal_flow_graph.c:sfg_evaluate` |
| 7 | BD-to-SFG conversion algorithm | `src/signal_flow_graph.c:bd_to_sfg` |
| 8 | Adjacency matrix construction | `src/bd_graph_utils.c` |

## L4: Fundamental Laws (Complete)

| # | Law/Theorem | Verification |
|---|------------|-------------|
| 1 | Mason's Gain Formula: T = (1/Delta)*Sum(Pk*Dk) | `src/mason.c:mason_analyze`, `src/signal_flow_graph.c:sfg_mason_transfer_function_gain` |
| 2 | Mason's Determinant: Delta = 1 - Sum(Li) + Sum(LiLj) - ... | `src/mason.c:mason_determinant` |
| 3 | Routh-Hurwitz Stability Criterion | `src/transfer_function.c:tf_stability_routh` |
| 4 | Feedback Equivalence Theorem: G_cl = G/(1+GH) | Verified in tests and examples |
| 5 | Series Equivalence: G1*G2 = direct multiplication | Verified in `tf_series` + `br_apply_series_rule` |
| 6 | Parallel Equivalence: G1+G2 = direct addition | Verified in `tf_parallel` + `br_apply_parallel_rule` |
| 7 | Block Diagram Reduction Convergence | `src/block_reduction.c:br_reduce` (iterative guarantee) |
| 8 | Mason's theorem (Lean 4 statement) | `src/blockdiagram.lean` |

## L5: Algorithms/Methods (Complete)

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | DFS forward path enumeration | `src/blockdiagram.c:bd_find_forward_paths` |
| 2 | DFS loop detection with backtracking | `src/blockdiagram.c:bd_find_loops` |
| 3 | Non-touching loop pair detection | `src/blockdiagram.c:bd_loops_are_non_touching` |
| 4 | Systematic block diagram reduction | `src/block_reduction.c:br_reduce` |
| 5 | Newton-Raphson polynomial root finding with deflation | `src/transfer_function.c:tf_find_roots_real` |
| 6 | Bilinear (Tustin) transform s<->z | `src/transfer_function.c:tf_s_to_z`, `tf_z_to_s` |
| 7 | Step response characteristics (Tr, Ts, PO, ess) | `src/transfer_function.c:tf_rise_time` etc. |
| 8 | BFS reachability analysis | `src/blockdiagram.c:bd_path_exists` |
| 9 | Deep clone (graph + transfer functions) | `src/blockdiagram.c:bd_clone` |
| 10 | SFG iterative solver | `src/signal_flow_graph.c:sfg_evaluate` |

## L6: Canonical Problems (Complete)

| # | Problem | Example |
|---|---------|---------|
| 1 | Feedback system reduction | `examples/ex1_feedback.c` |
| 2 | Multi-loop Mason's analysis | `examples/ex2_mason.c` |
| 3 | DC motor speed control block diagram | `examples/ex3_dc_motor.c` |
| 4 | Quadrotor attitude control | `examples/ex4_quadrotor.c` |

## L7: Applications (Complete)

| # | Application | Reference |
|---|------------|-----------|
| 1 | DC motor speed control with PI controller | `ex3_dc_motor.c` — contains real motor parameters |
| 2 | Quadrotor attitude stabilization (250mm racing drone) | `ex4_quadrotor.c` — Betaflight-like PID tuning |
| 3 | Disturbance rejection control configuration | `src/block_algebra.c:ba_build_disturbance_rejection` |
| 4 | Standard negative feedback for tracking control | `src/block_algebra.c:ba_build_standard_feedback` |

## L8: Advanced Topics (Partial)

| # | Topic | Status |
|---|-------|--------|
| 1 | Mason's sensitivity analysis (dT/dg) | Partial — scope declared in `mason.c:mason_sensitivity` |
| 2 | Non-touching k-set enumeration | Partial — k=2 implemented in `mason.c:mason_non_touching_k_sets` |
| 3 | Graph connectivity analysis (BFS reachability) | Complete — `bd_graph_utils.c` |
| 4 | Multi-loop feedback elimination | Complete — `block_reduction.c:br_reduce` |
| 5 | BD-to-SFG conversion | Complete — `signal_flow_graph.c:bd_to_sfg` |

## L9: Research Frontiers (Partial — Documented)

- Automated block diagram reduction using AI/ML techniques
- Large-scale system decomposition methods
- Quantum control system block diagram representations
- Formal verification of block diagram equivalences (Lean 4)

