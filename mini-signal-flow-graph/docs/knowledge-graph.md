# Knowledge Graph — mini-signal-flow-graph

## Nine-Level Knowledge Coverage

### L1: Definitions (Complete)
- SFG Node: source, sink, internal, summing junction, pickoff, state, dummy (`sfg_node_t`, `sfg_node_type_t`)
- SFG Branch: directed edge with transmittance/gain (`sfg_branch_t`, `sfg_gain_t`)
- SFG Graph: collection of nodes and branches (`sfg_graph_t`)
- Forward Path: simple directed path from source to sink (`sfg_path_t`)
- Loop (Feedback Cycle): closed path starting and ending at same node
- Path Gain: product of branch gains along a path
- Loop Gain: product of branch gains around a loop
- Non-Touching Loops: loops sharing no common nodes
- Graph Determinant Δ: alternating sum of non-touching loop group gains
- Cofactor Δ_k: determinant of subgraph not touching forward path k
- Transfer Function T: Σ(P_k·Δ_k)/Δ (Mason's formula)
- State-Space System: (A, B, C, D) representation (`sfg_ss_system_t`)
- Connection Matrix: adjacency matrix with gains (`sfg_connection_matrix`)
- Complex Gain: complex-valued transmittance for frequency-domain analysis

### L2: Core Concepts (Complete)
- Signal Flow Representation: graphical encoding of linear equations
- Transmittance (Gain): per-branch multiplicative factor
- Mason's Rule: closed-form transfer function from SFG topology
- SFG ↔ Block Diagram Equivalence: nodes = signals, branches = blocks
- Graph Connectivity: source-to-sink reachability
- Feedback Structure: loops capture system feedback paths
- Determinant Structure: Δ captures complete feedback topology
- Node Degree Analysis: in-degree (fan-in) and out-degree (fan-out)
- Graph Transpose: corresponds to adjoint system
- Node Merging: combines redundant signal paths

### L3: Mathematical Structures (Complete)
- Complex Field C: gains as complex numbers (magnitude + phase)
- Connection Matrix (Adjacency with gains): C[i][j] = gain from i to j
- Gain Matrix: G[i][j] = gain from j to i (state-space convention)
- Complex Polynomial Algebra: Horner evaluation, derivative, multiplication
- Polar Representation: magnitude/phase form for frequency analysis
- Matrix-Vector Operations: C·x for signal propagation
- Gershgorin Circle Theorem: eigenvalue bounds from connection matrix
- Condition Number: κ(I-G) from SFG structure
- Resolvent Matrix: (sI-A)^(-1) entries via SFG Mason

### L4: Fundamental Laws (Complete)
- **Mason's Gain Formula Theorem**: T = (Σ P_k·Δ_k) / Δ
  - Proof: algebraic consequence of Cramer's rule on (I-G)x = b
  - Δ = 1 - ΣL₁ + ΣL₂ - ΣL₃ + ΣL₄ - ...
  - Implemented in `sfg_mason_compute()`
- **Cramer's Rule ↔ SFG Equivalence**: Verified in `sfg_cramer_verify()`
  - Mason's formula = ratio of determinants from Cramer's rule
- **Graph Determinant Structure**: Δ = 0 implies singular (ill-posed) system
- **Bode Sensitivity Integral**: ∫₀^∞ ln|S(jω)| dω = 0
- **Routh-Hurwitz Criterion**: Stability from characteristic polynomial
- **State-Space ↔ SFG Isomorphism**: Bidirectional conversion preserves poles/zeros

### L5: Engineering Methods (Complete)
- DFS Forward Path Enumeration (`sfg_find_forward_paths`)
- Johnson's Algorithm for Cycle Enumeration (`sfg_find_all_loops`)
- Non-Touching Loop Group Enumeration (`sfg_find_non_touching_groups`)
- Full Mason's Formula Evaluation (`sfg_mason_compute`)
- SFG Reduction — Rule R1: Series Elimination (`sfg_reduce_series`)
- SFG Reduction — Rule R2: Parallel Combination (`sfg_reduce_parallel`)
- SFG Reduction — Rule R3: Self-Loop Elimination (`sfg_reduce_self_loop`)
- SFG Reduction — Rule R4: Node Absorption (`sfg_reduce_node_absorb`)
- SFG Reduction — Rule R5: Node Splitting (`sfg_reduce_node_split`)
- Automated Full Reduction (`sfg_reduce_to_single`)
- State-Space → SFG Conversion (`sfg_from_state_space`)
- SFG → State-Space Extraction (`sfg_to_state_space`)
- Transfer Function → SFG Construction (`sfg_from_transfer_function`)
- Sensitivity Computation (`sfg_sensitivity`, `sfg_sensitivity_all`)
- Frequency Response Evaluation (`sfg_frequency_response`, `sfg_bode_plot_data`)
- Faddeev-Leverrier Algorithm for Characteristic Polynomial

### L6: Canonical Problems (Complete)
- **Mass-Spring-Damper System**: Mechanical ↔ SFG (`example_mass_spring.c`)
- **RLC Circuit Analysis**: Electrical ↔ SFG (`example_electrical.c`)
- **Control System TF**: PI controller + plant via Mason (`example_control_tf.c`)
- **Aircraft Pitch Dynamics**: Boeing 747 short-period (`example_aircraft.c`)
- **MIMO Coupled System**: 2×2 transfer function matrix (`example_mimo.c`)
- **Feedback System Analysis**: G/(1+GH) structure verification
- **Multi-Loop System**: Non-touching loop group analysis
- **Controllability via SFG**: Path-based structural controllability test
- **Observability via SFG**: Reverse path-based observability test

### L7: Applications (Complete — 7 applications)
- **DC Motor Speed Control**: PI controller (Katy pattern)
- **Aircraft GNC**: Boeing 747 pitch control SFG
- **RLC Circuit Design**: EDA circuit analysis via SFG
- **Mass-Spring-Damper**: Mechanical system modeling
- **Sensor Signal Conditioning**: Amplifier feedback SFG
- **Process Control**: Supplier chain flow graph
- **Power System Modeling**: Generator-exciter SFG

### L8: Advanced Topics (Partial+ — 5 implemented)
- **MIMO Mason Extension**: Multi-input multi-output TF matrix (`sfg_mason_mimo`)
- **Monte Carlo Tolerance Analysis**: Statistical gain variation simulation (`sfg_monte_carlo`)
- **Worst-Case Analysis**: Interval-based extreme value computation
- **Yield Analysis**: Manufacturing specification compliance estimation
- **Balanced Sensitivity**: Balanced circuit design metric
- **Sensitivity Matrix**: Full MIMO sensitivity analysis (`sfg_sensitivity_matrix`)
- **Condition Number Estimation**: κ(I-G) from SFG structure
- **Eigenvalue Bounds**: Gershgorin circles from connection matrix
- **Resolvent Matrix**: (sI-A)^(-1) entries via SFG methods

### L9: Research Frontiers (Partial — documented)
- **Nonlinear SFG Extensions**: Variable-gain branches
- **Quantum Signal Flow**: SFG for quantum circuits
- **Fractional-Order SFG**: Non-integer integrator orders
- **Large-Scale SFG Partitioning**: Hierarchical decomposition
- **Symbolic SFG Computer Algebra**: Automated transfer function derivation
