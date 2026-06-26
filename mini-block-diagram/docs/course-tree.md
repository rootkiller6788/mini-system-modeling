# Course Tree — mini-block-diagram

## Prerequisites (What this module depends on)

```
mini-transfer-function
├── Polynomial arithmetic
├── Laplace transform basics
└── Complex numbers

mini-signal-flow-graph
├── Directed graph theory
├── DFS/BFS algorithms
└── Mason's formula theory in signal_flow_graph.h
```

## Dependents (What depends on this module)

```
mini-block-diagram (this module)
│
├── 2. mini-time-domain-analysis
│   └── Uses block diagrams for system interconnection
│
├── 3. mini-root-locus-method
│   └── Root locus requires closed-loop TF from block diagram
│
├── 4. mini-frequency-domain
│   └── Bode/Nyquist analysis uses TFs from block diagrams
│
├── 5. mini-classical-compensator
│   └── Compensator design modifies block diagram structure
│
├── 6. mini-state-space-theory
│   └── State-space to TF conversion bridges to block diagrams
│
└── 7. mini-pole-placement-observer
    └── Observer-based control uses block diagram feedback
```

## Knowledge Dependency Graph

```
L1: Definitions
├── L2: Core Concepts (depends on L1 type definitions)
│   ├── L3: Math Structures (depends on L2 for operations)
│   │   ├── L4: Fundamental Laws (depends on L3 for formal statements)
│   │   │   ├── L5: Algorithms (depends on L4 for correctness guarantees)
│   │   │   │   ├── L6: Canonical Problems (depends on L5 for solution methods)
│   │   │   │   │   ├── L7: Applications (depends on L6 for validated models)
│   │   │   │   │   │   ├── L8: Advanced Topics
│   │   │   │   │   │   │   └── L9: Research Frontiers
```

## Module Internal Dependencies

```
blockdiagram.h (core types: TransferFunction, BlockDiagram)
├── transfer_function.h (TF operations, stability)
│   └── transfer_function.c
├── signal_flow_graph.h (SFG types and operations)
│   ├── signal_flow_graph.c (includes blockdiagram.h)
│   └── mason.h (Mason's formula, includes signal_flow_graph.h)
│       └── mason.c
├── block_algebra.h (interconnection operations)
│   └── block_algebra.c (includes block_reduction.h)
├── block_reduction.h (reduction rules and trace)
│   └── block_reduction.c
└── bd_graph_utils.c (auxiliary graph utilities)
```

