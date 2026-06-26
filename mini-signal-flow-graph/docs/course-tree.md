# Course Tree — mini-signal-flow-graph

## Prerequisite Knowledge

```
Linear Algebra (matrices, determinants)
    │
    ├── Complex Numbers (C field, magnitude/phase)
    │
    ├── Graph Theory (directed graphs, cycles, paths)
    │
    └── Control Theory Basics (transfer functions, feedback)
            │
            └── mini-signal-flow-graph
                    │
                    ├── Graph construction (nodes, branches, gains)
                    ├── Path/loop enumeration (DFS, Johnson's algorithm)
                    ├── Mason's Gain Formula
                    ├── SFG reduction rules (R1-R5)
                    ├── State-space ↔ SFG conversion
                    ├── Sensitivity & stability analysis
                    └── Applications (mechanical, electrical, aerospace)
                            │
                            ├── mini-block-diagram (next)
                            ├── mini-transfer-function
                            └── mini-linearization
```

## Dependency Graph

```
mini-linear-algebra-control ──→ mini-signal-flow-graph
                                        │
                    ┌───────────────────┼───────────────────┐
                    ↓                   ↓                   ↓
          mini-block-diagram    mini-transfer-function   mini-linearization
```

## This Module Provides

1. **Core SFG theory**: Mason's formula, loop enumeration, determinant computation
2. **Graph reduction**: Systematic simplification rules
3. **State-space bridge**: Bidirectional conversion
4. **Analysis tools**: Sensitivity, stability, Monte Carlo
5. **Application examples**: Mechanical, electrical, aerospace, MIMO
