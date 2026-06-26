# Course Tree — mini-linearization

## Prerequisite Knowledge

```
Linear Algebra (matrix operations, eigenvalues, SVD)
  └── Calculus (Taylor series, partial derivatives)
        └── ODE Theory (stability, phase portraits)
              └── Linear Control Theory (state-space, pole placement, LQR)
                    └── mini-linearization
                          ├── mini-transfer-function
                          ├── mini-state-space
                          ├── mini-block-diagram
                          └── mini-signal-flow-graph
```

## Internal Dependencies

```
linearization_core.h  ← (no internal deps)
  ├── taylor_expansion.h  ← linearization_core.h
  ├── feedback_linearization.h  ← linearization_core.h
  ├── koopman_operator.h  ← (standalone, complex.h)
  ├── describing_function.h  ← (standalone, complex.h)
  ├── jacobian_methods.h  ← linearization_core.h
  └── gain_scheduling.h  ← linearization_core.h
        └── carleman_linearization.c  ← koopman_operator.h, linearization_core.h
```

## Knowledge Progression

1. **L1-L3**: Definitions and math structures (struct types, API)
2. **L4**: Fundamental theorems (proof validation)
3. **L5**: Methods implementation (algorithms)
4. **L6**: Problem solving (examples)
5. **L7-L9**: Applications, advanced topics, frontiers
