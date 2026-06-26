# Course Dependency Tree — mini-fluid-thermal-modeling

## Prerequisites (External)

```
Calculus (ODE/PDE)
  ├── Vector Calculus (grad, div, curl)
  ├── Differential Equations (analytical solutions)
  └── Numerical Methods (FD, iterative solvers)

Physics — Mechanics & Thermodynamics
  ├── Newtonʼs Laws (force balance → conservation)
  ├── Energy Conservation (1st Law)
  └── Kinetic Theory (transport properties)

Linear Algebra
  ├── Matrix Operations (tridiagonal solvers)
  └── Iterative Methods (Gauss-Seidel, SOR)
```

## Internal Dependency Tree

```
L1: Definitions
  └── fluid_thermal_core.h structs/enums ← all other modules depend on this

L2: Core Concepts (dimensionless groups, classifications)
  └── fluid_thermal_core.c ← builds on L1 definitions

L3: Engineering Quantities
  └── Typical value functions ← depends on L1 types

L4: Conservation Laws
  ├── fluid_mechanics.c (continuity, Bernoulli, NS solutions) ← L1+L2
  ├── heat_transfer.c (Fourier, Newton, Stefan-Boltzmann) ← L1+L2
  └── fluid_thermal.lean (formal theorems) ← L1 types

L5: Engineering Methods
  ├── fluid_mechanics.c (pipe flow, BL, dimensional analysis) ← L4
  ├── heat_transfer.c (Nu correlations, fins, transient) ← L4
  └── numerical_fluid_thermal.c (FD, Thomas, GS) ← L4

L6: Engineering Problems
  └── thermal_fluid_systems.c (pipe heating, pumps, natural circ) ← L4+L5

L7: Applications
  └── thermal_fluid_systems.c (electronics, HVAC, automotive) ← L6

L8: Advanced Methods
  └── numerical_fluid_thermal.c (SIMPLE), thermal_fluid_systems.c (multiphase) ← L5+L6

L9: Research Frontiers
  └── docs/knowledge-graph.md (documentation only) ← L8
```

## Cross-Module Dependencies

```
fluid_thermal_core.h ← NO internal dependencies (base module)
fluid_mechanics.h    ← fluid_thermal_core.h
heat_transfer.h      ← fluid_thermal_core.h
thermal_fluid_systems.h ← fluid_thermal_core.h
numerical_fluid_thermal.h ← fluid_thermal_core.h

fluid_thermal_core.c     ← fluid_thermal_core.h
fluid_mechanics.c        ← fluid_thermal_core.h, fluid_mechanics.h
heat_transfer.c          ← fluid_thermal_core.h, heat_transfer.h
thermal_fluid_systems.c  ← fluid_thermal_core.h, thermal_fluid_systems.h
numerical_fluid_thermal.c ← fluid_thermal_core.h, numerical_fluid_thermal.h
```

## Build Order

1. `fluid_thermal_core.o` (base — no internal deps)
2. `fluid_mechanics.o` → depends on (1)
3. `heat_transfer.o` → depends on (1)
4. `thermal_fluid_systems.o` → depends on (1)
5. `numerical_fluid_thermal.o` → depends on (1)
6. Link all .o into test binary / examples
