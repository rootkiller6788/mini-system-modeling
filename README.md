# Mini System Modeling

A collection of **from-scratch, zero-dependency C implementations** of university-level system modeling and control theory. Each module maps to MIT (and other top-tier university) courses, bridging theory and practice by translating textbook equations into runnable C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-block-diagram](mini-block-diagram/) | Block diagram data structures, algebra, systematic reduction, Mason's gain formula, signal flow graphs | MIT 6.302, Stanford ENGR105, ETH 151-0591 |
| [mini-electrical-modeling](mini-electrical-modeling/) | Circuit topology (KVL/KCL), RLC elements, state-space models, electromechanical systems, two-port networks | MIT 6.302, Stanford ENGR105, Berkeley EE 219 |
| [mini-fluid-thermal-modeling](mini-fluid-thermal-modeling/) | Fluid mechanics (Bernoulli, Navier-Stokes), heat transfer (conduction, convection, radiation), numerical methods, coupled thermal-fluid systems | MIT 2.006, MIT 2.51, Stanford ME 351B |
| [mini-linearization](mini-linearization/) | Taylor expansion, Jacobian linearization, describing functions, feedback linearization, gain scheduling, Koopman operator | MIT 6.241J, MIT 2.154, Stanford ENGR207B |
| [mini-mechanical-modeling](mini-mechanical-modeling/) | Mass-spring-damper elements, transfer functions, state-space, multi-DOF modal analysis, rotational systems, suspension and vibration isolation | MIT 2.003, MIT 2.004, Stanford ENGR105 |
| [mini-mechatronic-modeling](mini-mechatronic-modeling/) | Actuator selection and sizing, electromechanical coupling, sensor modeling, state-space for mechatronics, integrated system design | MIT 2.154, Stanford ENGR105, ETH 151-0591 |
| [mini-signal-flow-graph](mini-signal-flow-graph/) | SFG core structures, Mason's gain formula, path/loop enumeration, topological reduction, matrix methods, state-space conversion, sensitivity analysis | MIT 6.302, Stanford ENGR105, Cambridge 3F2 |
| [mini-transfer-function](mini-transfer-function/) | Polynomial algebra, s-domain rational functions, interconnection algebra, Routh-Hurwitz, Nyquist, Bode, SS↔TF conversion, model order reduction | MIT 6.302, Stanford ENGR105, Berkeley ME132, ETH 151-0591 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `tests/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes
- **Practical demos** — block diagram reduction, circuit simulation, fluid-thermal analysis, linearization workflows, mechanical system identification, mechatronic design, signal flow graph analysis, transfer function evaluation

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-block-diagram
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-system-modeling/
├── mini-block-diagram/           # Block diagram algebra, reduction, Mason's formula
├── mini-electrical-modeling/     # Circuit topology, RLC elements, two-port networks
├── mini-fluid-thermal-modeling/  # Fluid mechanics, heat transfer, coupled systems
├── mini-linearization/           # Taylor, Jacobian, describing functions, gain scheduling
├── mini-mechanical-modeling/     # Mass-spring-damper, MDOF, rotational systems
├── mini-mechatronic-modeling/    # Actuators, sensors, electromechanical coupling, integration
├── mini-signal-flow-graph/       # SFG core, Mason's formula, path enumeration, reduction
└── mini-transfer-function/       # Polynomial algebra, stability, conversion, order reduction
```

## License

MIT
