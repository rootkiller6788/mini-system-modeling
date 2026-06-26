# Course Tree — mini-electrical-modeling

## Prerequisite Dependency Graph

```
mini-electrical-modeling
│
├── Prerequisites (external modules)
│   ├── mini-control-mathematics (complex numbers, linear algebra)
│   ├── mini-system-modeling/mini-block-diagram (system representation)
│   └── mini-system-modeling/mini-transfer-function (TF concepts)
│
├── Core Dependencies (within module)
│   │
│   ├── electrical_elements (L1-L2) [NO dependencies]
│   │   └── Defines: R, L, C, sources, impedance, Ohm's Law
│   │
│   ├── circuit_topology (L4-L5) [depends: electrical_elements]
│   │   └── Depends on: element types for nodal/MNA analysis
│   │
│   ├── transfer_function_electrical (L3-L5) [depends: electrical_elements]
│   │   └── Depends on: RLC element values for TF creation
│   │
│   ├── electrical_state_space (L2-L5) [depends: electrical_elements]
│   │   └── Depends on: RLC values for state-space construction
│   │
│   ├── two_port_network (L1-L8) [depends: electrical_elements]
│   │   └── Depends on: complex impedance for Z-parameter formulation
│   │
│   └── electromechanical_systems (L1-L7) [depends: electrical_elements]
│       └── Depends on: resistance, inductance for motor electrical model
│
└── Knowledge Dependencies (cross-cutting)
    │
    ├── Complex numbers → phasors, impedance → all AC analysis
    ├── Linear algebra → state-space matrices, eigenvalue → SS analysis
    ├── Polynomial theory → transfer functions, root-finding → TF analysis
    ├── Matrix 2x2 operations → two-port parameter conversions
    ├── Numerical methods → Gaussian elimination, QR, RK4, Pade
    └── Physics → Maxwell/Faraday/Ohm/Joule → foundational laws
```

## Learning Path (Linearized)

1. **Start**: `electrical_elements` — Ohm's Law, RLC, impedance, power
2. **Then**: `circuit_topology` — KVL/KCL, Thevenin, MNA, Gaussian elimination
3. **Parallel**: `transfer_function_electrical` — s-domain, filters, Bode, step response
4. **Parallel**: `electrical_state_space` — eigenvalues, controllability, simulation
5. **Then**: `two_port_network` — Z/Y/S/ABCD, cascades, stability
6. **Finally**: `electromechanical_systems` — motors, generators, PWM

## Research Frontiers (L9)

| Topic | Dependency Chain |
|-------|-----------------|
| Nanoelectronic modeling | SS → quantum transport (NEGF) → ballistic MOSFET |
| Memristor modeling | electrical_elements → nonlinear R → Chua memristor |
| Non-Foster matching | two_port → negative impedance → NIC circuits |
| Wireless power transfer | electromechanical → resonant coupling → Qi/A4WP |
