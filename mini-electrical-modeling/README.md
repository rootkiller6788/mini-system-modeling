# mini-electrical-modeling

Electrical system modeling for control and automation — RLC circuits, transfer functions, state-space, two-port networks, electromechanical systems.

## Module Status: COMPLETE ✅

- **L1-L6: Complete** (all core definitions, concepts, math, laws, methods, problems)
- **L7: Partial** (3 applications: PWM drive, motor sizing, RF amplifier)
- **L8: Partial** (1/5 advanced: mixed-mode S-parameters)
- **L9: Partial** (documented: nanoelectronics, memristor, non-Foster, WPT)

**Line Count**: `include/` + `src/` >= 4500 lines ✅ (exceeds 3000 minimum)

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| **L1** | Definitions | ✅ Complete | Resistor/Cap/Inductor/Source/Phasor/TF/SS/Z/Y/H/ABCD/S/Motor params |
| **L2** | Core Concepts | ✅ Complete | Ohm/KVL/KCL/Impedance/AC power/Thevenin/MaxPower/Reciprocity/DC motor |
| **L3** | Math Structures | ✅ Complete | Complex Z/τ/ω0/ζ/Q/Matrix/QR eigenvalue/Polynomial TF |
| **L4** | Fundamental Laws | ✅ Complete | Ohm/Joule/Faraday/Kirchhoff/Tellegen/Hurwitz/Rollett |
| **L5** | Engineering Methods | ✅ Complete | MNA/Gaussian/QR/RK4/Pade/Durand-Kerner/Bode/Filters/ZOH |
| **L6** | Canonical Problems | ✅ Complete | RC filter/DC motor/RF cascade/RLC resonance (3 examples) |
| **L7** | Applications | ⚠️ Partial | PWM drive, motor sizing, 2.4GHz amplifier |
| **L8** | Advanced Topics | ⚠️ Partial | Mixed-mode S-parameters (differential circuits) |
| **L9** | Research Frontiers | ⚠️ Partial | Nanoelectronics, memristor, non-Foster, WPT (documented) |

## Core Definitions

- **Resistance R** [Ω] — `Resistor` struct, `ohms_law_voltage/current`
- **Capacitance C** [F] — `Capacitor` struct, `energy_capacitor`, `capacitor_current`
- **Inductance L** [H] — `Inductor` struct, `energy_inductor`, `inductor_voltage`
- **Impedance Z(s)** — `impedance_rlc_series/parallel`, complex s-domain
- **Phasor** — `Phasor` struct, sinusoidal steady-state representation
- **Transfer Function H(s)** — `TransferFunction` struct, rational polynomial
- **State-Space (A,B,C,D)** — `StateSpace` struct, canonical forms
- **Z-Parameters** — `ZParameters` struct, impedance matrix
- **Y-Parameters** — `YParameters` struct, admittance matrix
- **S-Parameters** — `SParameters` struct, scattering matrix with Z0 reference
- **ABCD Parameters** — `ABCDParameters` struct, transmission/cascade
- **DC Motor Params** — `DCMotorParams` struct, Ra/La/Kt/Kb/J/B

## Core Theorems (with formulas)

| Theorem | Formula | Implementation |
|---------|---------|---------------|
| Ohm's Law | V = I·R | `ohms_law_voltage`, `ohms_law_current` |
| Joule's Law | P = I²R = V²/R | `power_dc_i2r`, `power_dc_v2r` |
| Series RLC Impedance | Z(s) = R + sL + 1/(sC) | `impedance_rlc_series` |
| RC Time Constant | τ = RC | `time_constant_rc` |
| Resonant Frequency | ω₀ = 1/√(LC) | `resonant_frequency_rad_s` |
| Damping Factor | ζ = (R/2)·√(C/L) | `damping_factor_series_rlc` |
| Quality Factor | Q = ω₀L/R | `quality_factor_series_rlc` |
| KVL | Σv_i = 0 around loop | `kvl_verify` |
| KCL | Σi_i = 0 at node | `kcl_verify` |
| Thevenin | V_th = V_oc, R_th = V_th/I_sc | `thevenin_resistive_divider` |
| Max Power Transfer | P_max = V_th²/(4R_th) at R_L=R_th | `maximum_power_transfer` |
| Tellegen | Σv_k·i_k = 0 | `tellegen_power_sum` |
| Black's Feedback | H(s) = H₁/(1+H₁·H₂) | `tf_negative_feedback` |
| Matrix Exponential | Φ(t) = e^{At} (Pade [8/8]) | `ss_state_transition` |
| Rollett Stability | K = (1-|S₁₁|²-|S₂₂|²+|Δ|²)/(2|S₁₂S₂₁|) | `sp_rollett_stability` |
| DC Motor Speed | ω = (Kt·Va-Ra·TL)/(Kt·Kb+Ra·B) | `dc_motor_steady_state_speed` |

## Core Algorithms

| Algorithm | Complexity | Implementation |
|-----------|-----------|---------------|
| Nodal analysis (G matrix) | O(B + N²) | `nodal_analysis_build` |
| Modified Nodal Analysis (MNA) | O(B + (N+M)²) | `mna_build_matrix` |
| Gaussian elimination (pivoted) | O(n³) | `gaussian_elimination` |
| Durand-Kerner root-finding | O(n²) per iter | `tf_find_poles/zeros` |
| Binary search cutoff | O(log(range/ε)) | `tf_find_cutoff_freq` |
| Step response (partial fractions) | O(n) | `tf_step_response` |
| QR eigenvalue (Hessenberg+shift) | O(n³) | `ss_eigenvalues` |
| Pade matrix exponential [8/8] | O(n³) | `ss_state_transition` |
| RK4 state simulation | O(n²·steps) | `ss_simulate_rk4` |
| ZOH discretization (Van Loan) | O((n+m)³) | `ss_discretize_zoh` |
| Polynomial multiply (convolution) | O(n·m) | TF cascade/parallel/feedback |
| Filter design (E12/E6 lookup) | O(|E12|·|E6|) | `design_rc_lowpass_to_cutoff` |

## Classic Problems (with examples)

1. **RC Filter Design** — `examples/example_rc_filter.c`
   Design first-order RC low-pass for specific cutoff. Analyze frequency and time response.

2. **DC Motor Modeling** — `examples/example_dc_motor.c`
   Model 12V brushed DC motor: steady-state, speed-torque curve, startup transient.

3. **RF Amplifier Cascade** — `examples/example_two_port.c`
   Two-stage 2.4 GHz amplifier: S-parameter cascade, stability, gain analysis.

## Nine-School Course Mapping

| School | Key Course | Module Topics Covered |
|--------|-----------|----------------------|
| **MIT** | 6.302 Feedback Systems | Electrical networks, DC motor, two-port |
| **Stanford** | ENGR105 Feedback Control | s-domain elements, servomotor TF |
| **Berkeley** | ME 132/232 Dynamic Systems | RLC state-space, eigenvalue stability |
| **Caltech** | CDS 110/212 Robust | Controllability Gramians, matrix exp |
| **ETH** | 151-0591 Control I | Circuit TF, canonical forms, ZOH |
| **Cambridge** | 3F2/4F2 Systems | Two-port theory, frequency response |
| **Georgia Tech** | ECE 6550/ME 6401 | Electromechanical, PWM |
| **Purdue** | ECE 602/ME 675 | MNA, Gaussian elimination |
| **Tsinghua** | 自动控制原理 | TF analysis, SS fundamentals, canonical |

## Build and Test

```bash
# Build and run all tests
make test

# Build and run all examples
make examples
make run-examples

# Clean
make clean
```

## File Structure

```
mini-electrical-modeling/
├── Makefile
├── README.md
├── include/
│   ├── electrical_elements.h          (178 lines, L1-L4 definitions)
│   ├── circuit_topology.h             (273 lines, L4-L5 analysis)
│   ├── transfer_function_electrical.h (311 lines, L3-L5 TF)
│   ├── electrical_state_space.h       (331 lines, L2-L5 SS)
│   ├── two_port_network.h             (371 lines, L1-L8 two-port)
│   └── electromechanical_systems.h    (321 lines, L1-L7 motors)
├── src/
│   ├── electrical_elements.c          (509 lines, 30+ functions)
│   ├── circuit_topology.c             (470 lines, KVL/MNA/Gaussian)
│   ├── transfer_function_electrical.c (734 lines, filters/step/Bode)
│   ├── electrical_state_space.c       (843 lines, QR/RK4/Pade/ZOH)
│   ├── two_port_network.c             (300+ lines, Z/Y/S all conversions)
│   └── electromechanical_systems.c    (300+ lines, motor/generator/PWM)
├── tests/
│   └── test_electrical_modeling.c     (40+ assertion tests)
├── examples/
│   ├── example_rc_filter.c            (end-to-end filter design)
│   ├── example_dc_motor.c             (motor modeling and simulation)
│   └── example_two_port.c             (RF amplifier cascade)
└── docs/
    ├── knowledge-graph.md             (L1-L9 complete mapping)
    ├── coverage-report.md             (coverage assessment)
    ├── gap-report.md                  (missing items + priorities)
    ├── course-alignment.md            (9-school curriculum mapping)
    └── course-tree.md                 (prerequisite dependency graph)
```

## Reference Textbooks

- Dorf & Svoboda, *Introduction to Electric Circuits* (2018)
- Chua, Desoer, Kuh, *Linear and Nonlinear Circuits* (1987)
- Ogata, *Modern Control Engineering* (2010)
- Pozar, *Microwave Engineering* (2012)
- Krause et al., *Analysis of Electric Machinery* (2002)
- Golub & Van Loan, *Matrix Computations* (2013)
- Higham, *Functions of Matrices* (2008)
