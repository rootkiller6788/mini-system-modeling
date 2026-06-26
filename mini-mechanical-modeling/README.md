# mini-mechanical-modeling

Mechanical system modeling for control and automation — mass-spring-damper, transfer functions, state-space, rotational systems, multi-DOF modal analysis, and real-world applications.

## Module Status: COMPLETE

- **L1-L6: Complete** (all core definitions, concepts, math, laws, methods, problems)
- **L7: Complete** (3+ applications: automotive suspension Detroit/Tesla, vibration isolation, robot joints ISO 9283, UAV, spacecraft)
- **L8: Partial** (3/6 advanced: Guyan reduction, MAC criterion, Rayleigh damping, skyhook, Newmark-beta)
- **L9: Partial** (documented: nonlinear friction, active suspension, nanoscale)

**Line Count**: `include/` + `src/` = **4644 lines** (exceeds 3000 minimum)

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| **L1** | Definitions | Complete | Mass/Spring/Damper/Friction structs, wn/zeta/wd/Q/TF/SS/MDOF/Rotor/Gear/QuarterCar/RobotJoint |
| **L2** | Core Concepts | Complete | Newton F=ma, Hooke F=kx, viscous damping, Coulomb/Stribeck, energy, impedance, resonance, transmissibility |
| **L3** | Math Structures | Complete | 2nd-order ODE, complex Z_m, TF polynomials, SS matrices, eigenvalue K*phi=lambda*M*phi, FRF |
| **L4** | Fundamental Laws | Complete | Newton's 2nd law, Hooke's law, viscous law, energy conservation, D'Alembert, orthogonality |
| **L5** | Engineering Methods | Complete | Equivalent combinations, transmissibility, RK4, ZOH/FOH, Pade, modal extraction, Rayleigh, Guyan, Newmark |
| **L6** | Canonical Problems | Complete | SDOF vibration, 2-DOF modal, quarter-car (3 complete examples >=30 lines) |
| **L7** | Applications | Complete | Automotive suspension (Detroit, Tesla), vibration isolation, robot joints (ISO 9283), UAV quadrotor, spacecraft |
| **L8** | Advanced Topics | Partial | Guyan reduction, MAC, Rayleigh damping, skyhook, Newmark-beta, two-stage isolation |
| **L9** | Research Frontiers | Partial | Nonlinear friction, active/semi-active suspension, nanomechanics (documented) |

## Core Definitions

- **Mass m [kg]** — `MassEl` struct, `newton_second_law_force/acceleration`
- **Spring k [N/m]** — `SpringEl` struct, `hookes_law_force`, `spring_potential_energy`
- **Damper b [N*s/m]** — `DamperEl` struct, `viscous_damping_force`
- **Friction Model** — `FrictionModel` struct, `coulomb_friction`, `stribeck_friction`
- **Natural Frequency omega_n** — `natural_frequency_undamped/hz`
- **Damping Ratio zeta** — `damping_ratio`, `critical_damping_coefficient`
- **Logarithmic Decrement delta** — `logarithmic_decrement`, `damping_ratio_from_logdec`
- **Mechanical Impedance Z(jw)** — `mechanical_impedance_sdof`, `mechanical_mobility`
- **Transfer Function H(s)** — `MechanicalTF` struct, rational polynomial
- **State-Space (A,B,C,D)** — `MechSS` struct, `MechDSS` discrete
- **Moment of Inertia I** — 8 rigid body formulas, parallel axis theorem
- **Torsional Spring k_t** — `TorsionSpring`, `torsional_spring_torque`
- **Gear Ratio N** — `GearPair`, `GearTrain`, reflected parameters
- **MDOF System** — `MDOFSystem` with M/K/B/F matrices
- **Quarter-Car Model** — `QuarterCarModel`, body bounce, wheel hop, ride comfort
- **Vibration Isolator** — `IsolatorReq/Design`, transmissibility-based design
- **Robot Joint** — `RobotJoint`, effective inertia, PD gains, gravity compensation

## Core Theorems (with formulas)

| Theorem | Formula | Implementation |
|---------|---------|---------------|
| Newton's 2nd (translational) | F = m*a | `newton_second_law_force/acceleration` |
| Newton's 2nd (rotational) | T = J*alpha | `rotational_torque/accel` |
| Hooke's Law | F = -k*x | `hookes_law_force/displacement` |
| Torsional Hooke | T = -k_t*theta | `torsional_spring_torque` |
| Viscous Damping | F_b = -b*v | `viscous_damping_force` |
| Coulomb Friction | F_c = mu_k*N*sgn(v) | `coulomb_friction_force` |
| Spring Energy | U = 1/2*k*x^2 | `spring_potential_energy` |
| Kinetic Energy | K = 1/2*m*v^2 | `kinetic_energy_translational` |
| Natural Frequency | omega_n = sqrt(k/m) | `natural_frequency_undamped` |
| Damping Ratio | zeta = b/(2*sqrt(k*m)) | `damping_ratio` |
| Log Decrement | delta = (1/n)*ln(x0/xn) | `logarithmic_decrement` |
| Force Transmissibility | TR = sqrt((1+(2*zeta*r)^2)/((1-r^2)^2+(2*zeta*r)^2)) | `force_transmissibility` |
| Parallel Axis Theorem | I = I_cm + m*d^2 | `parallel_axis_theorem` |
| Percent Overshoot | %OS = 100*exp(-pi*zeta/sqrt(1-zeta^2)) | `percent_overshoot` |
| Half-Power Bandwidth | Delta_omega = 2*zeta*omega_n | `half_power_bandwidth` |
| SDOF TF | H(s) = 1/(m*s^2+b*s+k) | `tf_sdof_displacement` |
| SDOF SS | A=[[0,1],[-k/m,-b/m]], B=[[0],[1/m]] | `ss_sdof` |
| 2-DOF Eigenvalue | det(K - lambda*M) = 0 | `mdof_eigen_solve` |
| Rayleigh Damping | zeta_i = alpha/(2*omega_i)+beta*omega_i/2 | `rayleigh_damping_coeffs` |
| Gear Reflection | J_ref = J_load/N^2 | `reflected_inertia` |
| Isolator Design | k = m*(w_exc/r)^2, r = sqrt(1+1/TR) | `design_vibration_isolator` |
| PD Gains | Kp = J*wn^2, Kd = 2*zeta*wn*J | `robot_joint_pd_gains` |

## Core Algorithms

| Algorithm | Complexity | Implementation |
|-----------|-----------|---------------|
| Equivalent spring/damper (series/parallel) | O(n) | `spring_series_equivalent`, etc. |
| SDOF step response (analytic) | O(1) | `sdof_step_response_displacement` |
| Polynomial evaluation (Horner) | O(n) | `poly_eval` |
| TF cascade (polynomial convolution) | O(n1*n2) | `tf_cascade` |
| TF feedback | O(n1*n2) | `tf_neg_feedback`, `tf_pos_feedback` |
| RK4 state simulation | O(n^2*steps) | `ss_simulate_rk4` |
| Matrix determinant (LU+pivot) | O(n^3) | `mat_determinant` |
| Matrix inverse (Gauss-Jordan) | O(n^3) | `mat_inverse` |
| Linear system solve (Gaussian+pivot) | O(n^3) | `mat_solve` |
| Controllability/observability rank | O(n^3) | `ss_ctrb_rank`, `ss_obsv_rank` |
| QR eigenvalue iteration | O(n^3) | `ss_eigenvalues` |
| State transition matrix (Pade) | O(n^3) | `ss_state_transition` |
| ZOH discretization | O(n^3) | `ss_discretize_zoh` |
| Gramian computation (iterative) | O(k*n^3) | `ss_gramian_ctrb/obsv` |
| 2x2 generalized eigenvalue | O(1) | `mdof_eigen_solve` |
| Modal parameter extraction | O(n^2) | `mdof_modal_params` |
| FRF via modal superposition | O(n_modes) | `mdof_receptance_frf` |
| Rayleigh damping solver | O(1) | `rayleigh_damping_coeffs` |
| MAC criterion | O(n) | `mdof_mac` |
| Newmark-beta step | O(n^2) | `mdof_newmark_step` |

## Classic Problems (with examples)

1. **SDOF Vibration Analysis** — `examples/example_sdof_vibration.c`
   Free vibration parameters, step response, frequency response, transmissibility.

2. **2-DOF Modal Analysis** — `examples/example_2dof_modal.c`
   Natural frequencies, mode shapes, orthogonality, MAC, FRF, Rayleigh damping.

3. **Quarter-Car Suspension Design** — `examples/example_quarter_car.c`
   Ride comfort vs handling tradeoff, spring/damper selection, Detroit/Tesla application.

## Nine-School Course Mapping

| School | Key Course | Module Topics Covered |
|--------|-----------|----------------------|
| **MIT** | 2.003 Dynamics and Control I | SDOF/MDOF, TF, SS modeling |
| **MIT** | 2.004 Dynamics and Control II | Modal analysis, FRF, isolation |
| **Stanford** | ENGR105 Feedback Control | Mechanical ODEs, TF, SS |
| **Berkeley** | ME 132/232 Dynamic Systems | Newton-Euler, SDOF, eigenvalues |
| **Caltech** | CDS 110/212 Robust | Controllability, observability |
| **ETH** | 151-0591 Control I | Mechanical TF, SS canonical forms |
| **Cambridge** | 3F2/4F2 Systems & Control | Vibration, frequency response |
| **Georgia Tech** | AE 6530/ME 6401 | Rotor dynamics, geared systems |
| **Purdue** | ME 575/675 Multivariable | MDOF, modal superposition |
| **Tsinghua** | 自动控制原理 | Mechanical SS, Newton-Euler |

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
mini-mechanical-modeling/
  Makefile
  README.md
  include/
    mechanical_elements.h              (L1-L5: core physical elements)
    mechanical_transfer_function.h     (L3-L5: mechanical TFs)
    mechanical_state_space.h           (L2-L5: state-space models)
    rotational_systems.h              (L1-L6: rotation, gears, unbalance)
    multi_dof_systems.h               (L5-L8: MDOF, modal, Rayleigh, Guyan)
    mechanical_applications.h         (L6-L7: suspension, isolation, robotics, aero)
  src/
    mechanical_elements.c              (891 lines, 45+ functions)
    mechanical_transfer_function.c     (587 lines, 25+ functions)
    mechanical_state_space.c           (798 lines, 30+ functions)
    rotational_systems.c              (424 lines, 35+ functions)
    multi_dof_systems.c               (605 lines, 30+ functions)
    mechanical_applications.c          (482 lines, 35+ functions)
    mechanical_modeling.lean           (265 lines, 15+ theorems)
  tests/
    test_mechanical_modeling.c         (50+ assertion tests)
  examples/
    example_sdof_vibration.c          (SDOF vibration analysis)
    example_2dof_modal.c              (2-DOF modal analysis + FRF)
    example_quarter_car.c             (quarter-car suspension design)
  docs/
    knowledge-graph.md                 (L1-L9 complete mapping)
    coverage-report.md                 (coverage assessment)
    gap-report.md                      (missing items + priorities)
    course-alignment.md                (9-school curriculum mapping)
    course-tree.md                     (prerequisite dependency tree)
```

## Reference Textbooks

- Thomson & Dahleh, *Theory of Vibration with Applications* (1998)
- Rao, *Mechanical Vibrations* (2017)
- Goldstein, Poole & Safko, *Classical Mechanics* (2002)
- Ewins, *Modal Testing: Theory and Application* (2000)
- Craig & Kurdila, *Fundamentals of Structural Dynamics* (2006)
- Gillespie, *Fundamentals of Vehicle Dynamics* (1992)
- Rivin, *Passive Vibration Isolation* (2003)
- Siciliano et al., *Robotics: Modelling, Planning and Control* (2009)
- Ogata, *Modern Control Engineering* (2010)
- Brogan, *Modern Control Theory* (1991)
