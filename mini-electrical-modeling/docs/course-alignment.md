# Course Alignment — mini-electrical-modeling

## Nine-School Curriculum Mapping

### MIT — 6.302 Feedback Systems
| Topic | Module Coverage |
|-------|----------------|
| Electrical network modeling | `electrical_elements.c`, `transfer_function_electrical.c` |
| State-space from circuits | `electrical_state_space.c` |
| DC motor as canonical plant | `electromechanical_systems.c`, `example_dc_motor.c` |
| Two-port network analysis | `two_port_network.c` |

### Stanford — ENGR105 Feedback Control
| Topic | Module Coverage |
|-------|----------------|
| Circuit element s-domain models | `electrical_elements.c` (impedance functions) |
| Transfer functions of electrical networks | `transfer_function_electrical.c` |
| Servomotor transfer functions | `electromechanical_systems.c` |

### Berkeley — ME 132/232 Dynamic Systems
| Topic | Module Coverage |
|-------|----------------|
| RLC circuit state-space | `ss_from_series_rlc`, `ss_from_parallel_rlc` |
| Eigenvalue stability analysis | `ss_eigenvalues`, `ss_is_stable` |
| Electromechanical systems | `electromechanical_systems.c` |

### Caltech — CDS 110/212
| Topic | Module Coverage |
|-------|----------------|
| Controllability/observability of circuits | `ss_controllability_rank`, `ss_observability_rank` |
| Gramian computation | `ss_controllability_gramian`, `ss_observability_gramian` |
| Matrix exponential (Pade) | `ss_state_transition` |

### ETH — 151-0591 Control I
| Topic | Module Coverage |
|-------|----------------|
| Circuit transfer functions | `transfer_function_electrical.c` |
| State-space realization | `ss_controllable_canonical`, `ss_observable_canonical` |
| Discretization (ZOH) | `ss_discretize_zoh` |

### Cambridge — 3F2/4F2 Systems & Control
| Topic | Module Coverage |
|-------|----------------|
| Two-port parameter theory | `two_port_network.c` |
| Stability (Routh-Hurwitz via eigenvalues) | `ss_is_stable` |
| Frequency response methods | `tf_magnitude_db`, `tf_bode_magnitude_asymptotic` |

### Georgia Tech — ECE 6550/ME 6401
| Topic | Module Coverage |
|-------|----------------|
| Electromechanical actuators | `electromechanical_systems.c` |
| PWM drive modeling | `pwm_average_voltage`, `pwm_current_ripple` |
| Motor sizing and selection | `dc_motor_sizing`, `optimal_gear_ratio` |

### Purdue — ECE 602/ME 675
| Topic | Module Coverage |
|-------|----------------|
| Modified Nodal Analysis | `mna_build_matrix` |
| Gaussian elimination for circuits | `gaussian_elimination` |
| Industrial motor control | `electromechanical_systems.c` |

### Tsinghua — 自动控制原理 / 现代控制理论
| Topic | Module Coverage |
|-------|----------------|
| Transfer function analysis | `transfer_function_electrical.c` |
| State-space fundamentals | `electrical_state_space.c` |
| Canonical realizations | `ss_controllable_canonical`, `ss_observable_canonical` |
| Two-port network theory | `two_port_network.c` |

## Reference Textbooks
| Topic | Textbook | Coverage |
|-------|----------|----------|
| Electric Circuits | Dorf & Svoboda (2018) | Ch.1-14 — All passive elements, laws, analysis |
| Engineering EM | Hayt & Buck (2018) | Ch.13 — Transmission lines, S-parameters |
| Microwave Engineering | Pozar (2012) | Ch.4,12 — Two-port, amplifier design |
| Electric Machinery | Fitzgerald et al. | Ch.7-8 — DC machines |
| Matrix Computations | Golub & Van Loan (2013) | Ch.7 — QR algorithm |
| SPICE2 | Nagel (1975) | MNA formulation |
