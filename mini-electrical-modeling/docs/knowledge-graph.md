# Knowledge Graph — mini-electrical-modeling

## L1: Definitions (Complete)

| # | Concept | C Implementation | Header |
|---|---------|-----------------|--------|
| 1 | Resistance R [ohm] | `Resistor` struct | `electrical_elements.h` |
| 2 | Capacitance C [F] | `Capacitor` struct | `electrical_elements.h` |
| 3 | Inductance L [H] | `Inductor` struct | `electrical_elements.h` |
| 4 | Voltage/Current Sources | `ElectricalSource` struct | `electrical_elements.h` |
| 5 | Phasor representation | `Phasor` struct | `electrical_elements.h` |
| 6 | Transfer function H(s) | `TransferFunction` struct | `transfer_function_electrical.h` |
| 7 | State-space model (A,B,C,D) | `StateSpace` struct | `electrical_state_space.h` |
| 8 | Z-parameters | `ZParameters` struct | `two_port_network.h` |
| 9 | Y-parameters | `YParameters` struct | `two_port_network.h` |
| 10 | H-parameters (hybrid) | `HParameters` struct | `two_port_network.h` |
| 11 | ABCD (transmission) parameters | `ABCDParameters` struct | `two_port_network.h` |
| 12 | S-parameters (scattering) | `SParameters` struct | `two_port_network.h` |
| 13 | DC motor parameters | `DCMotorParams` struct | `electromechanical_systems.h` |
| 14 | Node, branch, loop, mesh | `CircuitBranch`, `Circuit` | `circuit_topology.h` |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Ohm's Law (V=IR) | `ohms_law_voltage`, `ohms_law_current` |
| 2 | Joule's Law (P=I^2R) | `power_dc_i2r`, `power_dc_v_i`, `power_dc_v2r` |
| 3 | Impedance Z(s) = R + sL + 1/(sC) | `impedance_rlc_series`, `_parallel` |
| 4 | AC power (S=P+jQ) | `ac_real_power`, `ac_reactive_power`, `ac_complex_power` |
| 5 | KVL / KCL | `kvl_verify`, `kcl_verify` |
| 6 | Thevenin/Norton equivalence | `thevenin_resistive_divider`, `norton_from_thevenin` |
| 7 | Maximum power transfer | `maximum_power_transfer` |
| 8 | Superposition principle | `superposition_two_voltage_sources` |
| 9 | Tellegen's theorem | `tellegen_power_sum` |
| 10 | Reciprocity (Z12=Z21) | `tp_is_reciprocal_z` |
| 11 | Passivity (Re(Z)>=0) | `tp_is_passive_z` |
| 12 | DC motor steady-state | `dc_motor_steady_state_speed` |
| 13 | Generator model | `dc_generator_terminal_voltage` |
| 14 | Frequency response | `tf_magnitude_db`, `tf_phase_rad` |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Complex impedance | `double complex` in phasors, Z(s) evaluation |
| 2 | Time constants (RC, RL) | `time_constant_rc`, `time_constant_rl` |
| 3 | Resonant frequency ω0=1/sqrt(LC) | `resonant_frequency_rad_s` |
| 4 | Damping factor ζ | `damping_factor_series_rlc` |
| 5 | Quality factor Q | `quality_factor_series_rlc` |
| 6 | Series/parallel equivalents | `equivalent_resistance_series`, etc. |
| 7 | Polynomial representation (TF) | `TransferFunction` numerator/denominator |
| 8 | Matrix operations (SS) | A,B,C,D matrix accessors, state transition |
| 9 | Eigenvalues (QR algorithm) | `ss_eigenvalues` |
| 10 | Modal parameters | `ss_modal_from_eigenvalue` |
| 11 | 2x2 complex matrices (two-port) | All Z/Y/S/ABCD parameter structs |

## L4: Fundamental Laws (Complete)

| # | Law/Theorem | Implementation |
|---|-------------|---------------|
| 1 | Ohm's Law | `ohms_law_voltage`, `ohms_law_current` |
| 2 | Kirchhoff Voltage Law (KVL) | `kvl_verify` |
| 3 | Kirchhoff Current Law (KCL) | `kcl_verify` |
| 4 | Joule's Law (power dissipation) | `power_dc_i2r` |
| 5 | Faraday's Law (induction) | `inductor_voltage` |
| 6 | Energy conservation in RLC | `energy_capacitor`, `energy_inductor` |
| 7 | Tellegen's Theorem | `tellegen_power_sum` |
| 8 | Maximum Power Transfer Theorem | `maximum_power_transfer` |
| 9 | Reciprocity Theorem | `tp_is_reciprocal_z` |
| 10 | Hurwitz Stability (eigenvalues) | `ss_is_stable` |
| 11 | Rollett Stability Criterion | `sp_rollett_stability` |

## L5: Engineering Methods (Complete)

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Nodal analysis (G matrix) | `nodal_analysis_build` |
| 2 | Modified Nodal Analysis (MNA) | `mna_build_matrix` |
| 3 | Gaussian elimination (pivoting) | `gaussian_elimination` |
| 4 | Durand-Kerner root-finding | `tf_find_poles`, `tf_find_zeros` |
| 5 | Binary search for cutoff freq | `tf_find_cutoff_freq` |
| 6 | Asymptotic Bode approximation | `tf_bode_magnitude_asymptotic` |
| 7 | Partial fraction step response | `tf_step_response` |
| 8 | Filter topology synthesis | `tf_sallen_key_lowpass`, `tf_mfb_lowpass` |
| 9 | TF arithmetic (cascade/parallel/feedback) | `tf_cascade`, `tf_parallel`, `tf_negative_feedback` |
| 10 | QR eigenvalue algorithm | `ss_eigenvalues` |
| 11 | Pade matrix exponential (Higham 2009) | `ss_state_transition` |
| 12 | RK4 state simulation | `ss_simulate_rk4` |
| 13 | RLC→SS (systematic method) | `ss_from_series_rlc`, etc. |
| 14 | Controllable/observable canonical forms | `ss_controllable_canonical`, `ss_observable_canonical` |
| 15 | ZOH discretization (Van Loan 1978) | `ss_discretize_zoh` |
| 16 | Z↔Y↔H↔ABCD↔S conversions | All `zp_to_*` and `*_to_zp` functions |
| 17 | Two-port interconnection | `tp_series/parallel/cascade_interconnect` |
| 18 | Pi/T network synthesis | `tp_abcd_pi_network`, `tp_abcd_t_network` |
| 19 | Standard filter design with E12/E6 | `design_rc_lowpass_to_cutoff` |
| 20 | Group delay computation | `tf_group_delay` |
| 21 | Motor parameter identification | `dc_motor_identify_parameters` |

## L6: Canonical Problems (Complete)

| # | Problem | Implementation |
|---|---------|---------------|
| 1 | RC low-pass/high-pass filter | `tf_rc_lowpass`, `tf_rc_highpass` + example |
| 2 | RLC resonance analysis | `tf_rlc_bandpass/lowpass/highpass` + example |
| 3 | DC motor speed-torque curve | `dc_motor_steady_state_*` + example |
| 4 | Two-stage RF amplifier cascade | Example: S-parameter cascade analysis |
| 5 | RLC state-space realization | `ss_from_series_rlc`, `ss_from_parallel_rlc` |
| 6 | Pi/T network analysis | `tp_abcd_pi_network`, `tp_abcd_t_network` |
| 7 | Filter design to specification | `design_rc_lowpass_to_cutoff` |
| 8 | Motor startup transient | Example: RK4 simulation |

## L7: Applications (Partial+ — 3 implemented)

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | PWM motor drive (Tesla/SpaceX relevant) | `pwm_average_voltage`, `pwm_current_ripple` |
| 2 | Motor sizing for robotics (DC motor selection) | `dc_motor_sizing`, `optimal_gear_ratio` |
| 3 | RF amplifier design (2.4 GHz WiFi) | Two-port S-parameter analysis example |

## L8: Advanced Topics (Partial — 1 implemented)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Mixed-mode S-parameters (differential) | `sp_to_mixed_mode` |
| 2 | Nonlinear circuit analysis | Documented, not yet implemented |
| 3 | Time-varying systems | Documented, not yet implemented |
| 4 | Electromagnetic compatibility (EMC) | Documented, not yet implemented |
| 5 | Multiphysics coupling (thermal-electrical) | `self_heating_temp_rise`, documented |

## L9: Research Frontiers (Partial — documented)

| # | Topic | Status |
|---|-------|--------|
| 1 | Nanoelectronic device modeling (quantum transport) | Documented |
| 2 | Memristor modeling (Chua 1971) | Documented |
| 3 | Non-Foster matching networks | Documented |
| 4 | Wireless power transfer optimization | Documented |
