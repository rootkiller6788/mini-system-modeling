# Knowledge Graph — mini-mechanical-modeling

## Nine-Layer Knowledge Coverage

### L1: Definitions (Complete ✅)
- Mass, Spring, Damper structs with physical parameters
- FrictionModel: Coulomb, viscous, Stribeck, static
- Natural frequency, damping ratio, logarithmic decrement
- Mass moments of inertia (point, rod, cylinder, sphere, plate)
- Mechanical transfer function (rational polynomial)
- State-space (A,B,C,D matrices), discrete-time SS
- Rotor, TorsionSpring, RotDamper, GearPair structs
- MDOFSystem with M,K,B,F matrices
- ModalParams: natural freq, damping, modal mass/stiffness, participation
- QuarterCarModel, IsolatorReq/Design, RobotJoint structs

### L2: Core Concepts (Complete ✅)
- Newton's 2nd law: F=ma (translational and rotational)
- Hooke's law: F=kx (linear), T=k_t*theta (torsional)
- Viscous damping: F=bv, T=b_t*omega
- Coulomb/Stribeck friction models
- Kinetic/potential energy, energy dissipation
- Mechanical impedance, mobility
- Resonance, antiresonance, vibration isolation
- Gear ratio, reflected inertia/damping/stiffness
- Modal superposition, orthogonality
- Rayleigh proportional damping

### L3: Mathematical Structures (Complete ✅)
- 2nd-order ODE parameters: omega_n, zeta, omega_d
- Complex mechanical impedance Z_m(jw)
- Rational polynomial transfer functions H(s)=N(s)/D(s)
- State-space matrices (n x n, n x m, p x n)
- Generalized eigenvalue problem K*phi = lambda*M*phi
- Mass/stiffness/damping matrix algebra
- Frequency response functions (receptance, mobility, accelerance)
- Modal Assurance Criterion (MAC)
- Polynomial convolution for TF operations

### L4: Fundamental Laws (Complete ✅)
- Newton's 2nd law (translational F=ma, rotational T=J*alpha)
- Hooke's law for springs (linear F=kx, torsional T=k_t*theta)
- Viscous damping law F=bv
- Conservation of mechanical energy
- D'Alembert principle for MDOF systems
- Lagrange equations (implicit in MDOF formulation)
- Orthogonality of mode shapes (mass and stiffness weighted)

### L5: Engineering Methods (Complete ✅)
- Equivalent spring/damper combinations (series, parallel)
- Vibration transmissibility analysis
- Time-domain response specification computation
- Half-power bandwidth damping estimation
- Transfer function derivation from ODE (SDOF, 2-DOF, base excitation)
- TF cascading, parallel, feedback operations
- RK4 state-space simulation
- State transition matrix (Pade approximation)
- ZOH and FOH discretization
- Controllability/observability rank test
- Modal parameter extraction (mass, stiffness, participation)
- Rayleigh damping coefficient computation
- Guyan static condensation
- Modal superposition for forced response
- Duhamel integral for transient analysis
- Newmark-beta integration for MDOF

### L6: Canonical Problems (Complete ✅)
- SDOF free/forced vibration analysis (example_sdof_vibration.c)
- 2-DOF modal analysis and FRF (example_2dof_modal.c)
- Quarter-car suspension design (example_quarter_car.c)
- Inverted pendulum state-space model
- Rotating unbalance response
- Rotor dynamics (Jeffcott model)
- Vibration isolator design
- Dynamic vibration absorber tuning

### L7: Applications (Partial+ ✅ — 3 applications)
- Automotive suspension design (quarter-car, half-car, full-car heave/pitch) — Detroit/Tesla
- Precision vibration isolation (semiconductor equipment)
- Industrial robot joint mechanics and PD tuning (ISO 9283)
- Quadrotor arm vibration (UAV dynamics)
- Spacecraft reaction wheel and momentum wheel sizing
- Satellite solar panel deployment dynamics
- Launch vehicle Pogo stability
- Machine tool spindle critical speed

### L8: Advanced Topics (Partial+ ✅ — 3 topics)
- Guyan static condensation for model reduction
- Modal Assurance Criterion for mode shape correlation
- Rayleigh damping model fitting
- Two-stage vibration isolation
- Skyhook damping concept
- Newmark-beta transient integration
- Harmonic drive stiffness modeling

### L9: Research Frontiers (Partial ⚠️ — documented)
- Nonlinear friction modeling (Stribeck curve)
- Active/semi-active suspension control
- Morphing wing actuator dynamics
- Nanoscale mechanical resonators
