# Knowledge Graph — mini-mechatronic-modeling

## L1: Definitions ✅ Complete
- DC Motor parameters (R_a, L_a, K_e, K_t, J_m, B_m)
- BLDC motor parameters (R_ph, L_ph, K_e_ll)
- Stepper motor parameters (holding torque, step angle)
- Gear train (ratio, efficiency, backlash, stiffness)
- Lead screw (pitch, efficiency, J_screw)
- Belt drive (ratio, belt stiffness, preload)
- Rotary sensors (encoder incremental/absolute, resolver, tachometer, Hall)
- Linear sensors (LVDT, strain gauge, laser, capacitive, eddy current)
- IMU parameters (accel/gyro range, noise, bias)
- Mechatronic load model (friction types, elastic coupling)
- Motor selection result (margins, feasibility)
- Mechatronic subsystem (reflected dynamics, ω_n, ζ)
- Servo axis (cascade gains, bandwidths, limits)
- Multi-axis manipulator (DH parameters)
- PWM drive (V_dc, f_pwm, dead time)
- Voice coil actuator (K_f, stroke, M_moving)
- Piezoelectric actuator (d33, stiffness, capacitance)

## L2: Core Concepts ✅ Complete
- Electromechanical analogy (F-V, F-I mappings)
- Lorentz force and its motor application
- Back-EMF as motional EMF
- K_t = K_e equivalence in SI
- DC motor torque-speed characteristic
- Motor efficiency and power flow
- Inertia matching principle
- Cascade control structure (current→velocity→position)
- Sensor measurement principles (encoder, resolver, Hall, LVDT, strain gauge)
- Friction regimes (Coulomb, viscous, Stribeck, stiction)
- Backlash dead-zone behavior
- Regenerative braking
- Bandwidth allocation in cascaded loops
- Thermal time constants (electrical vs mechanical)

## L3: Mathematical Structures ✅ Complete
- State-space matrices (A, B, C, D) in row-major storage
- Matrix operations (multiply, inverse, determinant, rank)
- Complex eigenvalues and modal analysis
- Natural frequency ω_n and damping ratio ζ extraction
- Transfer function from state-space (Cramer's rule)
- DC gain matrix computation
- Gaussian elimination with partial pivoting
- QR algorithm (Hessenberg reduction + Wilkinson shifts)
- Householder reflections for Hessenberg form
- Discretization (ZOH, Pade approximation)
- Quadrature encoder Gray code state machine
- Sensor noise models (ARW, bias instability, Allan variance)

## L4: Fundamental Laws ✅ Complete
- Lorentz force law: F = I·L×B
- Faraday's law: motional EMF = B·L·v
- DC motor back-EMF: V_bemf = K_e·ω
- DC motor torque: T_em = K_t·I_a
- K_t = K_e equivalence (energy conservation proof)
- Motor dynamic equations: V_a = R·I + L·dI/dt + K_e·ω
- Newton's 2nd law (rotation): J·dω/dt = ΣT
- Kalman controllability condition: rank(C) = n
- Kalman observability condition: rank(O) = n
- Lyapunov stability: Re(λ_i) < 0 for all i
- Power conservation in gear train: P_out = η·P_in
- Reciprocal electromechanical coupling

## L5: Engineering Methods ✅ Complete
- RMS torque calculation for motor sizing
- Peak torque and average speed from duty cycle
- Load reflection through gear train
- Optimal gear ratio (inertia matching): N_opt = √(J_L/J_m)
- Motor selection algorithm with feasibility checks
- First-order motor thermal model
- Resistance temperature correction
- Stepper motor torque-speed curve
- Stepper start-stop and slew rate limits
- Encoder velocity: M-method, T-method, M/T hybrid
- Encoder direction detection (Gray code)
- Resolver angle extraction (atan2)
- Sensor low-pass filtering (IIR)
- Moving average filter (FIR with circular buffer)
- Complementary filter for IMU
- 1D Kalman filter for position-velocity
- Wheatstone bridge (quarter, half, full)
- Load cell force calibration
- PI controller with anti-windup
- PID with filtered derivative
- Cascade bandwidth allocation
- Current PI design (pole-zero cancellation)
- Velocity PI design (symmetrical optimum)
- Position P design
- Pade [1/1] matrix exponential approximation
- ZOH discretization

## L6: Canonical Problems ✅ Complete
- DC servo motor design (motor + gear + load + controller)
- Ball screw CNC axis (sizing, torque, state-space)
- Industrial robot joint (shoulder joint, harmonic drive)
- Servo axis step response analysis
- Motor-gear-load inertia matching
- Ball screw linear-to-rotary conversion
- Robot joint gravity and inertial torque
- HDD seek time (bang-bang acceleration)
- Pick-and-place cycle time (trapezoidal profile)
- Antenna scan pattern time (raster scan)

## L7: Applications ⚠️ Partial
- HDD head positioning servo (hard disk drive)
- Pick-and-place Cartesian robot
- CNC ball screw axis (manufacturing automation)
- Industrial robot arm joint
- Antenna pointing system (radar/satellite)
- Printer carriage analysis
- Regenerative braking energy recovery
- Motion profile energy analysis

## L8: Advanced Topics ⚠️ Partial
- Flexible joint dynamics (Spong model, 4-state)
- Piezoelectric actuators (stroke under load, power dissipation)
- BLDC FOC torque model
- Kalman filtering for velocity estimation
- Allan variance noise characterization

## L9: Research Frontiers ⚠️ Partial (Documented)
- MEMS mechatronics (micro-scale sensors and actuators)
- Soft robotics actuation
- Smart materials (SMA, EAP, magnetostrictive)
- Multi-axis coordination and vibration suppression
- Digital twin for mechatronic systems
