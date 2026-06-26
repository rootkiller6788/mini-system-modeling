# Course Tree — mini-mechatronic-modeling

## Prerequisite Dependency Graph

```
                    Physics Fundamentals
                    ├── Mechanics (F=ma, τ=Jα)
                    ├── Electromagnetism (Lorentz, Faraday)
                    └── Circuit Theory (Ohm, KVL, KCL)
                              │
                    ┌─────────┼─────────┐
                    │         │         │
              Motor Theory  Sensors  Transmissions
              (DC/BLDC/    (Encoder/  (Gear/Screw/
               Stepper)     Resolver)  Belt)
                    │         │         │
                    └─────────┼─────────┘
                              │
                   Mechatronic Modeling
                   (this module)
                    ├── Electromechanical Coupling
                    ├── State-Space Representation
                    ├── Actuator Selection & Sizing
                    ├── Sensor Signal Processing
                    └── System Integration & Design
                              │
                    ┌─────────┼─────────┐
                    │         │         │
              Control Theory  Robotics  Manufacturing
              (PID/Cascade/  (Forward/   (CNC/Pick-
               State-Space)   Inverse)    and-Place)
```

## Linear Learning Path

1. **Foundations (L1-L2)**
   - Read: mechatronic_definitions.h, electromechanical_coupling.h
   - Understand: motor parameters, K_t=K_e equivalence, energy flow

2. **Dynamic Modeling (L3-L4)**
   - Read: mechatronic_state_space.h, electromechanical_coupling.c
   - Master: state-space construction, DC motor ODEs, eigenvalues

3. **Component Design (L5)**
   - Read: actuator_modeling.h, sensor_modeling.h
   - Master: motor selection, encoder methods, sensor filtering

4. **System Integration (L5-L6)**
   - Read: mechatronic_system_design.h
   - Master: cascade controller design, servo axis integration

5. **Applications (L6-L7)**
   - Run: example_dc_servo, example_ball_screw, example_robot_joint
   - Understand: end-to-end mechatronic design workflow

## Cross-Module Dependencies within mini-automation-theory

```
mini-mechatronic-modeling
  ├── depends on: mini-mechanical-modeling (gear, load dynamics)
  ├── depends on: mini-electrical-modeling (motor, circuit models)
  ├── feeds into: mini-state-space-theory (controllability, observability)
  ├── feeds into: mini-pole-placement-observer (cascade control)
  └── feeds into: mini-automation-application-systems (robot, CNC)
```

## Knowledge Level Progression

| Level | Prerequisites | This Module | Leads To |
|-------|--------------|-------------|----------|
| L1 | Basic physics | Motor/gear/sensor types | Component datasheet reading |
| L2 | Calculus | Electromechanical analogy | Bond graph modeling |
| L3 | Linear algebra | State-space matrices | Optimal control |
| L4 | ODEs | Newton+Lorentz+Faraday laws | Multi-physics FEM |
| L5 | Numerical methods | RK4, QR, Gaussian elimination | Real-time simulation |
| L6 | System design | Servo axis, robot joint | Industrial automation |
| L7 | Domain knowledge | HDD, CNC, pick-and-place | Product development |
| L8 | Graduate control | Flexible joint, piezo, Kalman | Research |
| L9 | PhD-level | MEMS, smart materials, digital twin | Frontier research |
