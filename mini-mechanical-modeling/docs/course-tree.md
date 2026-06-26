# Course Tree — mini-mechanical-modeling

## Prerequisite Dependencies

```
Newton's Laws (F=ma, Hooke's)
  +-- L1 Definitions: mass, spring, damper
  +-- L2 Concepts: energy, damping, friction
  +-- L3 Math: ODE parameters (wn, zeta, wd)
  +-- L4 Conservation laws
      |
      +-- L5 SDOF Transfer Functions
      |     +-- L5 Frequency/Step Response
      |     +-- L5 Transmissibility Analysis
      |     +-- L5 State-Space Construction
      |           +-- L5 Eigenvalue Analysis
      |           +-- L5 RK4 Simulation
      |           +-- L5 Discretization
      |
      +-- L5 Rotational Systems
      |     +-- L5 Gear Ratio & Reflected Parameters
      |     +-- L5 Rotational-Translational Conversion
      |     +-- L6 Rotating Unbalance
      |
      +-- L5 MDOF Systems
            +-- L5 Eigenvalue Problem (K*phi=lambda*M*phi)
            +-- L5 Modal Parameters & Orthogonality
            +-- L5 Frequency Response Functions
            +-- L8 Rayleigh Damping
            +-- L8 Guyan Reduction
            +-- L8 MAC Criterion

Applications depend on all above:
  +-- L6/L7 Quarter-Car Suspension
  +-- L7 Vibration Isolation
  +-- L7 Robot Joint Mechanics
  +-- L7 Aerospace/UAV Dynamics
  +-- L7 Precision Machine Dynamics

L9 Frontiers (documented):
  +-- Nonlinear friction
  +-- Active vibration control
  +-- Nanoscale mechanics

## Knowledge Dependencies Graph

SDOF Modeling → 2-DOF/MDOF → Modal Analysis → FRF → Applications
Mechanical Elements → TF → SS → Simulation → Control Design
Rigid Body Dynamics → Rotational → Geared Systems → Robotics
```

## Cross-Module Dependencies
This module provides mechanical plant models for:
- `mini-transfer-function` (TF from Newton-Euler ODE)
- `mini-block-diagram` (mechanical subsystem blocks)
- `mini-state-space-theory` (state matrices)
- `mini-pole-placement-observer` (plant model for control)
- `mini-automation-application-systems` (robot, vehicle models)
