# mini-linearization

Nonlinear system linearization for control and automation — Jacobian, Taylor expansion, feedback linearization, Koopman operator, describing functions, and gain scheduling.

## Module Status: COMPLETE ✅

- **L1-L6: Complete** (all core definitions, concepts, math, laws, methods, problems)
- **L7: Partial** (3 applications: maglev, DC motor saturation, quadrotor gain scheduling)
- **L8: Partial** (Koopman DMD, Carleman embedding, Monte Carlo Jacobian, balanced truncation)
- **L9: Partial** (Koopman operator theory, Carleman linearization, neural feedback linearization)

**Line Count**: `include/` + `src/` >= 3000 lines ✅

## Nine-Layer Knowledge Coverage

| Level | Name | Status | Key Items |
|-------|------|--------|-----------|
| **L1** | Definitions | ✅ Complete | Equilibrium point, Jacobian, Lie derivative/bracket, diffeomorphism, relative degree, describing function, Koopman operator, gain schedule |
| **L2** | Core Concepts | ✅ Complete | Local validity, small-signal model, input-state/input-output linearization, zero dynamics, Koopman spectral theory, LPV systems |
| **L3** | Math Structures | ✅ Complete | State-space (A,B,C,D), Jacobian/Hessian matrices, Lie algebra, Taylor coefficients, Gramians, DMD matrix |
| **L4** | Fundamental Laws | ✅ Complete | Lyapunov indirect method, Hartman-Grobman, Taylor's theorem, Frobenius theorem, Brockett's condition |
| **L5** | Engineering Methods | ✅ Complete | Numerical/analytical Jacobian, feedback linearization, DMD/EDMD, gain scheduling, DF analysis, Newton-Raphson |
| **L6** | Canonical Problems | ✅ Complete | Inverted pendulum, magnetic levitation, DC motor saturation, quadrotor (4 examples) |
| **L7** | Applications | ⚠️ Partial | Toyota DC motor, maglev train, quadrotor (NASA/Tesla/SpaceX), aerospace scheduling (F-35) |
| **L8** | Advanced Topics | ⚠️ Partial | Koopman DMD/EDMD, Carleman, ensemble Jacobian, Monte Carlo, balanced truncation |
| **L9** | Research Frontiers | ⚠️ Partial | Koopman operator, Carleman linearization, deep Koopman, neural feedback linearization |

## Core Definitions

- **Equilibrium Point** — `EquilibriumPoint` struct, `create_equilibrium`
- **Jacobian Matrix** J = df/dx — `JacobianMatrix` struct, `numerical_jacobian_central`
- **Linearized System** (A,B,C,D) — `LinearizedSystem` struct, `linearize_system`
- **Lie Derivative** L_f h(x) — `LieDerivative` struct, `compute_lie_derivative`
- **Lie Bracket** [f,g] — `LieBracket` struct, `compute_lie_bracket`
- **Relative Degree** — `RelativeDegree` struct, `compute_relative_degree`
- **Describing Function** N(A) — `DescribingFunction` struct, `df_ideal_relay/saturation/dead_zone/backlash/cubic/quantizer`
- **Koopman Operator** — `KoopmanObservable`, `KoopmanEigenfunction`, `KoopmanMode` structs
- **DMD** — `DMDResult` struct, `compute_dmd`, `dmd_predict`
- **Gain Schedule** — `GainSchedule` struct, `create_gain_schedule`, `interpolate_gain_schedule`

## Core Theorems

| Theorem | Implementation |
|---------|---------------|
| Lyapunov Indirect Method (1892) | `apply_lyapunov_indirect` |
| Hartman-Grobman Theorem | `check_equilibrium_hyperbolicity` |
| Taylor's Theorem with Remainder | `taylor_remainder_lagrange_1d` |
| Frobenius Integrability | `check_frobenius_involutivity` |
| Brockett's Necessary Condition | `brockett_necessary_condition` |
| Kalman Controllability/Observability | `compute_controllability_matrix` / `compute_observability_matrix` |
| Koopman Spectral Theorem | `compute_dmd`, `identify_koopman_eigenfunctions` |
| DF Harmonic Balance | `predict_limit_cycle_df` |

## Core Algorithms

| Algorithm | Implementation |
|-----------|---------------|
| Numerical Jacobian (central FD O(eps^2)) | `numerical_jacobian_central` |
| Complex-step Jacobian | `compute_jacobian_complex_step` |
| QR Eigenvalue (Hessenberg + Francis) | `compute_eigenvalues_qr` |
| Feedback Linearization Controller | `input_output_linearization_controller` |
| DMD (SVD-based) | `compute_dmd` |
| EDMD with basis dictionary | `compute_edmd` |
| Gain Schedule Interpolation | `interpolate_gain_schedule` |
| Newton-Raphson Equilibrium Finding | `find_equilibrium_newton` |
| LQR via Iterative ARE | `design_lqr_schedule` |
| Matrix Exponential (Pade + scaling) | `matrix_exponential_pade` |

## Classic Problems (examples/)

1. **Inverted Pendulum** — `examples/example_pendulum.c`
2. **Magnetic Levitation** — `examples/example_maglev.c`
3. **DC Motor Saturation** (Toyota) — `examples/example_dc_motor_saturation.c`
4. **Quadrotor Gain Scheduling** (NASA/SpaceX) — `examples/example_quadrotor.c`

## Nine-School Course Mapping

| School | Key Course | Topics |
|--------|-----------|--------|
| **MIT** | 6.302 Feedback Systems | Linearization, gain scheduling |
| **Stanford** | ENGR209/AA271 | Geometric nonlinear control, aerospace scheduling |
| **Berkeley** | ME 237 Nonlinear Control | Lie derivatives, normal forms |
| **Caltech** | CDS 240 | LPV gain scheduling, robust linearization |
| **ETH** | 151-0591 Control I/II | Jacobian linearization, Lyapunov |
| **Cambridge** | 4F2 Nonlinear Systems | Describing functions, limit cycles |
| **Michigan** | AERO 550 | Flight dynamics linearization |
| **Purdue** | ME 675 | DMD/Koopman identification |
| **Tsinghua** | 自动控制原理 | Linearization fundamentals |

## Build and Test

```bash
make test       # Run 15 tests
make examples   # Build all examples
make run-examples
make clean
```

## File Structure

```
mini-linearization/
├── Makefile, README.md
├── include/ (7 files) — linearization_core, taylor, feedback_linearization,
│        koopman, describing_function, jacobian_methods, gain_scheduling
├── src/ (8 files) — all core implementations
├── tests/ — test_linearization.c (15 tests)
├── examples/ — pendulum, maglev, dc_motor, quadrotor
└── docs/ — knowledge-graph, coverage-report, gap-report,
          course-alignment, course-tree
```

## References

- Khalil, *Nonlinear Systems* (3rd ed, 2002)
- Isidori, *Nonlinear Control Systems* (3rd ed, 1995)
- Slotine & Li, *Applied Nonlinear Control* (1991)
- Ogata, *Modern Control Engineering* (2010)
- Brunton & Kutz, *Data-Driven Science and Engineering* (2019)
- Golub & Van Loan, *Matrix Computations* (2013)
