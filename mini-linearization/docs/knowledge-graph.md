# Knowledge Graph — mini-linearization

## L1 Definitions
EquilibriumPoint, JacobianMatrix, LinearizedSystem, SmallSignalVariables, LinearizationSpectrum,
LieDerivative, LieBracket, RelativeDegree, Diffeomorphism, NormalForm, MIMORelativeDegree,
TaylorExpansion, MultiIndex, TaylorRemainderForm, DescribingFunction, LimitCyclePrediction,
DualInputDF, KoopmanObservable, KoopmanEigenfunction, KoopmanMode, DMDResult, EDMDConfig,
E DMDBasisType, KoopmanDecomposition, GainSchedule, GainScheduleEntry, GainScheduledController,
GSInterpMethod, NeuralFeedbackLinearizer, DeepKoopmanModel

## L2 Core Concepts
Local validity of linearization, Hartman-Grobman topological conjugacy, small-signal vs large-signal,
Input-state vs input-output linearization, zero dynamics, internal dynamics, minimum-phase systems,
Koopman linear representation of nonlinear dynamics, LPV gain scheduling, limit cycle existence,
Bumpless transfer, hyperbolic equilibrium

## L3 Mathematical Structures
State-space matrices (A,B,C,D), Jacobian/Hessian, Taylor coefficients (multi-index), Lie algebra,
Controllability/Observability matrices, Gramians, DMD/EDMD matrices, Hankel matrix, decoupling matrix,
Carleman embedding, matrix exponential, condition number, SVD

## L4 Fundamental Theorems
Lyapunov Indirect Method (1892), Hartman-Grobman Theorem, Taylor's Theorem with remainder,
Frobenius Theorem, Brockett's Necessary Condition (1983), Center Manifold Theorem,
Koopman Spectral Theorem, Kalman Controllability/Observability

## L5 Engineering Methods
Numerical Jacobian (forward/central O2/O4), complex-step, analytical Jacobian, finite differences,
Feedback linearization (IO/IS), DMD/EDMD/Kernel DMD, describing function analysis,
Gain scheduling (nearest/linear/blending/LPV), LQR via iterative ARE, pole placement,
Newton-Raphson equilibrium finding, continuation method, balanced truncation

## L6 Canonical Problems
Inverted pendulum linearization, magnetic levitation, DC motor saturation analysis,
Quadrotor gain scheduling across flight envelope

## L7 Applications (Partial)
- Toyota Prius MG2 DC motor saturation modeling
- Magnetic levitation (maglev train reference)
- Quadrotor attitude control (NASA/SpaceX/Tesla reference)
- F-35 flight control gain scheduling (documented)

## L8 Advanced Methods (Partial)
- Koopman DMD/EDMD with basis dictionaries
- Carleman linearization embedding
- Ensemble (Monte Carlo) Jacobian
- Balanced truncation model reduction
- Robust feedback linearization (sliding mode)
- Adaptive feedback linearization
- Sparse Jacobian via graph coloring

## L9 Research Frontiers (Partial)
- Koopman operator theory for nonlinear control (Koopman 1931, Mezic 2005)
- Carleman linearization (Carleman 1932, Forets-Pouly 2017)
- Deep Koopman autoencoders (Lusch-Kutz-Brunton 2018)
- Neural network-based feedback linearization
- Hankel-DMD with delay embedding (Arbabi-Mezic 2017)
- Online/streaming DMD for real-time analysis
