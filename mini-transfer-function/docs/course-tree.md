# Course Tree — mini-transfer-function

## Prerequisites
- Complex variables (Laplace transform)
- Linear algebra (matrix operations)
- Differential equations
- Basic control concepts (feedback)

## Internal Dependencies

```
Polynomial Algebra (tf_polynomial)
    |
    v
Transfer Function Core (transfer_function)
    |
    +---> TF Analysis (tf_analysis) -- stability, frequency, time
    |
    +---> TF Interconnections (tf_interconnections) -- algebra, reduction
    |
    +---> TF Conversion (tf_conversion) -- SS, ZPK, discrete
    |
    +---> TF Advanced (tf_advanced) -- reduction, apps, fractional
```

## External Dependencies
This module is a dependency for:
- mini-block-diagram (uses transfer function)
- mini-time-domain-analysis (response analysis)
- mini-frequency-domain (Bode, Nyquist)
- mini-classical-compensator (lead/lag design)
- mini-state-space-theory (conversion, analysis)

## L9: Research Frontiers
- Fractional-order control (CRONE, fractional PID)
- Nanoscale/non-Fourier heat transfer TF models
- Quantum control transfer functions
