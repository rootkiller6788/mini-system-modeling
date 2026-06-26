/**
 * @file two_port_network.h
 * @brief Two-port network parameters and transformations
 *
 * Knowledge Coverage:
 *   L1 - Definitions: Z-parameters, Y-parameters, H-parameters (hybrid),
 *        ABCD (transmission) parameters, S-parameters (scattering)
 *   L2 - Core Concepts: impedance matching, insertion loss, return loss,
 *        VSWR, network reciprocity, symmetry, passivity conditions
 *   L3 - Mathematical Structures: complex 2x2 matrix algebra,
 *        parameter conversion formulas
 *   L4 - Fundamental Laws: reciprocity (Z12 = Z21 for passive networks),
 *        passivity (positive real matrix conditions)
 *   L5 - Engineering Methods: parameter measurement techniques,
 *        network interconnection (series/parallel/cascade)
 *
 * Reference:
 *   MIT 6.630 Electromagnetic Theory — S-parameter fundamentals
 *   Berkeley EE 117 — Transmission lines and two-port theory
 *   Pozar, Microwave Engineering (2012) §4
 *   Carson, Electromagnetic Theory (section on network parameters)
 */

#ifndef TWO_PORT_NETWORK_H
#define TWO_PORT_NETWORK_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define TP_EPSILON 1e-12

/* L1: Z-parameters (impedance matrix) */
/* V1 = Z11*I1 + Z12*I2, V2 = Z21*I1 + Z22*I2 */
typedef struct {
    double complex z11, z12, z21, z22;
} ZParameters;

/* L1: Y-parameters (admittance matrix) */
/* I1 = Y11*V1 + Y12*V2, I2 = Y21*V1 + Y22*V2 */
typedef struct {
    double complex y11, y12, y21, y22;
} YParameters;

/* L1: H-parameters (hybrid — common for BJT) */
/* V1 = H11*I1 + H12*V2, I2 = H21*I1 + H22*V2 */
typedef struct {
    double complex h11, h12, h21, h22;
} HParameters;

/* L1: G-parameters (inverse hybrid) */
/* I1 = G11*V1 + G12*I2, V2 = G21*V1 + G22*I2 */
typedef struct {
    double complex g11, g12, g21, g22;
} GParameters;

/* L1: ABCD transmission parameters */
/* V1 = A*V2 - B*I2, I1 = C*V2 - D*I2 (IEEE convention: output current into network) */
/* Alternative: V1 = A*V2 + B*I2, I1 = C*V2 + D*I2 (output current out of network) */
typedef struct {
    double complex a, b, c, d;
} ABCDParameters;

/* L1: S-parameters (scattering — RF/microwave) */
/* b1 = S11*a1 + S12*a2, b2 = S21*a1 + S22*a2 */
/* Where ai, bi are normalized incident/reflected power waves */
typedef struct {
    double complex s11, s12, s21, s22;
    double z0; /* Reference impedance typically 50 ohm */
} SParameters;

/* L2: Parameter conversion functions */

/**
 * @brief Convert Z-parameters to Y-parameters
 *
 * Theorem: Y = Z^{-1} for a 2x2 matrix:
 *   Y11 = Z22/Δ, Y12 = -Z12/Δ, Y21 = -Z21/Δ, Y22 = Z11/Δ
 *   where Δ = Z11·Z22 - Z12·Z21 (determinant of Z)
 *
 * @return 0 on success, -1 if singular (Δ = 0)
 * Complexity: O(1)
 */
int zp_to_yp(const ZParameters *zp, YParameters *yp);

/**
 * @brief Convert Y to Z — inverse of above
 * Complexity: O(1)
 */
int yp_to_zp(const YParameters *yp, ZParameters *zp);

/**
 * @brief Convert Z to H parameters
 *
 * Formulas:
 *   h11 = Δ/Z22, h12 = Z12/Z22, h21 = -Z21/Z22, h22 = 1/Z22
 *   where Δ = Z11·Z22 - Z12·Z21
 *
 * Reference: Chua, Desoer, and Kuh (1987) Table 5.1
 * Complexity: O(1)
 */
int zp_to_hp(const ZParameters *zp, HParameters *hp);

/**
 * @brief Convert H to Z parameters
 * Complexity: O(1)
 */
int hp_to_zp(const HParameters *hp, ZParameters *zp);

/**
 * @brief Convert Z to ABCD (transmission) parameters
 *
 * Formulas:
 *   A = Z11/Z21, B = Δ/Z21, C = 1/Z21, D = Z22/Z21
 *   where Δ = Z11·Z22 - Z12·Z21
 *
 * @return 0 on success, -1 if Z21 = 0
 * Complexity: O(1)
 */
int zp_to_abcd(const ZParameters *zp, ABCDParameters *abcd);

/**
 * @brief Convert ABCD to Z parameters
 * Complexity: O(1)
 */
int abcd_to_zp(const ABCDParameters *abcd, ZParameters *zp);

/**
 * @brief Convert S to Z parameters
 *
 * Z = Z0 * (I + S) * (I - S)^{-1}
 * For 2x2:
 *   Z11 = Z0 * ((1+S11)*(1-S22) + S12*S21) / Δs
 *   Z12 = Z0 * (2*S12) / Δs
 *   Z21 = Z0 * (2*S21) / Δs
 *   Z22 = Z0 * ((1-S11)*(1+S22) + S12*S21) / Δs
 *   where Δs = (1-S11)*(1-S22) - S12*S21
 *
 * Reference: Pozar (2012) §4.3
 * Complexity: O(1)
 */
int sp_to_zp(const SParameters *sp, ZParameters *zp);

/**
 * @brief Convert Z to S parameters
 *
 * S = (Z/Z0 - I) * (Z/Z0 + I)^{-1}
 * Complexity: O(1)
 */
int zp_to_sp(const ZParameters *zp, double z0, SParameters *sp);

/* L2: Network classification */

/**
 * @brief Check reciprocity: Z12 == Z21 for Z-parameters
 *
 * Theorem (Rayleigh 1873): A network composed only of passive elements
 *   (R, L, C) is reciprocal. Reciprocity means that the transfer
 *   function is the same in both directions: Z12 = Z21, Y12 = Y21.
 *   For S-parameters: S12 = S21.
 *
 * In Z-parameters: Z12 == Z21 within tolerance.
 * Complexity: O(1)
 */
int tp_is_reciprocal_z(const ZParameters *zp, double tolerance);

/**
 * @brief Check passivity: Re(Z) matrix must be positive semi-definite
 *
 * Theorem: A passive network cannot generate energy. For Z-parameters:
 *   (a) Z(jω) is analytic in Re(s) > 0
 *   (b) Z*(jω) = Z(-jω) (Hermitian)
 *   (c) Re[Z(jω)] = (Z + Z*)/2 >= 0 (positive semi-definite)
 *
 * For 2x2 at a given frequency: requires
 *   Re(Z11) >= 0, Re(Z22) >= 0
 *   4·Re(Z11)·Re(Z22) >= |Z12 + Z21*|^2
 *
 * Complexity: O(1)
 */
int tp_is_passive_z(const ZParameters *zp, double tolerance);

/* L5: Network interconnection */

/**
 * @brief Series-series interconnection: Z = Z_a + Z_b
 *
 * Two two-ports connected in series at both ports.
 * Validity condition: port currents must be equal (needs ideal transformer test).
 * Complexity: O(1)
 */
ZParameters tp_series_interconnect(const ZParameters *za, const ZParameters *zb);

/**
 * @brief Parallel-parallel interconnection: Y = Y_a + Y_b
 *
 * Two two-ports connected in parallel at both ports.
 * Validity condition: port voltages must be equal (Brune test).
 * Complexity: O(1)
 */
YParameters tp_parallel_interconnect(const YParameters *ya, const YParameters *yb);

/**
 * @brief Cascade interconnection: ABCD = ABCD_a * ABCD_b
 *
 * Theorem: For two networks in cascade (output of first → input of second),
 *   the overall ABCD matrix is the product of individual ABCD matrices:
 *   [A B; C D] = [A1 B1; C1 D1] * [A2 B2; C2 D2]
 *   = [A1*A2+B1*C2  A1*B2+B1*D2; C1*A2+D1*C2  C1*B2+D1*D2]
 *
 * This is the key advantage of ABCD parameters for cascaded systems.
 * Complexity: O(1)
 */
ABCDParameters tp_cascade_interconnect(const ABCDParameters *abcda,
                                         const ABCDParameters *abcdb);

/* L2: Common two-port topologies */

/**
 * @brief Z-parameters for a series impedance (shunt-series configuration)
 *
 * Circuit: single impedance Z in series with both ports
 *   Z = [Z  Z; Z  Z]
 *
 * Complexity: O(1)
 */
ZParameters tp_series_impedance(double complex z);

/**
 * @brief Y-parameters for a shunt admittance
 *
 * Circuit: single admittance Y from each port line to ground
 *   Y = [Y  -Y; -Y  Y]
 *
 * Complexity: O(1)
 */
YParameters tp_shunt_admittance(double complex y);

/**
 * @brief ABCD parameters for series impedance
 *   A=1, B=Z, C=0, D=1
 * Complexity: O(1)
 */
ABCDParameters tp_abcd_series_impedance(double complex z);

/**
 * @brief ABCD parameters for shunt admittance
 *   A=1, B=0, C=Y, D=1
 * Complexity: O(1)
 */
ABCDParameters tp_abcd_shunt_admittance(double complex y);

/**
 * @brief ABCD parameters for an ideal transformer (turns ratio n:1)
 *   A=n, B=0, C=0, D=1/n
 * Complexity: O(1)
 */
ABCDParameters tp_abcd_ideal_transformer(double n);

/**
 * @brief ABCD parameters for a π-network
 *
 * Circuit:
 *        Y2
 *   o----[ ]----o
 *   |          |
 *   Y1        Y3
 *   |          |
 *   o----------o
 *
 *   A = 1 + Y2/Y3, B = 1/Y3, C = Y1 + Y2 + Y1*Y2/Y3, D = 1 + Y2/Y1
 * Complexity: O(1)
 */
ABCDParameters tp_abcd_pi_network(double complex y1, double complex y2, double complex y3);

/**
 * @brief ABCD parameters for a T-network
 *
 * Circuit:
 *   o---Z1---o---Z2---o
 *            |
 *           Z3
 *            |
 *   o--------o--------o
 *
 *   A = 1 + Z1/Z3, B = Z1 + Z2 + Z1*Z2/Z3, C = 1/Z3, D = 1 + Z2/Z3
 * Complexity: O(1)
 */
ABCDParameters tp_abcd_t_network(double complex z1, double complex z2, double complex z3);

/* L5: S-parameter derived quantities (RF engineering) */

/**
 * @brief Input reflection coefficient with load Z_L
 *
 * Γ_in = S11 + S12·S21·Γ_L / (1 - S22·Γ_L)
 * where Γ_L = (Z_L - Z0)/(Z_L + Z0)
 *
 * Complexity: O(1)
 */
double complex sp_input_reflection(const SParameters *sp, double complex z_load);

/**
 * @brief Output reflection coefficient with source Z_S
 *
 * Γ_out = S22 + S12·S21·Γ_S / (1 - S11·Γ_S)
 * Complexity: O(1)
 */
double complex sp_output_reflection(const SParameters *sp, double complex z_source);

/**
 * @brief Transducer power gain: G_T = |S21|^2 * (1-|Γ_S|^2)*(1-|Γ_L|^2) / |1-Γ_S·Γ_in|^2 / |1-S22·Γ_L|^2
 *
 * Reference: Pozar (2012) §12.1
 * Complexity: O(1)
 */
double sp_transducer_gain(const SParameters *sp, double complex z_source, double complex z_load);

/**
 * @brief Maximum available gain (MAG) when unconditionally stable
 *
 * MAG = |S21/S12| * (K - sqrt(K^2 - 1)) where K = Rollett stability factor
 *
 * Reference: Rollett (1962) IRE Trans. CT-9
 * Complexity: O(1)
 */
double sp_max_available_gain(const SParameters *sp);

/**
 * @brief Rollett stability factor K
 *
 * K = (1 - |S11|^2 - |S22|^2 + |Δ|^2) / (2·|S12·S21|)
 * where Δ = S11·S22 - S12·S21
 * K > 1: unconditionally stable
 *
 * Reference: Rollett, IRE Trans. Circuit Theory (1962)
 * Complexity: O(1)
 */
double sp_rollett_stability(const SParameters *sp);

/**
 * @brief Voltage standing wave ratio: VSWR = (1 + |Γ|)/(1 - |Γ|)
 *
 * Complexity: O(1)
 */
double sp_vswr(double complex reflection_coef);

/* L8: Balanced two-port (differential mode) */

/**
 * @brief Convert single-ended S-parameters to mixed-mode S-parameters
 *
 * For differential circuit analysis, converts standard 2-port S-params
 * to mixed-mode representation with differential- and common-mode ports.
 *
 * S_mm = [S_dd  S_dc;  S_cd  S_cc] where each is 2x2
 *
 * Reference: Bockelman and Eisenstadt, IEEE T-MTT (1995)
 * Complexity: O(1)
 */
typedef struct {
    double complex sdd11, sdd12, sdd21, sdd22; /* Differential-mode */
    double complex sdc11, sdc12, sdc21, sdc22; /* Differential-to-common */
    double complex scd11, scd12, scd21, scd22; /* Common-to-differential */
    double complex scc11, scc12, scc21, scc22; /* Common-mode */
    double z0;
} MixedModeSParameters;

int sp_to_mixed_mode(const SParameters *sp, MixedModeSParameters *mm);

#endif /* TWO_PORT_NETWORK_H */
