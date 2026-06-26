/**
 * @file circuit_topology.h
 * @brief Circuit topology analysis KVL, KCL, nodal analysis, mesh analysis, MNA
 *
 * Knowledge Coverage:
 *   L1 - Definitions: node, branch, mesh, loop, cut-set, incidence matrix,
 *        conductance matrix, node-voltage method
 *   L2 - Core Concepts: KVL (Kirchhoff Voltage Law), KCL (Kirchhoff Current Law),
 *        superposition, Thevenin/Norton equivalents, maximum power transfer
 *   L3 - Mathematical Structures: sparse matrices, LU decomposition for MNA,
 *        graph theory for circuit topology
 *   L4 - Fundamental Laws: KVL (sum of voltages around loop = 0),
 *        KCL (sum of currents at node = 0), Tellegen theorem
 *   L5 - Engineering Methods: nodal analysis, mesh analysis,
 *        modified nodal analysis (MNA), Gaussian elimination
 *
 * Reference:
 *   MIT 6.302 Feedback Systems - SPICE internals
 *   Berkeley EE 219 - Circuit Simulation (Nagel, SPICE2 author)
 *   Chua, Desoer, and Kuh, Linear and Nonlinear Circuits (1987)
 *   McCalla, Fundamentals of Computer-Aided Circuit Simulation (1988)
 */

#ifndef CIRCUIT_TOPOLOGY_H
#define CIRCUIT_TOPOLOGY_H

#include <stddef.h>
#include <stdlib.h>
#include <math.h>

/* Maximum circuit size for static allocation methods */
#define MAX_NODES 64
#define MAX_BRANCHES 128

/* L1: Element type enumeration */
typedef enum {
    ELEM_RESISTOR = 0,
    ELEM_CAPACITOR = 1,
    ELEM_INDUCTOR = 2,
    ELEM_VOLTAGE_SOURCE = 3,
    ELEM_CURRENT_SOURCE = 4,
    ELEM_VCVS = 5,   /* Voltage-controlled voltage source */
    ELEM_VCCS = 6,   /* Voltage-controlled current source */
    ELEM_CCVS = 7,   /* Current-controlled voltage source */
    ELEM_CCCS = 8,   /* Current-controlled current source */
    ELEM_OPAMP = 9    /* Ideal operational amplifier */
} ElementType;

/* L1: Circuit branch descriptor */
typedef struct {
    ElementType type;
    int node_plus;    /* Positive terminal node (0 = ground) */
    int node_minus;   /* Negative terminal node (0 = ground) */
    double value;     /* R[ohm], C[F], L[H], gain[V/V or A/A] */
    double value2;    /* Secondary value used for transfer (gm, beta, etc.) */
    int ctrl_node_p;  /* Controlling node + for dependent sources */
    int ctrl_node_n;  /* Controlling node - for dependent sources */
    int ctrl_branch;  /* Controlling branch for CCVS/CCCS */
} CircuitBranch;

/* L1: Circuit descriptor */
typedef struct {
    CircuitBranch *branches;
    size_t num_branches;
    size_t capacity_branches;
    size_t num_nodes;       /* Excluding ground (node 0) */
    size_t num_voltage_sources; /* Number of independent + controlled V-sources */
} Circuit;

/* L2: Circuit matrix (for nodal/MNA) dense representation */
typedef struct {
    double *data;
    size_t n;          /* Matrix dimension NxN */
    size_t capacity;
} CircuitMatrix;

/* L2: Solution vector */
typedef struct {
    double *data;
    size_t n;
} CircuitVector;

/* KVL and KCL verification */

/**
 * @brief Verify KVL for a loop sum of voltages around closed loop = 0
 *
 * Theorem (Kirchhoff 1845): In any closed loop of a circuit, the algebraic sum
 *   of branch voltages is zero: sum_i v_i = 0.
 *   This follows from conservation of energy the voltage is a potential function.
 *
 * @param voltages Array of branch voltages in the loop [V]
 * @param polarity Signs (+1 or -1) indicating reference direction
 * @param n Number of branches in the loop
 * @return Sum of signed voltages should be 0 for valid KVL
 * Complexity: O(n)
 */
double kvl_verify(const double *voltages, const int *polarity, size_t n);

/**
 * @brief Verify KCL for a node sum of currents entering = sum leaving
 *
 * Theorem (Kirchhoff 1845): At any node, the algebraic sum of currents is zero:
 *   sum_i i_i = 0. Currents entering = currents leaving.
 *   This follows from conservation of charge (charge cannot accumulate at a node).
 *
 * @param currents Array of branch currents at the node [A]
 * @param direction +1 for entering, -1 for leaving
 * @param n Number of branches connected to the node
 * @return Sum of signed currents should be 0 for valid KCL
 * Complexity: O(n)
 */
double kcl_verify(const double *currents, const int *direction, size_t n);

/* L5: Nodal analysis */

/**
 * @brief Build conductance matrix G for nodal analysis (DC resistive circuits)
 *
 * Theorem: Nodal analysis writes KCL at each node, expressing branch currents
 *   in terms of node voltages using Ohms law.
 *   G * V = I where:
 *     G[i][i] = sum of conductances connected to node i
 *     G[i][j] = -sum of conductances between node i and node j (i != j)
 *     I[i] = sum of current sources entering node i
 *
 * The circuit must not contain voltage sources (use MNA instead).
 *
 * @param circuit Circuit descriptor
 * @param G Output: NxN conductance matrix (caller allocates N*N doubles)
 * @param I Output: N current source vector (caller allocates N doubles)
 * @param N Number of nodes (excluding ground)
 * @return 0 on success, -1 if circuit contains voltage sources
 * Complexity: O(B + N^2) where B = branches, N = nodes
 */
int nodal_analysis_build(const Circuit *circuit, double *G, double *Ivec, size_t N);

/* L5: Modified Nodal Analysis (MNA) */

/**
 * @brief Build MNA matrix for circuits with voltage sources
 *
 * Theorem (Ho, Ruehli, Brennan 1975): MNA extends nodal analysis by adding
 *   unknown currents through voltage sources as extra variables.
 *   [G  E] [v]   [i_s]
 *   [E' 0] [i_v] = [v_s]
 *
 * where G = conductance matrix, E = voltage source incidence matrix,
 * v = node voltages, i_v = currents through voltage sources,
 * i_s = independent current sources, v_s = independent voltage sources.
 *
 * @return Size of MNA matrix on success, -1 on error
 * Complexity: O(B + (N+M)^2) where M = number of voltage sources
 */
int mna_build_matrix(const Circuit *circuit, CircuitMatrix *A, CircuitVector *b);

/* L5: Gaussian elimination for circuit solution */

/**
 * @brief Solve linear system Ax = b using Gaussian elimination with partial pivoting
 *
 * Theorem: Gaussian elimination (Gauss 1810) with partial pivoting provides
 *   O(n^3) solution for dense linear systems. For circuit matrices which
 *   are typically sparse and symmetric positive definite, specialized
 *   methods exist (Cholesky, LU with Markowitz ordering).
 *
 * @param A Input: NxN matrix; Output: destroyed (contains LU factors)
 * @param x Output: Solution vector
 * @param b Input: Right-hand side vector
 * @param n Dimension
 * @return 0 on success, -1 if singular matrix detected
 * Complexity: O(n^3)
 */
int gaussian_elimination(double *A, double *x, const double *b, size_t n);

/* L2: Thevenin and Norton equivalent circuits */

/**
 * @brief Compute Thevenin equivalent: V_th = V_open_circuit, R_th = V_th / I_short
 *
 * Theorem (Thevenin 1883, Helmholtz 1853): Any linear two-terminal network
 *   of voltage sources, current sources, and resistances is equivalent to
 *   a single voltage source V_th in series with a single resistance R_th.
 *
 * For a simple resistive divider from V_s through R1 to output, with R2 to ground:
 *   V_th = V_s * R2/(R1 + R2)
 *   R_th = R1 || R2 = (R1*R2)/(R1+R2)
 *
 * Complexity: O(1)
 */
typedef struct {
    double v_th;   /* Thevenin equivalent voltage [V] */
    double r_th;   /* Thevenin equivalent resistance [ohm] */
} TheveninEquiv;

TheveninEquiv thevenin_resistive_divider(double vs, double r1, double r2);

/**
 * @brief Norton equivalent: I_n = I_short_circuit, R_n = R_th
 *
 * Theorem (Norton 1926, Mayer 1926): Any linear two-terminal network is
 *   equivalent to a single current source I_n in parallel with R_n.
 *   Norton and Thevenin are duals: I_n = V_th/R_th, R_n = R_th.
 *
 * Complexity: O(1)
 */
typedef struct {
    double i_n;   /* Norton equivalent current [A] */
    double r_n;   /* Norton equivalent resistance [ohm] */
} NortonEquiv;

NortonEquiv norton_from_thevenin(const TheveninEquiv *th);
TheveninEquiv thevenin_from_norton(const NortonEquiv *no);

/* L2: Maximum power transfer */

/**
 * @brief Maximum power to load occurs when R_load = R_th
 *
 * Theorem (Jacobi 1840): For a Thevenin source with fixed V_th and R_th,
 *   maximum power is delivered to a resistive load when R_load = R_th.
 *   P_max = V_th^2 / (4 * R_th)
 *
 * Complexity: O(1)
 */
double maximum_power_transfer(double v_th, double r_th);

/* L2: Superposition principle */

/**
 * @brief Superposition for two-source DC circuit
 *
 * Theorem: In a linear circuit with multiple independent sources, the response
 *   (voltage/current at any point) equals the sum of responses due to
 *   each source acting alone, with all other independent sources deactivated
 *   (voltage sources shorted, current sources opened).
 *
 * For a simple case: Vx = a*Vs1 + b*Vs2 + c*Is1
 * Complexity: O(1)
 */
double superposition_two_voltage_sources(
    double vs1, double coef1,
    double vs2, double coef2);

/* Tellegen theorem verification */

/**
 * @brief Tellegen theorem: sum of v_k * i_k over all branches = 0
 *
 * Theorem (Tellegen 1952): For any lumped circuit satisfying KVL and KCL,
 *   sum_{all branches k} v_k * i_k = 0.
 *   This is a topological theorem independent of branch constitutive relations.
 *   Special case: conservation of instantaneous power.
 *
 * Complexity: O(n)
 */
double tellegen_power_sum(const double *voltages, const double *currents, size_t n);

/* L5: Numerical accuracy and conditioning for circuit matrices */

/**
 * @brief Estimate condition number via 1-norm approximation
 *
 * Theorem: For circuits with very small/large values (ratio > 1e12),
 *   the MNA matrix becomes ill-conditioned. Condition number k(A) = ||A||·||A^{-1}||.
 *   Rule of thumb: log10(k) digits of accuracy lost.
 *
 * @return Approximate log10 of condition number
 * Complexity: O(n^2)
 */
double circuit_matrix_condition_estimate(const CircuitMatrix *A);

#endif /* CIRCUIT_TOPOLOGY_H */
