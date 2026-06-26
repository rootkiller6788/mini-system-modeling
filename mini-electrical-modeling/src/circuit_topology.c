/**
 * @file circuit_topology.c
 * @brief Circuit topology analysis — KVL, KCL, nodal analysis, MNA, Gaussian elimination
 *
 * Independent knowledge points:
 *   - KVL verification (Kirchhoff 1845)
 *   - KCL verification (Kirchhoff 1845)
 *   - Nodal analysis matrix construction (Maxwell 1873 method)
 *   - MNA matrix construction (Ho, Ruehli, Brennan 1975)
 *   - Gaussian elimination with partial pivoting (Gauss 1810)
 *   - Thevenin/Norton equivalent computation
 *   - Maximum power transfer theorem (Jacobi 1840)
 *   - Superposition principle
 *   - Tellegen theorem (1952)
 *   - Matrix condition number estimation for circuit simulation
 */

#include "circuit_topology.h"
#include <math.h>
#include <string.h>
#include <float.h>

/* ===== KVL Verification ===== */

double kvl_verify(const double *voltages, const int *polarity, size_t n)
{
    /* Σ v_i * polarity_i = 0 for any closed loop */
    /* polarity: +1 if branch direction matches loop direction, -1 otherwise */
    if (!voltages || !polarity || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += voltages[i] * (double)polarity[i];
    }
    return sum;
}

/* ===== KCL Verification ===== */

double kcl_verify(const double *currents, const int *direction, size_t n)
{
    /* Σ i_i = 0 at any node */
    /* direction: +1 for entering, -1 for leaving (or vice versa, sign convention) */
    if (!currents || !direction || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += currents[i] * (double)direction[i];
    }
    return sum;
}

/* ===== Nodal Analysis: Build G and I matrices ===== */

int nodal_analysis_build(const Circuit *circuit, double *G, double *Ivec, size_t N)
{
    /* Build conductance matrix G and current source vector I */
    /* G[i][i] = sum of conductances at node i */
    /* G[i][j] = -sum of conductances between node i and node j (i != j) */
    /* Ivec[i] = net current source entering node i */

    if (!circuit || !G || !Ivec || N == 0) return -1;
    if (N > MAX_NODES) return -2;

    size_t NN = N * N;

    /* Initialize to zero */
    for (size_t k = 0; k < NN; k++) G[k] = 0.0;
    for (size_t k = 0; k < N;  k++) Ivec[k] = 0.0;

    for (size_t b = 0; b < circuit->num_branches; b++) {
        const CircuitBranch *br = &circuit->branches[b];
        int np = br->node_plus;
        int nm = br->node_minus;

        /* Nodes are indexed 1..N, 0 = ground. We store 0..N-1 internally. */
        int ip = np - 1;
        int im = nm - 1;

        switch (br->type) {
        case ELEM_RESISTOR: {
            /* Conductance g = 1/R */
            if (fabs(br->value) < 1e-30) continue; /* short circuit */
            double g = 1.0 / br->value;

            if (np != 0) { /* Not connected to ground */
                G[ip * N + ip] += g;
                if (nm != 0) G[ip * N + im] -= g;
            }
            if (nm != 0) {
                G[im * N + im] += g;
                if (np != 0) G[im * N + ip] -= g;
            }
            break;
        }
        case ELEM_CURRENT_SOURCE: {
            /* Current source I flows from node_plus to node_minus */
            if (np != 0) Ivec[ip] += br->value;
            if (nm != 0) Ivec[im] -= br->value;
            break;
        }
        case ELEM_VOLTAGE_SOURCE:
        case ELEM_VCVS:
        case ELEM_CCVS:
            /* Nodal analysis cannot handle voltage sources directly */
            /* Caller should use MNA */
            return -1;
        case ELEM_VCCS: {
            /* VCCS: I_out = gm * V_ctrl */
            /* gm = br->value, control voltage: V(ctrl_p) - V(ctrl_n) */
            int cp = br->ctrl_node_p - 1;
            int cn = br->ctrl_node_n - 1;
            double gm = br->value;

            if (np != 0) {
                if (br->ctrl_node_p != 0) G[ip * N + cp] += gm;
                if (br->ctrl_node_n != 0) G[ip * N + cn] -= gm;
            }
            if (nm != 0) {
                if (br->ctrl_node_p != 0) G[im * N + cp] -= gm;
                if (br->ctrl_node_n != 0) G[im * N + cn] += gm;
            }
            break;
        }
        default:
            /* Skip unknown element types */
            break;
        }
    }
    return 0;
}

/* ===== Modified Nodal Analysis (MNA) ===== */

int mna_build_matrix(const Circuit *circuit, CircuitMatrix *A, CircuitVector *bvec)
{
    /* MNA: [G  E] [v]   [i_s]
     *      [E' 0] [i_v] = [v_s]
     *
     * G: conductance matrix (N x N)
     * E: voltage source incidence (N x M)
     * v: node voltages (N x 1)
     * i_v: voltage source currents (M x 1)
     * i_s: independent current sources (N x 1)
     * v_s: independent voltage sources (M x 1)
     *
     * Total dimension: N + M
     */

    if (!circuit || !A || !bvec) return -1;

    size_t N = circuit->num_nodes;
    size_t M = circuit->num_voltage_sources;
    size_t dim = N + M;

    if (dim > MAX_NODES + MAX_BRANCHES) return -2;

    A->n = dim;
    A->data = (double*)calloc(dim * dim, sizeof(double));
    bvec->n = dim;
    bvec->data = (double*)calloc(dim, sizeof(double));

    if (!A->data || !bvec->data) {
        free(A->data);
        free(bvec->data);
        return -3;
    }

    /* Index tracking for voltage source unknowns */
    size_t vs_idx = 0;
    int *vs_branch_idx = (int*)malloc(circuit->num_branches * sizeof(int));
    if (!vs_branch_idx) {
        free(A->data); free(bvec->data);
        return -3;
    }
    for (size_t i = 0; i < circuit->num_branches; i++) vs_branch_idx[i] = -1;

    /* First pass: identify voltage source branches and assign MNA indices */
    for (size_t b = 0; b < circuit->num_branches; b++) {
        const CircuitBranch *br = &circuit->branches[b];
        if (br->type == ELEM_VOLTAGE_SOURCE || br->type == ELEM_VCVS || br->type == ELEM_CCVS) {
            vs_branch_idx[b] = (int)vs_idx++;
        }
    }

    /* Second pass: stamp element contributions */
    for (size_t b = 0; b < circuit->num_branches; b++) {
        const CircuitBranch *br = &circuit->branches[b];
        int np = br->node_plus;
        int nm = br->node_minus;

        switch (br->type) {
        case ELEM_RESISTOR: {
            if (fabs(br->value) < 1e-30) continue;
            double g = 1.0 / br->value;
            /* Stamp into G block */
            if (np != 0) {
                size_t ip = (size_t)(np - 1);
                A->data[ip * dim + ip] += g;
                if (nm != 0) {
                    size_t im = (size_t)(nm - 1);
                    A->data[ip * dim + im] -= g;
                    A->data[im * dim + ip] -= g;
                }
            }
            if (nm != 0) {
                size_t im = (size_t)(nm - 1);
                A->data[im * dim + im] += g;
            }
            break;
        }
        case ELEM_VOLTAGE_SOURCE: {
            int idx = vs_branch_idx[b];
            if (idx < 0) continue;
            size_t row = N + (size_t)idx;

            /* E matrix: +1 at node_plus, -1 at node_minus */
            if (np != 0) {
                size_t ip = (size_t)(np - 1);
                A->data[row * dim + ip] = 1.0;   /* E' block: lower-left */
                A->data[ip * dim + row] = 1.0;   /* E block: upper-right */
            }
            if (nm != 0) {
                size_t im = (size_t)(nm - 1);
                A->data[row * dim + im] -= 1.0;  /* E' block */
                A->data[im * dim + row] -= 1.0;  /* E block */
            }

            /* RHS: V_source */
            bvec->data[row] = br->value;
            break;
        }
        case ELEM_CURRENT_SOURCE: {
            if (np != 0) bvec->data[np - 1] += br->value;
            if (nm != 0) bvec->data[nm - 1] -= br->value;
            break;
        }
        case ELEM_VCCS: {
            /* gm * (V_ctrl_p - V_ctrl_n) stamped into G */
            int cp = br->ctrl_node_p;
            int cn = br->ctrl_node_n;
            double gm = br->value;

            if (np != 0 && cp != 0)
                A->data[(np-1) * dim + (cp-1)] += gm;
            if (np != 0 && cn != 0)
                A->data[(np-1) * dim + (cn-1)] -= gm;
            if (nm != 0 && cp != 0)
                A->data[(nm-1) * dim + (cp-1)] -= gm;
            if (nm != 0 && cn != 0)
                A->data[(nm-1) * dim + (cn-1)] += gm;
            break;
        }
        default:
            break;
        }
    }

    free(vs_branch_idx);
    return (int)dim;
}

/* ===== Gaussian Elimination with Partial Pivoting ===== */

int gaussian_elimination(double *A, double *x, const double *b, size_t n)
{
    /* Solve Ax = b using Gaussian elimination with partial pivoting */
    /* A is modified in-place to hold LU factors */
    /* Returns 0 on success, -1 on singular matrix */

    if (!A || !x || !b || n == 0) return -1;

    /* Create augmented matrix [A|b] as working copy */
    double *aug = (double*)malloc(n * (n + 1) * sizeof(double));
    if (!aug) return -1;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            aug[i * (n + 1) + j] = A[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }

    /* Forward elimination with partial pivoting */
    for (size_t k = 0; k < n; k++) {
        /* Find pivot row (maximum absolute value in column k) */
        size_t max_row = k;
        double max_val = fabs(aug[k * (n + 1) + k]);
        for (size_t i = k + 1; i < n; i++) {
            double val = fabs(aug[i * (n + 1) + k]);
            if (val > max_val) {
                max_val = val;
                max_row = i;
            }
        }

        /* Check for singular matrix */
        if (max_val < 1e-15) {
            free(aug);
            return -1; /* Singular */
        }

        /* Swap rows if needed */
        if (max_row != k) {
            for (size_t j = k; j <= n; j++) {
                double tmp = aug[k * (n + 1) + j];
                aug[k * (n + 1) + j] = aug[max_row * (n + 1) + j];
                aug[max_row * (n + 1) + j] = tmp;
            }
        }

        /* Store pivot */
        double pivot = aug[k * (n + 1) + k];

        /* Eliminate rows below */
        for (size_t i = k + 1; i < n; i++) {
            double factor = aug[i * (n + 1) + k] / pivot;
            aug[i * (n + 1) + k] = factor; /* Store L factor */
            for (size_t j = k + 1; j <= n; j++) {
                aug[i * (n + 1) + j] -= factor * aug[k * (n + 1) + j];
            }
        }
    }

    /* Back substitution */
    for (size_t i_int = n; i_int > 0; i_int--) {
        size_t i = i_int - 1;
        double sum = aug[i * (n + 1) + n];
        for (size_t j = i + 1; j < n; j++) {
            sum -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] = sum / aug[i * (n + 1) + i];
    }

    /* Copy LU factors back to A (upper triangular + implicit L) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            A[i * n + j] = aug[i * (n + 1) + j];
        }
    }

    free(aug);
    return 0;
}

/* ===== Thevenin Equivalent ===== */

TheveninEquiv thevenin_resistive_divider(double vs, double r1, double r2)
{
    /* For a divider: V_s through R1 to output, R2 from output to ground */
    /* V_th = V_s * R2/(R1+R2) */
    /* R_th = R1 || R2 */
    TheveninEquiv th;
    double sum = r1 + r2;
    if (fabs(sum) < 1e-30) {
        th.v_th = 0.0;
        th.r_th = 0.0;
        return th;
    }
    th.v_th = vs * r2 / sum;
    th.r_th = (r1 * r2) / sum;
    return th;
}

/* ===== Norton / Thevenin Duality ===== */

NortonEquiv norton_from_thevenin(const TheveninEquiv *th)
{
    NortonEquiv no;
    if (!th || fabs(th->r_th) < 1e-30) {
        no.i_n = DBL_MAX;
        no.r_n = 0.0;
        return no;
    }
    no.i_n = th->v_th / th->r_th;
    no.r_n = th->r_th;
    return no;
}

TheveninEquiv thevenin_from_norton(const NortonEquiv *no)
{
    TheveninEquiv th;
    if (!no) {
        th.v_th = 0.0;
        th.r_th = 0.0;
        return th;
    }
    th.v_th = no->i_n * no->r_n;
    th.r_th = no->r_n;
    return th;
}

/* ===== Maximum Power Transfer Theorem ===== */

double maximum_power_transfer(double v_th, double r_th)
{
    /* P_max = V_th^2 / (4 * R_th) when R_load = R_th */
    if (fabs(r_th) < 1e-30) return DBL_MAX;
    return (v_th * v_th) / (4.0 * r_th);
}

/* ===== Superposition Principle ===== */

double superposition_two_voltage_sources(
    double vs1, double coef1,
    double vs2, double coef2)
{
    /* Vx = a*Vs1 + b*Vs2 (coefficients from circuit topology) */
    return coef1 * vs1 + coef2 * vs2;
}

/* ===== Tellegen Theorem ===== */

double tellegen_power_sum(const double *voltages, const double *currents, size_t n)
{
    /* Σ v_k * i_k = 0 over all branches */
    /* This is a topological identity, not dependent on element types */
    if (!voltages || !currents || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t k = 0; k < n; k++) {
        sum += voltages[k] * currents[k];
    }
    return sum;
}

/* ===== Matrix Condition Number Estimate ===== */

double circuit_matrix_condition_estimate(const CircuitMatrix *A)
{
    /* Estimate log10 of condition number using 1-norm */
    /* For ill-conditioned matrices in circuit simulation:
     *   - Ratio of largest to smallest time constant too large
     *   - Nearly disconnected nodes
     *   - Extreme impedance ratios (>1e12)
     */
    if (!A || !A->data || A->n == 0) return 0.0;

    size_t n = A->n;
    double norm_a = 0.0;
    double norm_inv_est = 0.0;

    /* 1-norm of A: max column sum */
    for (size_t j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (size_t i = 0; i < n; i++) {
            col_sum += fabs(A->data[i * n + j]);
        }
        if (col_sum > norm_a) norm_a = col_sum;
    }

    /* Quick estimate of ||A^{-1}|| using power method on A^T A */
    /* For circuit matrices, diagonal dominance gives reasonable estimates */
    double min_pivot = DBL_MAX;
    for (size_t i = 0; i < n; i++) {
        double diag = fabs(A->data[i * n + i]);
        double off_diag = 0.0;
        for (size_t j = 0; j < n; j++) {
            if (i != j) off_diag += fabs(A->data[i * n + j]);
        }
        /* Gershgorin radius estimate: ||A^{-1}||_inf >= 1/min(|a_ii| - R_i) */
        double gersh = diag - off_diag;
        if (gersh > 0.0 && gersh < min_pivot) min_pivot = gersh;
    }

    if (min_pivot < 1e-30 || min_pivot > DBL_MAX / 2.0) norm_inv_est = 1e15;
    else norm_inv_est = 1.0 / min_pivot;

    double cond = norm_a * norm_inv_est;
    if (cond < 1.0) cond = 1.0;

    return log10(cond);
}
