/**
 * sfg_state_space.h — State-Space ↔ Signal Flow Graph Conversion
 *
 * State-space representations and signal flow graphs are dual
 * representations of linear dynamical systems. This module provides
 * bidirectional conversion between the two forms, enabling SFG-based
 * analysis of state-space systems and vice versa.
 *
 * L4 Conservation: The conversion preserves all system properties
 * (poles, zeros, controllability, observability) because it is an
 * isomorphic mapping between equivalent algebraic representations.
 *
 * State-Space Form:
 *   ẋ(t) = A x(t) + B u(t)
 *   y(t) = C x(t) + D u(t)
 *
 * SFG Representation:
 *   - Each state variable x_i → integrator (1/s) node pair
 *   - Each input u_j → source node
 *   - Each output y_k → sink node
 *   - Branches from x_j to x_i with gain A[i][j]
 *   - Branches from u_j to x_i with gain B[i][j]
 *   - Branches from x_j to y_i with gain C[i][j]
 *   - Branches from u_j to y_i with gain D[i][j]
 *
 * References:
 *  - Ogata, K. "Modern Control Engineering" Prentice-Hall, 5th Ed.,
 *    2010. Chapter 9 (Control Systems Analysis in State Space).
 *  - Kuo, B.C. "Automatic Control Systems" Prentice-Hall, 9th Ed.,
 *    2014. Chapter 5 (State-Variable Analysis).
 *  - D'Azzo & Houpis, "Linear Control System Analysis and Design",
 *    1995. Chapter 11.
 */

#ifndef SFG_STATE_SPACE_H
#define SFG_STATE_SPACE_H

#include "sfg_core.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SFG_SS_MAX_ORDER  64  /* Maximum state-space order */

/* ================================================================
 * L1: State-Space Data Types
 * ================================================================ */

/**
 * sfg_ss_system_t — Linear time-invariant state-space system
 *
 * L1 Definition: A state-space system (A, B, C, D) with:
 *   - n: number of state variables (system order)
 *   - m: number of inputs
 *   - p: number of outputs
 */
typedef struct {
    int         n;  /* Number of states (system order) */
    int         m;  /* Number of inputs */
    int         p;  /* Number of outputs */

    sfg_real_t* A;  /* System matrix (n×n, row-major) */
    sfg_real_t* B;  /* Input matrix (n×m, row-major) */
    sfg_real_t* C;  /* Output matrix (p×n, row-major) */
    sfg_real_t* D;  /* Feedthrough matrix (p×m, row-major) */

    /* Derived quantities (computed on demand) */
    int         A_alloc;  /* 1 if A was internally allocated */
    int         B_alloc;
    int         C_alloc;
    int         D_alloc;

    char        name[64]; /* System name */
} sfg_ss_system_t;

/* ================================================================
 * L5: State-Space → SFG Conversion
 * ================================================================ */

/**
 * sfg_from_state_space — Convert state-space system to signal flow graph
 *
 * L5 Method: Creates an SFG with:
 *   - m source nodes (inputs u_0 ... u_{m-1})
 *   - n integrator node pairs (state x_i and its derivative ẋ_i)
 *   - p sink nodes (outputs y_0 ... y_{p-1})
 *   - Branches encoding A, B, C, D matrices
 *
 * The integrator 1/s is represented as a branch with gain 1.0
 * from ẋ_i to x_i. This is equivalent to the Laplace-domain
 * relationship: x_i = (1/s) · ẋ_i.
 *
 * Complexity: O(n² + n·m + p·n + p·m) = O(n²) dominant
 *
 * @param ss       Source state-space system
 * @param g        Output SFG (pre-initialized by caller)
 * @return         0 on success, -1 if system too large
 */
int sfg_from_state_space(const sfg_ss_system_t* ss, sfg_graph_t* g);

/**
 * sfg_from_state_space_controllable — Convert to SFG in controllable
 * canonical form
 *
 * L5 Method: For SISO systems, the controllable canonical form places
 * the system matrix A in companion form, with the characteristic
 * polynomial coefficients appearing directly in the top row.
 *
 * The resulting SFG has a chain structure with feedback from each
 * state to the first state's input summing junction.
 *
 * Complexity: O(n)
 */
int sfg_from_state_space_controllable(const sfg_ss_system_t* ss,
                                       sfg_graph_t* g);

/**
 * sfg_from_state_space_observable — Convert to SFG in observable
 * canonical form
 *
 * L5 Method: The observable canonical form is the dual of the
 * controllable form. The SFG has a chain structure with feedforward
 * from each state.
 *
 * Complexity: O(n)
 */
int sfg_from_state_space_observable(const sfg_ss_system_t* ss,
                                     sfg_graph_t* g);

/**
 * sfg_from_transfer_function — Construct SFG from transfer function
 * coefficients (SISO)
 *
 * L5 Method: Given transfer function
 *   T(s) = (b_0 + b_1 s + ... + b_m s^m) / (a_0 + a_1 s + ... + a_n s^n)
 * construct the SFG in controllable canonical form.
 *
 * Complexity: O(n)
 */
int sfg_from_transfer_function(sfg_graph_t* g,
                                const sfg_real_t* num, int num_order,
                                const sfg_real_t* den, int den_order);

/* ================================================================
 * L5: SFG → State-Space Conversion
 * ================================================================ */

/**
 * sfg_to_state_space — Extract state-space representation from SFG
 *
 * L5 Method: Identify integrator (1/s) nodes as state variables and
 * extract the A, B, C, D matrices from the branch gains between
 * state nodes, input nodes, and output nodes.
 *
 * The caller must identify which nodes correspond to:
 *   - state variables (integrator outputs)
 *   - inputs (sources)
 *   - outputs (sinks or state combinations)
 *
 * Complexity: O(n² · B) where n = number of state nodes
 *
 * @param g           Source SFG
 * @param state_ids   Array of n state node IDs
 * @param n_states    Number of state variables
 * @param input_ids   Array of m input node IDs
 * @param n_inputs    Number of inputs
 * @param output_ids  Array of p output node IDs
 * @param n_outputs   Number of outputs
 * @param ss          Output state-space system (caller must free)
 * @return            0 on success, -1 on error
 */
int sfg_to_state_space(const sfg_graph_t* g,
                        const sfg_node_id_t* state_ids, int n_states,
                        const sfg_node_id_t* input_ids, int n_inputs,
                        const sfg_node_id_t* output_ids, int n_outputs,
                        sfg_ss_system_t* ss);

/* ================================================================
 * L5: State-Space via SFG Analysis
 * ================================================================ */

/**
 * sfg_ss_transfer_function — Compute transfer function via SFG reduction
 *
 * Alternative to standard (sI - A)^(-1) formula. Uses SFG reduction
 * to compute the transfer function matrix.
 *
 * Complexity: O(n³) for standard, O(n²) for sparse SFG
 */
int sfg_ss_transfer_function(const sfg_ss_system_t* ss,
                              sfg_gain_t* T, int* T_rows, int* T_cols);

/**
 * sfg_ss_controllability — Check controllability via SFG path analysis
 *
 * L6: A system is controllable iff there exists a path (possibly through
 * integrator nodes) from every input to every state variable.
 *
 * This is the SFG-based controllability test, complementing the
 * standard Kalman rank test.
 *
 * Returns: 1 if controllable, 0 if not
 */
int sfg_ss_controllability(const sfg_ss_system_t* ss);

/**
 * sfg_ss_observability — Check observability via SFG path analysis
 *
 * L6: A system is observable iff there exists a path from every state
 * variable to every output.
 *
 * Returns: 1 if observable, 0 if not
 */
int sfg_ss_observability(const sfg_ss_system_t* ss);

/* ================================================================
 * L3: State-Space Utility Functions
 * ================================================================ */

/**
 * sfg_ss_init — Initialize a state-space system structure
 */
void sfg_ss_init(sfg_ss_system_t* ss);

/**
 * sfg_ss_alloc — Allocate matrices for a state-space system
 *
 * Returns: 0 on success, -1 on allocation failure
 */
int sfg_ss_alloc(sfg_ss_system_t* ss, int n, int m, int p);

/**
 * sfg_ss_free — Free internal matrices of state-space system
 */
void sfg_ss_free(sfg_ss_system_t* ss);

/**
 * sfg_ss_set — Set a specific matrix element (row, col)
 */
void sfg_ss_set_A(sfg_ss_system_t* ss, int row, int col, sfg_real_t val);
void sfg_ss_set_B(sfg_ss_system_t* ss, int row, int col, sfg_real_t val);
void sfg_ss_set_C(sfg_ss_system_t* ss, int row, int col, sfg_real_t val);
void sfg_ss_set_D(sfg_ss_system_t* ss, int row, int col, sfg_real_t val);

/**
 * sfg_ss_get — Get a specific matrix element (row, col)
 */
sfg_real_t sfg_ss_get_A(const sfg_ss_system_t* ss, int row, int col);
sfg_real_t sfg_ss_get_B(const sfg_ss_system_t* ss, int row, int col);
sfg_real_t sfg_ss_get_C(const sfg_ss_system_t* ss, int row, int col);
sfg_real_t sfg_ss_get_D(const sfg_ss_system_t* ss, int row, int col);

/**
 * sfg_ss_print — Print state-space matrices
 */
void sfg_ss_print(const sfg_ss_system_t* ss);

/**
 * sfg_ss_characteristic_polynomial — Compute coefficients of |sI - A|
 *
 * L4: The characteristic polynomial determines stability.
 * Computed from the SFG by Mason's formula: Δ(s) = |sI - A|.
 *
 * @param poly  Output: coefficients [a_n ... a_0] of s^n + ... + a_0
 * @param n     Order of system (= size of poly - 1)
 */
int sfg_ss_characteristic_polynomial(const sfg_ss_system_t* ss,
                                      sfg_real_t* poly, int n);

#ifdef __cplusplus
}
#endif

#endif /* SFG_STATE_SPACE_H */
