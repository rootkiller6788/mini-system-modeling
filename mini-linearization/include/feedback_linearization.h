/**
 * @file feedback_linearization.h
 * @brief Feedback linearization
 */
#ifndef FEEDBACK_LINEARIZATION_H
#define FEEDBACK_LINEARIZATION_H
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "linearization_core.h"

typedef struct { double value; double *gradient; int n_states; } LieDerivative;
typedef struct { double *value; int n_states; } LieBracket;

typedef struct {
    int relative_degree;
    bool is_well_defined;
    double *lie_derivatives;
    double *input_coeffs;
} RelativeDegree;

typedef struct {
    double *z; double *Phi; double *Phi_inv;
    int n_states; bool is_diffeomorphism;
} Diffeomorphism;

typedef struct {
    double *z; double *eta; double *a_coeff; double *b_coeff;
    int relative_degree; int internal_dim;
} NormalForm;

void compute_lie_derivative(
    void (*f)(const double *x, double *fx, int n, void *data),
    double (*h)(const double *x, int n, void *data),
    const double *x, int n, LieDerivative *result, void *user_data);

void compute_lie_bracket(
    void (*f)(const double *x, double *fx, int n, void *data),
    void (*g)(const double *x, double *gx, int n, void *data),
    const double *x, int n, LieBracket *result, void *user_data);

void compute_iterated_lie_bracket(
    void (*f)(const double *x, double *fx, int n, void *data),
    void (*g)(const double *x, double *gx, int n, void *data),
    const double *x, int n, int k, double *result, void *user_data);

RelativeDegree *compute_relative_degree(
    void (*f)(const double *x, double *fx, int n, void *data),
    void (*g)(const double *x, double *gx, int n, void *data),
    double (*h)(const double *x, int n, void *data),
    const double *x0, int n, int max_r, void *user_data);

void free_relative_degree(RelativeDegree *rd);

bool brockett_necessary_condition(
    void (*f_map)(const double *x, const double *u, double *fx, int n, int m, void *data),
    int n, int m, void *user_data);

bool check_frobenius_involutivity(
    void (*distribution)(const double *x, double *v, int n, int d, void *data),
    int n, int d, const double *x, void *user_data);

double input_output_linearization_controller(
    void (*f)(const double *x, double *fx, int n, void *data),
    void (*g)(const double *x, double *gx, int n, void *data),
    double (*h)(const double *x, int n, void *data),
    const double *x, double v, int n, int r, void *user_data);

double input_state_linearization_controller(
    void (*f)(const double *x, double *fx, int n, void *data),
    void (*g)(const double *x, double *gx, int n, void *data),
    const double *x, double v, int n,
    double *A_c, double *B_c, void *user_data);

void transform_to_normal_form(
    void (*f)(const double *x, double *fx, int n, void *data),
    void (*g)(const double *x, double *gx, int n, void *data),
    double (*h)(const double *x, int n, void *data),
    const double *x, int n, int r, NormalForm *nf, void *user_data);

double analyze_zero_dynamics(
    const NormalForm *nf, const double *eta0,
    double t_end, int n_steps, double *eta_out);

 /* FEEDBACK_LINEARIZATION_H */

/**
 * L2: Multiple-input multiple-output (MIMO) extension definitions.
 */

/**
 * L3: Vector relative degree for MIMO systems.
 *
 * For a system with m inputs and m outputs, the vector relative degree
 * (r_1, ..., r_m) defines the number of differentiations of each output
 * until at least one input appears.
 *
 * The decoupling matrix D(x) is the m x m matrix with entries:
 *   D_{ij}(x) = L_{g_j} L_f^{r_i-1} h_i(x)
 *
 * The system has a well-defined vector relative degree if D(x) is
 * nonsingular in a neighborhood of the operating point.
 *
 * Reference: Isidori (1995), section 5.1
 */
typedef struct {
    int *relative_degrees;  /**< (r_1, ..., r_m) for each output */
    double *decoupling_matrix; /**< D(x), m x m, row-major */
    int m_outputs;          /**< Number of outputs/inputs */
    bool is_nonsingular;    /**< True if det(D) != 0 */
} MIMORelativeDegree;

/**
 * L3: Compute MIMO vector relative degree and decoupling matrix.
 *
 * Complexity: O(m * sum(r_i) * cost(Lie derivatives))
 * Reference: Isidori (1995), section 5.2
 */
MIMORelativeDegree *compute_mimo_relative_degree(
    void (*f)(const double *x, double *fx, int n, void *data),
    void (*g)(const double *x, double *Gx, int n, int m, void *data),
    void (*h)(const double *x, double *hx, int m, int n, void *data),
    const double *x0, int n, int m, int max_r, void *user_data);

/**
 * Free MIMORelativeDegree. O(1).
 */
void free_mimo_relative_degree(MIMORelativeDegree *rd);

/**
 * L5: MIMO input-output linearization controller.
 *
 * For a square MIMO system with nonsingular decoupling matrix D(x):
 *   u = D(x)^{-1} * (v - b(x))
 *
 * where b_i(x) = L_f^{r_i} h_i(x)
 *
 * @param f         Drift vector field
 * @param g         Input matrix field (n x m)
 * @param h         Output function
 * @param x         Current state
 * @param v         Synthetic control vector (length m)
 * @param n         State dimension
 * @param m         Input/output dimension
 * @param rd        Pre-computed MIMO relative degree
 * @param u_out     Output physical control (length m)
 * @param user_data Arbitrary data
 *
 * Complexity: O(m * max(r_i) * cost(Lie) + m^3)
 * Reference: Isidori (1995), section 5.3
 */
void mimo_io_linearization_controller(
    void (*f)(const double *x, double *fx, int n, void *data),
    void (*g)(const double *x, double *Gx, int n, int m, void *data),
    void (*h)(const double *x, double *hx, int m_out, int n, void *data),
    const double *x, const double *v, int n, int m,
    const MIMORelativeDegree *rd, double *u_out, void *user_data);

/**
 * L8: Adaptive feedback linearization with parameter uncertainty.
 *
 * When system parameters are unknown, combine feedback linearization
 * with parameter estimation (certainty equivalence).
 *
 * Complexity: O(max_r * (cost + param_update))
 * Reference: Sastry & Isidori (1989), IEEE TAC
 */
double adaptive_feedback_linearization(
    void (*f)(const double *x, double *fx, const double *theta, int n, void *data),
    void (*g)(const double *x, double *gx, int n, void *data),
    double (*h)(const double *x, int n, void *data),
    const double *x, double v, int n, int r,
    double *theta_hat, double adapt_gain, void *user_data);

/**
 * L9: Approximate feedback linearization via neural networks.
 *
 * Uses a neural network to learn the feedback linearizing transformation
 * when the analytical model is unavailable.
 *
 * Reference: He et al. (2017), IEEE TNNLS
 */
typedef struct {
    double *weights_input_hidden;  /**< W1 (n x n_hidden) */
    double *bias_hidden;           /**< b1 (n_hidden) */
    double *weights_hidden_output; /**< W2 (n_hidden x 1) */
    int n_input, n_hidden;
} NeuralFeedbackLinearizer;

/**
 * L8: Robust feedback linearization with sliding mode component.
 *
 * Combines nominal feedback linearization with a sliding mode term
 * to handle bounded uncertainties:
 *   u = u_nominal - k * sign(s)
 *
 * where s is the sliding surface.
 *
 * @param u_nominal  Nominal feedback linearization control
 * @param s          Sliding surface value
 * @param k          Sliding gain (must exceed uncertainty bound)
 * @return           Robust control input
 *
 * Reference: Slotine & Li (1991), section 7.3
 */
double robust_feedback_linearization_control(
    double u_nominal, double s, double k, double boundary_layer);




#endif /* FEEDBACK_LINEARIZATION_H */
