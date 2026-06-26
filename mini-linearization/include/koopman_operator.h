/**
 * @file koopman_operator.h
 * @brief Koopman operator theory for linearization of nonlinear dynamics.
 */
#ifndef KOOPMAN_OPERATOR_H
#define KOOPMAN_OPERATOR_H
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>

typedef struct {
    double complex value;
    double *state;
    int n_states;
    int observable_id;
} KoopmanObservable;

typedef struct {
    double complex eigenvalue;
    double complex *eigenfunction;
    int n_samples;
} KoopmanEigenfunction;

typedef struct {
    double complex *mode;
    double complex eigenvalue;
    double amplitude;
    double frequency_hz;
    double growth_rate;
} KoopmanMode;

typedef struct {
    double complex *eigenvalues;
    double complex *eigenvectors;
    double complex *modes;
    double *singular_values;
    int n_state_dim, n_snapshots, r;
    double reconstruction_error;
} DMDResult;

typedef enum {
    EDMD_BASIS_MONOMIAL,
    EDMD_BASIS_THIN_PLATE,
    EDMD_BASIS_GAUSSIAN_RBF,
    EDMD_BASIS_FOURIER,
    EDMD_BASIS_CUSTOM
} EDMDBasisType;

typedef struct {
    EDMDBasisType basis_type;
    int n_basis_functions;
    int max_poly_order;
    double rbf_sigma;
    void *basis_params;
} EDMDConfig;

typedef struct {
    KoopmanMode *modes;
    int n_modes;
    double complex *amplitudes;
    double prediction_error;
} KoopmanDecomposition;

void compute_dmd(const double *X, const double *Y, int n, int m, int r, double dt, DMDResult *result);
void dmd_predict(const DMDResult *dmd, const double *x0, double t, double *x_pred);
void dmd_mode_amplitudes(const DMDResult *dmd, const double *x0, double complex *amplitudes);
void compute_exact_dmd(const double *X, const double *Y, int n, int m, int r, double dt, DMDResult *result);
void compute_edmd(const double *X, const double *Y, int n, int m, const EDMDConfig *config, DMDResult *result);
void evaluate_basis_functions(const double *x, int n, const EDMDConfig *config, double *psi);
void compute_kernel_edmd(const double *X, const double *Y, int n, int m, double (*kernel)(const double *x, const double *y, int n, double gamma), double gamma, DMDResult *result);
void identify_koopman_eigenfunctions(const DMDResult *dmd, const EDMDConfig *config, const double *X, int n, int m, KoopmanEigenfunction *eigfuncs);
void koopman_mode_decomposition(const DMDResult *dmd, const double *x0, int n, KoopmanDecomposition *decomp);
void koopman_predict(const KoopmanDecomposition *decomp, double t, int n, double *x_pred);
int carleman_embedding_dimension(int n, int max_k);
void carleman_embedding_matrix(void (*f_poly)(const double *x, double *dx, int n, void *data), int n, int poly_deg, int trunc_k, double *A, int N, void *user_data);
void free_dmd_result(DMDResult *result);
void free_koopman_decomposition(KoopmanDecomposition *decomp);
void free_koopman_eigenfunctions(KoopmanEigenfunction *eigfuncs);

 /* KOOPMAN_OPERATOR_H */

/**
 * L9: Deep Koopman -- neural network-based observable discovery.
 *
 * Uses autoencoder architecture to learn optimal observables
 * for Koopman operator approximation from data.
 *
 * Reference: Lusch, Kutz, Brunton (2018), Nature Communications
 */
typedef struct {
    double *encoder_weights;   /**< Encoding layer weights */
    double *decoder_weights;   /**< Decoding layer weights */
    double *koopman_matrix;    /**< Learned Koopman operator */
    int state_dim, latent_dim, hidden_dim;
    double training_loss;
} DeepKoopmanModel;

/**
 * L9: Train a deep Koopman model from trajectory data.
 *
 * Minimizes: ||x - decode(encode(x))||^2 + ||encode(y) - K*encode(x)||^2
 *
 * @param X, Y        Snapshot pairs
 * @param n, m        Dimensions
 * @param latent_dim  Dimension of Koopman-invariant subspace
 * @param model       Output trained model
 *
 * Complexity: O(epochs * m * hidden_dim^2)
 * Reference: Lusch, Kutz, Brunton (2018)
 */
void train_deep_koopman(const double *X, const double *Y, int n, int m,
    int latent_dim, DeepKoopmanModel *model);

/**
 * L9: Time-delay embedding for Koopman analysis (Hankel-DMD).
 *
 * Uses Hankel matrix of time-delayed observations to enrich
 * the observable space for Koopman approximation, following
 * Takens' embedding theorem.
 *
 * @param X           Time series (n x m, col-major)
 * @param n           State dimension
 * @param m           Number of snapshots
 * @param delay       Embedding delay
 * @param result      Output DMDResult
 *
 * Complexity: O(n*delay * m^2)
 * Reference: Arbabi & Mezic (2017), SIAM J. Appl. Dyn. Syst.
 */
void hankel_dmd(const double *X, int n, int m, int delay, DMDResult *result);

/**
 * L8: Sparsity-promoting DMD (SP-DMD) for selecting dominant modes.
 *
 * Solves a regularized least-squares problem to select a sparse
 * set of DMD modes that best represent the data.
 *
 * @param dmd     Full DMD result (modified in place to sparse subset)
 * @param gamma   Sparsity regularization parameter
 * @param n_sparse Output: number of retained modes
 *
 * Complexity: O(r^3 * iterations)
 * Reference: Jovanovic, Schmid, Nichols (2014), Phys. Fluids
 */
void sparsity_promoting_dmd(DMDResult *dmd, double gamma, int *n_sparse);

/**
 * L8: Online/streaming DMD for real-time Koopman analysis.
 *
 * Updates DMD incrementally as new data arrives, without storing
 * the full data matrix.
 *
 * @param x_new        New state snapshot
 * @param n            State dimension
 * @param dmd          Current DMD result (updated in place)
 * @param forgetting   Forgetting factor (0 < lambda <= 1)
 *
 * Complexity: O(r^2 * n) per update
 * Reference: Zhang et al. (2019), SIAM J. Appl. Dyn. Syst.
 */
void online_dmd_update(const double *x_new, int n, DMDResult *dmd, double forgetting);

/**
 * L9: Koopman-based model predictive control (K-MPC).
 *
 * Uses Koopman operator to formulate a linear MPC problem
 * in the lifted observable space.
 *
 * @param dmd          Koopman/DMD model
 * @param x0           Current state
 * @param horizon      Prediction horizon
 * @param ref          Reference trajectory (horizon x n)
 * @param dt           Time step
 * @param u_seq        Output: optimal input sequence (horizon x m)
 * @param n, m         State and input dimensions
 *
 * Complexity: O(horizon * r^3)
 * Reference: Korda & Mezic (2018), Automatica
 */
void koopman_mpc(const DMDResult *dmd, const double *x0,
    int horizon, const double *ref, double dt,
    double *u_seq, int n, int m);




#endif /* KOOPMAN_OPERATOR_H */
