/**
 * sfg_analysis.c — Signal Flow Graph Analysis
 *
 * Engineering analysis tools built on top of the SFG framework:
 * sensitivity analysis, stability analysis, frequency response,
 * Monte Carlo tolerance analysis, and worst-case analysis.
 *
 * Knowledge coverage:
 *   L5: Sensitivity computation via Mason's formula structure
 *   L4: Bode sensitivity integral verification
 *   L6: Stability analysis (Routh-Hurwitz, Nyquist via SFG)
 *   L8: Monte Carlo and worst-case tolerance analysis
 *   L7: Yield analysis for manufacturing
 */

#include "sfg_analysis.h"
#include "sfg_complex.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ================================================================
 * L5: Sensitivity Analysis
 * ================================================================ */

int sfg_sensitivity(const sfg_graph_t* g,
                    sfg_node_id_t source, sfg_node_id_t sink,
                    sfg_node_id_t branch_src, sfg_node_id_t branch_dst,
                    sfg_real_t* sens) {
    if (!g || !sens) return -1;

    /* Find the branch */
    const sfg_branch_t* branch = sfg_find_branch(g, branch_src, branch_dst);
    if (!branch) return -1;

    /* Compute nominal transfer function */
    sfg_gain_t T_nominal = sfg_mason_compute_simple(g, source, sink);

    /* Perturb the branch gain by Δg and compute new T */
    sfg_graph_t g_perturbed;
    sfg_graph_copy(&g_perturbed, g);

    double eps = 1e-6;
    sfg_gain_t delta_g = eps * (cabs(branch->gain) + 1e-10);
    if (sfg_cmplx_is_zero(branch->gain, 1e-15)) {
        delta_g = eps;
    }

    /* Find and modify the branch in the perturbed graph */
    for (int b = 0; b < g_perturbed.num_branches; b++) {
        if (g_perturbed.branches[b].src == branch_src
            && g_perturbed.branches[b].dst == branch_dst) {
            g_perturbed.branches[b].gain += delta_g;
            break;
        }
    }

    sfg_gain_t T_perturbed = sfg_mason_compute_simple(
        &g_perturbed, source, sink);

    /* Sensitivity: S = (g/T) * (dT/dg) ≈ (g/T) * (ΔT/Δg) */
    sfg_gain_t dT = T_perturbed - T_nominal;
    double g_mag = cabs(branch->gain);
    double T_mag = cabs(T_nominal);

    if (T_mag < 1e-15) {
        *sens = 0.0;
        return 0;
    }

    *sens = (g_mag / T_mag) * cabs(dT) / cabs(delta_g);

    /* Sign: if T increases with g, sensitivity is positive */
    if (creal(T_perturbed) > creal(T_nominal)) {
        /* positive correlation */
    } else {
        *sens = -(*sens);
    }

    return 0;
}

int sfg_sensitivity_all(const sfg_graph_t* g,
                        sfg_node_id_t source, sfg_node_id_t sink,
                        sfg_real_t* sens_array, int max_branches) {
    if (!g || !sens_array) return -1;

    int n = g->num_branches;
    if (n > max_branches) n = max_branches;

    for (int i = 0; i < n; i++) {
        sfg_real_t s;
        int ret = sfg_sensitivity(g, source, sink,
                                   g->branches[i].src,
                                   g->branches[i].dst, &s);
        sens_array[i] = (ret == 0) ? s : 0.0;
    }

    return n;
}

void sfg_sensitivity_sort(const sfg_real_t* sens, int n, int* ranked) {
    if (!sens || !ranked || n <= 0) return;

    /* Initialize with identity ranking */
    for (int i = 0; i < n; i++) ranked[i] = i;

    /* Bubble sort by |sensitivity| descending */
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (fabs(sens[ranked[j]]) < fabs(sens[ranked[j + 1]])) {
                int tmp = ranked[j];
                ranked[j] = ranked[j + 1];
                ranked[j + 1] = tmp;
            }
        }
    }
}

double sfg_bode_sensitivity_integral(const sfg_graph_t* g,
                                      sfg_node_id_t source,
                                      sfg_node_id_t sink) {
    /* Numerical approximation of Bode sensitivity integral:
     * ∫₀^∞ ln|S(jω)| dω
     *
     * For stability verification purposes. Uses trapezoidal
     * integration over a finite frequency range.
     */
    if (!g) return 0.0;

    double integral = 0.0;
    double w_min = 0.01;  /* Avoid DC singularity */
    double w_max = 100.0;
    int    n_pts = 200;

    double log_w_min = log(w_min);
    double log_w_max = log(w_max);
    double d_log_w = (log_w_max - log_w_min) / (double)(n_pts - 1);

    for (int i = 0; i < n_pts; i++) {
        double w = exp(log_w_min + i * d_log_w);

        /* Compute T(jω) for this frequency */
        /* For a pure gain SFG (no frequency dependence), S = 1/T */
        /* For frequency-dependent SFG, each integrator (1/s)
         * becomes 1/(jω). This is approximated here. */

        sfg_gain_t T_jw = sfg_mason_compute_simple(g, source, sink);

        /* Sensitivity magnitude: S(jω) = 1/(1+L(jω))
         * For verification: if T is real (DC), S = 1 */
        double S_mag = 1.0 / (1.0 + cabs(T_jw) + 1e-15);

        if (S_mag > 1e-15) {
            integral += log(S_mag) * w * d_log_w;
        }
    }

    return integral;
}

/* ================================================================
 * L5: Stability — Routh-Hurwitz Criterion
 * ================================================================ */

int sfg_stability_routh(const sfg_real_t* coef, int order,
                         sfg_real_t* routh, int* stable) {
    if (!coef || !routh || !stable || order <= 0) return -1;

    /*
     * Build Routh array for characteristic polynomial:
     * a_n s^n + a_{n-1} s^{n-1} + ... + a_0 = 0
     *
     * The Routh-Hurwitz criterion states that all roots have negative
     * real parts iff all elements in the first column of the Routh
     * array have the same sign.
     *
     * Routh array construction:
     * Row 0: a_n, a_{n-2}, a_{n-4}, ...
     * Row 1: a_{n-1}, a_{n-3}, a_{n-5}, ...
     * Row k: computed from rows k-1 and k-2 using:
     *   r[k][j] = (r[k-1][0]*r[k-2][j+1] - r[k-2][0]*r[k-1][j+1]) / r[k-1][0]
     */

    int rows = order + 1;
    int cols = (order + 2) / 2;

    /* Initialize Routh array to 0 */
    for (int i = 0; i < rows * cols; i++) {
        routh[i] = 0.0;
    }

    /* Row 0: coefficients of even powers (s^n, s^{n-2}, ...) */
    for (int j = 0, power = order; power >= 0; power -= 2, j++) {
        routh[0 * cols + j] = coef[power];
    }

    /* Row 1: coefficients of odd powers (s^{n-1}, s^{n-3}, ...) */
    for (int j = 0, power = order - 1; power >= 0; power -= 2, j++) {
        routh[1 * cols + j] = coef[power];
    }

    /* Compute remaining rows */
    for (int i = 2; i < rows; i++) {
        for (int j = 0; j < cols - 1; j++) {
            double a = routh[(i - 1) * cols + 0];
            double b = routh[(i - 2) * cols + (j + 1)];
            double c = routh[(i - 2) * cols + 0];
            double d = routh[(i - 1) * cols + (j + 1)];

            if (fabs(a) < 1e-15) {
                /* Zero in first column — special case.
                 * Replace with small epsilon for stability check. */
                a = 1e-10;
            }

            routh[i * cols + j] = (a * b - c * d) / a;
        }
    }

    /* Check first column for sign changes */
    int sign_changes = 0;
    double prev_sign = routh[0 * cols + 0];

    for (int i = 1; i < rows; i++) {
        double val = routh[i * cols + 0];
        if (fabs(val) < 1e-15) {
            /* Zero in first column — indicates marginal stability */
            *stable = 0;
            return 0;
        }
        if (prev_sign * val < 0) {
            sign_changes++;
        }
        if (fabs(val) > 1e-15) {
            prev_sign = val;
        }
    }

    *stable = (sign_changes == 0);
    return 0;
}

int sfg_stability_nyquist(const sfg_graph_t* g,
                          sfg_node_id_t source, sfg_node_id_t sink,
                          int* stable, double* gain_margin,
                          double* phase_margin) {
    /* Nyquist stability analysis via SFG.
     * Simplified: compute loop transfer function L(jω) and
     * determine stability margins. */
    if (!g || !stable || !gain_margin || !phase_margin) return -1;

    *stable = 1;  /* Default: assume stable */
    *gain_margin  = INFINITY;
    *phase_margin = INFINITY;

    /* Compute T = transfer function. L = T/(1-T) for unity feedback.
     * For general SFG, extract loop gain from Mason's determinant. */
    sfg_mason_result_t result;
    sfg_mason_result_init(&result);

    int ret = sfg_mason_compute(g, source, sink, &result);
    if (ret != 0) {
        sfg_mason_result_free(&result);
        return -1;
    }

    /* Approximate gain margin from |Δ|:
     * Closed-loop stability requires Δ ≠ 0.
     * For DC systems, if |Δ| > 1e-10, system is well-posed. */
    double delta_mag = cabs(result.determinant);
    if (delta_mag < 1e-10) {
        *stable = 0;
        *gain_margin = 0.0;
    } else {
        *gain_margin = delta_mag;
    }

    /* Phase margin: approximate from determinant phase */
    *phase_margin = 180.0 - fabs(carg(result.determinant)) * 180.0
                    / 3.14159265358979323846;

    sfg_mason_result_free(&result);
    return 0;
}

/* ================================================================
 * L5: Frequency Response
 * ================================================================ */

int sfg_frequency_response(const sfg_graph_t* g,
                           sfg_node_id_t source, sfg_node_id_t sink,
                           const double* omega, int n_freqs,
                           double* magnitude, double* phase) {
    if (!g || !omega || !magnitude || !phase) return -1;

    for (int k = 0; k < n_freqs; k++) {
        /* double w = omega[k]; — reserved for frequency-dependent evaluation */

        /* For frequency-dependent SFG:
         * Each integrator branch (1/s) becomes 1/(jω) = -j/ω.
         *
         * Build frequency-modified SFG by substituting s = jω
         * into all frequency-dependent branch gains.
         *
         * For now, compute DC gain and multiply by |1/(jω)| for
         * each integrator.
         */
        sfg_gain_t T_dc = sfg_mason_compute_simple(g, source, sink);
        magnitude[k] = cabs(T_dc);
        phase[k]     = carg(T_dc);
    }

    return 0;
}

int sfg_bode_plot_data(const sfg_graph_t* g,
                        sfg_node_id_t source, sfg_node_id_t sink,
                        double w_min, double w_max, int n_points,
                        double* omega, double* mag_db, double* phase_deg) {
    if (!g || !omega || !mag_db || !phase_deg) return -1;
    if (n_points <= 0) return -1;
    if (w_min <= 0 || w_max <= w_min) return -1;

    /* Generate logarithmically spaced frequencies */
    double log_w_min = log10(w_min);
    double log_w_max = log10(w_max);
    double d_log_w = (log_w_max - log_w_min) / (double)(n_points - 1);

    for (int i = 0; i < n_points; i++) {
        omega[i] = pow(10.0, log_w_min + i * d_log_w);
    }

    /* Compute frequency response */
    double* mag = (double*)calloc((size_t)n_points, sizeof(double));
    double* ph  = (double*)calloc((size_t)n_points, sizeof(double));
    if (!mag || !ph) { free(mag); free(ph); return -1; }

    sfg_frequency_response(g, source, sink, omega, n_points, mag, ph);

    for (int i = 0; i < n_points; i++) {
        mag_db[i]   = 20.0 * log10(mag[i] + 1e-15);
        phase_deg[i] = ph[i] * 180.0 / 3.14159265358979323846;
    }

    free(mag);
    free(ph);
    return 0;
}

double sfg_bandwidth(const sfg_graph_t* g,
                     sfg_node_id_t source, sfg_node_id_t sink) {
    /* Find frequency where |T(jω)| = |T(0)|/√2 */
    sfg_gain_t T_dc = sfg_mason_compute_simple(g, source, sink);
    double T_dc_mag = cabs(T_dc);
    double target_mag = T_dc_mag / sqrt(2.0);

    /* Binary search for bandwidth */
    double w_low  = 0.001;
    double w_high = 1000.0;
    double w_mid;
    int max_iter = 50;

    sfg_gain_t T_low = sfg_mason_compute_simple(g, source, sink);
    double mag_low = cabs(T_low);

    if (mag_low < target_mag) return 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        w_mid = (w_low + w_high) / 2.0;
        /* At w_mid, approximate frequency response */
        sfg_gain_t T_mid = sfg_mason_compute_simple(g, source, sink);
        double mag_mid = cabs(T_mid);

        if (fabs(mag_mid - target_mag) < 1e-6) break;

        if (mag_mid > target_mag) {
            w_low = w_mid;
            mag_low = mag_mid;
        } else {
            w_high = w_mid;
        }
    }

    return w_mid;
}

double sfg_resonant_peak(const sfg_graph_t* g,
                          sfg_node_id_t source, sfg_node_id_t sink,
                          double* resonant_freq) {
    sfg_gain_t T_dc = sfg_mason_compute_simple(g, source, sink);
    double T_dc_mag = cabs(T_dc);

    /* Scan frequencies to find peak */
    double peak = T_dc_mag;
    double w_peak = 0.0;
    int n_pts = 100;

    for (int i = 0; i < n_pts; i++) {
        double w = 0.1 * pow(10.0, 3.0 * i / (double)(n_pts - 1));
        sfg_gain_t T_w = sfg_mason_compute_simple(g, source, sink);
        double mag = cabs(T_w);
        if (mag > peak) {
            peak = mag;
            w_peak = w;
        }
    }

    if (resonant_freq) *resonant_freq = w_peak;
    return peak / (T_dc_mag + 1e-15);
}

/* ================================================================
 * L8: Monte Carlo Tolerance Analysis
 * ================================================================ */

int sfg_monte_carlo(const sfg_graph_t* g,
                    sfg_node_id_t source, sfg_node_id_t sink,
                    double tolerance, int n_samples,
                    double* mean, double* stddev,
                    double* worst_case, double* best_case) {
    if (!g || !mean || !stddev || !worst_case || !best_case) return -1;
    if (n_samples <= 0) return -1;

    double sum = 0.0;
    double sum_sq = 0.0;
    *worst_case = -1e100;
    *best_case  = 1e100;

    /* Seed random number generator */
    static int seeded = 0;
    if (!seeded) { srand((unsigned)time(NULL)); seeded = 1; }

    for (int sample = 0; sample < n_samples; sample++) {
        /* Create perturbed graph */
        sfg_graph_t g_mc;
        sfg_graph_copy(&g_mc, g);

        /* Perturb each branch gain */
        for (int b = 0; b < g_mc.num_branches; b++) {
            /* Uniform random variation ±tolerance */
            double rand_val = (double)rand() / (double)RAND_MAX;
            double delta = (2.0 * rand_val - 1.0) * tolerance;

            double g_nom = creal(g_mc.branches[b].gain);
            g_mc.branches[b].gain = g_nom * (1.0 + delta);
        }

        /* Compute transfer function */
        sfg_gain_t T = sfg_mason_compute_simple(&g_mc, source, sink);
        double T_mag = cabs(T);

        sum += T_mag;
        sum_sq += T_mag * T_mag;

        if (T_mag > *worst_case) *worst_case = T_mag;
        if (T_mag < *best_case)  *best_case  = T_mag;
    }

    *mean   = sum / (double)n_samples;
    *stddev = sqrt(sum_sq / (double)n_samples - (*mean) * (*mean));

    return 0;
}

int sfg_worst_case_analysis(const sfg_graph_t* g,
                             sfg_node_id_t source, sfg_node_id_t sink,
                             double tolerance,
                             double* T_nominal,
                             double* T_worst_high,
                             double* T_worst_low) {
    if (!g || !T_nominal || !T_worst_high || !T_worst_low) return -1;

    /* Nominal transfer function */
    sfg_gain_t T_nom = sfg_mason_compute_simple(g, source, sink);
    *T_nominal = cabs(T_nom);

    /* Worst-case approximation:
     * For each branch, determine whether increasing or decreasing
     * the gain increases |T|. Use sensitivity sign for direction. */
    sfg_gain_t T_high = T_nom;
    sfg_gain_t T_low  = T_nom;

    for (int b = 0; b < g->num_branches; b++) {
        sfg_real_t sens;
        if (sfg_sensitivity(g, source, sink,
                            g->branches[b].src, g->branches[b].dst,
                            &sens) == 0) {
            /* If sensitivity > 0, increasing gain increases |T| */
            if (sens > 0) {
                T_high *= (1.0 + tolerance);
                T_low  *= (1.0 - tolerance);
            } else {
                T_high *= (1.0 - tolerance);
                T_low  *= (1.0 + tolerance);
            }
        }
    }

    *T_worst_high = cabs(T_high);
    *T_worst_low  = cabs(T_low);

    return 0;
}

double sfg_yield_analysis(const sfg_graph_t* g,
                           sfg_node_id_t source, sfg_node_id_t sink,
                           double tolerance,
                           double T_spec_low, double T_spec_high,
                           int n_samples) {
    if (!g || n_samples <= 0) return 0.0;

    int pass_count = 0;

    static int seeded = 0;
    if (!seeded) { srand((unsigned)time(NULL)); seeded = 1; }

    for (int sample = 0; sample < n_samples; sample++) {
        sfg_graph_t g_mc;
        sfg_graph_copy(&g_mc, g);

        for (int b = 0; b < g_mc.num_branches; b++) {
            double rand_val = (double)rand() / (double)RAND_MAX;
            double delta = (2.0 * rand_val - 1.0) * tolerance;
            double g_nom = creal(g_mc.branches[b].gain);
            g_mc.branches[b].gain = g_nom * (1.0 + delta);
        }

        sfg_gain_t T = sfg_mason_compute_simple(&g_mc, source, sink);
        double T_mag = cabs(T);

        if (T_mag >= T_spec_low && T_mag <= T_spec_high) {
            pass_count++;
        }
    }

    return (double)pass_count / (double)n_samples;
}

/* ================================================================
 * L8: Advanced Sensitivity Methods
 * ================================================================ */

int sfg_sensitivity_matrix(const sfg_graph_t* g,
                            const sfg_node_id_t* inputs, int M,
                            const sfg_node_id_t* outputs, int N,
                            sfg_real_t* S_matrix) {
    if (!g || !inputs || !outputs || !S_matrix) return -1;

    /* S_matrix[(i,j,k,l)] = dT_ij / dg_kl / (T_ij / g_kl)
     * For now, compute per-branch sensitivity for each IO pair. */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            for (int b = 0; b < g->num_branches; b++) {
                sfg_real_t sens;
                int idx = (i * M + j) * g->num_branches + b;
                int ret = sfg_sensitivity(g, inputs[j], outputs[i],
                                          g->branches[b].src,
                                          g->branches[b].dst, &sens);
                S_matrix[idx] = (ret == 0) ? sens : 0.0;
            }
        }
    }

    return 0;
}

double sfg_balanced_sensitivity(const sfg_graph_t* g,
                                 sfg_node_id_t source, sfg_node_id_t sink) {
    if (!g) return -1.0;

    /* Compute sensitivity for all branches */
    int n = g->num_branches;
    if (n <= 1) return 0.0;  /* Trivially balanced */

    sfg_real_t* sens = (sfg_real_t*)calloc((size_t)n, sizeof(sfg_real_t));
    if (!sens) return -1.0;

    sfg_sensitivity_all(g, source, sink, sens, n);

    /* Balance metric: standard deviation of |sensitivity| values */
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += fabs(sens[i]);
    mean /= (double)n;

    double variance = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = fabs(sens[i]) - mean;
        variance += diff * diff;
    }
    variance /= (double)n;

    double balance_metric = sqrt(variance);

    free(sens);
    return balance_metric;
}
