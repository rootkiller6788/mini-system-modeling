/**
 * @file gain_scheduling.h
 * @brief Gain scheduling for nonlinear control across operating points.
 *
 * Knowledge Coverage:
 *   L1 - Gain scheduling, scheduling variable, operating point interpolation
 *   L2 - Linear parameter-varying (LPV) systems, bumpless transfer
 *   L5 - Scheduling surface, interpolation methods, hysteresis switching
 *   L6 - Gain-scheduled controller implementation
 *   L7 - Aerospace/missile/engine control (F-35, Toyota, Tesla)
 *
 * Reference:
 *   Rugh & Shamma (2000), "Research on gain scheduling", Automatica
 *   Apkarian & Gahinet (1995), "LPV control"
 *   Khalil, Nonlinear Systems (2002), Ch.12
 */

#ifndef GAIN_SCHEDULING_H
#define GAIN_SCHEDULING_H

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "linearization_core.h"

/** L1: Single entry in a gain schedule. */
typedef struct {
    double scheduling_value[4];
    int n_scheduling_vars;
    LinearizedSystem *sys;
    double *K;     /**< Feedback gain (m x n) */
    double *L;     /**< Observer gain (n x p) */
    int n_inputs, n_states, n_outputs;
} GainScheduleEntry;

/** L1: Full gain schedule spanning the operating envelope. */
typedef struct {
    GainScheduleEntry *entries;
    int n_entries;
    int n_scheduling_vars;
    int n_states, n_inputs, n_outputs;
    bool use_hysteresis;
    double hysteresis_band;
} GainSchedule;

/** L5: Interpolation methods for scheduling. */
typedef enum {
    GS_NEAREST_NEIGHBOR,
    GS_LINEAR_INTERP,
    GS_BLENDING_FUNCTION,
    GS_LPV_LINEAR
} GSInterpMethod;

/** L6: Runtime gain-scheduled controller. */
typedef struct {
    GainSchedule *schedule;
    GSInterpMethod method;
    int current_entry_index;
    double *current_K;
    double *current_L;
    double *observer_state;
    double prev_u;
} GainScheduledController;

/** L5: Create empty gain schedule. O(1). */
GainSchedule *create_gain_schedule(int n_entries, int n_sched_vars);

/** L5: Add entry at specific operating point. O(1) + memcpy. */
void add_gain_schedule_entry(GainSchedule *gs, int index,
    const double *sched_val, const LinearizedSystem *sys,
    const double *K, const double *L,
    int n_states, int n_inputs, int n_outputs);

/** L5: Interpolate controller gains at scheduling point. O(2^d * n^2). */
void interpolate_gain_schedule(const GainSchedule *gs,
    const double *sched_val, double *K_out, double *L_out,
    GSInterpMethod method);

/** L6: Compute control via gain-scheduled controller. O(n*m + obs_update). */
double compute_gain_scheduled_control(
    GainScheduledController *ctrl, const double *sched_val,
    const double *y, double dt,
    const double *x_eq, const double *u_eq);

/** L6: Bumpless transfer first-order filter. O(1). */
double bumpless_transfer_filter(double u_new, double u_old, double dt, double tau);

/** L6: Hysteresis-based scheduling variable filtering. O(n_entries). */
int hysteresis_scheduling(const double *sched_val, const double *prev_val,
    int current_entry, const GainSchedule *gs, double band);

/** L7: Design LQR gains for all operating points. O(n_entries * n^3). */
void design_lqr_schedule(GainSchedule *gs, const double *Q, const double *R);

/** L2: Extract LPV model at given scheduling parameter. O(2^d * n^2). */
void gain_schedule_to_lpv(const GainSchedule *gs, const double *sched_val,
    double *A_out, double *B_out, double *C_out, double *D_out);

/** Initialize gain-scheduled controller. */
GainScheduledController *create_gain_scheduled_controller(
    GainSchedule *gs, GSInterpMethod method, int n_states);

/** Free GainSchedule. O(n_entries). */
void free_gain_schedule(GainSchedule *gs);

/** Free GainScheduledController. O(1). */
void free_gain_scheduled_controller(GainScheduledController *ctrl);

/**
 * L7: Design pole-placement gains for each operating point.
 *
 * Uses Ackermann's formula for SISO, or robust pole placement for MIMO.
 * The desired poles should be chosen faster than the plant dynamics
 * and should respect actuator limits.
 *
 * @param gs       Gain schedule
 * @param poles    Desired closed-loop poles for each entry
 *                 (n_entries x n_states, row-major)
 *
 * Complexity: O(n_entries * n_states^3)
 * Reference: Ogata (2010), section 12.6; Kautsky et al. (1985)
 */
void design_pole_placement_schedule(GainSchedule *gs, const double *poles);

/**
 * L8: Design robust gain schedule with Lyapunov stability guarantee.
 *
 * Ensures stability across the entire operating envelope by
 * solving coupled Lyapunov inequalities for adjacent operating points.
 *
 * @param gs       Gain schedule
 * @param alpha    Decay rate (negative real part bound)
 *
 * Complexity: O(n_entries * n_states^3)
 * Reference: Apkarian & Gahinet (1995), section 4
 */
void design_robust_gain_schedule(GainSchedule *gs, double alpha);

#endif /* GAIN_SCHEDULING_H */
