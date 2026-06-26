#ifndef SENSOR_MODELING_H
#define SENSOR_MODELING_H

/**
 * @file sensor_modeling.h
 * @brief L1-L5 Sensor Modeling and Signal Processing
 *
 * Models for common mechatronic sensors: encoders, resolvers, tachometers,
 * Hall sensors, accelerometers, gyroscopes, LVDT, strain gauges.
 *
 * Key References:
 *  - Fraden, J. "Handbook of Modern Sensors" (2016)
 *  - Nyce, D.S. "Linear Position Sensors" (2004)
 *  - Titterton & Weston "Strapdown Inertial Navigation Technology" (2004)
 *
 * 9-Layer Coverage:
 *   L1 — Sensor type definitions, parameter structs
 *   L2 — Measurement principles, noise models
 *   L3 — Signal processing (filtering, interpolation, quadrature decode)
 *   L5 — Sensor fusion methods, calibration
 */

#include "mechatronic_definitions.h"
#include <stddef.h>

/*===========================================================================*/
/* L1-L2: Encoder Functions                                                  */
/*===========================================================================*/

/** Convert encoder counts to angular position
 *  θ = (2π · counts) / PPR
 *  @param counts   Encoder count value
 *  @param ppr      Pulses per revolution
 *  @return Angular position [rad] */
double encoder_counts_to_angle(int32_t counts, int ppr);

/** Convert encoder counts to angular position with 4x quadrature decode
 *  Full quadrature yields 4x resolution improvement.
 *  @param counts   Raw count value (after 4x decode)
 *  @param ppr      Pulses per revolution (before 4x)
 *  @return Angular position [rad] */
double encoder_quadrature_to_angle(int32_t counts, int ppr);

/** Incremental encoder velocity estimation (M-method: pulse count over time)
 *  ω = (2π · Δcounts) / (PPR · ΔT)
 *  @param delta_counts  Counts accumulated in interval
 *  @param ppr           Pulses per revolution
 *  @param delta_t       Time interval [s]
 *  @return Angular velocity [rad/s] */
double encoder_velocity_m_method(int delta_counts, int ppr, double delta_t);

/** Encoder velocity estimation (T-method: time between pulses)
 *  ω = (2π) / (PPR · T_between_pulses)
 *  Suitable for low speeds.
 *  @param time_between  Time between consecutive edges [s]
 *  @param ppr           Pulses per revolution
 *  @param quad_encode   1 if quadrature (4 edges/pulse), 0 if single
 *  @return Angular velocity [rad/s] */
double encoder_velocity_t_method(double time_between, int ppr, int quad_encode);

/** Encoder velocity via M/T hybrid method (best of both)
 *  Uses fixed time window with simultaneous pulse counting and timing.
 *  @param pulse_count   Pulses counted in window
 *  @param window_time   Measurement window [s]
 *  @param ppr           Pulses per revolution
 *  @param clock_freq    High-speed counter clock [Hz]
 *  @return Angular velocity [rad/s]
 *
 *  Ref: Lorenz et al. "High-Resolution Velocity Estimation" (1995) */
double encoder_velocity_mt_method(int pulse_count, double window_time,
                                   int ppr, double clock_freq);

/** Encoder resolution in radians
 *  @param ppr   Pulses per revolution
 *  @param quad  1 if quadrature decode, 0 otherwise
 *  @return Resolution [rad] */
double encoder_resolution(int ppr, int quad);

/** Maximum trackable speed for an encoder
 *  ω_max = (2π · f_counter_limit) / PPR
 *  @param ppr            Pulses per revolution
 *  @param counter_freq   Max counter frequency [Hz]
 *  @return Max speed [rad/s] */
double encoder_max_speed(int ppr, double counter_freq);

/*===========================================================================*/
/* L2: Resolver Functions                                                    */
/*===========================================================================*/

/** Resolver output signals (amplitude-modulated sine/cosine)
 *
 *  V_sin = k · V_exc · sin(θ) · sin(ω_exc · t)
 *  V_cos = k · V_exc · cos(θ) · sin(ω_exc · t)
 *
 *  @param theta        Mechanical angle [rad]
 *  @param V_exc        Excitation amplitude [V]
 *  @param k            Transformation ratio
 *  @param omega_exc    Excitation frequency [rad/s]
 *  @param t            Time [s]
 *  @param V_sin        Output: sine winding voltage [V]
 *  @param V_cos        Output: cosine winding voltage [V] */
void resolver_outputs(double theta, double V_exc, double k,
                       double omega_exc, double t,
                       double *V_sin, double *V_cos);

/** Resolver-to-digital (RDC) angle extraction via arctangent
 *  θ = atan2(V_sin_envelope, V_cos_envelope)
 *
 *  @param V_sin_demod  Demodulated sine amplitude [V]
 *  @param V_cos_demod  Demodulated cosine amplitude [V]
 *  @return Angle [rad] */
double resolver_angle_from_demod(double V_sin_demod, double V_cos_demod);

/** Resolver maximum tracking rate
 *  ω_max_track = ω_exc / (2π · pole_pairs)
 *  @param excitation_freq  Excitation frequency [Hz]
 *  @param pole_pairs       Resolver pole pairs
 *  @return Max tracking speed [rad/s] */
double resolver_max_tracking_rate(double excitation_freq, int pole_pairs);

/** Resolver accuracy expressed as angular error
 *  @param accuracy_arcmin  Accuracy in arc-minutes
 *  @return Angular error [rad] */
double resolver_accuracy_rad(double accuracy_arcmin);

/*===========================================================================*/
/* L2: Tachometer Functions                                                  */
/*===========================================================================*/

/** DC tachometer output voltage
 *  V_out = K_tach · ω
 *  @param K_tach   Tachometer constant [V·s/rad]
 *  @param omega    Angular velocity [rad/s]
 *  @return Output voltage [V] */
double tachometer_voltage(double K_tach, double omega);

/** Tachometer ripple frequency
 *  f_ripple = (ω · n_commutator_segments) / (2π)
 *  @param omega      Speed [rad/s]
 *  @param n_segments Number of commutator segments
 *  @return Ripple frequency [Hz] */
double tachometer_ripple_frequency(double omega, int n_segments);

/*===========================================================================*/
/* L2: Hall Effect Sensor                                                    */
/*===========================================================================*/

/** Hall sensor voltage (linear mode)
 *  V_H = S_H · B · I_bias
 *  @param S_H      Hall sensitivity [V/(A·T)]
 *  @param B        Magnetic flux density perpendicular to sensor [T]
 *  @param I_bias   Bias current [A]
 *  @return Hall voltage [V] */
double hall_voltage(double S_H, double B, double I_bias);

/** Hall sensor for BLDC commutation (binary mode)
 *  Returns 1 if B exceeds threshold B_th, 0 otherwise.
 *  @param B        Magnetic flux density [T]
 *  @param B_th     Switching threshold [T]
 *  @param hysteresis Hysteresis band [T]
 *  @param prev_state Previous output state (for hysteresis)
 *  @return 0 or 1 */
int hall_digital_output(double B, double B_th, double hysteresis, int prev_state);

/** Estimate BLDC electrical angle from 3 Hall sensors
 *  Returns electrical angle in [0, 2π) based on Hall pattern.
 *  Pattern: HA, HB, HC (1=high, 0=low) → 6-step commutation
 *
 *  @param HA  Hall sensor A state (0 or 1)
 *  @param HB  Hall sensor B state (0 or 1)
 *  @param HC  Hall sensor C state (0 or 1)
 *  @return Electrical angle [rad] (center of sector) */
double bldc_angle_from_halls(int HA, int HB, int HC);

/*===========================================================================*/
/* L2: Accelerometer and Gyroscope Models                                    */
/*===========================================================================*/

/** Accelerometer measurement model
 *  a_measured = a_true + bias + scale_factor_error · a_true + noise
 *
 *  @param a_true            True acceleration [m/s²]
 *  @param bias              Constant bias [m/s²]
 *  @param scale_factor_err  Scale factor error (fraction)
 *  @param noise_std         Noise standard deviation [m/s²]
 *  @param noise_sample      Random noise sample (e.g., from Gaussian RNG)
 *  @return Measured acceleration [m/s²] */
double accelerometer_measurement(double a_true, double bias,
                                  double scale_factor_err,
                                  double noise_std, double noise_sample);

/** Gyroscope measurement model
 *  ω_measured = ω_true + bias + scale_factor_error · ω_true + noise
 *
 *  @param omega_true        True angular rate [rad/s]
 *  @param bias              Constant bias [rad/s]
 *  @param scale_factor_err  Scale factor error (fraction)
 *  @param noise_std         Noise standard deviation [rad/s]
 *  @param noise_sample      Random noise sample
 *  @return Measured angular rate [rad/s] */
double gyroscope_measurement(double omega_true, double bias,
                              double scale_factor_err,
                              double noise_std, double noise_sample);

/** Allan variance noise characterization for IMU
 *  Computes angle random walk and bias instability from noise density.
 *
 *  @param noise_density  Noise spectral density [(m/s²)/√Hz or (rad/s)/√Hz)]
 *  @param bandwidth      Sensor bandwidth [Hz]
 *  @param arw            Output: Angle/Velocity random walk
 *  @param bias_instab    Output: Bias instability estimate */
void imu_allan_variance_parameters(double noise_density, double bandwidth,
                                    double *arw, double *bias_instab);

/*===========================================================================*/
/* L3: Sensor Signal Processing                                              */
/*===========================================================================*/

/** First-order low-pass filter for sensor signals (IIR)
 *  y[k] = α · x[k] + (1-α) · y[k-1]
 *  α = dt / (τ + dt)  where τ = 1/(2π·f_c)
 *
 *  @param input       Current raw measurement
 *  @param prev_output  Previous filtered output
 *  @param cutoff_freq  Cutoff frequency [Hz]
 *  @param dt           Sampling period [s]
 *  @return Filtered output */
double sensor_lowpass_filter(double input, double prev_output,
                              double cutoff_freq, double dt);

/** Moving average filter for sensor signals
 *  y[k] = (1/N) · Σ_{i=0}^{N-1} x[k-i]
 *
 *  Uses buffer array and circular indexing.
 *
 *  @param input    New sample
 *  @param buffer   Circular buffer of N previous samples
 *  @param N        Window size
 *  @param index    Current buffer index (updated in-place to (index+1)%N)
 *  @param sum      Running sum (updated in-place: sum + new - old)
 *  @return Moving average output */
double sensor_moving_average(double input, double *buffer, int N,
                              int *index, double *sum);

/** Sensor signal derivative estimation (backward difference + lowpass)
 *  v[k] = (x[k] - x[k-1]) / dt
 *
 *  Due to noise amplification, typically followed by lowpass filter.
 *
 *  @param current   Current position
 *  @param previous  Previous position
 *  @param dt        Time step [s]
 *  @return Velocity estimate [unit/s] */
double sensor_derivative(double current, double previous, double dt);

/** Quadrature encoder direction detection
 *  Given current and previous A/B states, determine direction.
 *
 *  @param A_prev, B_prev  Previous channel states (0 or 1)
 *  @param A_curr, B_curr  Current channel states (0 or 1)
 *  @return +1 for forward, -1 for reverse, 0 for no change (error) */
int encoder_direction(int A_prev, int B_prev, int A_curr, int B_curr);

/*===========================================================================*/
/* L3: LVDT (Linear Variable Differential Transformer)                       */
/*===========================================================================*/

/** LVDT output voltage
 *  V_out = sensitivity · position · V_exc
 *
 *  @param position     Displacement from null [m]
 *  @param sensitivity  LVDT sensitivity [V/(V_exc·m)]
 *  @param V_exc        Primary excitation voltage [V]
 *  @return Output voltage [V] */
double lvdt_output(double position, double sensitivity, double V_exc);

/** LVDT null voltage (residual at zero position)
 *  @param quadrature_factor  Residual quadrature component fraction
 *  @param V_exc             Excitation voltage [V]
 *  @return Null voltage [V] */
double lvdt_null_voltage(double quadrature_factor, double V_exc);

/*===========================================================================*/
/* L3: Strain Gauge                                                          */
/*===========================================================================*/

/** Strain gauge resistance change
 *  ΔR/R = GF · ε
 *  where GF = gauge factor, ε = strain = ΔL/L
 *
 *  @param R_nominal  Nominal resistance [Ω]
 *  @param GF         Gauge factor (typically ~2.0 for metal foil)
 *  @param strain     Mechanical strain [m/m]
 *  @return Resistance change ΔR [Ω] */
double strain_gauge_resistance_change(double R_nominal, double GF,
                                       double strain);

/** Wheatstone quarter-bridge output voltage
 *  V_out = V_exc · (GF · ε) / 4  (linearized for small ε)
 *
 *  @param V_exc   Bridge excitation voltage [V]
 *  @param GF      Gauge factor
 *  @param strain  Mechanical strain [m/m]
 *  @return Bridge output voltage [V] */
double wheatstone_quarter_bridge(double V_exc, double GF, double strain);

/** Wheatstone half-bridge output (two active gauges, opposite strain)
 *  V_out = V_exc · (GF · ε) / 2
 *
 *  @param V_exc   Bridge excitation voltage [V]
 *  @param GF      Gauge factor
 *  @param strain  Mechanical strain [m/m]
 *  @return Bridge output voltage [V] */
double wheatstone_half_bridge(double V_exc, double GF, double strain);

/** Wheatstone full-bridge output (four active gauges)
 *  V_out = V_exc · GF · ε
 *
 *  @param V_exc   Bridge excitation voltage [V]
 *  @param GF      Gauge factor
 *  @param strain  Mechanical strain [m/m]
 *  @return Bridge output voltage [V] */
double wheatstone_full_bridge(double V_exc, double GF, double strain);

/** Load cell calibration: Force from bridge voltage
 *  F = (V_out / V_exc) · (4 / GF) · (A · E) for full bridge
 *  where A = cross-sectional area, E = Young's modulus
 *
 *  @param V_out     Bridge output [V]
 *  @param V_exc     Excitation voltage [V]
 *  @param GF        Gauge factor
 *  @param A         Cross-sectional area of load cell element [m²]
 *  @param E_modulus  Young's modulus of load cell material [Pa]
 *  @return Force [N] */
double load_cell_force(double V_out, double V_exc, double GF,
                        double A, double E_modulus);

/*===========================================================================*/
/* L5: Sensor Fusion                                                         */
/*===========================================================================*/

/** Complementary filter for IMU attitude (simplified 1-axis)
 *  angle[k] = α · (angle[k-1] + gyro·dt) + (1-α) · accel_angle
 *  where α = τ / (τ + dt), τ = filter time constant
 *
 *  Gyro: good for short-term (integrates drift)
 *  Accel: good for long-term (noisy, but no drift)
 *
 *  @param gyro_rate      Gyro angular rate measurement [rad/s]
 *  @param accel_angle    Angle estimate from accelerometer [rad]
 *  @param prev_angle     Previous filtered angle [rad]
 *  @param dt             Time step [s]
 *  @param tau            Filter time constant [s] (higher = trust gyro more)
 *  @return Fused angle estimate [rad]
 *
 *  Ref: Mahony et al. "Complementary Filter Design" (2008) */
double complementary_filter(double gyro_rate, double accel_angle,
                             double prev_angle, double dt, double tau);

/** Simplified 1D Kalman filter for position-velocity estimation
 *
 *  State: [x, v]^T
 *  Measurement: x_measured (with noise)
 *
 *  Predict:  x = x + v·dt + 0.5·a·dt²
 *            v = v + a·dt
 *            P = F·P·F^T + Q
 *
 *  Update:   K = P·H^T · (H·P·H^T + R)^(-1)
 *            x = x + K·(z - H·x)
 *            P = (I - K·H)·P
 *
 *  @param pos_measured   Measured position [m]
 *  @param acceleration   Known acceleration input [m/s²]
 *  @param dt             Time step [s]
 *  @param pos_est        Position estimate (input/output) [m]
 *  @param vel_est        Velocity estimate (input/output) [m/s]
 *  @param P              Covariance matrix [2x2] (input/output)
 *  @param Q_pos          Process noise variance for position
 *  @param Q_vel          Process noise variance for velocity
 *  @param R_meas         Measurement noise variance */
void kalman_filter_1d_position(double pos_measured, double acceleration,
                                double dt,
                                double *pos_est, double *vel_est,
                                double *P,
                                double Q_pos, double Q_vel, double R_meas);

#endif /* SENSOR_MODELING_H */
