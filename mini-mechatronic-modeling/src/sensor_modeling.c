/**
 * @file sensor_modeling.c
 * @brief L1-L5 Sensor modeling and signal processing.
 *
 * Implements encoder (incremental/absolute), resolver, tachometer,
 * Hall sensor, IMU, LVDT, strain gauge, and sensor fusion algorithms.
 *
 * 9-Layer Coverage:
 *   L1 — Sensor definitions and parameter structs
 *   L2 — Measurement principles and noise models
 *   L3 — Signal processing (filtering, quadrature, M/T method)
 *   L5 — Sensor fusion (complementary filter, Kalman filter)
 *
 * References:
 *   Fraden "Handbook of Modern Sensors" (2016)
 *   Nyce "Linear Position Sensors" (2004)
 *   Titterton & Weston "Strapdown Inertial Navigation Technology" (2004)
 *   Mahony et al. (2008)
 */

#include "sensor_modeling.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================*/
/* L1-L2: Encoder Functions                                                  */
/*===========================================================================*/

double encoder_counts_to_angle(int32_t counts, int ppr)
{
    /*
     * Incremental encoder: convert pulse count to angle.
     *
     *   θ = (2π / PPR) · counts
     *
     * For incremental encoder: counts is relative.
     * For absolute encoder: counts encodes absolute position within one rev.
     *
     * Ref: Fraden (2016) §5.3.1
     */
    if (ppr <= 0) return 0.0;
    return (2.0 * M_PI * (double)counts) / (double)ppr;
}

double encoder_quadrature_to_angle(int32_t counts, int ppr)
{
    /*
     * Quadrature encoder: 4x decode yields 4·PPR counts per revolution.
     *
     *   θ = (2π / (4·PPR)) · counts
     *
     * Quadrature encoding uses two channels (A, B) 90° out of phase,
     * detecting all four edges (A↑, B↑, A↓, B↓) per pulse for 4x resolution.
     *
     * Ref: Fraden (2016) §5.3.2
     */
    if (ppr <= 0) return 0.0;
    return (2.0 * M_PI * (double)counts) / (4.0 * (double)ppr);
}

double encoder_velocity_m_method(int delta_counts, int ppr, double delta_t)
{
    /*
     * M-Method velocity estimation (pulse counting over fixed time):
     *
     *   ω = (2π · Δcounts) / (PPR · Ts)
     *
     * Good for high speeds (many counts per sampling period).
     * Poor accuracy at low speeds (few counts → quantization error).
     *
     * Accuracy: ±(2π/(PPR·Ts)) at 1 count uncertainty.
     *
     * Ref: Lorenz et al. "High-Resolution Velocity Estimation for
     *      All-Digital AC Servo Drives" (1991) IEEE Trans. IA
     */
    if (ppr <= 0 || delta_t <= 0.0) return 0.0;
    return (2.0 * M_PI * (double)delta_counts) / ((double)ppr * delta_t);
}

double encoder_velocity_t_method(double time_between, int ppr, int quad_encode)
{
    /*
     * T-Method velocity estimation (pulse timing):
     *
     *   ω = 2π / (PPR · T_between)
     *
     * Measures time between successive encoder pulses using a high-
     * speed counter. Good for low speeds.
     *
     * With quadrature (quad_encode=1), resolution is 4x:
     *   ω = 2π / (4·PPR · T_between)
     *
     * Accuracy degrades at high speeds (timer quantization).
     *
     * Ref: Lorenz et al. (1991)
     */
    if (ppr <= 0 || time_between <= 0.0) return 1e30;
    int edges_per_rev = quad_encode ? (4 * ppr) : ppr;
    return (2.0 * M_PI) / ((double)edges_per_rev * time_between);
}

double encoder_velocity_mt_method(int pulse_count, double window_time,
                                   int ppr, double clock_freq)
{
    /*
     * M/T Hybrid Method:
     *
     * Combines M-method (pulse count) and T-method (time measurement):
     *
     *   ω = (2π · m1) / (PPR · m2 · T_clk)
     *
     * where:
     *   m1 = number of encoder pulses in measurement window
     *   m2 = number of high-speed clock cycles in same window
     *   T_clk = 1/f_clk
     *
     * Advantages:
     *   - Good accuracy over wide speed range
     *   - No switching between M and T methods
     *
     * Ref: Ohmae et al. "A Microprocessor-Controlled High-Accuracy
     *      Wide-Range Speed Regulator for Motor Drives" (1982)
     *      IEEE Trans. IE
     */
    if (ppr <= 0 || window_time <= 0.0 || clock_freq <= 0.0) return 0.0;

    /* m2 = clock cycles in measurement window */
    double m2 = window_time * clock_freq;
    if (m2 < 1.0) return 0.0;

    return (2.0 * M_PI * (double)pulse_count)
           / ((double)ppr * m2 / clock_freq);
}

double encoder_resolution(int ppr, int quad)
{
    /*
     * Angular resolution of an incremental encoder:
     *
     *   θ_res = 2π / (PPR · multiplier)
     *
     * where multiplier = 4 for quadrature, 1 for single-channel.
     *
     * Example: 1000 PPR with quadrature → resolution = 2π/4000 ≈ 0.00157 rad
     */
    if (ppr <= 0) return 0.0;
    int mult = quad ? 4 : 1;
    return (2.0 * M_PI) / ((double)ppr * (double)mult);
}

double encoder_max_speed(int ppr, double counter_freq)
{
    /*
     * Maximum trackable speed limited by counter bandwidth:
     *
     *   ω_max = (2π · f_counter) / PPR
     *
     * Counter frequency must be at least:
     *   f_counter ≥ PPR · ω_max / (2π)
     *
     * For 1000 PPR encoder at 3000 RPM (314 rad/s):
     *   f_counter ≥ 1000 · 314 / (2π) ≈ 50 kHz
     *
     * Ref: Fraden (2016) §5.3
     */
    if (ppr <= 0) return 0.0;
    return (2.0 * M_PI * counter_freq) / (double)ppr;
}

/*===========================================================================*/
/* L2: Resolver Functions                                                    */
/*===========================================================================*/

void resolver_outputs(double theta, double V_exc, double k,
                       double omega_exc, double t,
                       double *V_sin, double *V_cos)
{
    /*
     * Resolver output signals:
     *
     *   V_sin(t) = k · V_exc · sin(θ) · sin(ω_exc·t)
     *   V_cos(t) = k · V_exc · cos(θ) · sin(ω_exc·t)
     *
     * The resolver is a rotary transformer that amplitude-modulates
     * the excitation carrier with sin(θ) and cos(θ).
     *
     * The excitation frequency is typically 400 Hz – 10 kHz,
     * much higher than the mechanical bandwidth.
     *
     * Ref: Fraden (2016) §5.5
     */
    double carrier = sin(omega_exc * t);
    double mod_amp = k * V_exc * carrier;

    if (V_sin) *V_sin = mod_amp * sin(theta);
    if (V_cos) *V_cos = mod_amp * cos(theta);
}

double resolver_angle_from_demod(double V_sin_demod, double V_cos_demod)
{
    /*
     * Angle extraction from demodulated resolver signals:
     *
     *   θ = atan2(V_sin_demod, V_cos_demod)
     *
     * Demodulation removes the carrier, leaving the sin(θ) and cos(θ)
     * envelopes. The tracking RDC (Resolver-to-Digital Converter)
     * typically uses a Type-II servo loop for continuous tracking.
     *
     * Ref: Analog Devices "Resolver-to-Digital Converters" App Note
     */
    return atan2(V_sin_demod, V_cos_demod);
}

double resolver_max_tracking_rate(double excitation_freq, int pole_pairs)
{
    /*
     * Maximum tracking rate for resolver:
     *
     *   ω_max = ω_exc / (2 · pole_pairs)
     *
     * The RDC must sample at least twice per carrier cycle per pole pair.
     *
     * For 10 kHz excitation, 2-pole resolver:
     *   ω_max = 2π·10000 / 4 ≈ 15708 rad/s (150,000 RPM!)
     *
     * Ref: Fraden (2016) §5.5
     */
    if (pole_pairs <= 0) pole_pairs = 1;
    return (2.0 * M_PI * excitation_freq) / (2.0 * (double)pole_pairs);
}

double resolver_accuracy_rad(double accuracy_arcmin)
{
    /*
     * Convert resolver accuracy from arc-minutes to radians.
     *
     * 1 arcmin = 1/60 degree = π / (180·60) radians
     *
     * Typical resolver accuracy: 5–20 arcmin
     */
    return accuracy_arcmin * M_PI / (180.0 * 60.0);
}

/*===========================================================================*/
/* L2: Tachometer                                                            */
/*===========================================================================*/

double tachometer_voltage(double K_tach, double omega)
{
    /*
     * DC tachometer: permanent-magnet DC generator.
     *
     *   V_out = K_tach · ω
     *
     * K_tach is the tachometer constant in V·s/rad.
     * Typical values: 1–20 V / 1000 RPM
     *
     * Ref: Fraden (2016) §5.2
     */
    return K_tach * omega;
}

double tachometer_ripple_frequency(double omega, int n_segments)
{
    /*
     * DC tachometer output ripple:
     *
     *   f_ripple = (ω · n_commutator_segments) / (2π)
     *
     * The ripple comes from commutation between segments.
     * Typically n_segments = 5–20.
     *
     * Ripple amplitude: ~2-5% of output for good tachometers,
     * up to 10% for cheap ones.
     */
    if (n_segments <= 0) return 0.0;
    return omega * (double)n_segments / (2.0 * M_PI);
}

/*===========================================================================*/
/* L2: Hall Effect Sensor                                                    */
/*===========================================================================*/

double hall_voltage(double S_H, double B, double I_bias)
{
    /*
     * Hall effect voltage:
     *
     *   V_H = S_H · B · I_bias
     *
     * where S_H is the Hall coefficient [V/(A·T)].
     * This is the linear (analog) Hall sensor mode.
     *
     * For indium arsenide (InAs): S_H ≈ 100–300 V/(A·T)
     *
     * Ref: Ramsden "Hall-Effect Sensors" (2006)
     */
    return S_H * B * I_bias;
}

int hall_digital_output(double B, double B_th, double hysteresis, int prev_state)
{
    /*
     * Digital Hall sensor with hysteresis (Schmitt trigger behavior):
     *
     *   ON threshold:  B > B_th + hysteresis/2
     *   OFF threshold: B < B_th - hysteresis/2
     *
     * Hysteresis prevents oscillation near the switching threshold.
     * This is used in BLDC motor commutation and position detection.
     *
     * Ref: Allegro Microsystems "Hall-Effect Sensor ICs" App Note
     */
    double B_on = B_th + hysteresis * 0.5;
    double B_off = B_th - hysteresis * 0.5;

    if (B > B_on) return 1;
    if (B < B_off) return 0;
    return prev_state;  /* in hysteresis band, maintain previous state */
}

double bldc_angle_from_halls(int HA, int HB, int HC)
{
    /*
     * Estimate BLDC rotor electrical angle from three Hall sensors
     * spaced 120° electrically (or 60°).
     *
     * Six-step commutation sectors:
     *   HA HB HC | Sector | Center angle
     *    1  0  1 |    1   |   60° (π/3)
     *    1  0  0 |    2   |  120° (2π/3)
     *    1  1  0 |    3   |  180° (π)
     *    0  1  0 |    4   |  240° (4π/3)
     *    0  1  1 |    5   |  300° (5π/3)
     *    0  0  1 |    6   |    0° (0)
     *
     * Returns the center angle of the detected sector.
     *
     * Ref: Krishnan "Electric Motor Drives" (2001) §9.4
     */
    int code = (HA ? 4 : 0) | (HB ? 2 : 0) | (HC ? 1 : 0);

    switch (code) {
    case 5: return M_PI / 3.0;        /* 101: 60° */
    case 4: return 2.0 * M_PI / 3.0;  /* 100: 120° */
    case 6: return M_PI;              /* 110: 180° */
    case 2: return 4.0 * M_PI / 3.0;  /* 010: 240° */
    case 3: return 5.0 * M_PI / 3.0;  /* 011: 300° */
    case 1: return 0.0;               /* 001: 0° / 360° */
    default: return 0.0;              /* 000 or 111: error */
    }
}

/*===========================================================================*/
/* L2: Accelerometer and Gyroscope Measurement Models                        */
/*===========================================================================*/

double accelerometer_measurement(double a_true, double bias,
                                  double scale_factor_err,
                                  double noise_std, double noise_sample)
{
    /*
     * Accelerometer measurement model:
     *
     *   a_meas = a_true + b_a + s_a · a_true + w_a
     *
     * where:
     *   b_a = constant bias (offset error)
     *   s_a = scale factor error (sensitivity error)
     *   w_a = random noise ~ N(0, σ²)
     *
     * Typical MEMS accelerometer:
     *   Bias: 10–50 mg (0.1–0.5 m/s²)
     *   Scale factor error: 0.5–2%
     *   Noise density: 100–500 μg/√Hz
     *
     * Ref: Titterton & Weston (2004) §3.2
     */
    return a_true + bias + scale_factor_err * a_true + noise_sample * noise_std;
}

double gyroscope_measurement(double omega_true, double bias,
                              double scale_factor_err,
                              double noise_std, double noise_sample)
{
    /*
     * Gyroscope measurement model (same structure as accelerometer):
     *
     *   ω_meas = ω_true + b_g + s_g · ω_true + w_g
     *
     * Typical MEMS gyro:
     *   Bias: 0.1–10 °/s
     *   Scale factor error: 0.5–2%
     *   Noise density: 0.01–0.1 (°/s)/√Hz (ARW)
     *
     * Ref: Titterton & Weston (2004) §3.3
     */
    return omega_true + bias + scale_factor_err * omega_true
           + noise_sample * noise_std;
}

void imu_allan_variance_parameters(double noise_density, double bandwidth,
                                    double *arw, double *bias_instab)
{
    /*
     * Allan variance noise characterization:
     *
     * Angle Random Walk (ARW) from noise density:
     *   ARW [rad/√s] = noise_density / √(2 · BW)
     *
     * Bias instability from noise floor:
     *   Bias_instab ≈ noise_density · √(BW)  (rough estimate)
     *
     * The full Allan variance curve shows:
     *   - ARW at short integration times (slope -1/2)
     *   - Bias instability at the minimum (slope 0)
     *   - Rate random walk at long times (slope +1/2)
     *
     * Ref: IEEE Std 952-1997 "Specification Format for
     *      Single-Axis Interferometric Fiber Optic Gyros"
     */
    if (bandwidth < 1e-30) bandwidth = 1.0;

    if (arw) {
        *arw = noise_density / sqrt(2.0 * bandwidth);
    }
    if (bias_instab) {
        *bias_instab = noise_density * sqrt(bandwidth) * 0.1;
    }
}

/*===========================================================================*/
/* L3: Sensor Signal Processing                                              */
/*===========================================================================*/

double sensor_lowpass_filter(double input, double prev_output,
                              double cutoff_freq, double dt)
{
    /*
     * First-order IIR low-pass filter:
     *
     *   y[k] = α · x[k] + (1-α) · y[k-1]
     *
     * where α = dt / (τ + dt), τ = 1/(2π·f_c)
     *
     * Frequency response:
     *   H(jω) = 1 / (1 + jω/ω_c)
     *   |H(jω_c)| = 1/√2 ≈ 0.707 (-3 dB point)
     *
     * Ref: Smith "Digital Signal Processing" (1997) §19.2
     */
    if (cutoff_freq <= 0.0) return input;
    if (dt <= 0.0) return input;

    double tau = 1.0 / (2.0 * M_PI * cutoff_freq);
    double alpha = dt / (tau + dt);

    /* Clamp alpha to valid range */
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    return alpha * input + (1.0 - alpha) * prev_output;
}

double sensor_moving_average(double input, double *buffer, int N,
                              int *index, double *sum)
{
    /*
     * Moving average filter (FIR):
     *
     *   y[k] = (1/N) · Σ_{i=0}^{N-1} x[k-i]
     *
     * Efficient implementation using circular buffer and running sum:
     *   sum = sum - buffer[idx] + input
     *   buffer[idx] = input
     *   idx = (idx + 1) % N
     *   y = sum / N
     *
     * Complexity: O(1) per sample (vs. O(N) naive)
     *
     * Ref: Smith (1997) §15.1
     */
    if (!buffer || !index || !sum || N <= 0) return input;

    *sum = *sum - buffer[*index] + input;
    buffer[*index] = input;
    *index = (*index + 1) % N;

    return *sum / (double)N;
}

double sensor_derivative(double current, double previous, double dt)
{
    /*
     * Backward difference derivative estimate:
     *
     *   v[k] = (x[k] - x[k-1]) / Ts
     *
     * Noise amplification: The derivative amplifies high-frequency
     * noise with gain ∝ f. Therefore:
     *
     *   v_noise_rms = √2 · σ_x / Ts
     *
     * This is why velocity from position requires good resolution
     * or substantial filtering.
     *
     * Ref: Franklin et al. "Digital Control of Dynamic Systems" (1998)
     */
    if (dt <= 0.0) return 0.0;
    return (current - previous) / dt;
}

int encoder_direction(int A_prev, int B_prev, int A_curr, int B_curr)
{
    /*
     * Quadrature encoder direction detection via state transition.
     *
     * Quadrature states (Gray code):
     *   Forward: 00→01→11→10→00→...
     *   Reverse: 00→10→11→01→00→...
     *
     * State mapping:
     *   AB=00→0, AB=01→1, AB=11→2, AB=10→3
     *
     * Forward: state increases (mod 4)
     * Reverse: state decreases (mod 4)
     *
     * Ref: Fraden (2016) §5.3.2
     */
    int prev_state = (A_prev << 1) | B_prev;
    int curr_state = (A_curr << 1) | B_curr;

    if (prev_state == curr_state) return 0;  /* no change */

    /* Map states to Gray code positions */
    int gray_pos[4] = {0, 1, 3, 2};
    int gp = gray_pos[prev_state];
    int gc = gray_pos[curr_state];

    int diff = (gc - gp + 4) % 4;

    if (diff == 1) return 1;   /* forward */
    if (diff == 3) return -1;  /* reverse */
    return 0;
}

/*===========================================================================*/
/* L3: LVDT                                                                  */
/*===========================================================================*/

double lvdt_output(double position, double sensitivity, double V_exc)
{
    /*
     * LVDT output voltage:
     *
     *   V_out = sensitivity · position · V_exc
     *
     * LVDT sensitivity is typically specified as mV/V/mm or mV/V/0.001in.
     *
     * Example: sensitivity = 30 mV/V/mm, V_exc = 3V, position = 1mm
     *   → V_out = 30·10⁻³ · 1 · 3 = 90 mV
     *
     * LVDTs have theoretically infinite resolution and are frictionless,
     * making them ideal for precision servo applications.
     *
     * Ref: Nyce (2004) §4.2
     */
    return sensitivity * position * V_exc;
}

double lvdt_null_voltage(double quadrature_factor, double V_exc)
{
    /*
     * LVDT null (residual) voltage:
     *
     *   V_null = quadrature_factor · V_exc
     *
     * At zero position, an ideal LVDT outputs zero. In practice,
     * capacitive coupling and core imbalances produce a residual
     * quadrature voltage (90° out of phase with normal output).
     *
     * Typical null: 0.1–1% of full-scale output.
     *
     * Ref: Nyce (2004) §4.2.3
     */
    return quadrature_factor * V_exc;
}

/*===========================================================================*/
/* L3: Strain Gauge                                                          */
/*===========================================================================*/

double strain_gauge_resistance_change(double R_nominal, double GF,
                                       double strain)
{
    /*
     * Strain gauge fundamental equation:
     *
     *   ΔR / R = GF · ε
     *
     * where GF (Gauge Factor) = (ΔR/R) / (ΔL/L)
     *
     * For metal foil gauges: GF ≈ 2.0
     * For semiconductor gauges: GF ≈ 50–200
     *
     * Piezoresistive effect: resistance changes due to both geometric
     * deformation and strain-induced resistivity change.
     *
     * Ref: Window "Strain Gauge Technology" (1992)
     */
    return R_nominal * GF * strain;
}

double wheatstone_quarter_bridge(double V_exc, double GF, double strain)
{
    /*
     * Wheatstone bridge, quarter-bridge configuration:
     *
     *   V_out = V_exc · (ΔR / (4R + 2ΔR))
     *
     * For small strains (ΔR << R), linearized:
     *   V_out ≈ V_exc · GF · ε / 4
     *
     * Non-linearity error for quarter bridge: ~GF·ε/2 (≈0.1% for 1000 με)
     *
     * Ref: Window (1992) §3.2
     */
    double dR_over_R = GF * strain;  /* ΔR/R */
    double denom = 4.0 + 2.0 * dR_over_R;
    if (fabs(denom) < 1e-30) return 0.0;
    return V_exc * dR_over_R / denom;
}

double wheatstone_half_bridge(double V_exc, double GF, double strain)
{
    /*
     * Half-bridge (two active gauges, opposite strain):
     *
     *   V_out = V_exc · GF · ε / 2
     *
     * Advantages over quarter bridge:
     *   - Twice the sensitivity
     *   - Inherent temperature compensation (both gauges at same temp)
     *   - No non-linearity error (for equal-and-opposite strains)
     *
     * Ref: Window (1992) §3.3
     */
    return V_exc * GF * strain / 2.0;
}

double wheatstone_full_bridge(double V_exc, double GF, double strain)
{
    /*
     * Full bridge (four active gauges):
     *
     *   V_out = V_exc · GF · ε
     *
     * All four arms active: two in tension, two in compression.
     * Maximum sensitivity and best temperature compensation.
     *
     * This is the preferred configuration for load cells and
     * precision force measurement.
     *
     * Ref: Window (1992) §3.4
     */
    return V_exc * GF * strain;
}

double load_cell_force(double V_out, double V_exc, double GF,
                        double A, double E_modulus)
{
    /*
     * Load cell calibration: convert bridge voltage to force.
     *
     * From full bridge: V_out/V_exc = GF · ε
     * Stress-strain: σ = E · ε, σ = F/A
     *
     * → F = (V_out / V_exc) · (A · E / GF)
     *
     * For aluminum load cell (E ≈ 70 GPa), A = 1 cm², GF = 2:
     *   F[N] = V_out/V_exc · 3.5e6
     *
     * Ref: Window (1992) §7.1
     */
    if (fabs(V_exc) < 1e-30 || GF < 1e-30) return 0.0;
    double strain_from_v = V_out / (V_exc * GF);
    double stress = E_modulus * strain_from_v;
    return stress * A;
}

/*===========================================================================*/
/* L5: Sensor Fusion                                                         */
/*===========================================================================*/

double complementary_filter(double gyro_rate, double accel_angle,
                             double prev_angle, double dt, double tau)
{
    /*
     * Complementary filter for IMU orientation (1-axis):
     *
     *   θ[k] = α · (θ[k-1] + ω_gyro · dt) + (1-α) · θ_accel
     *
     * where α = τ / (τ + dt)
     *
     * Principle:
     *   - Integrate gyro for high-frequency response (no delay, but drifts)
     *   - Use accelerometer for low-frequency correction (no drift, but noisy)
     *
     * The crossover frequency is 1/(2πτ):
     *   Below crossover: trust accelerometer more
     *   Above crossover: trust gyroscope more
     *
     * Typical τ: 0.5–5 seconds for human motion tracking
     *
     * Ref: Mahony et al. "Nonlinear Complementary Filters on the
     *      Special Orthogonal Group" (2008) IEEE Trans. AC
     */
    if (dt <= 0.0) return prev_angle;

    double alpha = tau / (tau + dt);
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;

    double gyro_angle = prev_angle + gyro_rate * dt;

    return alpha * gyro_angle + (1.0 - alpha) * accel_angle;
}

void kalman_filter_1d_position(double pos_measured, double acceleration,
                                double dt,
                                double *pos_est, double *vel_est,
                                double *P,
                                double Q_pos, double Q_vel, double R_meas)
{
    /*
     * Simplified 1D Kalman filter for position-velocity estimation.
     *
     * State: x = [p, v]^T
     * Measurement: z = p_measured (position only)
     *
     * Predict step:
     *   p = p + v·dt + 0.5·a·dt²
     *   v = v + a·dt
     *   P = F·P·F^T + Q
     *   where F = [1  dt]
     *             [0   1]
     *         Q = [Q_pos·dt³/3  Q_pos·dt²/2]
     *             [Q_pos·dt²/2  Q_pos·dt    ]
     *
     * Update step (position measurement only):
     *   H = [1, 0]
     *   S = H·P·H^T + R
     *   K = P·H^T / S
     *   p = p + K[0]·(z - p)
     *   v = v + K[1]·(z - p)
     *   P = (I - K·H)·P
     *
     * Ref: Kalman "A New Approach to Linear Filtering and
     *      Prediction Problems" (1960) Trans. ASME
     */
    if (!pos_est || !vel_est || !P) return;

    /* ---- Predict ---- */
    double p_prior = *pos_est + (*vel_est) * dt + 0.5 * acceleration * dt * dt;
    double v_prior = *vel_est + acceleration * dt;

    /* Covariance predict: P = F·P·F^T + Q */
    /* F = [[1, dt], [0, 1]] */
    double F00 = 1.0, F01 = dt;
    double F10 = 0.0, F11 = 1.0;

    /* P_pred = F·P·F^T */
    double P00 = F00*P[0] + F01*P[2];  /* P[0]*F00 + P[2]*F01 */
    double P01_temp = F00*P[1] + F01*P[3];
    double P10_temp = F10*P[0] + F11*P[2];
    double P11_temp = F10*P[1] + F11*P[3];

    double P00_pred = P00*F00 + P01_temp*F01;
    double P01_pred = P00*F10 + P01_temp*F11;
    double P10_pred = P10_temp*F00 + P11_temp*F01;
    double P11_pred = P10_temp*F10 + P11_temp*F11;

    /* Add process noise Q */
    double Q00 = Q_pos * dt * dt * dt / 3.0 + Q_vel * dt * dt * dt / 3.0;
    double Q01 = Q_pos * dt * dt / 2.0;
    double Q11 = Q_pos * dt;

    P[0] = P00_pred + Q00;
    P[1] = P01_pred + Q01;
    P[2] = P10_pred + Q01;
    P[3] = P11_pred + Q11;

    /* ---- Update (position measurement) ---- */
    /* H = [1, 0] */
    /* S = H·P·H^T + R = P[0] + R_meas */
    double S = P[0] + R_meas;
    if (fabs(S) < 1e-30) return;

    /* K = P·H^T / S */
    double K0 = P[0] / S;
    double K1 = P[2] / S;

    /* Innovation: z - H·x_prior */
    double y_tilde = pos_measured - p_prior;

    /* State update */
    *pos_est = p_prior + K0 * y_tilde;
    *vel_est = v_prior + K1 * y_tilde;

    /* Covariance update: P = (I - K·H)·P */
    double P00_new = (1.0 - K0) * P[0];      /* (1-K0)*P00 + (-K0)*P20(0) */
    double P01_new = (1.0 - K0) * P[1];
    double P10_new = -K1 * P[0] + 1.0 * P[2];  /* (-K1)*P00 + 1*P20 */
    double P11_new = -K1 * P[1] + 1.0 * P[3];

    P[0] = P00_new;
    P[1] = P01_new;
    P[2] = P10_new;
    P[3] = P11_new;
}
