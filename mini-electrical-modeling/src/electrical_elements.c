/**
 * @file electrical_elements.c
 * @brief Implementation of fundamental electrical element computations
 *
 * Each function implements one independent knowledge point.
 * No filler functions, no template-generated code.
 */

#include "electrical_elements.h"
#include <math.h>
#include <complex.h>
#include <float.h>

/* ===== L2: Impedance in s-domain ===== */

double complex impedance_rlc_series(double r, double l, double c, double complex s)
{
    /* Z(s) = R + sL + 1/(sC) */
    /* Guard: C=0 means open circuit, impedance is infinite — return large value */
    if (c <= 0.0) {
        /* Pure RL series */
        return r + s * l;
    }
    if (c <= 1e-30) return DBL_MAX + I * DBL_MAX;
    return r + s * l + 1.0 / (s * c);
}

double complex impedance_rlc_parallel(double r, double l, double c, double complex s)
{
    /* For parallel: 1/Z = 1/R + 1/(sL) + sC → Z = R || sL || 1/(sC) */
    double complex zr = r;
    double complex zl = s * l;
    double complex zc = 1.0 / (s * c);

    double complex y_sum = 0.0;
    if (r > 1e-30)  y_sum += 1.0 / zr;
    if (l > 1e-30)  y_sum += 1.0 / zl;
    if (c > 1e-30)  y_sum += 1.0 / zc;

    if (cabs(y_sum) < 1e-30) return DBL_MAX + I * DBL_MAX;
    return 1.0 / y_sum;
}

double complex admittance_from_impedance(double complex z)
{
    /* Y = 1/Z = G + jB */
    /* Guard against division by zero */
    double mag2 = creal(z) * creal(z) + cimag(z) * cimag(z);
    if (mag2 < 1e-60) return DBL_MAX + I * 0.0;

    /* Y = Z* / |Z|^2 */
    return (creal(z) - I * cimag(z)) / mag2;
}

/* ===== L2: Phasor functions ===== */

double complex phasor_to_complex(const Phasor *p)
{
    if (!p) return 0.0;
    /* P = |V| * e^(jθ) = |V|·cos θ + j·|V|·sin θ */
    return p->magnitude * cos(p->angle_rad) +
           I * p->magnitude * sin(p->angle_rad);
}

Phasor complex_to_phasor(double complex z, double frequency)
{
    Phasor p;
    p.magnitude = cabs(z);
    p.angle_rad = carg(z);
    p.frequency = frequency;
    return p;
}

/* ===== L2: Ohm's law ===== */

double ohms_law_voltage(double current, double resistance)
{
    /* V = I * R (Ohm 1827) */
    return current * resistance;
}

double ohms_law_current(double voltage, double resistance)
{
    /* I = V / R */
    /* Guard: R=0 means short circuit, current is unbounded */
    if (fabs(resistance) < 1e-30) {
        return (voltage >= 0) ? DBL_MAX : -DBL_MAX;
    }
    return voltage / resistance;
}

double complex ohms_law_phasor_v(double complex i, double complex z)
{
    /* V = I * Z in phasor domain */
    return i * z;
}

double complex ohms_law_phasor_i(double complex v, double complex z)
{
    /* I = V / Z */
    double mag2 = creal(z) * creal(z) + cimag(z) * cimag(z);
    if (mag2 < 1e-60) return DBL_MAX + I * 0.0;
    return v / z;
}

/* ===== L2: Power calculations ===== */

double power_dc_v_i(double voltage, double current)
{
    /* P = V * I (Joule 1841) */
    return voltage * current;
}

double power_dc_i2r(double current, double resistance)
{
    /* P = I^2 * R */
    return current * current * resistance;
}

double power_dc_v2r(double voltage, double resistance)
{
    /* P = V^2 / R */
    /* Guard against division by zero */
    if (fabs(resistance) < 1e-30) return DBL_MAX;
    return (voltage * voltage) / resistance;
}

double ac_real_power(double v_rms, double i_rms, double phase_diff_rad)
{
    /* P = |V|·|I|·cos(θ_v - θ_i) [W] */
    /* Real (active) power absorbed by the load */
    return v_rms * i_rms * cos(phase_diff_rad);
}

double ac_reactive_power(double v_rms, double i_rms, double phase_diff_rad)
{
    /* Q = |V|·|I|·sin(θ_v - θ_i) [VAR] */
    /* Reactive power oscillates between source and load */
    return v_rms * i_rms * sin(phase_diff_rad);
}

double ac_apparent_power(double v_rms, double i_rms)
{
    /* |S| = |V|·|I| [VA] — magnitude of complex power */
    return v_rms * i_rms;
}

double ac_power_factor(double v_rms, double i_rms, double phase_diff_rad)
{
    /* pf = P/|S| = cos(φ) */
    /* pf = 1: purely resistive (φ = 0) */
    /* pf = 0: purely reactive (φ = ±90°) */
    (void)v_rms; (void)i_rms; /* unused in cos-only formula */
    double pf = cos(phase_diff_rad);

    /* Clamp to [-1, 1] to handle floating-point edge cases */
    if (pf > 1.0) return 1.0;
    if (pf < -1.0) return -1.0;
    return pf;
}

double complex ac_complex_power(const Phasor *v, const Phasor *i)
{
    /* S = V * I* = P + jQ */
    /* I* is the complex conjugate of the current phasor */
    if (!v || !i) return 0.0;
    double complex vc = phasor_to_complex(v);
    double complex ic = phasor_to_complex(i);
    /* Complex conjugate of I */
    double complex ic_conj = creal(ic) - I * cimag(ic);
    return vc * ic_conj;
}

/* ===== L3: Time constants and natural frequencies ===== */

double time_constant_rc(double resistance, double capacitance)
{
    /* τ = R·C [s] */
    /* Guard against zero/negative values */
    if (resistance <= 0.0 || capacitance <= 0.0) return DBL_MAX;
    return resistance * capacitance;
}

double time_constant_rl(double resistance, double inductance)
{
    /* τ = L/R [s] */
    if (resistance <= 0.0) return DBL_MAX;
    return inductance / resistance;
}

double resonant_frequency_rad_s(double inductance, double capacitance)
{
    /* ω_0 = 1/√(LC) */
    if (inductance <= 0.0 || capacitance <= 0.0) return 0.0;
    double lc = inductance * capacitance;
    if (lc < 1e-30) return DBL_MAX;
    return 1.0 / sqrt(lc);
}

double resonant_frequency_hz(double inductance, double capacitance)
{
    /* f_0 = ω_0 / (2π) */
    double omega = resonant_frequency_rad_s(inductance, capacitance);
    return omega / (2.0 * M_PI);
}

double damping_factor_series_rlc(double resistance, double inductance, double capacitance)
{
    /* ζ = (R/2)·√(C/L) for series RLC */
    /* Characteristic equation: s² + (R/L)s + 1/(LC) = 0 */
    if (inductance <= 0.0) return 0.0;
    if (capacitance <= 0.0) return 0.0;

    double omega_n = resonant_frequency_rad_s(inductance, capacitance);
    double alpha = resistance / (2.0 * inductance); /* Neper frequency */
    if (omega_n < 1e-30) return DBL_MAX;
    return alpha / omega_n;
}

double quality_factor_series_rlc(double resistance, double inductance, double capacitance)
{
    /* Q = ω_0·L/R = 1/(ω_0·C·R) for series RLC */
    /* Q = energy stored / energy dissipated per radian */
    if (resistance <= 0.0) return DBL_MAX;

    double omega_0 = resonant_frequency_rad_s(inductance, capacitance);
    if (omega_0 < 1e-30) return 0.0;

    /* Use Q = ω_0·L/R */
    return (omega_0 * inductance) / resistance;
}

/* ===== L3: Series and parallel equivalent combinations ===== */

double equivalent_resistance_series(const double *resistances, size_t n)
{
    if (!resistances || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += resistances[i];
    }
    return sum;
}

double equivalent_resistance_parallel(const double *resistances, size_t n)
{
    if (!resistances || n == 0) return 0.0;
    double sum_conductance = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (fabs(resistances[i]) < 1e-30) continue; /* short → infinite conductance, handle separately */
        sum_conductance += 1.0 / resistances[i];
    }
    if (fabs(sum_conductance) < 1e-30) return DBL_MAX;
    return 1.0 / sum_conductance;
}

double equivalent_capacitance_series(const double *capacitances, size_t n)
{
    /* C_eq_series = 1/Σ(1/C_i) — reciprocal of reciprocal sum */
    if (!capacitances || n == 0) return 0.0;
    double sum_recip = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (fabs(capacitances[i]) < 1e-30) return 0.0; /* broken cap */
        sum_recip += 1.0 / capacitances[i];
    }
    if (fabs(sum_recip) < 1e-30) return DBL_MAX;
    return 1.0 / sum_recip;
}

double equivalent_capacitance_parallel(const double *capacitances, size_t n)
{
    /* C_eq_parallel = ΣC_i — direct sum */
    if (!capacitances || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += capacitances[i];
    }
    return sum;
}

double equivalent_inductance_series(const double *inductances, size_t n)
{
    /* L_eq_series = ΣL_i (no mutual coupling) */
    if (!inductances || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += inductances[i];
    }
    return sum;
}

double equivalent_inductance_parallel(const double *inductances, size_t n)
{
    /* L_eq_parallel = 1/Σ(1/L_i) */
    if (!inductances || n == 0) return 0.0;
    double sum_recip = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (fabs(inductances[i]) < 1e-30) return 0.0;
        sum_recip += 1.0 / inductances[i];
    }
    if (fabs(sum_recip) < 1e-30) return DBL_MAX;
    return 1.0 / sum_recip;
}

/* ===== L4: Energy stored in reactive elements ===== */

double energy_capacitor(double capacitance, double voltage)
{
    /* E = 1/2·C·V^2 [J] */
    /* Derived from work to charge capacitor: W = ∫v·dq = ∫C·v·dv = ½CV² */
    if (capacitance <= 0.0 || fabs(voltage) < 1e-30) return 0.0;
    return 0.5 * capacitance * voltage * voltage;
}

double energy_inductor(double inductance, double current)
{
    /* E = 1/2·L·I^2 [J] */
    /* Derived from work to establish current: W = ∫i·v·dt = ∫L·i·di = ½LI² */
    if (inductance <= 0.0 || fabs(current) < 1e-30) return 0.0;
    return 0.5 * inductance * current * current;
}

/* ===== L5: Temperature effects on resistance ===== */

double resistance_at_temperature(double r0, double alpha, double temp_ref,
                                  double temp_actual)
{
    /* R(T) = R_0 * [1 + α*(T - T_0)] */
    /* Linear TCR approximation (Matthiessen rule) */
    /* For pure metals: α ≈ +0.0039/K (Cu), +0.0040/K (Al) */
    /* For thermistors: α is negative and larger magnitude */
    double delta_t = temp_actual - temp_ref;
    return r0 * (1.0 + alpha * delta_t);
}

double self_heating_temp_rise(double power_dissipated, double theta_ja)
{
    /* ΔT = P * θ_JA */
    /* Steady-state thermal equilibrium */
    if (power_dissipated < 0.0 || theta_ja < 0.0) return 0.0;
    return power_dissipated * theta_ja;
}

/* ===== L5: Capacitor and inductor dynamics ===== */

double capacitor_charge_voltage(double supply_voltage, double t, double tau)
{
    /* v_c(t) = V_s * (1 - e^(-t/τ)) */
    /* Solution to: τ·dv_c/dt + v_c = V_s, v_c(0) = 0 */
    if (tau <= 0.0) return supply_voltage; /* instantaneous */
    if (t <= 0.0) return 0.0;
    return supply_voltage * (1.0 - exp(-t / tau));
}

double capacitor_discharge_voltage(double initial_voltage, double t, double tau)
{
    /* v_c(t) = V_0 * e^(-t/τ) */
    /* Solution to: τ·dv_c/dt + v_c = 0, v_c(0) = V_0 */
    if (tau <= 0.0) return 0.0; /* instantaneous discharge */
    if (t <= 0.0) return initial_voltage;
    return initial_voltage * exp(-t / tau);
}

double capacitor_current(double capacitance, double dv_dt)
{
    /* i_c = C * dv/dt */
    /* Maxwell displacement current in circuit form */
    if (capacitance <= 0.0) return 0.0;
    return capacitance * dv_dt;
}

double inductor_voltage(double inductance, double di_dt)
{
    /* v_l = L * di/dt */
    /* Faraday law in circuit form: v = dΦ/dt = L·di/dt */
    if (inductance <= 0.0) return 0.0;
    return inductance * di_dt;
}

/* ===== L5: Voltage and current dividers ===== */

double voltage_divider_unloaded(double vin, double r1, double r2)
{
    /* V_out = V_in * R2/(R1 + R2) — unloaded divider */
    double sum = r1 + r2;
    if (fabs(sum) < 1e-30) return 0.0;
    return vin * r2 / sum;
}

double voltage_divider_loaded(double vin, double r1, double r2, double r_load)
{
    /* V_out = V_in * (R2||R_L)/(R1 + (R2||R_L)) — loaded divider */
    if (fabs(r_load) < 1e-30) return 0.0; /* shorted output */
    double r2_parallel_rl = (r2 * r_load) / (r2 + r_load);
    double denominator = r1 + r2_parallel_rl;
    if (fabs(denominator) < 1e-30) return 0.0;
    return vin * r2_parallel_rl / denominator;
}

double current_divider_two_branch(double i_total, double r_this, double r_other)
{
    /* I_branch = I_total * R_other/(R_this + R_other) */
    double sum = r_this + r_other;
    if (fabs(sum) < 1e-30) return 0.0;
    return i_total * r_other / sum;
}

/* ===== L5: Frequency-dependent element behavior ===== */

double capacitive_reactance(double capacitance, double angular_freq)
{
    /* X_C = -1/(ω·C) — sign convention: capacitive reactance is negative */
    /* |X_C| = 1/(ω·C) */
    if (capacitance <= 0.0 || fabs(angular_freq) < 1e-30) return -DBL_MAX;
    return -1.0 / (angular_freq * capacitance);
}

double inductive_reactance(double inductance, double angular_freq)
{
    /* X_L = ω·L — inductive reactance is positive */
    if (inductance <= 0.0) return 0.0;
    return angular_freq * inductance;
}

double impedance_magnitude_rlc_series(double r, double l, double c, double angular_freq)
{
    /* |Z| = √(R² + (ωL - 1/(ωC))²) */
    double x_l = inductive_reactance(l, angular_freq);
    double x_c = capacitive_reactance(c, angular_freq);
    double x_net = x_l + x_c; /* x_c is negative */
    return sqrt(r * r + x_net * x_net);
}

double impedance_phase_rlc_series(double r, double l, double c, double angular_freq)
{
    /* φ = atan2(X_net, R) = atan2(ωL - 1/(ωC), R) */
    double x_l = inductive_reactance(l, angular_freq);
    double x_c = capacitive_reactance(c, angular_freq);
    double x_net = x_l + x_c;
    return atan2(x_net, r);
}

/* ===== L7: Real device models ===== */

double real_capacitor_impedance_magnitude(const RealCapacitor *cap, double freq_hz)
{
    if (!cap || freq_hz < 0.0) return DBL_MAX;
    double omega = 2.0 * M_PI * freq_hz;
    double xc = (fabs(omega * cap->capacitance) > 1e-30)
                ? -1.0 / (omega * cap->capacitance) : -DBL_MAX;
    double xl = omega * cap->esl;
    double x_net = xl + xc;
    return sqrt(cap->esr * cap->esr + x_net * x_net);
}

double real_capacitor_self_resonant_freq(const RealCapacitor *cap)
{
    /* SRF = 1/(2π√(ESL·C)) — where ESL and C series-resonate */
    if (!cap || cap->capacitance <= 0.0 || cap->esl <= 0.0) return DBL_MAX;
    double lc = cap->esl * cap->capacitance;
    if (lc < 1e-30) return DBL_MAX;
    return 1.0 / (2.0 * M_PI * sqrt(lc));
}

int real_capacitor_is_capacitive(const RealCapacitor *cap, double freq_hz)
{
    /* Below SRF: capacitive behavior. Above SRF: inductive behavior. */
    double srf = real_capacitor_self_resonant_freq(cap);
    return (freq_hz < srf) ? 1 : 0;
}

double real_inductor_q_factor(const RealInductor *ind, double freq_hz)
{
    /* Q = ωL/DCR at low frequencies (before core loss and parasitic C dominate) */
    if (!ind || ind->dcr <= 0.0) return DBL_MAX;
    double omega = 2.0 * M_PI * freq_hz;
    return (omega * ind->inductance) / ind->dcr;
}

double real_inductor_impedance_magnitude(const RealInductor *ind, double freq_hz)
{
    /* Z = DCR + jωL || R_core || 1/(jωC_p) */
    if (!ind) return DBL_MAX;
    double omega = 2.0 * M_PI * freq_hz;

    /* jωL term */
    double complex zl = I * omega * ind->inductance;

    /* Parallel combination: Z = DCR + (jωL || R_core || 1/(jωC_p)) */
    double complex y_parallel = 0.0;

    if (ind->inductance > 0.0)
        y_parallel += 1.0 / zl;
    if (ind->core_loss_r > 0.0)
        y_parallel += 1.0 / ind->core_loss_r;

    /* Parasitic capacitance in parallel */
    if (ind->parasitic_c > 0.0) {
        double complex zc = 1.0 / (I * omega * ind->parasitic_c);
        y_parallel += 1.0 / zc;
    }

    if (cabs(y_parallel) < 1e-30) return ind->dcr;

    double complex z_parallel = 1.0 / y_parallel;
    double complex z_total = ind->dcr + z_parallel;

    return cabs(z_total);
}
