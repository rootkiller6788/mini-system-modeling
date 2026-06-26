/**
 * @file electrical_elements.h
 * @brief Fundamental electrical element definitions R, L, C, sources
 *
 * Knowledge Coverage:
 *   L1 - Definitions: resistance, capacitance, inductance, impedance, admittance,
 *        time constant, resonance frequency, Q-factor, phasor, complex power
 *   L2 - Core Concepts: Ohms law, energy storage in L/C, series/parallel combinations,
 *        AC steady-state, power factor, reactance, susceptance
 *   L3 - Mathematical Structures: complex impedance, phasor representation,
 *        admittance triangle, s-domain element models
 *   L4 - Fundamental Laws: Ohms law, Joules law, energy conservation in RLC
 *
 * Reference:
 *   MIT 6.302 Feedback Systems - electrical network modeling
 *   Stanford ENGR105 - circuit element s-domain models
 *   Hayt and Buck, Engineering Electromagnetics (2018)
 *   Dorf and Svoboda, Introduction to Electric Circuits (2018)
 *
 * All functions documented with physical theorem source and O(·) complexity.
 */

#ifndef ELECTRICAL_ELEMENTS_H
#define ELECTRICAL_ELEMENTS_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* L1: Core Data Types */

/** Resistance element [ohm] dissipates power via P = I^2 R (Joules law) */
typedef struct {
    double resistance;
    double power_rating;
    double temperature_coeff;
    double tolerance;
} Resistor;

/** Capacitance element [F = farad] stores energy via E = 1/2 C V^2 */
typedef struct {
    double capacitance;
    double voltage_rating;
    double esr;
    double leakage_current;
    double temperature_coeff;
} Capacitor;

/** Inductance element [H = henry] stores energy via E = 1/2 L I^2 */
typedef struct {
    double inductance;
    double current_rating;
    double dcr;
    double saturation_current;
    double core_material;
} Inductor;

/* L1: Source definitions */
typedef enum {
    SOURCE_DC = 0,
    SOURCE_AC = 1,
    SOURCE_PULSE = 2,
    SOURCE_PWL = 3,
    SOURCE_ARBITRARY = 4
} SourceType;

typedef struct {
    SourceType type;
    double dc_value;
    double amplitude;
    double frequency;
    double phase;
    double rise_time;
    double fall_time;
    double pulse_width;
    double period;
    double (*arbitrary_fn)(double t);
} ElectricalSource;

/* L2: Phasor representation */
typedef struct {
    double magnitude;
    double angle_rad;
    double frequency;
} Phasor;

/* L2: Impedance in s-domain */
/* Z(s) = R + sL + 1/(sC) for series RLC */
double complex impedance_rlc_series(double r, double l, double c, double complex s);
double complex impedance_rlc_parallel(double r, double l, double c, double complex s);
double complex admittance_from_impedance(double complex z);

/* L2: Phasor functions */
double complex phasor_to_complex(const Phasor *p);
Phasor complex_to_phasor(double complex z, double frequency);

/* L2: Ohms law */
double ohms_law_voltage(double current, double resistance);
double ohms_law_current(double voltage, double resistance);
double complex ohms_law_phasor_v(double complex i, double complex z);
double complex ohms_law_phasor_i(double complex v, double complex z);

/* L2: Power calculations */
double power_dc_v_i(double voltage, double current);
double power_dc_i2r(double current, double resistance);
double power_dc_v2r(double voltage, double resistance);
double ac_real_power(double v_rms, double i_rms, double phase_diff_rad);
double ac_reactive_power(double v_rms, double i_rms, double phase_diff_rad);
double ac_apparent_power(double v_rms, double i_rms);
double ac_power_factor(double v_rms, double i_rms, double phase_diff_rad);
double complex ac_complex_power(const Phasor *v, const Phasor *i);

/* L3: Time constants and natural frequency */
double time_constant_rc(double resistance, double capacitance);
double time_constant_rl(double resistance, double inductance);
double resonant_frequency_rad_s(double inductance, double capacitance);
double resonant_frequency_hz(double inductance, double capacitance);
double damping_factor_series_rlc(double resistance, double inductance, double capacitance);
double quality_factor_series_rlc(double resistance, double inductance, double capacitance);

/* L3: Series and parallel equivalents */
double equivalent_resistance_series(const double *resistances, size_t n);
double equivalent_resistance_parallel(const double *resistances, size_t n);
double equivalent_capacitance_series(const double *capacitances, size_t n);
double equivalent_capacitance_parallel(const double *capacitances, size_t n);
double equivalent_inductance_series(const double *inductances, size_t n);
double equivalent_inductance_parallel(const double *inductances, size_t n);

/* L4: Energy storage */
double energy_capacitor(double capacitance, double voltage);
double energy_inductor(double inductance, double current);

/* L5: Temperature effects */
double resistance_at_temperature(double r0, double alpha, double temp_ref, double temp_actual);
double self_heating_temp_rise(double power_dissipated, double theta_ja);

/* L5: Charge/discharge dynamics */
double capacitor_charge_voltage(double supply_voltage, double t, double tau);
double capacitor_discharge_voltage(double initial_voltage, double t, double tau);
double capacitor_current(double capacitance, double dv_dt);
double inductor_voltage(double inductance, double di_dt);

/* L5: Dividers */
double voltage_divider_unloaded(double vin, double r1, double r2);
double voltage_divider_loaded(double vin, double r1, double r2, double r_load);
double current_divider_two_branch(double i_total, double r_this, double r_other);

/* L5: Frequency-dependent behavior */
double capacitive_reactance(double capacitance, double angular_freq);
double inductive_reactance(double inductance, double angular_freq);
double impedance_magnitude_rlc_series(double r, double l, double c, double angular_freq);
double impedance_phase_rlc_series(double r, double l, double c, double angular_freq);

/* L7: Real device models */
typedef struct {
    double capacitance;
    double esr;
    double esl;
    double leakage_r;
    double dielectric_absorption;
} RealCapacitor;

double real_capacitor_impedance_magnitude(const RealCapacitor *cap, double freq_hz);
double real_capacitor_self_resonant_freq(const RealCapacitor *cap);
int    real_capacitor_is_capacitive(const RealCapacitor *cap, double freq_hz);

typedef struct {
    double inductance;
    double dcr;
    double core_loss_r;
    double parasitic_c;
    double saturation_current;
} RealInductor;

double real_inductor_q_factor(const RealInductor *ind, double freq_hz);
double real_inductor_impedance_magnitude(const RealInductor *ind, double freq_hz);

#endif /* ELECTRICAL_ELEMENTS_H */
