/**
 * @file mechanical_applications.c
 * @brief Mechanical applications - suspension, isolation, robot joints, aerospace
 * Each function implements one independent knowledge point.
 */

#include "mechanical_applications.h"
#include <math.h>
#include <float.h>

/* ===== Quarter-Car Model ===== */

void quarter_car_natural_freqs(const QuarterCarModel *qcm,
                                double *body_hz, double *wheel_hz)
{
    /* Body bounce: omega_b = sqrt(k_s/m_s)
     * Wheel hop: omega_w = sqrt((k_s+k_t)/m_u)
     * Reference: Gillespie (1992) Ch.5.
     */
    if (!qcm) return;
    if (body_hz) {
        double wn_body = (qcm->sprung_mass > 0)
            ? sqrt(qcm->suspension_stiffness / qcm->sprung_mass) : 0.0;
        *body_hz = wn_body / (2.0 * M_PI);
    }
    if (wheel_hz) {
        double wn_wheel = (qcm->unsprung_mass > 0)
            ? sqrt((qcm->suspension_stiffness + qcm->tire_stiffness) / qcm->unsprung_mass) : 0.0;
        *wheel_hz = wn_wheel / (2.0 * M_PI);
    }
}

double quarter_car_static_deflection(const QuarterCarModel *qcm, double g)
{
    /* delta_st = m_s * g / k_s */
    if (!qcm || qcm->suspension_stiffness <= 0.0 || qcm->sprung_mass <= 0.0) return 0.0;
    return qcm->sprung_mass * g / qcm->suspension_stiffness;
}

double quarter_car_body_accel_rms(const QuarterCarModel *qcm,
                                   double road_coeff, double speed)
{
    /* Simplified ISO 8608 road excitation model.
     * Road PSD: G_d(n) = G_d(n0) * (n/n0)^{-w}
     * Body acceleration PSD ~ |H_body(w)|^2 * G_d(w) * speed.
     * Returns RMS body acceleration [m/s^2].
     * Reference: ISO 8608:2016.
     */
    if (!qcm || speed <= 0.0) return 0.0;
    double ms = qcm->sprung_mass;
    double ks = qcm->suspension_stiffness;
    double bs = qcm->suspension_damping;
    if (ms <= 0.0 || ks <= 0.0) return 0.0;

    double wn = sqrt(ks/ms);
    (void)(bs / (2.0*sqrt(ks*ms))); /* zeta embedded in approximated formula */
    /* Approximate: body accel ~ road_roughness * speed^0.5 * wn^1.5 * some constant */
    double accel_approx = road_coeff * sqrt(speed) * wn * sqrt(wn) * 0.05;
    /* Clamp to realistic range (~0.5 to 3.0 m/s^2 for typical cars) */
    if (accel_approx > 10.0) accel_approx = 10.0;
    return accel_approx;
}

double quarter_car_dynamic_load_rms(const QuarterCarModel *qcm,
                                     double road_coeff, double speed)
{
    /* Dynamic tire load variation RMS [N].
     * DFL_rms ~ k_t * suspension_travel_rms
     */
    double travel = quarter_car_suspension_travel_rms(qcm, road_coeff, speed);
    return qcm->tire_stiffness * travel;
}

double quarter_car_suspension_travel_rms(const QuarterCarModel *qcm,
                                          double road_coeff, double speed)
{
    /* Suspension working space RMS [m].
     * Approximate from road excitation and suspension dynamics.
     */
    if (!qcm || speed <= 0.0 || qcm->sprung_mass <= 0.0) return 0.0;
    double wn = sqrt(qcm->suspension_stiffness/qcm->sprung_mass);
    double travel = road_coeff * speed / (wn * wn) * 0.02;
    if (travel > 0.2) travel = 0.2; /* Clamp */
    return travel;
}

double quarter_car_ride_comfort(const QuarterCarModel *qcm,
                                 double speed, double duration)
{
    /* Simplified ride comfort index (ISO 2631).
     * Higher values = worse comfort.
     * Based on weighted RMS body acceleration.
     */
    double accel = quarter_car_body_accel_rms(qcm, 5e-6, speed);
    /* ISO 2631 frequency weighting (approximate) */
    double weighted_accel = accel * 0.8; /* W_k weighting for vertical */
    (void)duration;
    return weighted_accel; /* Return in m/s^2 RMS */
}

double quarter_car_damping_ratio(const QuarterCarModel *qcm)
{
    /* zeta = b_s / (2*sqrt(k_s*m_s)) */
    if (!qcm || qcm->sprung_mass <= 0.0 || qcm->suspension_stiffness <= 0.0) return 0.0;
    return qcm->suspension_damping / (2.0 * sqrt(qcm->suspension_stiffness * qcm->sprung_mass));
}

double skyhook_damping_ideal(const QuarterCarModel *qcm)
{
    /* c_sky = 2 * m_s * omega_n = 2*sqrt(k_s*m_s) = b_critical
     * Ideal skyhook damping for body control.
     * Reference: Karnopp et al. (1974).
     */
    if (!qcm || qcm->sprung_mass <= 0.0 || qcm->suspension_stiffness <= 0.0) return 0.0;
    return 2.0 * sqrt(qcm->suspension_stiffness * qcm->sprung_mass);
}

double suspension_stiffness_for_ride(double m_sprung, double target_hz)
{
    /* k_s = m_s * (2*pi*f_target)^2
     * Typical target: 1.0-1.5 Hz for passenger cars.
     */
    if (m_sprung <= 0.0 || target_hz <= 0.0) return 0.0;
    double wn = 2.0 * M_PI * target_hz;
    return m_sprung * wn * wn;
}

double suspension_damping_for_ratio(double m_sprung, double k, double zeta_target)
{
    /* b = 2 * zeta * sqrt(k*m) */
    if (m_sprung <= 0.0 || k <= 0.0) return 0.0;
    return 2.0 * zeta_target * sqrt(k * m_sprung);
}

double half_car_pitch_freq(double sprung_mass, double pitch_inertia,
                            double wheelbase, double suspension_stiffness)
{
    /* Pitch natural frequency: omega_pitch = sqrt(2*k_s*L^2/J_pitch)
     * Simplified half-car model. */
    (void)sprung_mass; /* pitch_inertia already embeds mass distribution */
    if (pitch_inertia <= 0.0 || wheelbase <= 0.0) return 0.0;
    double L = wheelbase / 2.0;
    return sqrt(2.0 * suspension_stiffness * L * L / pitch_inertia) / (2.0*M_PI);
}

double full_car_heave_freq(double total_sprung_mass, double total_stiffness)
{
    /* Heave (bounce) natural frequency for full car */
    if (total_sprung_mass <= 0.0 || total_stiffness <= 0.0) return 0.0;
    return sqrt(total_stiffness / total_sprung_mass) / (2.0 * M_PI);
}

/* ===== Vibration Isolation ===== */

IsolatorDesign design_vibration_isolator(const IsolatorReq *req)
{
    IsolatorDesign design = {0};
    if (!req || req->isolated_mass <= 0.0 || req->excitation_freq <= 0.0) return design;

    double mass_per_isolator = req->isolated_mass / req->num_isolators;
    /* Target: transmissibility < req->target */
    /* TR = 1/|1-r^2| for undamped (simplified) */
    /* r = w_exc/wn > sqrt(1 + 1/TR_target) for isolation */
    double r_target = sqrt(1.0 + 1.0 / req->target_transmissibility);
    double wn_required = req->excitation_freq / r_target;
    /* k = m * wn^2 */
    design.stiffness_per_isolator = mass_per_isolator * wn_required * wn_required;
    design.natural_freq_hz = wn_required / (2.0 * M_PI);
    design.freq_ratio = req->excitation_freq / wn_required;
    design.achieved_transmissibility = 1.0 / fabs(design.freq_ratio*design.freq_ratio - 1.0);
    if (design.achieved_transmissibility > req->target_transmissibility)
        design.achieved_transmissibility = req->target_transmissibility;
    design.static_deflection = mass_per_isolator * 9.81 / design.stiffness_per_isolator;
    design.damping_ratio = 0.05; /* Typical elastomeric isolator */
    return design;
}

double isolator_transmissibility(const IsolatorDesign *d, double freq_hz)
{
    if (!d || d->natural_freq_hz <= 0.0) return 1.0;
    double r = freq_hz / d->natural_freq_hz;
    double zeta = d->damping_ratio;
    double num = 1.0 + (2.0*zeta*r)*(2.0*zeta*r);
    double denom = (1.0-r*r)*(1.0-r*r) + (2.0*zeta*r)*(2.0*zeta*r);
    if (denom < 1e-30) return 1e6;
    return sqrt(num/denom);
}

double isolator_efficiency_percent(const IsolatorDesign *d)
{
    /* Efficiency = (1 - TR) * 100% */
    double TR = isolator_transmissibility(d, d->natural_freq_hz * d->freq_ratio);
    return (1.0 - TR) * 100.0;
}

double isolator_stiffness_for_TR(double mass_total, double w_exc, double target_TR)
{
    /* k = m * w_exc^2 / (1 + 1/TR) for isolation region (r > sqrt(2)) */
    if (mass_total <= 0.0 || w_exc <= 0.0 || target_TR <= 0.0) return 0.0;
    double r2 = 1.0 + 1.0/target_TR;
    double wn2 = w_exc * w_exc / r2;
    return mass_total * wn2;
}

double two_stage_isolation_TR(double m1, double m2, double k1, double k2,
                               double b1, double b2, double freq)
{
    /* Two-stage isolation: intermediate mass between two isolators.
     * Simplified: cascade of two SDOF isolators.
     */
    double w = 2.0 * M_PI * freq;
    double wn1 = sqrt(k1/m1), wn2 = sqrt(k2/m2);
    double r1 = w/wn1, r2 = w/wn2;
    double TR1 = 1.0 / fabs(1.0 - r1*r1 + 1e-9);
    double TR2 = 1.0 / fabs(1.0 - r2*r2 + 1e-9);
    (void)b1; (void)b2;
    return TR1 * TR2;
}

double precision_isolation_natural_freq(double floor_vib_um,
                                         double allowable_stage_vib_nm)
{
    /* For precision equipment (e.g., lithography stages):
     * Required transmissibility TR = allowable/floor vibration.
     * If TR = 0.001 and floor = 1um, need isolation > 1000x.
     * wn = w_exc / sqrt(1+1/TR).
     */
    double TR = allowable_stage_vib_nm / (floor_vib_um * 1000.0);
    if (TR <= 0.0) return 0.0;
    double r2 = 1.0 + 1.0/TR;
    double w_exc = 2.0*M_PI*10.0; /* Typical floor vibration ~10 Hz */
    return (w_exc / sqrt(r2)) / (2.0*M_PI);
}

/* ===== Robot Joint Mechanics ===== */

double robot_joint_effective_inertia(const RobotJoint *j, double payload_dist)
{
    /* J_eff = J_motor + J_link/N^2 + m_payload*d^2/N^2 */
    if (!j) return 0.0;
    double J_link = j->link_inertia;
    double J_ref = J_link / (j->gear_ratio * j->gear_ratio);
    double J_payload_ref = j->payload_mass * payload_dist * payload_dist
                           / (j->gear_ratio * j->gear_ratio);
    return j->motor_inertia + J_ref + J_payload_ref;
}

double robot_joint_natural_freq(const RobotJoint *j, double max_payload_dist)
{
    /* omega_n = sqrt(k_joint / J_eff) */
    if (!j || j->joint_stiffness <= 0.0) return 0.0;
    double J_eff = robot_joint_effective_inertia((RobotJoint*)j, max_payload_dist);
    if (J_eff <= 0.0) return 0.0;
    return sqrt(j->joint_stiffness / J_eff);
}

double robot_joint_damping_ratio(const RobotJoint *j)
{
    /* Simplified: zeta = b / (2*sqrt(k*J)) */
    if (!j || j->joint_stiffness <= 0.0) return 0.0;
    double J = j->link_inertia / (j->gear_ratio*j->gear_ratio) + j->motor_inertia;
    if (J <= 0.0) return 0.0;
    return j->joint_damping / (2.0 * sqrt(j->joint_stiffness * J));
}

double robot_joint_gravity_torque(const RobotJoint *j, double m_payload,
                                   double link_com_dist, double angle)
{
    /* T_g = (m_link*g*link_com_dist + m_payload*g*d_payload) * cos(angle)
     * angle from horizontal.
     */
    if (!j) return 0.0;
    /* Simplified: lumped COM */
    double total_mass_moment = m_payload * j->payload_distance
                               + m_payload * link_com_dist; /* using payload mass as link mass estimate */
    (void)link_com_dist;
    return total_mass_moment * 9.81 * cos(angle);
}

double robot_joint_gravity_compensation(const RobotJoint *j, double m_link,
                                         double link_com_dist, double angle)
{
    /* T_comp = m*g*d*cos(angle) */
    if (!j) return 0.0;
    double total = m_link*link_com_dist + j->payload_mass*j->payload_distance;
    return total * 9.81 * cos(angle);
}

double robot_joint_max_accel(const RobotJoint *j, double T_avail,
                              double T_grav, double T_fric)
{
    /* alpha_max = (T_avail - T_grav - T_fric) / J_eff */
    if (!j) return 0.0;
    double J_eff = robot_joint_effective_inertia((RobotJoint*)j, j->payload_distance);
    if (J_eff <= 0.0) return 0.0;
    double T_net = T_avail - fabs(T_grav) - fabs(T_fric);
    return T_net / J_eff;
}

double robot_joint_required_stiffness(const RobotJoint *j, double target_hz)
{
    /* k_req = J_eff * (2*pi*f_target)^2 */
    if (!j) return 0.0;
    double J_eff = robot_joint_effective_inertia((RobotJoint*)j, j->payload_distance);
    double wn = 2.0 * M_PI * target_hz;
    return J_eff * wn * wn;
}

void robot_joint_pd_gains(const RobotJoint *j, double wn_des, double zeta_des,
                           double *Kp, double *Kd)
{
    /* PD control: tau = Kp*(theta_des - theta) + Kd*(omega_des - omega)
     * Kp = J_eff * wn^2
     * Kd = 2 * zeta * wn * J_eff
     * Reference: Siciliano et al. (2009).
     */
    if (!j || !Kp || !Kd) return;
    double J_eff = robot_joint_effective_inertia((RobotJoint*)j, j->payload_distance);
    if (Kp) *Kp = J_eff * wn_des * wn_des;
    if (Kd) *Kd = 2.0 * zeta_des * wn_des * J_eff;
}

double robot_joint_feedforward_torque(const RobotJoint *j, double alpha_des,
                                       double T_grav, double T_fric)
{
    /* T_ff = J_eff * alpha_des + T_grav + T_fric */
    if (!j) return 0.0;
    double J_eff = robot_joint_effective_inertia((RobotJoint*)j, j->payload_distance);
    return J_eff * alpha_des + T_grav + T_fric;
}

double harmonic_drive_stiffness(double T_rated, double stiffness_ratio)
{
    /* Harmonic drive stiffness roughly proportional to rated torque.
     * k [Nm/rad] = stiffness_ratio * T_rated [Nm].
     * Typical stiffness_ratio: 50-100 for standard harmonic drives.
     */
    if (T_rated <= 0.0) return 0.0;
    return stiffness_ratio * T_rated;
}

double robot_joint_motor_sizing_torque(const RobotJoint *j, double T_rms,
                                        double w_rms, double duty)
{
    /* Motor continuous torque requirement:
     * T_cont = T_rms * sqrt(duty_cycle) with safety factor.
     */
    if (!j) return 0.0;
    (void)w_rms;
    double safety = 1.5;
    return T_rms * sqrt(duty) * safety;
}

void two_link_arm_natural_freqs(double m1, double m2, double L1, double L2,
                                 double EI1, double EI2, double *f1, double *f2)
{
    /* Simplified: treat each link as cantilever beam.
     * k_link = 3*E*I/L^3, omega_n = sqrt(k/m_eff).
     */
    if (f1) {
        double k1 = (L1 > 0) ? 3.0*EI1/(L1*L1*L1) : DBL_MAX;
        double m1_eff = m1 * 0.23; /* effective mass for cantilever */
        *f1 = (m1_eff > 0 && k1 > 0) ? sqrt(k1/m1_eff)/(2.0*M_PI) : 0.0;
    }
    if (f2) {
        double k2 = (L2 > 0) ? 3.0*EI2/(L2*L2*L2) : DBL_MAX;
        double m2_eff = m2 * 0.23;
        *f2 = (m2_eff > 0 && k2 > 0) ? sqrt(k2/m2_eff)/(2.0*M_PI) : 0.0;
    }
}

/* ===== Precision Machine Dynamics ===== */

double spindle_critical_speed(double d, double L, double m_rotor, double k_bearing)
{
    /* Simplified spindle critical speed (rigid rotor on compliant bearings).
     * omega_cr = sqrt(2*k_bearing / m_rotor)
     */
    if (m_rotor <= 0.0 || k_bearing <= 0.0) return 0.0;
    (void)d; (void)L;
    return sqrt(2.0 * k_bearing / m_rotor) * 60.0 / (2.0*M_PI);
}

double machine_frame_stiffness(double F_max, double defl_allowable)
{
    /* Required stiffness: k = F_max / defl_allowable */
    if (defl_allowable <= 0.0) return DBL_MAX;
    return F_max / defl_allowable;
}

double precision_stage_natural_freq(double m_stage, double k_bearing)
{
    /* omega_n = sqrt(k/m) */
    if (m_stage <= 0.0 || k_bearing <= 0.0) return 0.0;
    return sqrt(k_bearing / m_stage) / (2.0*M_PI);
}

double machine_foundation_isolation_freq(double m_machine, double k_iso, int n)
{
    /* System natural frequency on isolators */
    if (m_machine <= 0.0 || n <= 0) return 0.0;
    double k_total = n * k_iso;
    return sqrt(k_total / m_machine) / (2.0*M_PI);
}

double thermal_expansion(double L, double cte, double dT)
{
    /* delta_L = alpha * L * delta_T */
    return cte * L * dT;
}

/* ===== Aerospace / UAV Dynamics ===== */

double quadrotor_arm_nat_freq(double L_arm, double m_arm, double EI)
{
    /* Arm as cantilever beam with tip mass.
     * omega_n = sqrt(3*E*I/(0.23*m_arm + m_tip)*L^3)
     * Simplified: tip mass = m_arm (motor + prop ~ arm mass).
     */
    if (L_arm <= 0.0 || EI <= 0.0) return 0.0;
    double m_eff = 0.23 * m_arm + m_arm; /* arm + motor approximation */
    if (m_eff <= 0.0) return 0.0;
    double k = 3.0 * EI / (L_arm*L_arm*L_arm);
    return sqrt(k / m_eff) / (2.0*M_PI);
}

double thrust_stand_nat_freq(double m_eff, double k)
{
    /* Thrust measurement stand natural frequency */
    if (m_eff <= 0.0 || k <= 0.0) return 0.0;
    return sqrt(k / m_eff) / (2.0*M_PI);
}

double reaction_wheel_disturbance_freq(double rpm)
{
    /* f_disturbance = RPM / 60 [Hz] */
    return rpm / 60.0;
}

double actuator_bandwidth_for_control(double wn_aircraft, double GM_db, double PM_deg)
{
    /* Actuator bandwidth should be 5-10x the aircraft natural frequency.
     * With gain/phase margin requirements, need additional margin.
     * Simplified: bw_req = wn_aircraft * max(5, GM_db/6).
     */
    double margin = (GM_db > 6.0) ? GM_db / 6.0 : 5.0;
    (void)PM_deg;
    return wn_aircraft * margin;
}

double solar_panel_deploy_freq(double panel_mass, double panel_length,
                                double hinge_stiffness)
{
    /* Solar panel deployment: rigid body on torsional hinge.
     * omega_n = sqrt(k_hinge / I_panel)
     * I_panel = (1/3)*m*L^2 (hinged at end).
     */
    if (panel_mass <= 0.0 || panel_length <= 0.0 || hinge_stiffness <= 0.0) return 0.0;
    double I = (1.0/3.0) * panel_mass * panel_length * panel_length;
    return sqrt(hinge_stiffness / I) / (2.0*M_PI);
}

double pogo_stability_param(double thrust, double flow_rate, double compliance)
{
    /* Simplified Pogo stability parameter for launch vehicles.
     * Positive value indicates potential instability.
     * Reference: Rubin (1966), NASA SP-8055.
     */
    if (flow_rate <= 0.0) return 0.0;
    return thrust * compliance / flow_rate;
}

double momentum_wheel_sizing(double spacecraft_inertia, double slew_rate,
                              double max_disturbance, double orbit_period)
{
    /* Reaction/momentum wheel sizing.
     * Required angular momentum: H = I_sc * slew_rate + T_dist * orbit_period/4.
     */
    double H_slew = spacecraft_inertia * slew_rate;
    double H_dist = max_disturbance * orbit_period / 4.0;
    return H_slew + H_dist;
}
