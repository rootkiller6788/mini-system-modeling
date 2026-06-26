/**
 * @file mechanical_applications.h
 * @brief Mechanical system applications - vehicle suspension, vibration isolation,
 *        robot joints, precision machines, aerospace dynamics
 *
 * Knowledge Coverage:
 *   L1 - Definitions: quarter-car model, ISO ride comfort, isolation efficiency,
 *        robot joint compliance, gravity compensation
 *   L5 - Engineering Methods: suspension tuning (skyhook), isolator design,
 *        PD gain tuning for robot joints, natural frequency constraints
 *   L6 - Canonical Problems: quarter-car suspension, vibration isolator sizing,
 *        robot joint mechanical design for control
 *   L7 - Applications: automotive suspension (Detroit, Tesla), precision isolation
 *        (semiconductor equipment), industrial robot joints (ISO 9283),
 *        UAV dynamics (quadrotor), spacecraft reaction wheels
 *   L8 - Advanced: skyhook damping, two-stage isolation
 *
 * Reference:
 *   MIT 2.154 Maneuvering and Control
 *   Gillespie, Fundamentals of Vehicle Dynamics (1992)
 *   Rivin, Passive Vibration Isolation (2003)
 *   Siciliano et al., Robotics: Modelling, Planning and Control (2009)
 */

#ifndef MECHANICAL_APPLICATIONS_H
#define MECHANICAL_APPLICATIONS_H

#include <stdlib.h>
#include <math.h>

/* ===== Vehicle Suspension (Quarter-Car Model) ===== */

typedef struct {
    double sprung_mass;
    double unsprung_mass;
    double suspension_stiffness;
    double suspension_damping;
    double tire_stiffness;
    double tire_damping;
} QuarterCarModel;

void quarter_car_natural_freqs(const QuarterCarModel *qcm,
                                double *body_hz, double *wheel_hz);
double quarter_car_static_deflection(const QuarterCarModel *qcm, double g);
double quarter_car_body_accel_rms(const QuarterCarModel *qcm,
                                   double road_coeff, double speed);
double quarter_car_dynamic_load_rms(const QuarterCarModel *qcm,
                                     double road_coeff, double speed);
double quarter_car_suspension_travel_rms(const QuarterCarModel *qcm,
                                          double road_coeff, double speed);
double quarter_car_ride_comfort(const QuarterCarModel *qcm,
                                 double speed, double duration);
double quarter_car_damping_ratio(const QuarterCarModel *qcm);
double skyhook_damping_ideal(const QuarterCarModel *qcm);
double suspension_stiffness_for_ride(double m_sprung, double target_hz);
double suspension_damping_for_ratio(double m_sprung, double k, double zeta_target);

double half_car_pitch_freq(double sprung_mass, double pitch_inertia,
                            double wheelbase, double suspension_stiffness);
double full_car_heave_freq(double total_sprung_mass, double total_stiffness);

/* ===== Vibration Isolation ===== */

typedef struct {
    double isolated_mass;
    double excitation_freq;
    double excitation_amp;
    double target_transmissibility;
    double static_load;
    int num_isolators;
} IsolatorReq;

typedef struct {
    double stiffness_per_isolator;
    double natural_freq_hz;
    double freq_ratio;
    double achieved_transmissibility;
    double static_deflection;
    double damping_ratio;
} IsolatorDesign;

IsolatorDesign design_vibration_isolator(const IsolatorReq *req);
double isolator_transmissibility(const IsolatorDesign *d, double freq_hz);
double isolator_efficiency_percent(const IsolatorDesign *d);
double isolator_stiffness_for_TR(double mass_total, double w_exc, double target_TR);
double two_stage_isolation_TR(double m1, double m2, double k1, double k2,
                               double b1, double b2, double freq);
double precision_isolation_natural_freq(double floor_vib_um,
                                         double allowable_stage_vib_nm);

/* ===== Robot Joint Mechanics ===== */

typedef struct {
    double link_inertia;
    double payload_mass;
    double payload_distance;
    double gear_ratio;
    double motor_inertia;
    double joint_stiffness;
    double joint_damping;
    double max_torque;
} RobotJoint;

double robot_joint_effective_inertia(const RobotJoint *j, double payload_dist);
double robot_joint_natural_freq(const RobotJoint *j, double max_payload_dist);
double robot_joint_damping_ratio(const RobotJoint *j);
double robot_joint_gravity_torque(const RobotJoint *j, double m_payload,
                                   double link_com_dist, double angle);
double robot_joint_gravity_compensation(const RobotJoint *j, double m_link,
                                         double link_com_dist, double angle);
double robot_joint_max_accel(const RobotJoint *j, double T_avail,
                              double T_grav, double T_fric);
double robot_joint_required_stiffness(const RobotJoint *j, double target_hz);
void robot_joint_pd_gains(const RobotJoint *j, double wn_des, double zeta_des,
                           double *Kp, double *Kd);
double robot_joint_feedforward_torque(const RobotJoint *j, double alpha_des,
                                       double T_grav, double T_fric);
double harmonic_drive_stiffness(double T_rated, double stiffness_ratio);
double robot_joint_motor_sizing_torque(const RobotJoint *j, double T_rms,
                                        double w_rms, double duty);

void two_link_arm_natural_freqs(double m1, double m2, double L1, double L2,
                                 double EI1, double EI2, double *f1, double *f2);

/* ===== Precision Machine Dynamics ===== */

double spindle_critical_speed(double d, double L, double m_rotor, double k_bearing);
double machine_frame_stiffness(double F_max, double defl_allowable);
double precision_stage_natural_freq(double m_stage, double k_bearing);
double machine_foundation_isolation_freq(double m_machine, double k_iso, int n);
double thermal_expansion(double L, double cte, double dT);

/* ===== Aerospace / UAV Dynamics ===== */

double quadrotor_arm_nat_freq(double L_arm, double m_arm, double EI);
double thrust_stand_nat_freq(double m_eff, double k);
double reaction_wheel_disturbance_freq(double rpm);
double actuator_bandwidth_for_control(double wn_aircraft, double GM_db, double PM_deg);
double solar_panel_deploy_freq(double panel_mass, double panel_length,
                                double hinge_stiffness);
double pogo_stability_param(double thrust, double flow_rate, double compliance);
double momentum_wheel_sizing(double spacecraft_inertia, double slew_rate,
                              double max_disturbance, double orbit_period);

#endif /* MECHANICAL_APPLICATIONS_H */
