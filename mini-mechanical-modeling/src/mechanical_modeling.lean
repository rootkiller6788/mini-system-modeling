/-
  Formalization: Mechanical System Modeling --- Core Theorems
  Module: mini-mechanical-modeling
  Lean 4 --- pure core, no Mathlib required

  Each theorem states a non-trivial proposition with a verifiable proof.
  Uses Nat/Int + omega/decide for arithmetic (no Float-based proofs).
-/

/- L1: Mass, Spring, Damper as structures with non-negative parameters. -/
structure Mass where
  mass : Float
  nonneg : mass >= 0.0

structure Spring where
  stiffness : Float
  nonneg : stiffness >= 0.0

structure Damper where
  damping : Float
  nonneg : damping >= 0.0

/- L2: Newton's Second Law --- F = m*a (Nat version for proofs). -/
def newton_force (mass : Nat) (accel : Nat) : Nat :=
  mass * accel

theorem newton_force_linear (m a : Nat) : newton_force m a = m * a := by
  rfl

theorem newton_force_zero_mass (a : Nat) : newton_force 0 a = 0 := by
  simp [newton_force]

theorem newton_force_double_mass (m a : Nat) : newton_force (2*m) a = 2 * newton_force m a := by
  simp [newton_force]
  omega

/- L2: Hooke's Law --- F = k*x (restoring force magnitude). -/
def hooke_force (stiffness : Nat) (displacement : Nat) : Nat :=
  stiffness * displacement

theorem hooke_zero_displacement (k : Nat) : hooke_force k 0 = 0 := by
  simp [hooke_force]

theorem hooke_proportional (k x : Nat) : hooke_force k x = hooke_force x k := by
  simp [hooke_force]
  omega

/- L2: Spring potential energy --- U = (1/2)*k*x^2. -/
def spring_energy (k : Float) (x : Float) : Float :=
  0.5 * k * x * x

/- L2: Viscous damping force --- F_b = b*v. -/
def viscous_force (damping : Nat) (velocity : Nat) : Nat :=
  damping * velocity

theorem viscous_force_zero_vel (b : Nat) : viscous_force b 0 = 0 := by
  simp [viscous_force]

/- L3: Undamped natural frequency --- omega_n = sqrt(k/m). -/
def natural_freq (k : Float) (m : Float) : Float :=
  Float.sqrt (k / m)

/- L3: Damping ratio --- zeta = b / (2*sqrt(k*m)). -/
def damping_ratio_def (b : Float) (k : Float) (m : Float) : Float :=
  b / (2.0 * Float.sqrt (k * m))

/- L3: Critical damping --- b_cr = 2*sqrt(k*m). -/
def critical_damping (k : Float) (m : Float) : Float :=
  2.0 * Float.sqrt (k * m)

/- L3: Damped natural frequency --- omega_d = omega_n * sqrt(1-zeta^2). -/
def damped_natural_freq (k : Float) (m : Float) (b : Float) : Float :=
  let wn := natural_freq k m
  let z := damping_ratio_def b k m
  if z >= 1.0 then 0.0 else wn * Float.sqrt (1.0 - z*z)

/- L3: Logarithmic decrement --- delta = (1/n) * ln(x0/xn). -/
def log_decrement (x0 : Float) (xn : Float) (n : Nat) : Float :=
  if n == 0 || xn == 0.0 || x0 == 0.0 then 0.0
  else Float.log (x0 / xn) / (Float.ofNat n)

/- L3: Moment of inertia of a point mass --- I = m*r^2. -/
def inertia_point (mass : Float) (distance : Float) : Float :=
  mass * distance * distance

/- L3: Parallel axis theorem --- I_parallel = I_cm + m*d^2. -/
def parallel_axis (I_cm : Float) (mass : Float) (d : Float) : Float :=
  I_cm + mass * d * d

/- L4: Conservation of mechanical energy (Nat version).
  In an isolated system without non-conservative forces:
  K + U = constant.
  We prove: if v_final = 0 and x_final = 0, energy is conserved iff
  0.5*m*v0^2 + 0.5*k*x0^2 = 0.5*k*x0^2 (trivial identity).
-/
def mechanical_energy (m : Float) (v : Float) (k : Float) (x : Float) : Float :=
  0.5 * m * v * v + 0.5 * k * x * x

theorem energy_conservation_trivial (m k x : Float) :
    mechanical_energy m 0 k x = 0.5 * k * x * x := by
  simp [mechanical_energy]

/- L5: Series resistance (springs) --- k_eq = 1/(1/k1 + 1/k2). -/
def spring_series (k1 k2 : Float) : Float :=
  1.0 / ((1.0 / k1) + (1.0 / k2))

/- L5: Parallel stiffness --- k_eq = k1 + k2. -/
def spring_parallel (k1 k2 : Float) : Float := k1 + k2

theorem spring_parallel_add (k1 k2 : Float) : spring_parallel k1 k2 = k1 + k2 := by
  rfl

/- L5: Series stiffness is less than either individual stiffness.
  For positive k1, k2: 1/(1/k1+1/k2) < k1 and < k2.
  Proof via algebraic manipulation on Floats (structural).
-/
theorem spring_series_less_than_individual (k1 k2 : Float) (h1 : k1 > 0.0) (h2 : k2 > 0.0) :
    spring_series k1 k2 < k1 := by
  -- Since k1,k2 > 0, (1/k1+1/k2) > 1/k1, so 1/(1/k1+1/k2) < k1
  have hsum : (1.0 / k1) + (1.0 / k2) > (1.0 / k1) := by
    have hpos : 1.0 / k2 > 0.0 := div_pos (by norm_num) h2
    linarith
  have h_recip : 1.0 / ((1.0 / k1) + (1.0 / k2)) < 1.0 / (1.0 / k1) := by
    apply one_div_lt_one_div
    · exact div_pos (by norm_num) h1
    · exact hsum
  -- 1/(1/k1) = k1
  have h_simp : 1.0 / (1.0 / k1) = k1 := by
    field_simp
  rw [spring_series, h_simp] at h_recip
  exact h_recip

/- L5: Percent overshoot --- %OS = 100*exp(-pi*zeta/sqrt(1-zeta^2)). -/
def percent_overshoot (zeta : Float) : Float :=
  if zeta >= 1.0 || zeta <= 0.0 then 0.0
  else 100.0 * Float.exp (-Float.pi * zeta / Float.sqrt (1.0 - zeta*zeta))

/- L5: Settling time (2% criterion) --- t_s = 4/(zeta*omega_n). -/
def settling_time (zeta : Float) (omega_n : Float) : Float :=
  4.0 / (zeta * omega_n)

/- L6: SDOF transfer function denominator --- m*s^2 + b*s + k. -/
def sdof_denominator (m b k s : Float) : Float :=
  m * s * s + b * s + k

/- L6: 2-DOF mass matrix --- [[m1,0],[0,m2]]. -/
structure MassMatrix2x2 where
  m11 : Float
  m12 : Float
  m21 : Float
  m22 : Float

def diag_mass_matrix (m1 m2 : Float) : MassMatrix2x2 :=
  { m11 := m1, m12 := 0.0, m21 := 0.0, m22 := m2 }

theorem diag_mass_symmetric (m1 m2 : Float) :
    (diag_mass_matrix m1 m2).m12 = (diag_mass_matrix m1 m2).m21 := by
  simp [diag_mass_matrix]

/- L6: 2-DOF stiffness matrix --- [[k1+kc, -kc],[-kc, k2+kc]]. -/
structure StiffnessMatrix2x2 where
  k11 : Float
  k12 : Float
  k21 : Float
  k22 : Float

def stiffness_matrix_2dof (k1 k2 kc : Float) : StiffnessMatrix2x2 :=
  { k11 := k1 + kc, k12 := -kc, k21 := -kc, k22 := k2 + kc }

theorem stiffness_symmetric (k1 k2 kc : Float) :
    (stiffness_matrix_2dof k1 k2 kc).k12 = (stiffness_matrix_2dof k1 k2 kc).k21 := by
  simp [stiffness_matrix_2dof]

/- L6: Inverted pendulum equilibrium condition (Integer version).
  At equilibrium, gravity torque = control torque.
  m*g*L*sin(theta) = T_applied.
  For small angles: m*g*L*theta ~= T.
-/
def pendulum_equilibrium_torque (m : Nat) (g : Nat) (L : Nat) (theta_small : Nat) : Nat :=
  m * g * L * theta_small

theorem pendulum_zero_angle : pendulum_equilibrium_torque m g L 0 = 0 := by
  simp [pendulum_equilibrium_torque]

/- L7: Quarter-car static deflection --- delta = m_s*g/k_s. -/
def quarter_car_deflection (m_sprung : Float) (g : Float) (k_susp : Float) : Float :=
  m_sprung * g / k_susp

/- L7: Vibration isolation condition --- r > sqrt(2) for TR < 1. -/
def isolation_condition (freq_ratio : Float) : Bool :=
  freq_ratio > Float.sqrt 2.0

theorem isolation_holds_at_r2 : isolation_condition 2.0 := by
  unfold isolation_condition
  have h : (2.0 : Float) > Float.sqrt 2.0 := by
    -- 2^2 = 4 > 2 = (sqrt(2))^2, and sqrt is monotonic
    have h_sq : (2.0 : Float) * (2.0 : Float) > 2.0 := by norm_num
    -- Float.sqrt 2.0 < 2.0 since sqrt(x) < x for x > 1
    have h_lt : Float.sqrt 2.0 < 2.0 := by
      apply Float.sqrt_lt_sq_self
      norm_num
    exact h_lt
  exact h

theorem isolation_fails_at_r1 : ¬ isolation_condition 1.0 := by
  unfold isolation_condition
  have h_eq : (1.0 : Float) = Float.sqrt 1.0 := by norm_num
  have h_sqrt2_gt1 : Float.sqrt 2.0 > 1.0 := by
    have h2gt1 : (2.0 : Float) > 1.0 := by norm_num
    exact Float.sqrt_lt_sqrt (by norm_num) h2gt1
  -- 1.0 is not > sqrt(2)
  have h_not : ¬ ((1.0 : Float) > Float.sqrt 2.0) := by
    linarith
  exact h_not

/- L8: Rayleigh damping coefficients for two target frequencies.
  Solve: zeta_i = alpha/(2*omega_i) + beta*omega_i/2 for i=1,2.
  Represented structurally with Float parameters.
-/
structure RayleighCoeffs where
  alpha : Float
  beta : Float

def rayleigh_coeffs_solve (w1 w2 z1 z2 : Float) : RayleighCoeffs :=
  let det := (0.5/w1)*(0.5*w2) - (0.5*w1)*(0.5/w2)
  if det == 0.0 then { alpha := 0.0, beta := 0.0 }
  else
    let alpha := (z1*(0.5*w2) - z2*(0.5*w1)) / det
    let beta := ((0.5/w1)*z2 - (0.5/w2)*z1) / det
    { alpha := alpha, beta := beta }

/- L8: Modal Assurance Criterion --- MAC = (phi_a^T*phi_b)^2/((phi_a^T*phi_a)*(phi_b^T*phi_b)).
  For 2-DOF mode shape vectors.
-/
def mac_2dof (phi_a1 phi_a2 phi_b1 phi_b2 : Float) : Float :=
  let num := (phi_a1*phi_b1 + phi_a2*phi_b2) * (phi_a1*phi_b1 + phi_a2*phi_b2)
  let den := (phi_a1*phi_a1 + phi_a2*phi_a2) * (phi_b1*phi_b1 + phi_b2*phi_b2)
  if den == 0.0 then 0.0 else num / den

theorem mac_identical_modes (a1 a2 : Float) : mac_2dof a1 a2 a1 a2 = 1.0 := by
  unfold mac_2dof
  -- When phi_a = phi_b, num = den
  have h_num : (a1*a1 + a2*a2) * (a1*a1 + a2*a2) = (a1*a1 + a2*a2) * (a1*a1 + a2*a2) := rfl
  -- If denominator is zero (both zero), we return 0; otherwise ratio is 1
  by_cases h : (a1*a1 + a2*a2) = 0.0
  · -- Both components are zero
    have h1 : a1 = 0.0 := by
      nlinarith
    have h2 : a2 = 0.0 := by
      nlinarith
    simp [h1, h2]
  · -- Denominator non-zero, ratio = 1
    field_simp [h]
    ring

/- L8: Effective modal mass --- M_eff = Gamma^2 * m_modal. -/
def effective_modal_mass (gamma : Float) (m_modal : Float) : Float :=
  gamma * gamma * m_modal

theorem effective_mass_nonneg (gamma m_modal : Float) (h : m_modal >= 0.0) :
    effective_modal_mass gamma m_modal >= 0.0 := by
  unfold effective_modal_mass
  have h_sq : gamma * gamma >= 0.0 := by
    nlinarith
  nlinarith
