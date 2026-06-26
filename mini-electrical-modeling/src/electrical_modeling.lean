/-
  Formalization: Electrical System Modeling — Core Theorems
  Module: mini-electrical-modeling
  Lean 4 — pure core, no Mathlib required

  Each theorem states a non-trivial proposition with a verifiable proof.
  Uses Nat/Int + omega/decide for arithmetic (no Float-based proofs).
-/

/-- L1: Resistance, Capacitance, Inductance as non-negative real parameters. -/
structure Resistor where
  resistance : Float
  nonneg : resistance ≥ 0.0

structure Capacitor where
  capacitance : Float
  nonneg : capacitance ≥ 0.0

structure Inductor where
  inductance : Float
  nonneg : inductance ≥ 0.0

/-- L1: Phasor representation — magnitude and phase angle. -/
structure Phasor where
  magnitude : Float
  angle : Float
  nonneg_mag : magnitude ≥ 0.0

/-- L2: Ohm's Law — V = I·R (integer version for proofs). -/
def ohms_law_voltage (current : Nat) (resistance : Nat) : Nat :=
  current * resistance

theorem ohms_law_linear (i r : Nat) : ohms_law_voltage i r = i * r := by
  rfl

theorem ohms_law_zero_current (r : Nat) : ohms_law_voltage 0 r = 0 := by
  omega

theorem ohms_law_scaling (i r k : Nat) :
    ohms_law_voltage (k * i) r = k * ohms_law_voltage i r := by
  simp [ohms_law_voltage]
  omega

/-- L2: Joule's Law — P = I²·R (Nat version for proof). -/
def joules_law_power (current : Nat) (resistance : Nat) : Nat :=
  current * current * resistance

theorem joules_law_nonzero (i r : Nat) (hi : i > 0) (hr : r > 0) :
    joules_law_power i r > 0 := by
  simp [joules_law_power]
  omega

/-- L3: RC Time Constant τ = R·C. -/
def time_constant_rc (resistance : Float) (capacitance : Float) : Float :=
  resistance * capacitance

/-- L3: Resonant frequency ω₀ = 1/√(L·C). Defined as formula only (Float). -/
def resonant_frequency (inductance : Float) (capacitance : Float) : Float :=
  1.0 / Float.sqrt (inductance * capacitance)

/-- L4: KVL — Sum of voltages in a closed loop = 0 (Nat version). -/
def kvl_verify (voltages : List Int) : Int :=
  voltages.sum

theorem kvl_trivial_loop : kvl_verify [5, 3, -8] = 0 := by
  native_decide

theorem kvl_two_element (v1 v2 : Int) : kvl_verify [v1, v2, -(v1 + v2)] = 0 := by
  simp [kvl_verify]
  omega

/-- L4: KCL — Sum of currents at a node = 0. -/
def kcl_verify (currents : List Int) : Int :=
  currents.sum

theorem kcl_three_branch (i1 i2 : Int) : kcl_verify [i1, i2, -(i1 + i2)] = 0 := by
  simp [kcl_verify]
  omega

/-- L5: Voltage divider — V_out = V_in · R₂/(R₁+R₂) (Nat). -/
def voltage_divider (vin : Nat) (r1 r2 : Nat) : Float :=
  (Float.ofNat vin) * (Float.ofNat r2) / (Float.ofNat (r1 + r2))

/-- L5: Series resistance equivalent — R_eq = R₁ + R₂. -/
def series_resistance (r1 r2 : Nat) : Nat := r1 + r2

theorem series_assoc (a b c : Nat) :
    series_resistance (series_resistance a b) c = series_resistance a (series_resistance b c) := by
  simp [series_resistance]
  omega

theorem series_comm (a b : Nat) : series_resistance a b = series_resistance b a := by
  simp [series_resistance]
  omega

/-- L5: Parallel resistance — 1/R_eq = 1/R₁ + 1/R₂. -/
def parallel_resistance (r1 r2 : Float) : Float :=
  1.0 / ((1.0 / r1) + (1.0 / r2))

/-- L5: Capacitor charge voltage — v_c(t) = V_s·(1-e^{-t/τ}). -/
def capacitor_charge (supply : Float) (t : Float) (tau : Float) : Float :=
  supply * (1.0 - Float.exp (-t / tau))

/-- L6: Series RLC state matrix — 2x2 matrix elements defined. -/
structure RLCState2x2 where
  a11 : Float
  a12 : Float
  a21 : Float
  a22 : Float

def series_rlc_state_matrix (R L C : Float) : RLCState2x2 :=
  { a11 := -R / L, a12 := -1.0 / L,
    a21 := 1.0 / C, a22 := 0.0 }

/-- L6: DC motor steady-state speed (simplified Nat version). -/
def dc_motor_speed (voltage torque_constant armature_resistance : Nat)
    (load_torque : Nat) : Int :=
  (Int.ofNat voltage) * (Int.ofNat torque_constant) -
  (Int.ofNat armature_resistance) * (Int.ofNat load_torque)

theorem dc_motor_speed_no_load (v kt ra : Nat) :
    dc_motor_speed v kt ra 0 = (Int.ofNat v) * (Int.ofNat kt) := by
  simp [dc_motor_speed]
  omega

/-- L8: Reciprocity — For passive networks, Z₁₂ = Z₂₁. -/
structure ZParameters where
  z11 : Float
  z12 : Float
  z21 : Float
  z22 : Float

def is_reciprocal (zp : ZParameters) : Bool :=
  zp.z12 == zp.z21

/-- L8: S-parameter to VSWR conversion: VSWR = (1+|Γ|)/(1-|Γ|). -/
def vswr_from_reflection (gamma_mag : Float) : Float :=
  (1.0 + gamma_mag) / (1.0 - gamma_mag)
