# Knowledge Graph — mini-fluid-thermal-modeling

## L1: Definitions ✅ Complete

| # | Concept | C Definition | Lean Definition |
|---|---------|-------------|-----------------|
| 1 | FluidProperties | `typedef struct { double density; ... }` | `structure FluidProperties where ...` |
| 2 | FlowCondition | `typedef struct { double velocity; ... }` | `structure FlowCondition where ...` |
| 3 | ThermalCondition | `typedef struct { double ambient_temperature; ... }` | `structure ThermalCondition where ...` |
| 4 | GeometryDescriptor | `typedef struct { double diameter; ... }` | `structure GeometryDescriptor where ...` |
| 5 | FlowRegime | `typedef enum { FLOW_CREEPING, ... }` | `inductive FlowRegime where ...` |
| 6 | HeatTransferMode | `typedef enum { HT_MODE_CONDUCTION, ... }` | `inductive HeatTransferMode where ...` |
| 7 | ConvectionType | `typedef enum { CONV_FORCED, ... }` | `inductive ConvectionType where ...` |
| 8 | HXFlowArrangement | `typedef enum { HX_PARALLEL_FLOW, ... }` | `inductive HXFlowArrangement where ...` |
| 9 | CompressibilityRegime | `typedef enum { COMPRESS_INCOMPRESSIBLE, ... }` | `inductive CompressibilityRegime where ...` |
| 10 | PrandtlCategory | N/A | `inductive PrandtlCategory where ...` |

## L2: Core Concepts ✅ Complete

| # | Concept | Implementation |
|---|---------|----------------|
| 1 | Reynolds number | `compute_reynolds()` in `fluid_thermal_core.c` |
| 2 | Prandtl number | `compute_prandtl()` in `fluid_thermal_core.c` |
| 3 | Nusselt number | `compute_nusselt()` in `fluid_thermal_core.c` |
| 4 | Grashof number | `compute_grashof()` in `fluid_thermal_core.c` |
| 5 | Rayleigh number | `compute_rayleigh()` in `fluid_thermal_core.c` |
| 6 | Peclet number | `compute_peclet()` in `fluid_thermal_core.c` |
| 7 | Biot number | `compute_biot()` in `fluid_thermal_core.c` |
| 8 | Fourier number | `compute_fourier()` in `fluid_thermal_core.c` |
| 9 | Stanton number | `compute_stanton()` in `fluid_thermal_core.c` |
| 10 | Mach number | `compute_mach()` in `fluid_thermal_core.c` |
| 11 | Flow regime classification | `classify_flow_regime()` |
| 12 | Compressibility classification | `classify_compressibility()` |
| 13 | Convection type classification | `classify_convection()` |
| 14 | Boundary layer theory | `blasius_boundary_layer_thickness()`, etc. in `fluid_mechanics.c` |
| 15 | Heat exchanger effectiveness | `heat_exchanger_effectiveness()` in `heat_transfer.c` |
| 16 | Natural circulation | `natural_circulation_driving_pressure()` in `thermal_fluid_systems.c` |

## L3: Engineering Quantities ✅ Complete

| # | Quantity | API |
|---|----------|-----|
| 1 | Typical Re ranges | `typical_reynolds_range()` |
| 2 | Typical Prandtl values | `typical_prandtl()` |
| 3 | Typical h [W/m²K] ranges | `typical_h_range()` |
| 4 | Thermal conductivity k | `typical_thermal_conductivity()` |
| 5 | Specific heat cp | `typical_specific_heat()` |
| 6 | Density ρ | `typical_density()` |
| 7 | Minor loss coefficients | `minor_loss_coefficient()` |
| 8 | Pipe roughness values | embedded in fluid property tables |

## L4: Conservation Laws ✅ Complete

| # | Law | C Implementation | Lean Theorem |
|---|-----|-----------------|--------------|
| 1 | Continuity (mass) | `continuity_equation()` | `continuity_determines_velocity` |
| 2 | Bernoulli | `bernoulli_solve()` | `bernoulli_pressure_velocity_relation` |
| 3 | Navier-Stokes (Poiseuille) | `hagen_poiseuille_velocity()` | N/A |
| 4 | Navier-Stokes (Couette) | `couette_velocity()` | N/A |
| 5 | Fourier's Law | `fourier_conduction()` | `fourier_heat_flow_positive` |
| 6 | Newton's Law of Cooling | `newton_cooling()` | `newton_cooling_positive_when_surface_hotter` |
| 7 | Stefan-Boltzmann | `stefan_boltzmann_emission()` | N/A |
| 8 | Mass conservation positivity | N/A | `mass_conservation_positivity` |

## L5: Engineering Methods ✅ Complete

| # | Method | Implementation |
|---|--------|----------------|
| 1 | Darcy-Weisbach equation | `darcy_weisbach_head_loss()`, `darcy_weisbach_pressure_drop()` |
| 2 | Colebrook-White equation | `colebrook_friction_factor()` (Newton-Raphson) |
| 3 | Swamee-Jain explicit | `swamee_jain_friction_factor()` |
| 4 | Laminar friction factor | `darcy_friction_laminar()` |
| 5 | Dittus-Boelter correlation | `dittus_boelter_nusselt()` |
| 6 | Sieder-Tate correlation | `sieder_tate_nusselt()` |
| 7 | Churchill-Chu correlation | `churchill_chu_vertical_plate()` |
| 8 | McAdams horizontal plate | `mcadams_horizontal_plate()` |
| 9 | Morgan horizontal cylinder | `morgan_horizontal_cylinder()` |
| 10 | Flat plate Nusselt (laminar) | `flat_plate_laminar_nusselt_average()` |
| 11 | Flat plate Nusselt (turbulent) | `flat_plate_turbulent_nusselt_average()` |
| 12 | Zhukauskas tube bank | `zhukauskas_tube_bank_nusselt()` |
| 13 | LMTD method | `lmtd()`, `heat_exchanger_duty_lmtd()` |
| 14 | ε-NTU method | `heat_exchanger_effectiveness()`, `ntu_from_effectiveness_*()` |
| 15 | Fin analysis | `fin_efficiency()`, `fin_heat_transfer_rate()` |
| 16 | Lumped capacitance | `lumped_capacitance_temp()` |
| 17 | One-term approximation | `one_term_approximation_slab()` |
| 18 | Thermal resistance networks | `thermal_resistance_plane_wall()`, etc. |
| 19 | Overall U coefficient | `overall_heat_transfer_coefficient_*()` |
| 20 | Buckingham Pi theorem | `buckingham_pi_count()`, `dimensionless_group()` |
| 21 | Dimensional analysis | `compute_froude()`, `compute_weber()`, etc. |
| 22 | Finite difference (1D) | `solve_1d_steady_conduction()` |
| 23 | Thomas algorithm | `thomas_algorithm()` |
| 24 | Gauss-Seidel (2D) | `gauss_seidel_2d_laplace()` |
| 25 | SOR (2D) | `sor_2d_laplace()` |
| 26 | Simpson integration | `simpson_integral()` |
| 27 | RK4 ODE solver | `rk4_solve()` |
| 28 | Hardy Cross method | `hardy_cross_solver()` |

## L6: Engineering Problems ✅ Complete

| # | Problem | Implementation |
|---|---------|----------------|
| 1 | Pipe sizing (head loss) | `example_pipe_flow.c` |
| 2 | Heat exchanger design | `example_heat_exchanger.c` |
| 3 | Electronics cooling | `example_electronics_cooling.c` |
| 4 | Pipe heating (CHF) | `pipe_heating_outlet_temp()` |
| 5 | Pipe heating (CWT) | `pipe_outlet_temp_cwt()` |
| 6 | Pump system curve | `system_curve_head()`, `pump_system_operating_point()` |
| 7 | Natural circulation design | `natural_circulation_flow_rate()` |
| 8 | NPSH analysis | `npsh_available()` |
| 9 | Critical insulation radius | `critical_insulation_radius()` |

## L7: Applications ✅ Complete (4 applications)

| # | Application | Implementation |
|---|-------------|----------------|
| 1 | CPU heat sink design | `example_electronics_cooling.c`, `heat_sink_thermal_resistance()`, `junction_temperature()` |
| 2 | 3D IC stacked die thermal | `stacked_die_temperatures()` |
| 3 | HVAC sensible/latent load | `sensible_heat_load()`, `latent_heat_load()`, `air_changes_per_hour()` |
| 4 | Automotive radiator | `radiator_heat_rejection()`, `pump_power_scaling()` |

## L8: Advanced Methods ✅ Partial (4/7)

| # | Topic | Implementation |
|---|-------|----------------|
| 1 | Conjugate heat transfer | `conjugate_nusselt_correction()` |
| 2 | Multi-phase flow basics | `void_fraction_homogeneous()`, `two_phase_multiplier()` |
| 3 | SIMPLE pressure correction | `simple_pressure_correction_1d()` |
| 4 | Crank-Nicolson implicit | `solve_1d_transient_conduction_implicit()` |
| 5 | CFD turbulence modeling | Not implemented |
| 6 | Non-Newtonian fluids | Not implemented |
| 7 | Multi-scale modeling | Not implemented |

## L9: Research Frontiers ✅ Partial (documented only)

| # | Topic | Status |
|---|-------|--------|
| 1 | Nanofluid heat transfer | Documented reference |
| 2 | Microchannel flow boiling | Documented reference |
| 3 | Non-Fourier heat conduction | Documented reference |
| 4 | Bio-heat transfer (Pennes) | Documented reference |
