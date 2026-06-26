/**
 * mini-fluid-thermal-modeling — Core Definitions
 *
 * Combined fluid mechanics and heat transfer modeling module.
 * Covers incompressible/compressible flow, conduction, convection,
 * radiation, and coupled fluid-thermal systems.
 *
 * Reference Texts:
 *   - Incropera & DeWitt, "Fundamentals of Heat and Mass Transfer" (2007)
 *   - White, "Fluid Mechanics" (2016)
 *   - Cengel & Boles, "Thermodynamics: An Engineering Approach" (2014)
 *   - Bird, Stewart, Lightfoot, "Transport Phenomena" (2007)
 *
 * Knowledge Levels Covered:
 *   L1: Core type definitions, fluid/thermal properties
 *   L2: Flow regimes, heat transfer modes, dimensionless groups
 *   L3: Engineering quantity ranges and typical values
 */

#ifndef FLUID_THERMAL_CORE_H
#define FLUID_THERMAL_CORE_H

#include <stddef.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * L1: Core Enumerations — Flow and Thermal Classification
 * ------------------------------------------------------------------------- */

/** Flow regime classification based on Reynolds number */
typedef enum {
    FLOW_CREEPING,       /** Re < 1, Stokes flow */
    FLOW_LAMINAR,        /** 1 <= Re < 2300, smooth layered flow */
    FLOW_TRANSITIONAL,   /** 2300 <= Re < 4000, intermittent turbulence */
    FLOW_TURBULENT,      /** Re >= 4000, fully turbulent */
    FLOW_HYPERSONIC      /** Ma > 5, compressibility dominant */
} FlowRegime;

/** Compressibility regime */
typedef enum {
    COMPRESS_INCOMPRESSIBLE,  /** Ma < 0.3, density ~ constant */
    COMPRESS_SUBSONIC,        /** 0.3 <= Ma < 1.0 */
    COMPRESS_TRANSONIC,       /** 0.8 <= Ma <= 1.2 */
    COMPRESS_SUPERSONIC,      /** 1.2 < Ma <= 5 */
    COMPRESS_HYPERSONIC       /** Ma > 5 */
} CompressibilityRegime;

/** Heat transfer mode */
typedef enum {
    HT_MODE_CONDUCTION,       /** Fourier's law, stationary media */
    HT_MODE_CONVECTION,       /** Newton's law, fluid motion */
    HT_MODE_RADIATION,        /** Stefan-Boltzmann, EM waves */
    HT_MODE_BOILING,          /** Phase change: liquid → vapor */
    HT_MODE_CONDENSATION,     /** Phase change: vapor → liquid */
    HT_MODE_COMBINED          /** Multiple modes simultaneously */
} HeatTransferMode;

/** Convection type */
typedef enum {
    CONV_FORCED,              /** External driver (pump, fan, wind) */
    CONV_NATURAL,             /** Buoyancy-driven (free convection) */
    CONV_MIXED                /** Combined forced + natural */
} ConvectionType;

/** Heat exchanger flow arrangement */
typedef enum {
    HX_PARALLEL_FLOW,         /** Co-current */
    HX_COUNTER_FLOW,          /** Counter-current (most efficient) */
    HX_CROSS_FLOW,            /** Perpendicular streams */
    HX_SHELL_AND_TUBE,        /** Shell-side + tube-side */
    HX_COMPACT                /** Plate-fin, microchannel */
} HXFlowArrangement;

/** Pipe flow model */
typedef enum {
    PIPE_DARCY_WEISBACH,      /** Universal friction factor method */
    PIPE_HAZEN_WILLIAMS,      /** Empirical, water systems */
    PIPE_MANNING,             /** Open channel flow */
    PIPE_POISEUILLE           /** Exact laminar solution */
} PipeFlowModel;

/* ---------------------------------------------------------------------------
 * L1: Core Data Structures — Fluid and Thermal Properties
 * ------------------------------------------------------------------------- */

/** Fluid thermodynamic and transport properties at a state point */
typedef struct {
    double density;           /** Density ρ [kg/m³] */
    double dynamic_viscosity; /** Dynamic viscosity μ [Pa·s] */
    double kinematic_viscosity; /** Kinematic viscosity ν = μ/ρ [m²/s] */
    double thermal_conductivity; /** Thermal conductivity k [W/(m·K)] */
    double specific_heat_cp;  /** Specific heat at constant pressure cp [J/(kg·K)] */
    double specific_heat_cv;  /** Specific heat at constant volume cv [J/(kg·K)] */
    double thermal_diffusivity; /** Thermal diffusivity α = k/(ρ·cp) [m²/s] */
    double prandtl_number;    /** Pr = ν/α = cp·μ/k [dimensionless] */
    double thermal_expansion; /** Thermal expansion coefficient β [1/K] */
    double speed_of_sound;    /** Speed of sound c [m/s] */
    double surface_tension;   /** Surface tension σ [N/m] */
    double bulk_modulus;      /** Bulk modulus K [Pa] */
} FluidProperties;

/** Flow condition at a cross-section or point */
typedef struct {
    double velocity;           /** Flow velocity u [m/s] */
    double mass_flow_rate;     /** Mass flow rate ṁ [kg/s] */
    double volumetric_flow_rate; /** Volumetric flow rate Q [m³/s] */
    double pressure;           /** Static pressure p [Pa] */
    double temperature;        /** Temperature T [K] */
    double reynolds_number;    /** Re = ρ·u·L/μ [dimensionless] */
    double mach_number;        /** Ma = u/c [dimensionless] */
    FlowRegime regime;         /** Flow regime classification */
    CompressibilityRegime compressibility; /** Compressibility regime */
} FlowCondition;

/** Geometry descriptor for characteristic length calculation */
typedef struct {
    double diameter;           /** Pipe/hole diameter D [m] */
    double hydraulic_diameter; /** Hydraulic diameter Dh = 4A/P [m] */
    double length;             /** Streamwise length L [m] */
    double cross_section_area; /** Flow area A [m²] */
    double wetted_perimeter;   /** Wetted perimeter P [m] */
    double surface_area;       /** Heat transfer surface area As [m²] */
    double volume;             /** Volume V [m³] */
    double roughness;          /** Surface roughness ε [m] */
} GeometryDescriptor;

/** Thermal boundary/initial condition */
typedef struct {
    double ambient_temperature;       /** Far-field temperature T∞ [K] */
    double surface_temperature;       /** Wall/surface temperature Ts [K] */
    double heat_flux;                 /** Surface heat flux q'' [W/m²] */
    double heat_generation;           /** Volumetric heat generation q̇ [W/m³] */
    double convective_coefficient;    /** Heat transfer coefficient h [W/(m²·K)] */
    HeatTransferMode mode;            /** Dominant heat transfer mode */
    ConvectionType convection_type;   /** Forced/natural/mixed convection */
} ThermalCondition;

/** Complete fluid-thermal state at a point */
typedef struct {
    FluidProperties fluid;        /** Fluid properties at local conditions */
    FlowCondition flow;           /** Flow state */
    ThermalCondition thermal;     /** Thermal state */
    GeometryDescriptor geometry;  /** Local geometry */
    double prandtl_number;        /** Pr [dimensionless] */
    double nusselt_number;        /** Nu = h·L/k [dimensionless] */
    double grashof_number;        /** Gr = g·β·ΔT·L³/ν² [dimensionless] */
    double rayleigh_number;       /** Ra = Gr·Pr [dimensionless] */
    double peclet_number;         /** Pe = Re·Pr [dimensionless] */
    double stanton_number;        /** St = Nu/(Re·Pr) [dimensionless] */
    double fourier_number;        /** Fo = α·t/L² [dimensionless] */
    double biot_number;           /** Bi = h·L/k [dimensionless] */
} FluidThermalState;

/* ---------------------------------------------------------------------------
 * L1: Dimensionless Group Computation
 * ------------------------------------------------------------------------- */

/**
 * Compute Reynolds number.
 * Re = ρ·u·L / μ = u·L / ν
 *
 * Source: Osborne Reynolds (1883), "An experimental investigation of the
 * circumstances which determine whether the motion of water shall be
 * direct or sinuous."
 */
double compute_reynolds(double density, double velocity,
                        double char_length, double dynamic_viscosity);

/**
 * Compute Prandtl number.
 * Pr = ν/α = cp·μ/k
 *
 * Ratio of momentum diffusivity to thermal diffusivity.
 * Source: Ludwig Prandtl (1910)
 */
double compute_prandtl(double cp, double mu, double k);

/**
 * Compute Nusselt number from heat transfer coefficient.
 * Nu = h·L / k
 *
 * Ratio of convective to conductive heat transfer.
 * Source: Wilhelm Nusselt (1915), "The basic laws of heat transfer"
 */
double compute_nusselt(double h, double char_length, double k);

/**
 * Compute Grashof number.
 * Gr = g·β·(Ts - T∞)·L³ / ν²
 *
 * Ratio of buoyancy to viscous forces in natural convection.
 * Source: Franz Grashof (1863)
 */
double compute_grashof(double g, double beta, double delta_t,
                       double char_length, double nu);

/**
 * Compute Rayleigh number.
 * Ra = Gr·Pr
 *
 * Product of Grashof and Prandtl; determines natural convection regime.
 */
double compute_rayleigh(double gr, double pr);

/**
 * Compute Mach number.
 * Ma = u / c
 *
 * Source: Ernst Mach (1887)
 */
double compute_mach(double velocity, double speed_of_sound);

/**
 * Compute Peclet number.
 * Pe = Re·Pr = u·L / α
 *
 * Ratio of advective to diffusive transport rate.
 */
double compute_peclet(double re, double pr);

/**
 * Compute Biot number.
 * Bi = h·Lc / k_solid
 *
 * Ratio of internal conduction resistance to surface convection resistance.
 * Bi < 0.1 → lumped capacitance valid.
 * Source: Jean-Baptiste Biot (1804)
 */
double compute_biot(double h, double char_length, double k_solid);

/**
 * Compute Fourier number.
 * Fo = α·t / L²
 *
 * Dimensionless time in transient conduction.
 * Source: Joseph Fourier (1822), "Théorie analytique de la chaleur"
 */
double compute_fourier(double alpha, double time, double char_length);

/**
 * Compute Stanton number.
 * St = Nu / (Re·Pr) = h / (ρ·u·cp)
 *
 * Modified Nusselt number; ratio of heat transferred to fluid thermal capacity.
 */
double compute_stanton(double nu, double re, double pr);

/* ---------------------------------------------------------------------------
 * L2: Flow Regime Classification
 * ------------------------------------------------------------------------- */

/**
 * Determine flow regime from Reynolds number.
 *
 * Re < 1       → Creeping (Stokes flow)
 * 1 ≤ Re < 2300 → Laminar
 * 2300 ≤ Re < 4000 → Transitional
 * Re ≥ 4000    → Turbulent
 *
 * Source: Reynolds (1883), confirmed by Moody (1944)
 */
FlowRegime classify_flow_regime(double reynolds_number);

/**
 * Determine compressibility regime from Mach number.
 */
CompressibilityRegime classify_compressibility(double mach_number);

/**
 * Determine convection type from Gr/Re² ratio.
 * Gr/Re² ≪ 1 → Forced convection dominant
 * Gr/Re² ≫ 1 → Natural convection dominant
 * Gr/Re² ≈ 1 → Mixed convection
 *
 * Source: Churchill (1977), "A comprehensive correlating equation"
 */
ConvectionType classify_convection(double grashof, double reynolds);

/* ---------------------------------------------------------------------------
 * L2: Fluid State Updates
 * ------------------------------------------------------------------------- */

/**
 * Initialize FluidProperties for common fluids.
 *
 * @param fluid_name  "air", "water", "oil_engine", "glycerin", "mercury"
 * @param temperature Temperature in Kelvin (for T-dependent properties)
 * @param props       Output: populated FluidProperties struct
 * @return 0 on success, -1 if fluid not found
 *
 * Properties sourced from Incropera Appendix A.
 */
int fluid_properties_init(const char *fluid_name, double temperature,
                          FluidProperties *props);

/**
 * Update temperature-dependent fluid properties using empirical correlations.
 * Modifies viscosity, thermal conductivity, Pr based on temperature.
 *
 * Sutherland's law for gases: μ/μ₀ = (T/T₀)^1.5 · (T₀+S)/(T+S)
 * where S is the Sutherland constant [K].
 *
 * Source: Sutherland (1893), "The viscosity of gases and molecular force"
 */
void fluid_properties_update_temperature(FluidProperties *props,
                                         double temperature);

/**
 * Apply Sutherland's law for gas viscosity.
 * μ(T) = μ_ref · (T/T_ref)^1.5 · (T_ref + S) / (T + S)
 *
 * @param T          Temperature [K]
 * @param mu_ref     Reference viscosity at T_ref [Pa·s]
 * @param T_ref      Reference temperature [K]
 * @param S          Sutherland constant [K]
 * @return Viscosity at T [Pa·s]
 */
double sutherland_viscosity(double T, double mu_ref, double T_ref, double S);

#endif /* FLUID_THERMAL_CORE_H */
