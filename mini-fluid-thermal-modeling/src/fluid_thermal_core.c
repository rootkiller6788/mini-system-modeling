/**
 * mini-fluid-thermal-modeling — Core Implementation
 *
 * Implementation of core dimensionless group computations,
 * flow regime classification, and fluid property models.
 *
 * Knowledge Levels: L1-L2
 */

#include "fluid_thermal_core.h"
#include <string.h>
#include <stdio.h>

/* =====================================================================
 * L1: Dimensionless Group Computation
 * ===================================================================== */

double compute_reynolds(double density, double velocity,
                        double char_length, double dynamic_viscosity)
{
    if (dynamic_viscosity == 0.0) return INFINITY;
    if (char_length == 0.0) return 0.0;
    return density * velocity * char_length / dynamic_viscosity;
}

double compute_prandtl(double cp, double mu, double k)
{
    if (k == 0.0) return INFINITY;
    if (cp == 0.0) return 0.0;
    return cp * mu / k;
}

double compute_nusselt(double h, double char_length, double k)
{
    if (k == 0.0) return INFINITY;
    return h * char_length / k;
}

double compute_grashof(double g, double beta, double delta_t,
                       double char_length, double nu)
{
    if (nu == 0.0) return INFINITY;
    if (char_length == 0.0) return 0.0;
    double nu_sq = nu * nu;
    if (nu_sq == 0.0) return INFINITY;
    return fabs(g) * beta * fabs(delta_t)
           * char_length * char_length * char_length / nu_sq;
}

double compute_rayleigh(double gr, double pr)
{
    return gr * pr;
}

double compute_mach(double velocity, double speed_of_sound)
{
    if (speed_of_sound == 0.0) return INFINITY;
    return velocity / speed_of_sound;
}

double compute_peclet(double re, double pr)
{
    return re * pr;
}

double compute_biot(double h, double char_length, double k_solid)
{
    if (k_solid == 0.0) return INFINITY;
    return h * char_length / k_solid;
}

double compute_fourier(double alpha, double time, double char_length)
{
    if (char_length == 0.0) return INFINITY;
    return alpha * time / (char_length * char_length);
}

double compute_stanton(double nu, double re, double pr)
{
    double denom = re * pr;
    if (denom == 0.0) return INFINITY;
    return nu / denom;
}

/* =====================================================================
 * L2: Flow Regime Classification
 * ===================================================================== */

FlowRegime classify_flow_regime(double reynolds_number)
{
    if (reynolds_number < 0.0) return FLOW_CREEPING; /* invalid, treat as creeping */
    if (reynolds_number < 1.0) return FLOW_CREEPING;
    if (reynolds_number < 2300.0) return FLOW_LAMINAR;
    if (reynolds_number < 4000.0) return FLOW_TRANSITIONAL;
    return FLOW_TURBULENT;
}

CompressibilityRegime classify_compressibility(double mach_number)
{
    if (mach_number < 0.0) return COMPRESS_INCOMPRESSIBLE; /* invalid */
    if (mach_number < 0.3) return COMPRESS_INCOMPRESSIBLE;
    if (mach_number < 0.8) return COMPRESS_SUBSONIC;
    if (mach_number <= 1.2) return COMPRESS_TRANSONIC;
    if (mach_number <= 5.0) return COMPRESS_SUPERSONIC;
    return COMPRESS_HYPERSONIC;
}

ConvectionType classify_convection(double grashof, double reynolds)
{
    if (reynolds == 0.0) {
        if (grashof == 0.0) return CONV_FORCED; /* default */
        return CONV_NATURAL;
    }
    double re_sq = reynolds * reynolds;
    if (re_sq == 0.0) return CONV_NATURAL;
    double ratio = grashof / re_sq;
    if (ratio < 0.1) return CONV_FORCED;
    if (ratio > 10.0) return CONV_NATURAL;
    return CONV_MIXED;
}

/* =====================================================================
 * L2: Fluid State Updates
 * ===================================================================== */

/**
 * Lookup table for common fluid properties at reference temperature.
 * All values at ~300K, 1 atm unless noted.
 *
 * Sources:
 *   Incropera & DeWitt (2007), Appendix A
 *   White (2016), Appendix A
 */
typedef struct {
    const char *name;
    double rho;       /* density [kg/m³] */
    double mu;        /* dynamic viscosity [Pa·s] */
    double k;         /* thermal conductivity [W/(m·K)] */
    double cp;        /* specific heat [J/(kg·K)] */
    double cv;        /* specific heat cv [J/(kg·K)] */
    double c;         /* speed of sound [m/s] */
    double beta;      /* expansion coefficient [1/K] */
    double sigma;     /* surface tension [N/m] */
} FluidData;

static const FluidData fluid_table[] = {
    {"air", 1.1614, 1.846e-5, 0.0263, 1007.0, 720.0, 347.0, 0.00333, 0.0},
    {"water", 997.0, 8.55e-4, 0.613, 4179.0, 4150.0, 1490.0, 0.000257, 0.0728},
    {"oil_engine", 884.0, 0.486, 0.145, 1900.0, 1750.0, 1350.0, 0.00070, 0.035},
    {"glycerin", 1260.0, 1.49, 0.286, 2430.0, 2200.0, 1920.0, 0.00050, 0.0633},
    {"mercury", 13540.0, 1.53e-3, 8.54, 139.0, 100.0, 1450.0, 0.000181, 0.487},
    {"ethanol", 789.0, 1.10e-3, 0.171, 2440.0, 2300.0, 1140.0, 0.00109, 0.0227},
    {"silicone_oil", 960.0, 0.048, 0.15, 1550.0, 1450.0, 1000.0, 0.00092, 0.020},
    {"hydrogen", 0.0819, 8.96e-6, 0.186, 14320.0, 10180.0, 1310.0, 0.00333, 0.0},
    {"helium", 0.1625, 1.99e-5, 0.156, 5193.0, 3116.0, 1020.0, 0.00333, 0.0},
    {"r134a", 4.25, 1.11e-5, 0.014, 1030.0, 860.0, 165.0, 0.0042, 0.006}
};

static const int fluid_table_size =
    sizeof(fluid_table) / sizeof(fluid_table[0]);

int fluid_properties_init(const char *fluid_name, double temperature,
                          FluidProperties *props)
{
    if (!fluid_name || !props) return -1;

    for (int i = 0; i < fluid_table_size; i++) {
        if (strcmp(fluid_name, fluid_table[i].name) == 0) {
            const FluidData *fd = &fluid_table[i];
            props->density = fd->rho;
            props->dynamic_viscosity = fd->mu;
            props->kinematic_viscosity = fd->mu / fd->rho;
            props->thermal_conductivity = fd->k;
            props->specific_heat_cp = fd->cp;
            props->specific_heat_cv = fd->cv;
            props->thermal_diffusivity = fd->k / (fd->rho * fd->cp);
            props->prandtl_number = fd->cp * fd->mu / fd->k;
            props->thermal_expansion = fd->beta;
            props->speed_of_sound = fd->c;
            props->surface_tension = fd->sigma;
            props->bulk_modulus = fd->rho * fd->c * fd->c;

            /* Apply temperature corrections if significantly different from 300K */
            if (fabs(temperature - 300.0) > 5.0) {
                fluid_properties_update_temperature(props, temperature);
            }
            return 0;
        }
    }
    return -1; /* fluid not found */
}

void fluid_properties_update_temperature(FluidProperties *props,
                                         double temperature)
{
    if (!props) return;

    /* Sutherland's law for viscosity (for gases; approximate for liquids) */
    /* μ(T) = μ_ref · (T/T_ref)^1.5 · (T_ref + S) / (T + S) */
    double T_ref = 300.0;
    double Sutherland_S = 110.4; /* air default; approximation for others */
    double mu_old = props->dynamic_viscosity;

    /* For liquid-like fluids (Pr > 2), use exponential T dependence:
     * μ ∝ exp(A + B/T) simplified for small ranges */
    if (props->prandtl_number > 2.0) {
        /* Approximate liquid viscosity T-dependence:
         * μ/μ₀ ≈ (T/T₀)^(-1.6) for common liquids near room temp */
        double ratio = temperature / T_ref;
        props->dynamic_viscosity = mu_old * pow(ratio, -1.6);
    } else {
        /* Sutherland for gases */
        double T = temperature;
        double ratio = T / T_ref;
        props->dynamic_viscosity = mu_old * pow(ratio, 1.5)
            * (T_ref + Sutherland_S) / (T + Sutherland_S);
    }

    /* Update derived properties */
    if (props->density > 0.0) {
        props->kinematic_viscosity =
            props->dynamic_viscosity / props->density;
    }

    /* Thermal conductivity: k ∝ T^0.8 for gases, k ∝ T^{-0.3} for liquids */
    if (props->prandtl_number > 2.0) {
        /* Liquid: slight decrease with T */
        props->thermal_conductivity *= pow(temperature / T_ref, -0.3);
    } else {
        /* Gas: increase with T */
        props->thermal_conductivity *= pow(temperature / T_ref, 0.8);
    }

    /* Specific heat: approximately constant for small T ranges */
    /* (Deviations are fluid-specific; using small correction) */

    /* Recompute derived non-dimensional properties */
    if (props->thermal_conductivity > 0.0) {
        props->prandtl_number = props->specific_heat_cp
            * props->dynamic_viscosity / props->thermal_conductivity;
    }
    if (props->density > 0.0 && props->specific_heat_cp > 0.0
        && props->thermal_conductivity > 0.0) {
        props->thermal_diffusivity = props->thermal_conductivity
            / (props->density * props->specific_heat_cp);
    }
}

double sutherland_viscosity(double T, double mu_ref, double T_ref, double S)
{
    if (T <= 0.0 || T_ref <= 0.0) return mu_ref;
    return mu_ref * pow(T / T_ref, 1.5) * (T_ref + S) / (T + S);
}
