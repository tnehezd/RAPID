#ifndef THERMAL_MODEL_H
#define THERMAL_MODEL_H

typedef enum {
    THERMAL_SUCCESS = 0,
    THERMAL_INVALID_INPUT,
    THERMAL_NUMERICAL_ERROR,
    THERMAL_NOT_BRACKETED,
    THERMAL_NO_CONVERGENCE
} ThermalStatus;

typedef enum {
    THERMODYNAMICS_LOCALLY_ISOTHERMAL = 0,
    THERMODYNAMICS_RADIATIVE
} ThermodynamicsModel;

/* All dimensional quantities in this API use CGS, except temperature (K). */
typedef struct {
    double surface_density;       /* g cm^-2, strictly positive */
    double orbital_frequency;     /* s^-1, strictly positive */
    double alpha;                 /* dimensionless, nonnegative */
    double mean_molecular_weight; /* in hydrogen-atom masses */
    double background_temperature;
    double shear_parameter;       /* |d ln Omega/d ln r|; Keplerian = 1.5 */
} ThermalParameters;

typedef struct {
    double minimum_temperature;
    double maximum_temperature;
    double relative_tolerance;
    int maximum_iterations;
} ThermalSolverOptions;

typedef struct {
    double temperature;
    double sound_speed;           /* cm s^-1, isothermal */
    double scale_height;          /* cm */
    double midplane_density;      /* g cm^-3 */
    double viscosity;             /* cm^2 s^-1 */
    double opacity;               /* cm^2 g^-1 */
    double optical_depth;         /* midplane to one surface */
    double effective_optical_depth;
    double heating;               /* erg cm^-2 s^-1, both surfaces */
    double cooling;               /* net of background, both surfaces */
    int iterations;
} ThermalState;

/* Smoothed Lin-Papaloizou fit, including the electron-scattering branch.
 * Requires finite rho > 0 and 0 < T < 1e9 K. Output is unchanged on failure. */
ThermalStatus calculateLinPapaloizouOpacity(double density, double temperature,
                                          double *opacity);

/* Effective vertical closure H = cs/Omega; not a vertical SG solution.
 * Heating is nu*Sigma*(shear_parameter*Omega)^2. */
ThermalStatus calculateThermalState(const ThermalParameters *parameters,
                                   double temperature, ThermalState *state);

/* Defaults: 1--1500 K, relative tolerance 1e-8, at most 128 bisections. */
ThermalSolverOptions defaultThermalSolverOptions(void);

/* Solve cooling = heating within an explicit sign-changing bracket.
 * The lower bound is raised to T_bg. Multiple roots are not enumerated and
 * thermal stability is not inferred. State is unchanged on any failure. */
ThermalStatus solveMidplaneTemperature(const ThermalParameters *parameters,
                                      const ThermalSolverOptions *options,
                                      ThermalState *state);

const char *thermalStatusMessage(ThermalStatus status);

#endif
