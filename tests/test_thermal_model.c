#include "thermal_model.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #condition); \
        ++failures; \
    } \
} while (0)

static int closeRelative(double actual, double expected, double tolerance)
{
    return isfinite(actual) && fabs(actual - expected) <= tolerance * fabs(expected);
}

static ThermalParameters referenceParameters(void)
{
    ThermalParameters parameters = {100.0, 1.991013777793097e-7, 1e-3, 2.34, 10.0, 1.5};
    return parameters;
}

static void testOpacity(void)
{
    /* Independent direct-power CGS evaluations, not log-space evaluations. */
    const double temperatures[] = {10, 100, 200, 500, 1500, 3000, 10000, 100000};
    const double expected[] = {
        0.02, 1.99990000751815, 1.57919617991672, 2.49997440026283,
        0.520242886327743, 0.00101736933124135, 8.32050294337844,
        0.348030026941942
    };
    for (size_t i = 0; i < sizeof(temperatures) / sizeof(temperatures[0]); ++i) {
        double opacity = -1.0;
        CHECK(calculateLinPapaloizouOpacity(1e-9, temperatures[i], &opacity) ==
              THERMAL_SUCCESS);
        CHECK(closeRelative(opacity, expected[i], 2e-12));
    }
    /* The low-density high-temperature branch intentionally omits scattering. */
    double opacity;
    CHECK(calculateLinPapaloizouOpacity(1e-12, 1e5, &opacity) == THERMAL_SUCCESS);
    CHECK(opacity < 0.348);

    for (int i = 0; i <= 12; ++i) {
        double density = pow(10.0, -16.0 + i);
        double boundaries[] = {
            1600.0 * pow(density, 0.0444444444),
            5700.0 * pow(density, 0.02381),
            2.28e6 * pow(density, 0.2267)
        };
        for (size_t j = 0; j < sizeof(boundaries) / sizeof(boundaries[0]); ++j) {
            CHECK(calculateLinPapaloizouOpacity(density, boundaries[j] * (1.0 - 1e-9),
                                               &opacity) == THERMAL_SUCCESS);
            CHECK(opacity > 0.0 && isfinite(opacity));
            CHECK(calculateLinPapaloizouOpacity(density, boundaries[j] * (1.0 + 1e-9),
                                               &opacity) == THERMAL_SUCCESS);
            CHECK(opacity > 0.0 && isfinite(opacity));
        }
    }
    opacity = -123.0;
    CHECK(calculateLinPapaloizouOpacity(0.0, 100.0, &opacity) == THERMAL_INVALID_INPUT);
    CHECK(calculateLinPapaloizouOpacity(NAN, 100.0, &opacity) == THERMAL_INVALID_INPUT);
    CHECK(calculateLinPapaloizouOpacity(1e-9, 0.0, &opacity) == THERMAL_INVALID_INPUT);
    CHECK(calculateLinPapaloizouOpacity(1e-9, INFINITY, &opacity) == THERMAL_INVALID_INPUT);
    CHECK(calculateLinPapaloizouOpacity(1e-9, 1e9, &opacity) == THERMAL_INVALID_INPUT);
    CHECK(calculateLinPapaloizouOpacity(1e-9, 100.0, NULL) == THERMAL_INVALID_INPUT);
    CHECK(opacity == -123.0);
}

static void testState(void)
{
    ThermalParameters parameters = referenceParameters();
    ThermalState state;
    CHECK(calculateThermalState(&parameters, 100.0, &state) == THERMAL_SUCCESS);
    double cs_squared = 1.380649e-16 * 100.0 / (2.34 * 1.6735575e-24);
    CHECK(closeRelative(state.sound_speed * state.sound_speed, cs_squared, 1e-14));
    CHECK(closeRelative(state.scale_height, sqrt(cs_squared) /
                        parameters.orbital_frequency, 1e-14));
    CHECK(closeRelative(state.midplane_density * state.scale_height *
                        sqrt(2.0 * acos(-1.0)), parameters.surface_density, 1e-14));
    CHECK(closeRelative(state.viscosity, parameters.alpha * cs_squared /
                        parameters.orbital_frequency, 1e-14));
    CHECK(closeRelative(state.heating, 2.25 * parameters.alpha * cs_squared *
                        parameters.surface_density * parameters.orbital_frequency, 1e-14));
    CHECK(closeRelative(state.optical_depth,
                        state.opacity * parameters.surface_density / 2.0, 1e-14));
    parameters.shear_parameter = 1.0;
    CHECK(calculateThermalState(&parameters, 100.0, &state) == THERMAL_SUCCESS);
    CHECK(closeRelative(state.heating, parameters.alpha * cs_squared *
                        parameters.surface_density * parameters.orbital_frequency, 1e-14));
    parameters.shear_parameter = 0.0;
    CHECK(calculateThermalState(&parameters, 100.0, &state) == THERMAL_SUCCESS);
    CHECK(state.heating == 0.0 && state.viscosity > 0.0);
    parameters.shear_parameter = 1.5;

    const double surface_densities[] = {1e-8, 1e8};
    for (size_t i = 0; i < 2; ++i) {
        parameters.surface_density = surface_densities[i];
        CHECK(calculateThermalState(&parameters, 100.0, &state) == THERMAL_SUCCESS);
        double tau = state.optical_depth;
        if (i == 0) {
            CHECK(tau < 1e-6);
            CHECK(closeRelative(state.effective_optical_depth, 0.25 / tau, 1e-6));
        } else {
            CHECK(tau > 1e6);
            CHECK(closeRelative(state.effective_optical_depth, 0.375 * tau, 1e-6));
        }
        double expected_cooling = 2.0 * 5.670374419e-5 *
            (pow(100.0, 4) - pow(parameters.background_temperature, 4)) /
            (0.375 * tau + sqrt(3.0) / 4.0 + 0.25 / tau);
        CHECK(closeRelative(state.cooling, expected_cooling, 1e-14));
    }
}

static void testSolver(void)
{
    ThermalParameters parameters = referenceParameters();
    ThermalSolverOptions options = defaultThermalSolverOptions();
    ThermalState state;
    CHECK(solveMidplaneTemperature(&parameters, &options, &state) == THERMAL_SUCCESS);
    CHECK(state.temperature > parameters.background_temperature);
    CHECK(state.iterations > 0 && state.iterations <= options.maximum_iterations);
    CHECK(closeRelative(state.cooling, state.heating, options.relative_tolerance));

    /* Manufacture known equilibria by deriving alpha from independently
     * evaluated cooling. Narrow brackets isolate each branch of the fit. */
    const double targets[] = {20, 100, 200, 500, 1200, 3000, 10000, 100000};
    for (size_t i = 0; i < sizeof(targets) / sizeof(targets[0]); ++i) {
        double temperature = targets[i];
        ThermalState target;
        parameters = referenceParameters();
        parameters.alpha = 1.0;
        CHECK(calculateThermalState(&parameters, temperature, &target) == THERMAL_SUCCESS);
        parameters.alpha = target.cooling / target.heating;
        options.minimum_temperature = 0.99 * temperature;
        options.maximum_temperature = 1.01 * temperature;
        CHECK(solveMidplaneTemperature(&parameters, &options, &state) == THERMAL_SUCCESS);
        CHECK(closeRelative(state.temperature, temperature, 1e-8));
        CHECK(closeRelative(state.cooling, state.heating, options.relative_tolerance));
    }

    parameters = referenceParameters();
    parameters.alpha = 0.0;
    options = defaultThermalSolverOptions();
    CHECK(solveMidplaneTemperature(&parameters, &options, &state) == THERMAL_SUCCESS);
    CHECK(state.temperature == parameters.background_temperature);
    CHECK(state.heating == 0.0 && state.cooling == 0.0 && state.iterations == 0);

    parameters = referenceParameters();
    ThermalState sentinel = {0};
    sentinel.temperature = -123.0;
    state = sentinel;
    options.maximum_temperature = 10.01;
    CHECK(solveMidplaneTemperature(&parameters, &options, &state) == THERMAL_NOT_BRACKETED);
    CHECK(memcmp(&state, &sentinel, sizeof(state)) == 0);
    options = defaultThermalSolverOptions();
    options.maximum_iterations = 1;
    CHECK(solveMidplaneTemperature(&parameters, &options, &state) == THERMAL_NO_CONVERGENCE);
    CHECK(memcmp(&state, &sentinel, sizeof(state)) == 0);
    options = defaultThermalSolverOptions();
    options.relative_tolerance = NAN;
    CHECK(solveMidplaneTemperature(&parameters, &options, &state) == THERMAL_INVALID_INPUT);
    options = defaultThermalSolverOptions();
    parameters.surface_density = -1.0;
    CHECK(solveMidplaneTemperature(&parameters, &options, &state) == THERMAL_INVALID_INPUT);
    CHECK(calculateThermalState(&parameters, 100.0, &state) == THERMAL_INVALID_INPUT);
    parameters = referenceParameters();
    parameters.alpha = DBL_MAX;
    CHECK(solveMidplaneTemperature(&parameters, &options, &state) == THERMAL_NUMERICAL_ERROR);
    CHECK(memcmp(&state, &sentinel, sizeof(state)) == 0);
    CHECK(solveMidplaneTemperature(NULL, &options, &state) == THERMAL_INVALID_INPUT);
    CHECK(solveMidplaneTemperature(&parameters, NULL, &state) == THERMAL_INVALID_INPUT);
    CHECK(solveMidplaneTemperature(&parameters, &options, NULL) == THERMAL_INVALID_INPUT);
}

static int printProfile(void)
{
    ThermalSolverOptions options = defaultThermalSolverOptions();
    printf("# Synthetic 1 Msun disk, Sigma=100*(r/AU)^-1 g/cm^2, alpha=1e-3\n");
    printf("# r_AU T_K H_over_r rho_g_cm3 nu_cm2_s kappa_cm2_g tau Qplus Qminus\n");
    for (int i = 0; i <= 100; ++i) {
        double radius_au = pow(100.0, i / 100.0);
        ThermalParameters parameters = referenceParameters();
        parameters.surface_density = 100.0 / radius_au;
        parameters.orbital_frequency = sqrt(6.67430e-8 * 1.98847e33 /
                                            pow(radius_au * 1.495978707e13, 3));
        ThermalState state;
        ThermalStatus status = solveMidplaneTemperature(&parameters, &options, &state);
        if (status != THERMAL_SUCCESS) {
            fprintf(stderr, "r=%.9g AU: %s\n", radius_au, thermalStatusMessage(status));
            return 1;
        }
        printf("%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
               radius_au, state.temperature, state.scale_height / (radius_au * 1.495978707e13),
               state.midplane_density, state.viscosity, state.opacity, state.optical_depth,
               state.heating, state.cooling);
    }
    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--profile") == 0) {
        return printProfile();
    }
    if (argc != 1) {
        fprintf(stderr, "Usage: %s [--profile]\n", argv[0]);
        return 1;
    }
    testOpacity();
    testState();
    testSolver();
    if (failures != 0) {
        fprintf(stderr, "%d thermal checks failed\n", failures);
        return 1;
    }
    printf("Thermal model: all checks passed\n");
    return 0;
}
