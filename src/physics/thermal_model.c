#include "thermal_model.h"
#include "config.h"

#include <math.h>
#include <stddef.h>

static double logSum(double a, double b)
{
    double larger = fmax(a, b);
    return larger + log1p(exp(fmin(a, b) - larger));
}

static double logSmoothMinimum(double a, double b, double power)
{
    return -logSum(-power * a, -power * b) / power;
}

static double logSmoothMaximum(double a, double b, double power)
{
    return logSum(power * a, power * b) / power;
}

ThermalStatus calculateLinPapaloizouOpacity(double density, double temperature,
                                          double *opacity)
{
    if (opacity == NULL || !isfinite(density) || density <= 0.0 ||
        !isfinite(temperature) || temperature <= 0.0 || temperature >= 1e9) {
        return THERMAL_INVALID_INPUT;
    }

    /* Evaluate the blended power laws in log space to avoid intermediate
     * overflow. The fit uses CGS density, not the simulation's code density. */
    double log_rho = log(density);
    double log_t = log(temperature);
    double dust = log(5e-3) + log_t;
    double molecules = log(2e-8) + (2.0 / 3.0) * log_rho + 3.0 * log_t;
    double kramers = log(1.5e20) + log_rho - 2.5 * log_t;
    double log_opacity;

    if (log_t <= log(1600.0) + 0.0444444444 * log_rho) {
        double ice = log(2e-4) + 2.0 * log_t;
        double ice_evaporation = log(2e16) - 7.0 * log_t;
        double suppressed_dust = dust - logSum(0.0, log(1e22) - 10.0 * log_t);
        log_opacity = logSmoothMaximum(
            logSmoothMinimum(ice, ice_evaporation, 2.0), suppressed_dust, 4.0);
    } else if (log_t <= log(5700.0) + 0.02381 * log_rho) {
        double grain_evaporation =
            log(2e34) + (2.0 / 3.0) * log_rho - 9.0 * log_t;
        double suppressed_molecules =
            molecules - logSum(0.0, 8.0 * (log(3000.0) - log_t));
        log_opacity = logSmoothMaximum(
            logSmoothMinimum(dust, grain_evaporation, 4.0),
            suppressed_molecules, 4.0);
    } else if (density <= 1e-10 ||
               log_t < log(2.28e6) + 0.2267 * log_rho) {
        double h_minus = log(1e-36) + log_rho / 3.0 + 10.0 * log_t;
        double suppressed_molecules = molecules - logSum(
            0.0, 10.0 * (log_t - log(1.1e4) - 0.04762 * log_rho));
        log_opacity = logSmoothMaximum(
            logSmoothMinimum(h_minus, kramers, 2.0),
            suppressed_molecules, 4.0);
    } else {
        log_opacity = logSmoothMaximum(kramers, log(0.348), 4.0);
    }

    double value = exp(log_opacity);
    if (!isfinite(value) || value <= 0.0) {
        return THERMAL_NUMERICAL_ERROR;
    }
    *opacity = value;
    return THERMAL_SUCCESS;
}

static int validParameters(const ThermalParameters *parameters)
{
    return parameters != NULL &&
           isfinite(parameters->surface_density) && parameters->surface_density > 0.0 &&
           isfinite(parameters->orbital_frequency) && parameters->orbital_frequency > 0.0 &&
           isfinite(parameters->alpha) && parameters->alpha >= 0.0 &&
           isfinite(parameters->mean_molecular_weight) && parameters->mean_molecular_weight > 0.0 &&
           isfinite(parameters->shear_parameter) && parameters->shear_parameter >= 0.0 &&
           isfinite(parameters->background_temperature) && parameters->background_temperature > 0.0 &&
           parameters->background_temperature < 1e9 &&
           isfinite(parameters->irradiation_temperature) &&
           parameters->irradiation_temperature >= 0.0 &&
           parameters->irradiation_temperature < 1e9;
}

static double radiationTemperature(const ThermalParameters *parameters)
{
    if (parameters->irradiation_temperature == 0.0) return parameters->background_temperature;
    return sqrt(hypot(parameters->background_temperature * parameters->background_temperature,
                      parameters->irradiation_temperature * parameters->irradiation_temperature));
}

ThermalStatus calculateThermalState(const ThermalParameters *parameters,
                                   double temperature, ThermalState *state)
{
    if (!validParameters(parameters) || state == NULL ||
        !isfinite(temperature) || temperature <= 0.0 || temperature >= 1e9) {
        return THERMAL_INVALID_INPUT;
    }

    ThermalState value = {0};
    value.temperature = temperature;
    value.sound_speed = sqrt(BOLTZMANN_CONSTANT_CGS * temperature /
                             (parameters->mean_molecular_weight * HYDROGEN_MASS_IN_GRAMS));
    value.scale_height = value.sound_speed / parameters->orbital_frequency;
    value.midplane_density = parameters->surface_density /
                             (SQRT_TWO_PI * value.scale_height);
    if (!isfinite(value.sound_speed) || value.sound_speed <= 0.0 ||
        !isfinite(value.scale_height) || value.scale_height <= 0.0 ||
        !isfinite(value.midplane_density) || value.midplane_density <= 0.0) {
        return THERMAL_NUMERICAL_ERROR;
    }
    ThermalStatus status = calculateLinPapaloizouOpacity(
        value.midplane_density, temperature, &value.opacity);
    if (status != THERMAL_SUCCESS) {
        return status;
    }
    value.viscosity = parameters->alpha * value.sound_speed * value.scale_height;
    value.optical_depth = 0.5 * value.opacity * parameters->surface_density;
    value.effective_optical_depth = 0.375 * value.optical_depth +
                                   sqrt(3.0) / 4.0 +
                                   0.25 / value.optical_depth;
    value.heating = parameters->shear_parameter * parameters->shear_parameter *
                    value.viscosity * parameters->surface_density *
                    parameters->orbital_frequency * parameters->orbital_frequency;
    double background = radiationTemperature(parameters);
    /* Factored T^4 - T_rad^4 avoids cancellation near the irradiation floor. */
    value.cooling = 2.0 * STEFAN_BOLTZMANN_CONSTANT_CGS * (temperature - background) *
                    (temperature + background) *
                    (temperature * temperature + background * background) /
                    value.effective_optical_depth;
    if (!isfinite(value.viscosity) || !isfinite(value.optical_depth) ||
        value.optical_depth <= 0.0 || !isfinite(value.effective_optical_depth) ||
        !isfinite(value.heating) || !isfinite(value.cooling) ||
        (parameters->alpha > 0.0 && value.viscosity <= 0.0) ||
        (parameters->alpha > 0.0 && parameters->shear_parameter > 0.0 &&
         value.heating <= 0.0)) {
        return THERMAL_NUMERICAL_ERROR;
    }
    *state = value;
    return THERMAL_SUCCESS;
}

ThermalSolverOptions defaultThermalSolverOptions(void)
{
    ThermalSolverOptions options = {1.0, 1500.0, 1e-8, 128};
    return options;
}

static int isBalanced(const ThermalState *state, double tolerance)
{
    return fabs(state->cooling - state->heating) <=
           tolerance * fmax(fabs(state->cooling), state->heating);
}

ThermalStatus solveMidplaneTemperature(const ThermalParameters *parameters,
                                      const ThermalSolverOptions *options,
                                      ThermalState *state)
{
    if (!validParameters(parameters) || options == NULL || state == NULL ||
        !isfinite(options->minimum_temperature) || options->minimum_temperature <= 0.0 ||
        !isfinite(options->maximum_temperature) || options->maximum_temperature >= 1e9 ||
        options->maximum_temperature <= options->minimum_temperature ||
        !isfinite(options->relative_tolerance) || options->relative_tolerance <= 0.0 ||
        options->relative_tolerance >= 1.0 || options->maximum_iterations <= 0) {
        return THERMAL_INVALID_INPUT;
    }

    double lower = fmax(options->minimum_temperature, radiationTemperature(parameters));
    double upper = options->maximum_temperature;
    if (lower > upper) {
        return THERMAL_NOT_BRACKETED;
    }
    ThermalState left, right;
    ThermalStatus status = calculateThermalState(parameters, lower, &left);
    if (status != THERMAL_SUCCESS) {
        return status;
    }
    status = calculateThermalState(parameters, upper, &right);
    if (status != THERMAL_SUCCESS) {
        return status;
    }
    if (isBalanced(&left, options->relative_tolerance)) {
        *state = left;
        return THERMAL_SUCCESS;
    }
    if (isBalanced(&right, options->relative_tolerance)) {
        *state = right;
        return THERMAL_SUCCESS;
    }
    double left_residual = left.cooling - left.heating;
    double right_residual = right.cooling - right.heating;
    if ((left_residual < 0.0) == (right_residual < 0.0)) {
        return THERMAL_NOT_BRACKETED;
    }

    for (int iteration = 1; iteration <= options->maximum_iterations; ++iteration) {
        double temperature = lower + 0.5 * (upper - lower);
        ThermalState middle;
        status = calculateThermalState(parameters, temperature, &middle);
        if (status != THERMAL_SUCCESS) {
            return status;
        }
        middle.iterations = iteration;
        if (isBalanced(&middle, options->relative_tolerance) &&
            upper - lower <= options->relative_tolerance * temperature) {
            *state = middle;
            return THERMAL_SUCCESS;
        }
        if (temperature == lower || temperature == upper) {
            return THERMAL_NO_CONVERGENCE;
        }
        double residual = middle.cooling - middle.heating;
        if ((residual < 0.0) == (left_residual < 0.0)) {
            lower = temperature;
            left_residual = residual;
        } else {
            upper = temperature;
        }
    }
    return THERMAL_NO_CONVERGENCE;
}

const char *thermalStatusMessage(ThermalStatus status)
{
    switch (status) {
    case THERMAL_SUCCESS:
        return "success";
    case THERMAL_INVALID_INPUT:
        return "invalid thermal input or solver options";
    case THERMAL_NUMERICAL_ERROR:
        return "non-finite, underflowed or overflowed thermal quantity";
    case THERMAL_NOT_BRACKETED:
        return "temperature bounds do not bracket an equilibrium";
    case THERMAL_NO_CONVERGENCE:
        return "thermal solver did not meet temperature and energy tolerances";
    }
    return "unknown thermal status";
}
