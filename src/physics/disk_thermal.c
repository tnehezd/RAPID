#include "disk_thermal.h"

#include "config.h"
#include "ascii_output.h"
#include "boundary_conditions.h"
#include "disk_model.h"
#include "gas_physics.h"
#include "self_gravity.h"
#include "logger.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double thermalTimeUnitSeconds(void)
{
    return YEARS_PER_DAY_CONVERSION_FACTOR * 86400.0 / (2.0 * M_PI);
}

static double interpolateThermalQuantity(double radius, const DiskParameters *disk,
                                         int sound_speed)
{
    if (disk->thermal_states == NULL || !isfinite(radius) || radius <= 0.0) {
        LOG_ERROR("Uninitialized radiative state or invalid radius %.9g AU.\n", radius);
        exit(EXIT_FAILURE);
    }
    int low = 0;
    int high = disk->grid_number + 1;
    const double *r = disk->radial_grid;
    if (radius <= r[low]) high = low;
    else if (radius >= r[high]) low = high;
    else {
        while (high - low > 1) {
            int middle = (low + high) / 2;
            if (r[middle] <= radius) low = middle;
            else high = middle;
        }
    }
    double left = sound_speed ? disk->thermal_states[low].sound_speed :
                                disk->thermal_states[low].scale_height;
    double right = sound_speed ? disk->thermal_states[high].sound_speed :
                                 disk->thermal_states[high].scale_height;
    double fraction = high == low ? 0.0 : (radius - r[low]) / (r[high] - r[low]);
    double result = left + fraction * (right - left);
    if (!isfinite(result) || result <= 0.0) {
        LOG_ERROR("Invalid cached radiative quantity at r=%.9g AU.\n", radius);
        exit(EXIT_FAILURE);
    }
    return result;
}

double interpolateThermalSoundSpeed(double radius, const DiskParameters *disk)
{
    return interpolateThermalQuantity(radius, disk, 1) * thermalTimeUnitSeconds() / AU_IN_CM;
}

double interpolateThermalScaleHeight(double radius, const DiskParameters *disk)
{
    return interpolateThermalQuantity(radius, disk, 0) / AU_IN_CM;
}

static double frequencySquared(int i, const DiskParameters *disk)
{
    double radius = disk->radial_grid[i];
    double value = G_DIMENSIONLESS * disk->stellar_mass / (radius * radius * radius);
    if (disk->enable_selfgravity) value += disk->selfgravity_omega2_extra[i];
    if (!isfinite(value) || value <= 0.0) {
        LOG_ERROR("Radiative model requires positive Omega^2 at cell %d, r=%.9g AU.\n",
                  i, radius);
        exit(EXIT_FAILURE);
    }
    return value;
}

static double shearParameter(int i, const DiskParameters *disk)
{
    if (!disk->enable_selfgravity) return 1.5;
    int n = disk->grid_number;
    int middle = i == 0 ? 1 : (i == n + 1 ? n : i);
    double x = disk->radial_grid[i];
    double derivative = 0.0;
    /* Differentiate only the SG correction; keep the stellar part exact. */
    for (int j = middle - 1; j <= middle + 1; ++j) {
        int k = middle - 1 + (j - (middle - 1) + 1) % 3;
        int l = middle - 1 + (j - (middle - 1) + 2) % 3;
        derivative += disk->selfgravity_omega2_extra[j] *
            (2.0 * x - disk->radial_grid[k] - disk->radial_grid[l]) /
            ((disk->radial_grid[j] - disk->radial_grid[k]) *
             (disk->radial_grid[j] - disk->radial_grid[l]));
    }
    double stellar_omega2 = G_DIMENSIONLESS * disk->stellar_mass / (x * x * x);
    return fabs((-3.0 * stellar_omega2 + x * derivative) /
                (2.0 * frequencySquared(i, disk)));
}

static ThermalParameters cellThermalParameters(int i, double sigma, const DiskParameters *disk)
{
    ThermalParameters parameters = {
        .surface_density = sigma * SOLAR_MASS_IN_GRAMS / (AU_IN_CM * AU_IN_CM),
        .orbital_frequency = sqrt(frequencySquared(i, disk)) / thermalTimeUnitSeconds(),
        .alpha = calculateTurbulentAlpha(disk->radial_grid[i], disk),
        .mean_molecular_weight = disk->thermal_mean_molecular_weight,
        .background_temperature = disk->thermal_background_temperature,
        .shear_parameter = shearParameter(i, disk)
    };
    if (disk->irradiation_model == IRRADIATION_GEOMETRY) {
        parameters.irradiation_temperature = disk->thermal_states[i].irradiation_temperature;
    } else if (disk->irradiation_luminosity_lsun > 0.0 && disk->irradiation_factor > 0.0) {
        double radius_cm = disk->radial_grid[i] * AU_IN_CM;
        parameters.irradiation_temperature = exp(0.25 * (
            log(disk->irradiation_luminosity_lsun) + log(SOLAR_LUMINOSITY_CGS) +
            log(disk->irradiation_factor) - log(4.0 * M_PI * STEFAN_BOLTZMANN_CONSTANT_CGS) -
            2.0 * log(radius_cm)));
    }
    return parameters;
}

static void writeThermalTimescales(const DiskParameters *disk, const SimulationOptions *options,
                                  double time_years, double output_time)
{
    char filename[MAX_PATH_LEN];
    int length = snprintf(filename, sizeof(filename), "%s/%s/thermal_timescales_%08d.dat",
                          options->output_dir_name, kLogFilesDirectory, (int)output_time);
    if (length < 0 || (size_t)length >= sizeof(filename)) {
        LOG_ERROR("Thermal timescale output path is too long.\n");
        exit(EXIT_FAILURE);
    }
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        LOG_ERROR("Could not open thermal timescales '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
    double cv = BOLTZMANN_CONSTANT_CGS /
        (disk->thermal_mean_molecular_weight * HYDROGEN_MASS_IN_GRAMS * (disk->thermal_gamma - 1.0));
    double year = YEARS_PER_DAY_CONVERSION_FACTOR * 86400.0;
    HeaderData header = {.current_time = time_years};
    printFileHeader(file, FILE_TYPE_THERMAL_TIMESCALES, &header);
    for (int i = 1; i <= disk->grid_number; ++i) {
        ThermalParameters parameters = cellThermalParameters(i, disk->gas_surface_density_vector[i], disk);
        const ThermalState *state = &disk->thermal_states[i];
        double capacity = parameters.surface_density * cv;
        double energy = capacity * state->temperature;
        double delta = 1e-4 * state->temperature;
        ThermalState lower, upper;
        ThermalStatus low = calculateThermalState(&parameters, state->temperature - delta, &lower);
        ThermalStatus high = calculateThermalState(&parameters, state->temperature + delta, &upper);
        if (low != THERMAL_SUCCESS || high != THERMAL_SUCCESS) {
            LOG_ERROR("Thermal relaxation derivative failed at cell %d: %s / %s.\n",
                      i, thermalStatusMessage(low), thermalStatusMessage(high));
            fclose(file);
            exit(EXIT_FAILURE);
        }
        double derivative = ((upper.cooling - upper.heating) -
                             (lower.cooling - lower.heating)) / (2.0 * delta);
        int sign = derivative > 0.0 ? 1 : (derivative < 0.0 ? -1 : 0);
        double cooling_time = state->cooling > 0.0 ? energy / state->cooling / year : 0.0;
        double relaxation_time = sign != 0 ? capacity / fabs(derivative) / year : 0.0;
        if (!isfinite(cooling_time) || !isfinite(relaxation_time) || !isfinite(derivative)) {
            LOG_ERROR("Non-finite thermal timescale at cell %d.\n", i);
            fclose(file);
            exit(EXIT_FAILURE);
        }
        fprintf(file, "  %-23.15e %-23.15e %-23.15e %-23d %-23.15e %-23d %-23.15e\n",
                disk->radial_grid[i], energy, cooling_time, state->cooling > 0.0,
                relaxation_time, sign, derivative);
    }
    int failed = ferror(file);
    if (fclose(file) != 0) failed = 1;
    if (failed) {
        LOG_ERROR("Failed to write thermal timescales '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
}

void refreshRadiativeDisk(DiskParameters *disk, SimulationOptions *options)
{
    if (!isRadiativeThermodynamics(disk->thermodynamics_model)) return;
    int initialize_temperature = disk->thermal_states == NULL;
    int n = disk->grid_number;
    if (n < 3) {
        LOG_ERROR("Radiative model needs at least three radial cells.\n");
        exit(EXIT_FAILURE);
    }
    if (disk->irradiation_model == IRRADIATION_GEOMETRY &&
        (!isfinite(disk->irradiation_stellar_radius_rsun) ||
         disk->irradiation_stellar_radius_rsun <= 0.0 ||
         disk->irradiation_stellar_radius_rsun * SOLAR_RADIUS_IN_CM >=
         disk->radial_grid[1] * AU_IN_CM)) {
        LOG_ERROR("Geometric irradiation requires a positive stellar radius smaller than the inner physical grid radius.\n");
        exit(EXIT_FAILURE);
    }
    if (disk->thermal_states == NULL) {
        disk->thermal_states = calloc((size_t)n + 2, sizeof(*disk->thermal_states));
        if (disk->thermal_states == NULL) {
            LOG_ERROR("Failed to allocate radiative disk state.\n");
            exit(EXIT_FAILURE);
        }
        /* Initial guess only, for the SG kernel's temperature-dependent softening. */
        for (int i = 0; i <= n + 1; ++i) {
            double radius = disk->radial_grid[i];
            double height = disk->h_aspect_ratio * pow(radius, 1.0 + disk->flaring_index);
            if (!isfinite(radius) || radius <= 0.0 || !isfinite(height) || height <= 0.0) {
                LOG_ERROR("Radiative initialization requires positive radii and initial H; use a finer or logarithmic grid.\n");
                exit(EXIT_FAILURE);
            }
            disk->thermal_states[i].scale_height = height * AU_IN_CM;
        }
    }

    int geometric = disk->irradiation_model == IRRADIATION_GEOMETRY;
    double radii[n + 2], heights[n + 2];
    int shadow_hull[n];
    for (int i = 0; i <= n + 1; ++i) radii[i] = disk->radial_grid[i] * AU_IN_CM;
    for (int iteration = 0; iteration < 64; ++iteration) {
        updateSelfGravity(disk);
        ThermalState next[n + 2];
        for (int i = 0; i <= n + 1; ++i) heights[i] = disk->thermal_states[i].scale_height;
        double maximum_change = 0.0;
        int shadow_hull_count = 0;
        for (int index = 0; index <= n + 1; ++index) {
            /* Geometry is causal in radius; solve physical cells before ghosts. */
            int i = geometric ? (index < n ? index + 1 : (index == n ? 0 : n + 1)) : index;
            double sigma = disk->gas_surface_density_vector[i];
            if (!isfinite(sigma) || sigma < 0.0) {
                LOG_ERROR("Invalid surface density %.9g at thermal cell %d, r=%.9g AU.\n",
                          sigma, i, disk->radial_grid[i]);
                exit(EXIT_FAILURE);
            }
            /* Absorbing ghosts may have zero Sigma. Their thermal closure uses
             * the positive floor; the actual mass/pressure boundary is unchanged. */
            if ((i == 0 || i == n + 1) && sigma == 0.0) sigma = disk->density_floor;
            ThermalParameters parameters = cellThermalParameters(i, sigma, disk);
            IrradiationGeometry geometry = {
                .radii = radii, .heights = heights, .index = i,
                .surface_height_factor = disk->irradiation_surface_height,
                .stellar_radius = disk->irradiation_stellar_radius_rsun * SOLAR_RADIUS_IN_CM,
                .shadow_hull = shadow_hull, .shadow_hull_count = shadow_hull_count
            };
            if (geometric && i >= 1 && i <= n) {
                parameters.irradiation_temperature = 0.0;
                parameters.irradiation_geometry = &geometry;
                parameters.irradiation_luminosity =
                    disk->irradiation_luminosity_lsun * SOLAR_LUMINOSITY_CGS;
            }
            ThermalStatus status;
            if (geometric && (i == 0 || i == n + 1)) {
                int adjacent = i == 0 ? 1 : n;
                parameters.irradiation_temperature = next[adjacent].irradiation_temperature;
                status = calculateThermalState(&parameters, next[adjacent].temperature, &next[i]);
            } else if (disk->thermodynamics_model == THERMODYNAMICS_RADIATIVE_ENERGY &&
                initialize_temperature && (i == 0 || i == n + 1)) {
                /* Ghost temperatures are replaced by adjacent physical T below.
                 * In particular, a massless absorbing ghost needs no equilibrium. */
                status = calculateThermalState(&parameters,
                    disk->thermal_background_temperature, &next[i]);
            } else if (disk->thermodynamics_model == THERMODYNAMICS_RADIATIVE_ENERGY &&
                       !initialize_temperature) {
                status = calculateThermalState(&parameters,
                    disk->thermal_states[i].temperature, &next[i]);
            } else {
                status = solveMidplaneTemperature(
                    &parameters, &disk->thermal_solver_options, &next[i]);
            }
            if (status != THERMAL_SUCCESS) {
                LOG_ERROR("Radiative %s failed at cell %d, r=%.9g AU, Sigma=%.9g g/cm^2: %s (T bracket %.9g--%.9g K).\n",
                          disk->thermodynamics_model == THERMODYNAMICS_RADIATIVE_ENERGY &&
                          !initialize_temperature ? "closure" : "equilibrium",
                          i, disk->radial_grid[i], parameters.surface_density,
                          thermalStatusMessage(status),
                          disk->thermal_solver_options.minimum_temperature,
                          disk->thermal_solver_options.maximum_temperature);
                exit(EXIT_FAILURE);
            }
            heights[i] = next[i].scale_height;
            if (geometric && i >= 1 && i <= n &&
                !appendIrradiationShadowHull(&geometry, shadow_hull, &shadow_hull_count)) {
                LOG_ERROR("Invalid irradiation shadow hull at cell %d.\n", i);
                exit(EXIT_FAILURE);
            }
            maximum_change = fmax(maximum_change, fabs(next[i].scale_height /
                                            disk->thermal_states[i].scale_height - 1.0));
        }
        int converged = !disk->enable_selfgravity || maximum_change < 1e-8;
        for (int i = 0; i <= n + 1; ++i) {
            double old_height = disk->thermal_states[i].scale_height;
            disk->thermal_states[i] = next[i];
            if (!converged) {
                disk->thermal_states[i].scale_height =
                    0.5 * (old_height + next[i].scale_height);
            }
        }
        if (converged) {
            if (initialize_temperature &&
                disk->thermodynamics_model == THERMODYNAMICS_RADIATIVE_ENERGY) {
                disk->thermal_states[0].temperature = disk->thermal_states[1].temperature;
                disk->thermal_states[n + 1].temperature = disk->thermal_states[n].temperature;
                refreshRadiativeDisk(disk, options);
                return;
            }
            createInitialGasPressure(disk, options);
            createInitialGasPressureGradient(disk, options);
            createInitialGasVelocity(disk, options);
            if (disk->thermodynamics_model == THERMODYNAMICS_RADIATIVE_ENERGY) {
                for (int i = 1; i <= n; ++i) {
                    disk->gas_velocity_vector[i] =
                        0.5 * (radiativeMassFlux(i - 1, disk, options) +
                               radiativeMassFlux(i, disk, options)) /
                        (disk->radial_grid[i] * disk->gas_surface_density_vector[i]);
                }
                options->current_bc_target = 0;
                applyBoundaryConditions(disk->gas_velocity_vector, disk, options);
            }
            return;
        }
    }
    LOG_ERROR("Coupled radiative/SG equilibrium did not converge in 64 iterations.\n");
    exit(EXIT_FAILURE);
}

static void writeIrradiationSnapshot(const DiskParameters *disk, const SimulationOptions *options,
                                    double time_years, double output_time)
{
    char filename[MAX_PATH_LEN];
    int length = snprintf(filename, sizeof(filename), "%s/%s/irradiation_profile_%08d.dat",
                          options->output_dir_name, kLogFilesDirectory, (int)output_time);
    if (length < 0 || (size_t)length >= sizeof(filename)) {
        LOG_ERROR("Irradiation output path is too long.\n");
        exit(EXIT_FAILURE);
    }
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        LOG_ERROR("Could not open irradiation snapshot '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
    HeaderData header = {.current_time = time_years, .disk_params = disk};
    printFileHeader(file, FILE_TYPE_IRRADIATION_PROFILE, &header);
    for (int i = 1; i <= disk->grid_number; ++i) {
        const ThermalState *state = &disk->thermal_states[i];
        fprintf(file, "  %-25.15e %-25.15e %-25.15e %-25.15e %-25.15e\n", disk->radial_grid[i],
                disk->irradiation_surface_height * state->scale_height / AU_IN_CM,
                state->irradiation_factor, state->stellar_visibility, state->irradiation_temperature);
    }
    int failed = ferror(file);
    if (fclose(file) != 0) failed = 1;
    if (failed) {
        LOG_ERROR("Failed to write irradiation snapshot '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
}

void writeThermalSnapshot(const DiskParameters *disk, const SimulationOptions *options,
                          double time_years, double output_time)
{
    if (!isRadiativeThermodynamics(disk->thermodynamics_model)) return;
    char filename[MAX_PATH_LEN];
    int length = snprintf(filename, sizeof(filename), "%s/%s/temperature_profile_%08d.dat",
                          options->output_dir_name, kLogFilesDirectory, (int)output_time);
    if (length < 0 || (size_t)length >= sizeof(filename)) {
        LOG_ERROR("Thermal output path is too long.\n");
        exit(EXIT_FAILURE);
    }
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        LOG_ERROR("Could not open thermal snapshot '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
    HeaderData header = {.current_time = time_years, .disk_params = disk};
    printFileHeader(file, FILE_TYPE_THERMAL_PROFILE, &header);
    for (int i = 1; i <= disk->grid_number; ++i) {
        const ThermalState *state = &disk->thermal_states[i];
        fprintf(file, "  %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e\n",
                disk->radial_grid[i], state->temperature,
                state->scale_height / (disk->radial_grid[i] * AU_IN_CM),
                state->midplane_density, state->viscosity, state->opacity,
                state->optical_depth, state->heating, state->cooling);
    }
    int failed = ferror(file);
    if (fclose(file) != 0) failed = 1;
    if (failed) {
        LOG_ERROR("Failed to write thermal snapshot '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
    if (disk->thermodynamics_model == THERMODYNAMICS_RADIATIVE_ENERGY) {
        writeEnergyBudget(disk, options, time_years, output_time);
        writeThermalTimescales(disk, options, time_years, output_time);
    }
    if (disk->irradiation_model == IRRADIATION_GEOMETRY) {
        writeIrradiationSnapshot(disk, options, time_years, output_time);
    }
}
