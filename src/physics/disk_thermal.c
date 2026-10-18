#include "disk_thermal.h"

#include "config.h"
#include "disk_model.h"
#include "gas_physics.h"
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

void refreshRadiativeDisk(DiskParameters *disk, SimulationOptions *options)
{
    if (disk->thermodynamics_model != THERMODYNAMICS_RADIATIVE) return;
    int n = disk->grid_number;
    double time_unit = thermalTimeUnitSeconds();
    if (n < 3) {
        LOG_ERROR("Radiative model needs at least three radial cells.\n");
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

    for (int iteration = 0; iteration < 64; ++iteration) {
        updateSelfGravity(disk);
        ThermalState next[n + 2];
        double maximum_change = 0.0;
        for (int i = 0; i <= n + 1; ++i) {
            double sigma = disk->gas_surface_density_vector[i];
            if (!isfinite(sigma) || sigma < 0.0) {
                LOG_ERROR("Invalid surface density %.9g at thermal cell %d, r=%.9g AU.\n",
                          sigma, i, disk->radial_grid[i]);
                exit(EXIT_FAILURE);
            }
            /* Absorbing ghosts may have zero Sigma. Their thermal closure uses
             * the positive floor; the actual mass/pressure boundary is unchanged. */
            if ((i == 0 || i == n + 1) && sigma == 0.0) sigma = disk->density_floor;
            ThermalParameters parameters = {
                .surface_density = sigma * SOLAR_MASS_IN_GRAMS / (AU_IN_CM * AU_IN_CM),
                .orbital_frequency = sqrt(frequencySquared(i, disk)) / time_unit,
                .alpha = calculateTurbulentAlpha(disk->radial_grid[i], disk),
                .mean_molecular_weight = disk->thermal_mean_molecular_weight,
                .background_temperature = disk->thermal_background_temperature,
                .shear_parameter = shearParameter(i, disk)
            };
            ThermalStatus status = solveMidplaneTemperature(
                &parameters, &disk->thermal_solver_options, &next[i]);
            if (status != THERMAL_SUCCESS) {
                LOG_ERROR("Radiative equilibrium failed at cell %d, r=%.9g AU, Sigma=%.9g g/cm^2: %s (T bracket %.9g--%.9g K).\n",
                          i, disk->radial_grid[i], parameters.surface_density,
                          thermalStatusMessage(status),
                          disk->thermal_solver_options.minimum_temperature,
                          disk->thermal_solver_options.maximum_temperature);
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
            createInitialGasPressure(disk, options);
            createInitialGasPressureGradient(disk, options);
            createInitialGasVelocity(disk, options);
            return;
        }
    }
    LOG_ERROR("Coupled radiative/SG equilibrium did not converge in 64 iterations.\n");
    exit(EXIT_FAILURE);
}

void writeThermalSnapshot(const DiskParameters *disk, const SimulationOptions *options,
                          double time_years)
{
    if (disk->thermodynamics_model != THERMODYNAMICS_RADIATIVE) return;
    char filename[MAX_PATH_LEN];
    int length = snprintf(filename, sizeof(filename), "%s/%s/temperature_profile_%.12g.dat",
                          options->output_dir_name, kLogFilesDirectory, time_years);
    if (length < 0 || (size_t)length >= sizeof(filename)) {
        LOG_ERROR("Thermal output path is too long.\n");
        exit(EXIT_FAILURE);
    }
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        LOG_ERROR("Could not open thermal snapshot '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
    fprintf(file, "# Radiative equilibrium at %.15g years; SG=%s\n",
            time_years, disk->enable_selfgravity ? "on" : "off");
    fprintf(file, "# r_AU T_K H_over_r rho_g_cm3 nu_cm2_s kappa_cm2_g tau Qplus Qminus\n");
    for (int i = 1; i <= disk->grid_number; ++i) {
        const ThermalState *state = &disk->thermal_states[i];
        fprintf(file, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
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
}
