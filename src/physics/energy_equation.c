#include "disk_thermal.h"
#include "ascii_output.h"

#include "boundary_conditions.h"
#include "config.h"
#include "gas_physics.h"
#include "logger.h"
#include "photoevaporation.h"

#include <math.h>
#include <stdlib.h>

double totalDiskInternalEnergy(const DiskParameters *disk)
{
    double total = 0.0;
    double unit = AU_IN_CM / thermalTimeUnitSeconds();
    double cv = BOLTZMANN_CONSTANT_CGS /
        (disk->thermal_mean_molecular_weight * HYDROGEN_MASS_IN_GRAMS *
         (disk->thermal_gamma - 1.0) * unit * unit);
    for (int i = 1; i <= disk->grid_number; ++i) {
        double left = 0.5 * (disk->radial_grid[i - 1] + disk->radial_grid[i]);
        double right = 0.5 * (disk->radial_grid[i] + disk->radial_grid[i + 1]);
        total += M_PI * (right * right - left * left) *
                 disk->gas_surface_density_vector[i] * cv * disk->thermal_states[i].temperature;
    }
    return total;
}

void writeEnergyBudget(const DiskParameters *disk, const SimulationOptions *options,
                       double time_years, double output_time)
{
    char filename[MAX_PATH_LEN];
    int length = snprintf(filename, sizeof(filename), "%s/%s/energy_budget_%08d.dat",
                          options->output_dir_name, kLogFilesDirectory, (int)output_time);
    if (length < 0 || (size_t)length >= sizeof(filename)) {
        LOG_ERROR("Energy budget output path is too long.\n");
        exit(EXIT_FAILURE);
    }
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        LOG_ERROR("Could not open energy budget '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
    double unit = AU_IN_CM / thermalTimeUnitSeconds();
    double conversion = SOLAR_MASS_IN_GRAMS * unit * unit;
    double total = totalDiskInternalEnergy(disk);
    double initial = disk->energy_budget.initialized ? disk->energy_budget.initial : total;
    double expected = disk->energy_budget.heating - disk->energy_budget.cooling +
        disk->energy_budget.work + disk->energy_budget.boundary -
        disk->energy_budget.wind + disk->energy_budget.floor;
    double residual = total - initial - expected;
    double scale = fmax(fabs(initial), fabs(total));
    HeaderData header = {.current_time = time_years};
    printFileHeader(file, FILE_TYPE_ENERGY_BUDGET, &header);
    fprintf(file, "  %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e %-23.15e\n",
            time_years, initial * conversion, total * conversion,
            disk->energy_budget.heating * conversion, disk->energy_budget.cooling * conversion,
            disk->energy_budget.work * conversion, disk->energy_budget.boundary * conversion,
            disk->energy_budget.wind * conversion, disk->energy_budget.floor * conversion,
            residual * conversion, scale > 0.0 ? residual / scale : 0.0);
    int failed = ferror(file);
    if (fclose(file) != 0) failed = 1;
    if (failed) {
        LOG_ERROR("Failed to write energy budget '%s'.\n", filename);
        exit(EXIT_FAILURE);
    }
}

double radiativeMassFlux(int j, const DiskParameters *disk, const SimulationOptions *options)
{
    if ((j == 0 && options->inner_boundary_condition_type == 4) ||
        (j == disk->grid_number && options->outer_boundary_condition_type == 4)) return 0.0;
    double left_r = disk->radial_grid[j];
    double right_r = disk->radial_grid[j + 1];
    double left = calculateKinematicViscosity(left_r, disk) *
                  disk->gas_surface_density_vector[j] * sqrt(left_r);
    double right = calculateKinematicViscosity(right_r, disk) *
                   disk->gas_surface_density_vector[j + 1] * sqrt(right_r);
    double flux = -3.0 * sqrt(0.5 * (left_r + right_r)) * (right - left) / (right_r - left_r);
    /* Outflow is relative to the domain: negative at the inner edge,
     * positive at the outer edge. Use the limited flux for energy and work too. */
    if (j == 0 && options->inner_boundary_condition_type == 5) return fmin(flux, 0.0);
    if (j == disk->grid_number && options->outer_boundary_condition_type == 5) return fmax(flux, 0.0);
    return flux;
}

void evolveRadiativeEnergy(DiskParameters *disk, SimulationOptions *options)
{
    int n = disk->grid_number;
    double remaining = options->user_defined_time_step;
    double time_unit = thermalTimeUnitSeconds();
    double velocity_unit = AU_IN_CM / time_unit;
    double cv = BOLTZMANN_CONSTANT_CGS /
        (disk->thermal_mean_molecular_weight * HYDROGEN_MASS_IN_GRAMS *
         (disk->thermal_gamma - 1.0) * velocity_unit * velocity_unit);
    double source_conversion = time_unit * time_unit * time_unit / SOLAR_MASS_IN_GRAMS;
    if (!isfinite(remaining) || remaining <= 0.0 || !isfinite(cv) || cv <= 0.0) {
        LOG_ERROR("Invalid energy timestep or ideal-gas heat capacity.\n");
        exit(EXIT_FAILURE);
    }
    if (!disk->energy_budget.initialized) {
        disk->energy_budget.initial = totalDiskInternalEnergy(disk);
        disk->energy_budget.initialized = true;
    }

    double face[n + 1], mass_flux[n + 1], energy_flux[n + 1], volume_flux[n + 1];
    double energy[n + 2], source[n + 2], sigma_next[n + 2], temperature_next[n + 2];
    static int warned_floor = 0;
    for (int step = 0; remaining > 0.0; ++step) {
        if (step >= 100000) {
            LOG_ERROR("Energy subcycling exceeded 100000 steps; reduce the outer timestep.\n");
            exit(EXIT_FAILURE);
        }
        computePhotoevaporationSink(disk);
        for (int i = 0; i <= n + 1; ++i) {
            energy[i] = disk->gas_surface_density_vector[i] * cv *
                        disk->thermal_states[i].temperature;
        }
        for (int j = 0; j <= n; ++j) {
            double left_r = disk->radial_grid[j];
            double right_r = disk->radial_grid[j + 1];
            face[j] = 0.5 * (left_r + right_r);
            /* Outward r*Sigma*v_r; the same face flux transports mass and energy. */
            mass_flux[j] = radiativeMassFlux(j, disk, options);
            int donor = mass_flux[j] >= 0.0 ? j : j + 1;
            energy_flux[j] = mass_flux[j] * cv * disk->thermal_states[donor].temperature;
            double face_sigma = 0.5 * (disk->gas_surface_density_vector[j] +
                                       disk->gas_surface_density_vector[j + 1]);
            volume_flux[j] = mass_flux[j] / fmax(face_sigma, disk->density_floor);
        }

        double dt = remaining;
        for (int i = 1; i <= n; ++i) {
            if (!isfinite(disk->sigma_dot_photoevap[i]) || disk->sigma_dot_photoevap[i] < 0.0) {
                LOG_ERROR("Invalid photoevaporation sink at energy cell %d.\n", i);
                exit(EXIT_FAILURE);
            }
            double area = 0.5 * (face[i] * face[i] - face[i - 1] * face[i - 1]);
            double divergence = (volume_flux[i] - volume_flux[i - 1]) / area;
            source[i] = (disk->thermal_states[i].heating - disk->thermal_states[i].cooling) *
                        source_conversion - (disk->thermal_gamma - 1.0) * energy[i] * divergence;
            double outgoing = fmax(mass_flux[i], 0.0) + fmax(-mass_flux[i - 1], 0.0);
            double sigma = disk->gas_surface_density_vector[i];
            if (outgoing > 0.0) dt = fmin(dt, 0.2 * area * sigma / outgoing);
            double temperature = disk->thermal_states[i].temperature;
            double emission = 2.0 * STEFAN_BOLTZMANN_CONSTANT_CGS * pow(temperature, 4) /
                              disk->thermal_states[i].effective_optical_depth;
            /* Net cooling vanishes near the radiation floor, but its thermal
             * relaxation rate does not. Retain emission and absorption separately. */
            double source_magnitude =
                (fabs(disk->thermal_states[i].heating) + emission +
                 fabs(emission - disk->thermal_states[i].cooling)) *
                source_conversion + fabs((disk->thermal_gamma - 1.0) * energy[i] * divergence);
            if (source_magnitude > 0.0) dt = fmin(dt, 0.1 * energy[i] / source_magnitude);
            if (disk->sigma_dot_photoevap[i] > 0.0 && sigma > 1.01 * disk->density_floor) {
                dt = fmin(dt, 0.1 * sigma / disk->sigma_dot_photoevap[i]);
            }
            double dr = fmin(disk->radial_grid[i] - disk->radial_grid[i - 1],
                             disk->radial_grid[i + 1] - disk->radial_grid[i]);
            double nu = calculateKinematicViscosity(disk->radial_grid[i], disk);
            if (nu > 0.0) dt = fmin(dt, 0.1 * dr * dr / (3.0 * nu));
        }
        if (!isfinite(dt) || dt <= 0.0 || remaining - dt == remaining) {
            LOG_ERROR("Energy substep cannot advance time (dt=%.9g).\n", dt);
            exit(EXIT_FAILURE);
        }
        for (int i = 1; i <= n; ++i) {
            double area = 0.5 * (face[i] * face[i] - face[i - 1] * face[i - 1]);
            double sigma = disk->gas_surface_density_vector[i] -
                           dt * (mass_flux[i] - mass_flux[i - 1]) / area;
            double updated_energy = energy[i] -
                dt * (energy_flux[i] - energy_flux[i - 1]) / area + dt * source[i];
            /* Sinks remove the local specific internal energy; a density floor
             * adds mass at the transported temperature, not at zero energy. */
            double temperature = updated_energy / (sigma * cv);
            sigma_next[i] = fmax(disk->density_floor, sigma - dt * disk->sigma_dot_photoevap[i]);
            if (sigma_next[i] == disk->density_floor && !warned_floor) {
                LOG_WARN("Energy mode reached the density floor; floor mass carries local specific internal energy. "
                         "Subsequent floor warnings are suppressed.\n");
                warned_floor = 1;
            }
            if (!isfinite(sigma) || sigma <= 0.0 || !isfinite(temperature) ||
                temperature < disk->thermal_solver_options.minimum_temperature ||
                temperature > disk->thermal_solver_options.maximum_temperature) {
                LOG_ERROR("Invalid evolved energy at cell %d: Sigma=%.9g, T=%.9g K. "
                          "Check thermal bounds and timestep.\n", i, sigma, temperature);
                exit(EXIT_FAILURE);
            }
            temperature_next[i] = temperature;
            double weight = 2.0 * M_PI * area;
            double wind_removed = dt * disk->sigma_dot_photoevap[i];
            double raw_sigma = sigma - wind_removed;
            disk->energy_budget.heating += weight * dt *
                disk->thermal_states[i].heating * source_conversion;
            disk->energy_budget.cooling += weight * dt *
                disk->thermal_states[i].cooling * source_conversion;
            disk->energy_budget.work -= weight * dt * (disk->thermal_gamma - 1.0) *
                energy[i] * (volume_flux[i] - volume_flux[i - 1]) / area;
            disk->energy_budget.wind += weight * wind_removed * cv * temperature;
            disk->energy_budget.floor += weight * (sigma_next[i] - raw_sigma) * cv * temperature;
        }
        disk->energy_budget.boundary -= 2.0 * M_PI * dt * (energy_flux[n] - energy_flux[0]);
        options->current_bc_target = 0;
        applySurfaceDensityBoundaryConditions(sigma_next, disk, options);
        /* Open boundaries use zero-gradient specific energy; absorbing ghosts
         * carry no mass. Reflecting boundaries have exactly zero face flux. */
        temperature_next[0] = temperature_next[1];
        temperature_next[n + 1] = temperature_next[n];
        for (int i = 0; i <= n + 1; ++i) {
            disk->gas_surface_density_vector[i] = sigma_next[i];
            disk->thermal_states[i].temperature = temperature_next[i];
        }
        refreshRadiativeDisk(disk, options);
        remaining -= dt;
    }
}
