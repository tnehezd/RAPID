#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "simulation_types.h"
#include "disk_model.h"
#include "gas_physics.h"
#include "simulation_core.h"
#include "logger.h"
#include "print_panels.h"
#include "print_terminal.h"

void runMassConservationTest(DiskParameters *disk_params, SimulationOptions *sim_opts) {

    // Construct path for mass summary inside the logs folder
    char *filepath = NULL;
    asprintf(&filepath, "%s/%s/mass_conservation_100k.dat", sim_opts->output_dir_name, kLogFilesDirectory);

    FILE *fp = fopen(filepath, "w");
    if (fp) {
        fprintf(fp, "# Physical_Time_yrs\tTotal_Mass_Msun\tRelative_Error\n");
    } else {
        LOG_ERROR("Could not open file %s for writing mass conservation report.", filepath);
    }

    double interval = sim_opts->output_frequency; // Output frequency in years
    double last_snapshot_time = 0.0; // Track the last snapshot time for progress calculation
    double initial_mass = 0.0;
    double output_time = interval; // Initialize output time for snapshot tracking
    for (int i = 1; i <= disk_params->grid_number; i++) {
        initial_mass += 2.0 * M_PI * disk_params->radial_grid[i] * disk_params->gas_surface_density_vector[i] * (disk_params->radial_grid[i + 1] - disk_params->radial_grid[i]);
    }

    // Save initial profile at t = 0 yrs in the logs folder
    char *init_prof_path = NULL;
    asprintf(&init_prof_path, "%s/%s/sigma_profile_t_000000.dat", sim_opts->output_dir_name, kLogFilesDirectory);
    FILE *init_prof_fp = fopen(init_prof_path, "w");
    if (init_prof_fp) {
        fprintf(init_prof_fp, "# Radius_AU\tGas_Surface_Density_Msun_per_AU2\n");
        for (int i = 1; i <= disk_params->grid_number; i++) {
            fprintf(init_prof_fp, "%.6e\t%.6e\n", disk_params->radial_grid[i], disk_params->gas_surface_density_vector[i]);
        }
        fclose(init_prof_fp);
    }
    free(init_prof_path);

    // Set timestep to 1.0 year (converted to internal code units: 1 yr = 2 * PI radians)
    double dt_years = 1.0; 
    sim_opts->user_defined_time_step = dt_years * 2.0 * M_PI;

    double target_time = sim_opts->maximum_simulation_time; // 100,000 years
    double current_time = 0.0;
    long step = 0;

    if (fp) {
        fprintf(fp, "%.1f\t%.10e\t0.00e+00\n", current_time, initial_mass);
    }


    while (current_time < target_time) {
        // Refresh gas surface density and pressure gradient for one yearly step
        refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);
        
        current_time += dt_years;
        step++;

        // Calculate total disk mass at the current step
        double current_mass = 0.0;
        for (int i = 1; i <= disk_params->grid_number; i++) {
            if (disk_params->gas_surface_density_vector[i] < 0.0) {
                disk_params->gas_surface_density_vector[i] = 0.0;
            }
            current_mass += 2.0 * M_PI * disk_params->radial_grid[i] * disk_params->gas_surface_density_vector[i] * (disk_params->radial_grid[i + 1] - disk_params->radial_grid[i]);
        }

        double relative_error = fabs(current_mass - initial_mass) / initial_mass;
        int was_snapshot = (fmod(current_time, interval) < dt_years || current_time == dt_years);
        printBenchmarkStatus("Mass Conservation Test",
                            current_time,
                            target_time,
                            current_mass,
                            initial_mass,
                            was_snapshot,
                            step,
                            dt_years,
                            output_time,
                            last_snapshot_time,
                            interval,
                            sim_opts);


        // Save data to file
        if (fp) {
            fprintf(fp, "%.1f\t%.10e\t%.10e\n", current_time, current_mass, relative_error);
        }

        // Save radial profiles periodically (every 20,000 years) in the logs folder
        if (fmod(current_time, interval) < dt_years || current_time == dt_years) {
            char *prof_path = NULL;
            asprintf(&prof_path, "%s/%s/sigma_profile_t_%06d.dat", sim_opts->output_dir_name, kLogFilesDirectory, (int)current_time);
            FILE *prof_fp = fopen(prof_path, "w");
            if (prof_fp) {
                fprintf(prof_fp, "# Radius_AU\tGas_Surface_Density_Msun_per_AU2\n");
                for (int i = 1; i <= disk_params->grid_number; i++) {
                    fprintf(prof_fp, "%.6e\t%.6e\n", disk_params->radial_grid[i], disk_params->gas_surface_density_vector[i]);
                }
                fclose(prof_fp);
            }
            free(prof_path);

            last_snapshot_time = current_time;
            output_time += interval;
        }
    }

    if (fp) {
        fclose(fp);
    }

    if (filepath) {
        free(filepath);
    }

}