#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "simulation_types.h"
#include "disk_model.h"
#include "gas_physics.h"
#include "simulation_core.h"
#include "logger.h"

void runMassConservationTest(DiskParameters *disk_params, SimulationOptions *sim_opts) {
    MSG("=== [BENCHMARK] Starting Pure Gas Mass Conservation Test (100k yrs) ===");

    // Construct path for mass summary inside the logs folder
    char *filepath = NULL;
    asprintf(&filepath, "%s/test_outputs/mass_conservation_100k.dat", sim_opts->output_dir_name);

    FILE *fp = fopen(filepath, "w");
    if (fp) {
        fprintf(fp, "# Physical_Time_yrs\tTotal_Mass_Msun\tRelative_Error\n");
    } else {
        LOG_ERROR("Could not open file %s for writing mass conservation report.", filepath);
    }

    double initial_mass = 0.0;
    for (int i = 1; i <= disk_params->grid_number; i++) {
        initial_mass += 2.0 * M_PI * disk_params->radial_grid[i] * disk_params->gas_surface_density_vector[i] * disk_params->delta_r;
    }
    MSG("[BENCHMARK] Initial total disk mass: %.10e M_sun", initial_mass);

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

    double target_time = 100000.0; // 100,000 years
    double current_time = 0.0;
    long step = 0;

    if (fp) {
        fprintf(fp, "%.1f\t%.10e\t0.00e+00\n", current_time, initial_mass);
    }

    MSG("[BENCHMARK] Running evolution up to %.1f years with dt = %.1f yr...", target_time, dt_years);

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
            current_mass += 2.0 * M_PI * disk_params->radial_grid[i] * disk_params->gas_surface_density_vector[i] * disk_params->delta_r;
        }

        double relative_error = fabs(current_mass - initial_mass) / initial_mass;

        // Console logging every 5000 years or at the first step
        if ((long)current_time % 5000 == 0 || current_time == dt_years) {
            MSG("[BENCHMARK] Time = %.0f yrs | Mass = %.10e | Rel. Error = %.2e", 
                current_time, current_mass, relative_error);
        }

        // Save data to file
        if (fp) {
            fprintf(fp, "%.1f\t%.10e\t%.10e\n", current_time, current_mass, relative_error);
        }

        // Save radial profiles periodically (every 20,000 years) in the logs folder
        if ((long)current_time % 20000 == 0) {
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
        }
    }

    if (fp) {
        fclose(fp);
        MSG("[BENCHMARK] Data saved to '%s'.", filepath);
    }

    if (filepath) {
        free(filepath);
    }

    MSG("=== [BENCHMARK] Long-Term Mass Conservation Test Completed ===");
}