#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "simulation_types.h"
#include "disk_model.h"
#include "gas_physics.h"
#include "simulation_core.h"
#include "logger.h"


void runRingViscosityTest(DiskParameters *disk_params, SimulationOptions *sim_opts) {
    MSG("=== [BENCHMARK] Starting Viscous Ring Spreading Test ===");

    // Construct path for the summary output inside the logs/test folder
    char *filepath = NULL;
    asprintf(&filepath, "%s/%s/ring_viscosity_summary.dat", sim_opts->output_dir_name, kLogFilesDirectory);

    FILE *fp = fopen(filepath, "w");
    if (fp) {
        fprintf(fp, "# Physical_Time_yrs\tTotal_Mass_Msun\tMax_Density_Msun_per_AU2\tPeak_Radius_AU\n");
    } else {
        LOG_ERROR("Could not open file %s for writing ring viscosity report.", filepath);
    }

    // --- INITIALIZE A NARROW GAUSSIAN RING PROFILE ---
    double r_center = 10.0;  // Ring center [AU]
    double ring_width = 0.5; // Width parameter [AU]
    
    // Overwrite the gas surface density vector with a localized ring profile
    for (int i = 1; i <= disk_params->grid_number; i++) {
        double r = disk_params->radial_grid[i];
        double dr_diff = r - r_center;
        // Gaussian ring profile normalized roughly to give a sensible total mass
        disk_params->gas_surface_density_vector[i] = 1.0 * exp(-(dr_diff * dr_diff) / (2.0 * ring_width * ring_width));
    }

    double initial_mass = 0.0;
    for (int i = 1; i <= disk_params->grid_number; i++) {
        initial_mass += 2.0 * M_PI * disk_params->radial_grid[i] * disk_params->gas_surface_density_vector[i] * disk_params->delta_r;
    }
    MSG("[BENCHMARK] Initial ring total mass: %.10e M_sun", initial_mass);

    // Save initial profile at t = 0 yrs
    char *init_prof_path = NULL;
    asprintf(&init_prof_path, "%s/%s/ring_profile_t_000000.dat", sim_opts->output_dir_name, kLogFilesDirectory);
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
    double max_drift_velocity = 0.0;  // gas-only benchmark → no dust
    double dt_years = calculateTimeStep(disk_params, max_drift_velocity);
    sim_opts->user_defined_time_step = dt_years;


    double target_time = sim_opts->maximum_simulation_time;  // Total simulation time in years
    double current_time = 0.0;
    long step = 0;
    double interval = sim_opts->output_frequency; // Output frequency in years

    MSG("[BENCHMARK] Running viscous ring evolution up to %.1f years with dt = %.1f yr...", target_time, dt_years);

    while (current_time < target_time) {
        // Evolve/refresh gas surface density profile for one yearly step
        refreshGasSurfaceDensityPressurePressureGradient(sim_opts, disk_params);
        
        current_time += dt_years;
        step++;

        // Calculate diagnostics: total mass and peak location
        double current_mass = 0.0;
        double max_density = 0.0;
        double peak_radius = r_center;

        for (int i = 1; i <= disk_params->grid_number; i++) {
            if (disk_params->gas_surface_density_vector[i] < 0.0) {
                disk_params->gas_surface_density_vector[i] = 0.0;
            }
            double sigma = disk_params->gas_surface_density_vector[i];
            current_mass += 2.0 * M_PI * disk_params->radial_grid[i] * sigma * disk_params->delta_r;
            
            if (sigma > max_density) {
                max_density = sigma;
                peak_radius = disk_params->radial_grid[i];
            }
        }

        // Console logging every 5000 years
        if (fmod(current_time, interval) < dt_years || current_time == dt_years) {
            MSG("[BENCHMARK] Time = %.0f yrs | Mass = %.10e | Max Sigma = %.4e at R = %.2f AU", 
                current_time, current_mass, max_density, peak_radius);
        }

        // Save summary data to file
        if (fp) {
            fprintf(fp, "%.1f\t%.10e\t%.10e\t%.6e\n", current_time, current_mass, max_density, peak_radius);
        }

        // Save radial profiles periodically (every 10,000 years)
        if (fmod(current_time, interval) < dt_years || current_time == dt_years) {
            char *prof_path = NULL;
            asprintf(&prof_path, "%s/%s/ring_profile_t_%06d.dat", sim_opts->output_dir_name, kLogFilesDirectory, (int)current_time);
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
        MSG("[BENCHMARK] Ring viscosity summary saved to '%s'.", filepath);
    }

    if (filepath) {
        free(filepath);
    }

    MSG("=== [BENCHMARK] Viscous Ring Spreading Test Completed ===");
}