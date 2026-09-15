#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "simulation_types.h"
#include "disk_model.h"
#include "gas_physics.h"
#include "logger.h"
#include "photoevap_flux_test.h"
#include "photoevaporation_wrapper.h"
#include "simulation_core.h"
#include "print_panels.h"
#include "print_terminal.h"

void runPhotoevapFluxTest(DiskParameters *dp, SimulationOptions *opt)
{
    // --- Create flux output file ---
    int step = 0;
    char *filepath = NULL;
    asprintf(&filepath, "%s/%s/photoevap_flux_test.dat",
             opt->output_dir_name, kLogFilesDirectory);

    FILE *fp = fopen(filepath, "w");
    if (fp) {
        fprintf(fp, "# time_yrs\tM_disk\tM_loss_actual\tM_loss_integrated\n");
    }

    double max_drift_velocity = 0.0;  // gas-only benchmark → no dust
    double dt = calculateTimeStep(dp, max_drift_velocity);

    double interval = opt->output_frequency;
    double output_time = 0.0;
    double last_snapshot_time = 0.0;
    double current_time = 0.0;
    double target_time = opt->maximum_simulation_time;
    // --- Initial disk mass ---
    double initial_mass = 0.0;
    for (int i = 1; i <= dp->grid_number; i++) {
        initial_mass += dp->gas_surface_density_vector[i] *
                        2.0 * M_PI * dp->radial_grid[i] * dp->delta_r;
    }

    double integrated_mass_loss = 0.0;

    while (current_time < target_time) {

        // --- Compute photoevap sink ---
        computePhotoevaporationSink(dp);

        // --- Integrate sigma_dot over radius ---
        double sigma_dot_integral = 0.0;
        for (int i = 1; i <= dp->grid_number; i++) {
            sigma_dot_integral += dp->sigma_dot_photoevap[i] *
                                  2.0 * M_PI * dp->radial_grid[i] * dp->delta_r;
        }

        // --- Update gas surface density ---
        for (int i = 1; i <= dp->grid_number; i++) {
            dp->gas_surface_density_vector[i] -= dt * dp->sigma_dot_photoevap[i];
            if (dp->gas_surface_density_vector[i] < dp->density_floor)
                dp->gas_surface_density_vector[i] = dp->density_floor;
        }

        // --- Compute new disk mass ---
        double current_mass = 0.0;
        for (int i = 1; i <= dp->grid_number; i++) {
            current_mass += dp->gas_surface_density_vector[i] *
                            2.0 * M_PI * dp->radial_grid[i] * dp->delta_r;
        }

        // --- Accumulate expected mass loss ---
        integrated_mass_loss += sigma_dot_integral * dt;


        // --- Snapshot logic identical to main solver ---
        int periodic_output_time = (fmod(current_time, interval) < dt);
        int initial_output_time  = (current_time == 0.0);
        int output_time_sync     = ((output_time - current_time) < dt);

        int was_snapshot = (periodic_output_time || initial_output_time) && output_time_sync;

        printBenchmarkStatus("Photoevaporation Flux Test",
                            current_time,
                            target_time,
                            current_mass,
                            initial_mass,
                            was_snapshot,
                            step,
                            dt,
                            output_time,
                            last_snapshot_time,
                            interval);

        if ((periodic_output_time || initial_output_time) && output_time_sync) {

            // --- Save sigma profile ---
            char *prof_path = NULL;
            asprintf(&prof_path, "%s/%s/sigma_profile_t_%06d.dat",
                     opt->output_dir_name, kLogFilesDirectory, (int)output_time);

            FILE *prof_fp = fopen(prof_path, "w");
            if (prof_fp) {
                fprintf(prof_fp, "# Radius_AU\tGas_Surface_Density_Msun_per_AU2\n");
                for (int i = 1; i <= dp->grid_number; i++) {
                    fprintf(prof_fp, "%.6e\t%.6e\n",
                            dp->radial_grid[i],
                            dp->gas_surface_density_vector[i]);
                }
                fclose(prof_fp);
            }
            free(prof_path);


            char *dot_path = NULL;
            asprintf(&dot_path, "%s/%s/sigma_dot_profile_t_%06d.dat",
                    opt->output_dir_name, kLogFilesDirectory, (int)current_time);

            FILE *dot_fp = fopen(dot_path, "w");
            if (dot_fp) {
                fprintf(dot_fp, "# Radius_AU   SigmaDot_Msun_per_AU2_per_day\n");
                for (int i = 1; i <= dp->grid_number; i++) {
                    fprintf(dot_fp, "%.6e %.6e\n",
                            dp->radial_grid[i],
                            dp->sigma_dot_photoevap[i]);
                }
                fclose(dot_fp);
            }
            free(dot_path);

            // --- Advance snapshot time ---
            output_time += interval;
        }

        // --- Write flux benchmark line ---
        if (fp) {
            fprintf(fp, "%.1f\t%.10e\t%.10e\t%.10e\n",
                    current_time,
                    current_mass,
                    initial_mass - current_mass,
                    integrated_mass_loss);
        }

        if(was_snapshot) {
            last_snapshot_time = current_time;
        }

        // --- Advance time ---
        current_time += dt;
        step++;
        dt = calculateTimeStep(dp, max_drift_velocity);  // recalc dt each step
    }

    if (fp) fclose(fp);
    if (filepath) free(filepath);

}
