#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "simulation_types.h"
#include "disk_model.h"
#include "gas_physics.h"
#include "logger.h"
#include "photoevap_flux_test.h"
#include "photoevaporation.h"
#include "simulation_core.h"
#include "print_panels.h"
#include "print_terminal.h"
#include "boundary_conditions.h"
#include "utils.h"

static void writeSigmaProfile(const DiskParameters *dp,
                              const SimulationOptions *opt,
                              double output_time)
{
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
}

static void writeSigmaDotProfile(const DiskParameters *dp,
                                 const SimulationOptions *opt,
                                 double output_time)
{
    char *dot_path = NULL;
    asprintf(&dot_path, "%s/%s/sigma_dot_profile_t_%06d.dat",
             opt->output_dir_name, kLogFilesDirectory, (int)output_time);

    FILE *dot_fp = fopen(dot_path, "w");
    if (dot_fp) {
        fprintf(dot_fp, "# Radius_AU   SigmaDot_Msun_per_AU2_per_year\n");
        for (int i = 1; i <= dp->grid_number; i++) {
            fprintf(dot_fp, "%.6e %.6e\n",
                    dp->radial_grid[i],
                    dp->sigma_dot_photoevap[i]);
        }
        fclose(dot_fp);
    }
    free(dot_path);
}

static void syncSigmaFromViscousState(DiskParameters *dp,
                                      const double *viscous_sigma)
{
    for (int i = 1; i <= dp->grid_number; i++) {
        double viscosity = calculateKinematicViscosity(dp->radial_grid[i], dp);
        double sigma = viscous_sigma[i] / viscosity;
        dp->gas_surface_density_vector[i] =
            (sigma < dp->density_floor) ? dp->density_floor : sigma;
    }
}

static double viscousSigmaOperator(const DiskParameters *dp,
                                   const double *viscous_sigma,
                                   int i)
{
    double left_dr = dp->radial_grid[i] - dp->radial_grid[i - 1];
    double right_dr = dp->radial_grid[i + 1] - dp->radial_grid[i];
    double span_dr = left_dr + right_dr;
    double first_derivative = (viscous_sigma[i + 1] - viscous_sigma[i - 1]) /
                              span_dr;
    double second_derivative = 2.0 *
        (viscous_sigma[i - 1] / (left_dr * span_dr) -
         viscous_sigma[i] / (left_dr * right_dr) +
         viscous_sigma[i + 1] / (right_dr * span_dr));

    return ftcsSecondDerivativeCoefficient(dp->radial_grid[i], dp) *
               second_derivative +
           ftcsFirstDerivativeCoefficient(dp->radial_grid[i], dp) *
               first_derivative;
}

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
                        2.0 * M_PI * dp->radial_grid[i] *
                        (dp->radial_grid[i + 1] - dp->radial_grid[i]);
    }

    double integrated_mass_loss = 0.0;

    double viscous_sigma[dp->grid_number + 2];
    double viscous_sigma_next[dp->grid_number + 2];
    for (int i = 0; i <= dp->grid_number + 1; i++) {
        viscous_sigma[i] = dp->gas_surface_density_vector[i] *
                           calculateKinematicViscosity(dp->radial_grid[i], dp);
    }
    opt->current_bc_target = 0;
    applyBoundaryConditions(viscous_sigma, dp, opt);
    syncSigmaFromViscousState(dp, viscous_sigma);

    // Save the untouched initial condition before the first evolution step.
    writeSigmaProfile(dp, opt, output_time);
    output_time += interval;

    while (current_time < target_time) {

        syncSigmaFromViscousState(dp, viscous_sigma);

        // --- Compute photoevap sink ---
        computePhotoevaporationSink(dp);

        if (current_time == 0.0) {
            writeSigmaDotProfile(dp, opt, 0.0);
        }

        // --- Integrate sigma_dot over radius ---
        double sigma_dot_integral = 0.0;
        for (int i = 1; i <= dp->grid_number; i++) {
            sigma_dot_integral += dp->sigma_dot_photoevap[i] *
                                  2.0 * M_PI * dp->radial_grid[i] *
                                  (dp->radial_grid[i + 1] - dp->radial_grid[i]);
        }

        // Evolve viscous_sigma = nu*Sigma, including the photoevaporation sink.
        for (int i = 1; i <= dp->grid_number; i++) {
            double viscosity = calculateKinematicViscosity(dp->radial_grid[i], dp);
            viscous_sigma_next[i] = viscous_sigma[i] + dt *
                (viscousSigmaOperator(dp, viscous_sigma, i) -
                 viscosity * dp->sigma_dot_photoevap[i]);
        }

        for (int i = 1; i <= dp->grid_number; i++) {
            viscous_sigma[i] = viscous_sigma_next[i];
        }

        opt->current_bc_target = 0;
        applyBoundaryConditions(viscous_sigma, dp, opt);
        syncSigmaFromViscousState(dp, viscous_sigma);

        // --- Compute new disk mass ---
        double current_mass = 0.0;
        for (int i = 1; i <= dp->grid_number; i++) {
            current_mass += dp->gas_surface_density_vector[i] *
                            2.0 * M_PI * dp->radial_grid[i] *
                            (dp->radial_grid[i + 1] - dp->radial_grid[i]);
        }

        // --- Accumulate expected mass loss ---
        integrated_mass_loss += sigma_dot_integral * dt;


        // --- Snapshot logic identical to main solver ---
        int periodic_output_time = (fmod(current_time, interval) < dt);
        int output_time_sync     = ((output_time - current_time) < dt);

        int was_snapshot = periodic_output_time && output_time_sync;

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
                            interval,
                            opt);

        if (periodic_output_time && output_time_sync) {

            writeSigmaProfile(dp, opt, output_time);


            writeSigmaDotProfile(dp, opt, output_time);

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
