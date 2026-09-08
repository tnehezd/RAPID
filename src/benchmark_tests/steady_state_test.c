#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "simulation_types.h"
#include "disk_model.h"
#include "gas_physics.h"
#include "simulation_core.h"
#include "logger.h"
#include "steady_state_test.h"
#include "boundary_conditions.h"

/************************************************************
 * STEADY-STATE INITIALIZATION
 * Sets Sigma(r) and v_r(r) such that Mdot(r) is constant.
 ************************************************************/
void initializeSteadyStateProfile(DiskParameters *dp,
                                  SimulationOptions *opt,
                                  double mdot0)
{
    for (int i = 1; i <= dp->grid_number; i++) {

        double r = dp->radial_grid[i];

        // Keplerian angular velocity
        double Omega = sqrt(dp->stellar_mass / (r*r*r));

        // Aspect ratio (no flaring)
        double H_over_R = dp->h_aspect_ratio;

        // Viscosity ν(r) = α (H/R)^2 r^2 Ω
        double nu = dp->alpha_parameter * H_over_R * H_over_R * r * r * Omega;

        // Steady-state Sigma(r)
        dp->gas_surface_density_vector[i] = mdot0 / (3.0 * M_PI * nu);

        // Steady-state radial velocity v_r(r)
        dp->gas_velocity_vector[i] =
            -mdot0 / (2.0 * M_PI * r * dp->gas_surface_density_vector[i]);
    }
}


/************************************************************
 * STEADY-STATE REFRESH
 * Does NOT evolve Sigma. Only applies BCs and updates pressure.
 ************************************************************/
void refreshSteadyState(SimulationOptions *sim_opts, DiskParameters *dp)
{
    // Update pressure from Sigma (but do NOT evolve Sigma)
    for (int i = 1; i <= dp->grid_number; i++) {
        dp->gas_pressure_vector[i] =
            calculateGasPressure(dp->gas_surface_density_vector[i],
                                 dp->radial_grid[i], dp);
    }

    // Update pressure gradient
    calculateGasPressureGradient(dp);

    // Apply BCs (Sigma, P, dP/dr)
    sim_opts->current_bc_target = 0;
    applyBoundaryConditions(dp->gas_surface_density_vector, dp, sim_opts);

    sim_opts->current_bc_target = 1;
    applyBoundaryConditions(dp->gas_pressure_vector, dp, sim_opts);

    sim_opts->current_bc_target = 2;
    applyBoundaryConditions(dp->gas_pressure_gradient_vector, dp, sim_opts);
}



void runSteadyStateTest(DiskParameters *dp, SimulationOptions *opt) {

    MSG("=== [BENCHMARK] Starting Steady-State Accretion Test ===");

    char *filepath = NULL;
    asprintf(&filepath, "%s/%s/steady_state_summary.dat",
             opt->output_dir_name, kLogFilesDirectory);

    FILE *fp = fopen(filepath, "w");
    if (fp) {
        fprintf(fp, "# Physical_Time_yrs\tMean_Mdot\tMdot_Inner\tMdot_Outer\n");
    }

    double dt_years = 1.0;
    opt->user_defined_time_step = dt_years * 2.0 * M_PI;

    double target_time = 100000.0;
    double current_time = 0.0;
    double interval = opt->output_frequency; // Output frequency in years

    MSG("[BENCHMARK] Initializing steady-state profile...");

    double mdot0 = 1e-7;
    initializeSteadyStateProfile(dp, opt, mdot0);

    MSG("[BENCHMARK] Running steady-state validation...");

    while (current_time < target_time) {

        // DO NOT EVOLVE SIGMA
        refreshSteadyState(opt, dp);

        current_time += dt_years;

        double mdot_inner = 0.0;
        double mdot_outer = 0.0;
        double mdot_sum = 0.0;

        for (int i = 1; i <= dp->grid_number; i++) {
            double r = dp->radial_grid[i];
            double sigma = dp->gas_surface_density_vector[i];
            double vr = dp->gas_velocity_vector[i];

            double mdot_local = -2.0 * M_PI * r * sigma * vr;

            mdot_sum += mdot_local;

            if (i == 2) mdot_inner = mdot_local;
            if (i == dp->grid_number - 1) mdot_outer = mdot_local;
        }

        double mean_mdot = mdot_sum / dp->grid_number;

        if (fmod(current_time, interval) < dt_years || current_time == dt_years) {
            MSG("[STEADY] t=%.0f yrs | mean=%.4e | inner=%.4e | outer=%.4e",
                current_time, mean_mdot, mdot_inner, mdot_outer);
        }

        if (fp) {
            fprintf(fp, "%.1f\t%.10e\t%.10e\t%.10e\n",
                    current_time, mean_mdot, mdot_inner, mdot_outer);
        }
    }

    if (fp) fclose(fp);
    if (filepath) free(filepath);

    MSG("=== [BENCHMARK] Steady-State Accretion Test Completed ===");
}
