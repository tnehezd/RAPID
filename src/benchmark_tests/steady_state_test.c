#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
    // Prevent unused parameter warning
    (void)opt;

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
 * REAL STEADY-STATE EVOLUTION TEST
 * Runs simulation steps with actual viscous transport.
 ************************************************************/
void runSteadyStateTest(DiskParameters *dp, SimulationOptions *opt) {

    MSG("=== [BENCHMARK] Starting Real Steady-State Accretion Test ===");

    char *filepath = NULL;
    asprintf(&filepath, "%s/%s/steady_state_summary.dat",
             opt->output_dir_name, kLogFilesDirectory);

    FILE *fp = fopen(filepath, "w");
    if (fp) {
        fprintf(fp, "# Physical_Time_yrs\tMean_Mdot\tMdot_Inner\tMdot_Outer\n");
    }

    // Timing parameters in internal units
    double t = 0.0;
    double target_time_years = opt->maximum_simulation_time;

    LOG_INFO("Maximum simulation time: %.1f years", target_time_years);
    double target_time_internal = target_time_years * 2.0 * M_PI;
    double interval_years = opt->output_frequency; 

    double last_logged_time_years = 0.0;
    double dt_old = 0.0;

    MSG("[BENCHMARK] Initializing steady-state profile...");
    double mdot0 = 1e-7; 
    initializeSteadyStateProfile(dp, opt, mdot0);

    // Explicitly update velocity fields before entering loop
    calculateGasRadialVelocity(dp);

    // =========================================================================
    // --- TARTOZÉK: GARANTÁLT T=0 INITIAL KIÍRÁS A HUROK ELŐTT ---
    // =========================================================================
    {
        double mdot_inner = 0.0;
        double mdot_outer = 0.0;
        double mdot_sum = 0.0;

        // Az i=1 és i=grid_number széleket kihagyjuk a szomszédok miatt
        for (int i = 2; i < dp->grid_number; i++) {
            double r_j   = dp->radial_grid[i];
            double r_p   = dp->radial_grid[i+1];
            double r_m   = dp->radial_grid[i-1];

            double nu_j  = calculateKinematicViscosity(r_j, dp);
            double nu_p  = calculateKinematicViscosity(r_p, dp);
            double nu_m  = calculateKinematicViscosity(r_m, dp);

            double sig_p = dp->gas_surface_density_vector[i+1];
            double sig_m = dp->gas_surface_density_vector[i-1];

            // Centrális differencia a d(Sigma * nu * sqrt(r)) / dr tagra
            double flux_p = sig_p * nu_p * sqrt(r_p);
            double flux_m = sig_m * nu_m * sqrt(r_m);
            double dflux_dr = (flux_p - flux_m) / (2.0 * dp->delta_r);

            double mdot_local = 6.0 * M_PI * sqrt(r_j) * dflux_dr;
            mdot_sum += mdot_local;

            if (i == 2) mdot_inner = mdot_local;
            if (i == dp->grid_number - 1) mdot_outer = mdot_local;
        }
        double mean_mdot = mdot_sum / (dp->grid_number - 2);

        MSG("[STEADY-REAL] t=0.0 yrs | mean_Mdot=%.4e | inner=%.4e | outer=%.4e",
            mean_mdot, mdot_inner, mdot_outer);
        
        if (fp) {
            fprintf(fp, "0.0\t%.10e\t%.10e\t%.10e\n", mean_mdot, mdot_inner, mdot_outer);
            fflush(fp); // Kikényszerítjük a lemezre írást, hogy ne maradjon pufferben
        }

        char *prof_path = NULL;
        asprintf(&prof_path, "%s/%s/steady_profile_t_0.dat", opt->output_dir_name, kLogFilesDirectory);
        FILE *fprof = fopen(prof_path, "w");
        if (fprof) {
            fprintf(fprof, "# r\tSigma\tv_r\tLocal_Mdot\n");
            for (int j = 1; j <= dp->grid_number; j++) {
                double r_j = dp->radial_grid[j];
                double sig_j = dp->gas_surface_density_vector[j];
                double vr_j = dp->gas_velocity_vector[j];
                double mdot_j = -2.0 * M_PI * r_j * sig_j * vr_j;
                fprintf(fprof, "%.4f\t%.6e\t%.6e\t%.6e\n", r_j, sig_j, vr_j, mdot_j);
            }
            fclose(fprof);
        }
        if (prof_path) free(prof_path);
    }
    // =========================================================================

    MSG("[BENCHMARK] Running viscous evolution loop for %.0f years...", target_time_years);

    do {
        // Dynamic time-step allocation via internal CFL logic
        double dt_new = calculateTimeStep(dp, 0.0) / 5.0; 
        if (dt_old == 0.0) dt_old = dt_new;
        
        double deltat = 0.7 * dt_old + 0.3 * dt_new;
        dt_old = deltat;
        opt->user_defined_time_step = deltat;

        // REAL VISCOUS EVOLUTION STEP 
        refreshGasSurfaceDensityPressurePressureGradient(opt, dp);
        
        // Advance internal time
        t += deltat;

        double current_time_years = t / (2.0 * M_PI);

        // Verification of local accretion rate
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

        // Diagnostics and file printing at scheduled frequencies
        if (current_time_years - last_logged_time_years >= interval_years || t >= target_time_internal) {
            MSG("[STEADY-REAL] t=%.1f yrs | mean_Mdot=%.4e | inner=%.4e | outer=%.4e",
                current_time_years, mean_mdot, mdot_inner, mdot_outer);
            
            if (fp) {
                fprintf(fp, "%.1f\t%.10e\t%.10e\t%.10e\n",
                        current_time_years, mean_mdot, mdot_inner, mdot_outer);
                fflush(fp); // Biztonsági mentés lemezre minden logoláskor
            }

            // --- EXPORT RADIAL PROFILES FOR PYTHON PLOTTER ---
            char *prof_path = NULL;
            asprintf(&prof_path, "%s/%s/steady_profile_t_%.0f.dat",
                     opt->output_dir_name, kLogFilesDirectory, current_time_years);
            
            FILE *fprof = fopen(prof_path, "w");
            if (fprof) {
                fprintf(fprof, "# r\tSigma\tv_r\tLocal_Mdot\n");
                for (int j = 1; j <= dp->grid_number; j++) {
                    double r_j = dp->radial_grid[j];
                    double sig_j = dp->gas_surface_density_vector[j];
                    double vr_j = dp->gas_velocity_vector[j];
                    double mdot_j = -2.0 * M_PI * r_j * sig_j * vr_j;

                    fprintf(fprof, "%.4f\t%.6e\t%.6e\t%.6e\n", r_j, sig_j, vr_j, mdot_j);
                }
                fclose(fprof);
            }
            if (prof_path) free(prof_path);

            last_logged_time_years = current_time_years;
        }

    } while (t < target_time_internal);


    if (fp) fclose(fp);
    if (filepath) free(filepath);

    MSG("=== [BENCHMARK] Steady-State Accretion Test Completed ===");
}
