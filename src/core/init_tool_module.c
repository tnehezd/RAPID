#include "logger.h"
#include "init_tool_module.h"
#include "config.h" 
#include "disk_model.h" 
#include "dust_physics.h" 
#include "utils.h" 
#include "ascii_output.h" 
#include "gas_physics.h"
#include "boundary_conditions.h"
#include "simulation_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>   

#include "print_panels.h"

void initializeDefaultOptions(InitializeDefaultOptions *def) {
    def->use_cutoff             = false; 
    def->n_grid_points          = 1000; 
    def->n_dust_particles       = 2000;
    def->r_inner                = 0.1;
    def->r_outer                = 5.0;
    def->sigma0_gas_au          = 0.01; 
    def->sigma_exponent         = 0.5;  
    def->alpha_viscosity        = 1.0e-2;
    def->star_mass              = 1.0; 
    def->aspect_ratio           = 5.0e-2;
    def->flaring_index          = 0.0;  
    def->density_floor          = 1e-12;

    def->deadzone_r_inner       = 0.0;
    def->deadzone_r_outer       = 0.0;
    def->deadzone_dr_inner      = 0.0;
    def->deadzone_dr_outer      = 0.0;
    def->deadzone_alpha_mod     = 0.01;

    def->dust_to_gas_ratio      = 0.01;
    def->total_disk_mass        = 0.01; // FIXED: Renamed from disk_mass_dust
    def->one_size_particle_cm   = 1.0;
    def->two_pop_ratio          = 0.85; 
    def->micro_size_cm          = 1e-4; 
    def->drift_factor           = 1.0; 
    def->fragmentation_factor   = 1.0;

    def->output_base_path[0]    = '\0';
    def->dust_density_g_cm3     = 1.6;

    def->r_cutoff               = 30.0; 
    def->n_for_cutoff           = 1.5;  
}


/**
 * @brief Computes the initial surface density constant (Sigma0 at 1 AU) from total mass,
 * or dynamically BACK-CALCULATES the integrated disk mass if Sigma0 was provided.
 */
static long double calculateSigm0FromDiskMass(InitializeDefaultOptions *init_opts) {
    double gamma = fabs(init_opts->sigma_exponent);

    // -------------------------------------------------------------------------
    // CASE 1: Pure Power-Law Profile
    // -------------------------------------------------------------------------
    if (!init_opts->use_cutoff) {
        double exponent_for_integral = -gamma + 2.0;
        double denominator;

        if (fabs(exponent_for_integral) < 1e-9) {
            denominator = log(init_opts->r_outer) - log(init_opts->r_inner);
        } else {
            denominator = (pow(init_opts->r_outer, exponent_for_integral) - pow(init_opts->r_inner, exponent_for_integral)) / exponent_for_integral;
        }

        // IF SIGMA0 IS THE INPUT: Back-calculate the true total gas mass via analytical integral
        if (init_opts->total_disk_mass <= 0.0) { 
            long double true_gas_mass = 2.0 * M_PI * init_opts->sigma0_gas_au * denominator;
            
            // The total disk mass is the GAS mass (the dust is just a 1% component inside it)
            init_opts->total_disk_mass = (double)true_gas_mass; 
            return (long double)init_opts->sigma0_gas_au;
        }

        // IF MASS IS THE INPUT: Compute Sigma0 from the provided total gas mass
        return (long double)init_opts->total_disk_mass / (2.0 * M_PI * denominator);
    } 
    // -------------------------------------------------------------------------
    // CASE 2: Exponential Cutoff Profile (Lynden-Bell & Pringle 1974)
    // -------------------------------------------------------------------------
    else {
        double r_cut = init_opts->r_cutoff;
        double n_cut = init_opts->n_for_cutoff;
        double conversion_exponent = (2.0 - gamma) / n_cut;

        // IF SIGMA0 IS THE INPUT: Back-calculate the cutoff-integrated total gas mass
        if (init_opts->total_disk_mass <= 0.0) { 
            long double true_gas_mass = (2.0 * M_PI * init_opts->sigma0_gas_au * pow(r_cut, 2.0) * tgamma(conversion_exponent)) / n_cut;
            
            init_opts->total_disk_mass = (double)true_gas_mass; 
            return (long double)init_opts->sigma0_gas_au;
        }

        // IF MASS IS THE INPUT: Compute Sigma0 from the provided total gas mass with cutoff normalization
        long double sigma_1AU = (init_opts->total_disk_mass * n_cut) / 
                               (2.0 * M_PI * pow(r_cut, 2.0) * tgamma(conversion_exponent));
        return sigma_1AU;
    }
}

/**
 * @brief Evaluates the local gas surface density at a given radius depending on the IC type.
 */
static long double calculateGasSurfaceDensityInitTool(double r_au, InitializeDefaultOptions *init_opts, long double current_sigma0) {
    double gamma = fabs(init_opts->sigma_exponent);

    if (!init_opts->use_cutoff) {
        return current_sigma0 * pow(r_au, -gamma);
    } else {
        double r_scaled = r_au / init_opts->r_cutoff;
        return current_sigma0 * pow(r_au, -gamma) * exp(-pow(r_scaled, init_opts->n_for_cutoff));
    }
}

static long double calculateDustSurfaceDensityInitTool(double r_au, InitializeDefaultOptions *init_opts, long double current_sigma0) {
    return calculateGasSurfaceDensityInitTool(r_au, init_opts, current_sigma0) * init_opts->dust_to_gas_ratio;
}

static double findMinimumForThreeNumbersInitTool(double value1, double value2, double value3) {
    return fmin(value1, fmin(value2, value3));
}

static int validateInitializationInputs(const InitializeDefaultOptions *opts) {
    if (opts->r_inner <= 0.0) {
        LOG_ERROR("Inner radius must be positive. Value: %lg\n", opts->r_inner);
        return 1;
    }
    if (opts->r_outer <= opts->r_inner) {
        LOG_ERROR("Outer radius must be greater than inner. R_in: %lg, R_out: %lg\n", opts->r_inner, opts->r_outer);
        return 1;
    }
    if (opts->n_grid_points <= 0) {
        LOG_ERROR("Grid points must be positive. Value: %d\n", opts->n_grid_points);
        return 1;
    }
    if (opts->n_dust_particles <= 0) {
        LOG_ERROR("Dust particles must be positive. Value: %d\n", opts->n_dust_particles);
        return 1;
    }
    return 0; 
}

int runInitialization(InitializeDefaultOptions *default_options, DiskParameters *disk_params, SimulationOptions *sim_opts) {
    FILE *dust_ouputput_file = NULL; 
    FILE *disk_parameters_output_file = NULL; 
    FILE *gas_parameters_output_file = NULL; 

    long double current_sigma0_gas;

    if (validateInitializationInputs(default_options) != 0) { 
        return 1; 
    }

    double drdze_inner_calculated = pow(default_options->deadzone_r_inner, 1.0 + default_options->flaring_index) * default_options->aspect_ratio * default_options->deadzone_dr_inner;
    double drdze_outer_calculated = pow(default_options->deadzone_r_outer, 1.0 + default_options->flaring_index) * default_options->aspect_ratio * default_options->deadzone_dr_outer;

    current_sigma0_gas = calculateSigm0FromDiskMass(default_options);
    

    char *full_init_dust_profile_path = NULL;
    char *full_disk_param_path = NULL;
    char *full_init_density_path = NULL;

    asprintf(&full_init_dust_profile_path, "%s/%s%s", default_options->output_base_path, kInitialDustProfileFileName, kFileNamesSuffix);
    asprintf(&full_disk_param_path, "%s/%s%s", default_options->output_base_path, kDiskConfigFile, kFileNamesSuffix);
    asprintf(&full_init_density_path, "%s/%s%s", default_options->output_base_path, kInitialGasProfileFileName, kFileNamesSuffix);

    dust_ouputput_file = fopen(full_init_dust_profile_path, "w");
    disk_parameters_output_file = fopen(full_disk_param_path, "w");
    gas_parameters_output_file = fopen(full_init_density_path, "w");

    if (!dust_ouputput_file || !disk_parameters_output_file || !gas_parameters_output_file) {
        perror("Error opening diagnostic output files");
        if (dust_ouputput_file) fclose(dust_ouputput_file);
        if (disk_parameters_output_file) fclose(disk_parameters_output_file);
        if (gas_parameters_output_file) fclose(gas_parameters_output_file);
        return 1;
    }

    HeaderData initial_header_data;
    initial_header_data.current_time = 0.0;
    initial_header_data.is_initial_data = 1;
    initial_header_data.R_in = default_options->r_inner;
    initial_header_data.R_out = default_options->r_outer;
    initial_header_data.sigma_exponent = -fabs(default_options->sigma_exponent);
    initial_header_data.sigma0_gas_au = current_sigma0_gas;
    initial_header_data.grav_const = G_DIMENSIONLESS;
    initial_header_data.dz_r_inner = default_options->deadzone_r_inner;
    initial_header_data.dz_r_outer = default_options->deadzone_r_outer;
    initial_header_data.dz_dr_inner_calc = drdze_inner_calculated;
    initial_header_data.dz_dr_outer_calc = drdze_outer_calculated;
    initial_header_data.dz_alpha_mod = default_options->deadzone_alpha_mod;
    initial_header_data.dust_density_g_cm3 = default_options->dust_density_g_cm3;
    initial_header_data.alpha_viscosity = default_options->alpha_viscosity;
    initial_header_data.star_mass = default_options->star_mass;
    initial_header_data.flaring_index = default_options->flaring_index;
    initial_header_data.n_grid_points = default_options->n_grid_points;

    printFileHeader(dust_ouputput_file, FILE_TYPE_INTIIAL_DUST_PROFILE, &initial_header_data);
    printFileHeader(gas_parameters_output_file, FILE_TYPE_GAS_DENSITY, &initial_header_data);
    printFileHeader(disk_parameters_output_file, FILE_TYPE_DISK_PARAM, &initial_header_data);

    disk_params->grid_number = default_options->n_grid_points; 
    disk_params->r_min = default_options->r_inner;
    disk_params->r_max = default_options->r_outer;
    disk_params->sigma_0 = current_sigma0_gas;
    disk_params->sigma_power_law_index = -fabs(default_options->sigma_exponent);
    disk_params->r_dze_i = default_options->deadzone_r_inner;
    disk_params->r_dze_o = default_options->deadzone_r_outer;
    disk_params->dr_dze_i = default_options->deadzone_dr_inner;
    disk_params->dr_dze_o = default_options->deadzone_dr_outer;
    disk_params->alpha_parameter = default_options->alpha_viscosity;
    disk_params->alpha_parameter_modification = default_options->deadzone_alpha_mod;
    disk_params->h_aspect_ratio = default_options->aspect_ratio;
    disk_params->flaring_index = default_options->flaring_index;
    disk_params->stellar_mass = default_options->star_mass;
    disk_params->particle_density = default_options->dust_density_g_cm3;
    disk_params->fragmentation_factor = default_options->fragmentation_factor;
    disk_params->drift_factor = default_options->drift_factor;

    if (disk_params->grid_number > 1) {
        disk_params->delta_r = (disk_params->r_max - disk_params->r_min) / ((double)disk_params->grid_number - 1.0);
    } else {
        disk_params->delta_r = 0.0;
    }

    disk_params->radial_grid = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));
    disk_params->gas_surface_density_vector = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));
    disk_params->gas_pressure_vector = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));
    disk_params->gas_pressure_gradient_vector = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));
    disk_params->gas_velocity_vector = (double *)malloc((disk_params->grid_number + 2) * sizeof(double));

    if (!disk_params->radial_grid || !disk_params->gas_surface_density_vector || !disk_params->gas_pressure_vector || 
        !disk_params->gas_pressure_gradient_vector || !disk_params->gas_velocity_vector) {
        return 1;
    }

    readDiskParameters(disk_params);
    createRadialGrid(disk_params);

    for (int i = 1; i <= disk_params->grid_number; i++) {
        disk_params->gas_surface_density_vector[i] = (double)calculateGasSurfaceDensityInitTool(disk_params->radial_grid[i], default_options, current_sigma0_gas);
    }
    sim_opts->current_bc_target = 0;
    applyBoundaryConditions(disk_params->gas_surface_density_vector, disk_params, sim_opts);

    createInitialGasPressure(disk_params,sim_opts);
    createInitialGasPressureGradient(disk_params,sim_opts);
    createInitialGasVelocity(disk_params,sim_opts);

    for (int i_loop = 0; i_loop < default_options->n_grid_points; i_loop++) {
        double r_gas_grid_au = disk_params->radial_grid[i_loop + 1];
        long double sigma_gas_local_val = disk_params->gas_surface_density_vector[i_loop + 1];
        double pressure_local_val = disk_params->gas_pressure_vector[i_loop + 1];
        double dPdr_local_val = disk_params->gas_pressure_gradient_vector[i_loop + 1];

        fprintf(gas_parameters_output_file, "%-5d %-15.6e %-15.6Le %-15.6e %15.6e\n",
            i_loop, r_gas_grid_au, sigma_gas_local_val, pressure_local_val, dPdr_local_val);
    }
    fflush(gas_parameters_output_file); 

    double fragmentation_velocity_cm_s = disk_params->fragmentation_velocity > 0 ? disk_params->fragmentation_velocity : 1000.0;
    double fragmentation_velocity_au_yr2pi = fragmentation_velocity_cm_s * CM_PER_SEC_TO_AU_PER_YEAR_2PI;
    double fragmentation_velocity_sq_au_yr2pi_sq = fragmentation_velocity_au_yr2pi * fragmentation_velocity_au_yr2pi;

    for (int i_loop = 0; i_loop < default_options->n_dust_particles; i_loop++) {
        double r_dust_particle_au;

        if (default_options->n_dust_particles > 1) {
            r_dust_particle_au = default_options->r_inner + (default_options->r_outer - default_options->r_inner) * i_loop / ((double)default_options->n_dust_particles - 1.0);
        } else { 
            r_dust_particle_au = default_options->r_inner; 
        }

        double temp_sigma, temp_pressure, temp_dPdr;
        linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, r_dust_particle_au, &temp_sigma, disk_params->delta_r, disk_params);
        long double sigma_gas_local = temp_sigma;
        
        linearInterpolation(disk_params->gas_pressure_vector, disk_params->radial_grid, r_dust_particle_au, &temp_pressure, disk_params->delta_r, disk_params);
        double pressure_local = temp_pressure;
        
        linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, r_dust_particle_au, &temp_dPdr, disk_params->delta_r, disk_params);
        double dPdr_local = temp_dPdr;
        
        double s_max_cm;

        if (default_options->one_size_particle_cm > 0.0) {
            s_max_cm = default_options->one_size_particle_cm;
        } else {
            double calculateKeplerianVelocity_au_yr2pi = calculateKeplerianVelocity(r_dust_particle_au, disk_params);
            double sound_speed_au_yr2pi = calculateLocalSoundSpeed(r_dust_particle_au, disk_params);
            double sound_speed_sq = sound_speed_au_yr2pi * sound_speed_au_yr2pi;
            
            long double sigma_dust_local = calculateDustSurfaceDensityInitTool(r_dust_particle_au, default_options, current_sigma0_gas);
            long double sigma_dust_local_cgs = sigma_dust_local / SURFACE_DENSITY_CONVERSION_FACTOR;
            double sigma_gas_local_cgs = (double)sigma_gas_local / SURFACE_DENSITY_CONVERSION_FACTOR;
            double dlnPdlnr_local;

            if (fabs(pressure_local) < 1e-12) {
                dlnPdlnr_local = 0.0;
            } else {
                dlnPdlnr_local = r_dust_particle_au / pressure_local * dPdr_local;
            }

            double s_drift = default_options->drift_factor * 2.0 / M_PI * sigma_dust_local_cgs / default_options->dust_density_g_cm3 *
                             (calculateKeplerianVelocity_au_yr2pi * calculateKeplerianVelocity_au_yr2pi) / sound_speed_sq * fabs(1.0 / dlnPdlnr_local);

            double s_frag = default_options->fragmentation_factor * 2.0 / (3.0 * M_PI) * sigma_gas_local_cgs /
                            (default_options->dust_density_g_cm3 * calculateTurbulentAlpha(r_dust_particle_au, disk_params)) *
                            fragmentation_velocity_sq_au_yr2pi_sq / sound_speed_sq;

            double dlnPdlnr_abs_cs2_half = fabs(dlnPdlnr_local * sound_speed_sq * 0.5);
            double s_df;

            if (dlnPdlnr_abs_cs2_half < 1e-12) {
                s_df = 1e99;
            } else {
                s_df = fragmentation_velocity_au_yr2pi * calculateKeplerianVelocity_au_yr2pi / dlnPdlnr_abs_cs2_half * 2.0 * sigma_gas_local_cgs / (M_PI * default_options->dust_density_g_cm3);
            }

            s_max_cm = findMinimumForThreeNumbersInitTool(s_drift, s_frag, s_df);
        }

        if (s_max_cm <= 0) s_max_cm = 1e-10;

        long double representative_mass_total_in_cell = 2.0 * M_PI * r_dust_particle_au *
                                                        ((default_options->r_outer - default_options->r_inner) / ((double)default_options->n_dust_particles - 1.0)) *
                                                        calculateDustSurfaceDensityInitTool(r_dust_particle_au, default_options, current_sigma0_gas);

        long double repr_mass_pop1 = representative_mass_total_in_cell * default_options->two_pop_ratio;
        long double repr_mass_pop2 = representative_mass_total_in_cell * (1.0 - default_options->two_pop_ratio);

        fprintf(dust_ouputput_file, "%-5d %-15.6e %-20.12Lg %-20.12Lg %-15.6e %-15.6e\n",
                i_loop, r_dust_particle_au, repr_mass_pop1, repr_mass_pop2, s_max_cm, default_options->micro_size_cm);
    }

    fflush(dust_ouputput_file);
    fclose(dust_ouputput_file);
    fclose(gas_parameters_output_file);

    fprintf(disk_parameters_output_file, "%-15.6e %-15.6e %-10d %-15.6e %-20.12Lg %-15.6e %-15.6e %-15.6e %-20.12e %-20.12e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
            default_options->r_inner, default_options->r_outer, default_options->n_grid_points, -fabs(default_options->sigma_exponent), current_sigma0_gas,
            G_DIMENSIONLESS, default_options->deadzone_r_inner, default_options->deadzone_r_outer,
            drdze_inner_calculated, drdze_outer_calculated,
            default_options->deadzone_alpha_mod, default_options->dust_density_g_cm3, default_options->alpha_viscosity, default_options->star_mass, default_options->flaring_index);

    fflush(disk_parameters_output_file);
    fclose(disk_parameters_output_file);


    printInitializationParameters(default_options, current_sigma0_gas);

    return 0;
}