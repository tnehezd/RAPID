#include "init_tool_module.h"
#include "config.h" 
#include "disk_model.h" 
#include "dust_physics.h" 
#include "utils.h" 
#include "ascii_output.h" 
#include "gas_physics.h"
#include "boundary_conditions.h"
#include "simulation_types.h"
#include "vertical_settling.h"
#include "vertical_physics.h"
#include "vertical_profile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>   


void initializeDefaultOptions(InitializeDefaultOptions *def) {

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

    def->deadzone_r_inner       = 0.0;
    def->deadzone_r_outer       = 0.0;
    def->deadzone_dr_inner      = 0.0;
    def->deadzone_dr_outer      = 0.0;
    def->deadzone_alpha_mod     = 0.01;

    def->dust_to_gas_ratio      = 0.01;
    def->disk_mass_dust         = 0.0100000001; 
    def->one_size_particle_cm   = 1.0;
    def->two_pop_ratio          = 0.85; 
    def->micro_size_cm          = 1e-4; 
    def->drift_factor           = 1.0; 
    def->fragmentation_factor   = 1.0;

    def->output_base_path[0]    = '\0';
    def->dust_density_g_cm3     = 1.6;

    def->vertical_grid_number   = 100;     
    def->vertical_grid_max_height      = 5.0;    
    def->vertical_grid          = NULL;   
    def->dust_scaleheight       = NULL;   

}

static long double calculateSigm0FromDiskMass(InitializeDefaultOptions *init_opts) {

    double exponent_for_integral = -init_opts->sigma_exponent + 2.0;
    double denominator;

    if (fabs(exponent_for_integral) < 1e-9) {
        denominator = log(init_opts->r_outer) - log(init_opts->r_inner);
    } else {
        denominator = (pow(init_opts->r_outer, exponent_for_integral) - pow(init_opts->r_inner, exponent_for_integral)) / exponent_for_integral;
    }

    if (fabs(denominator) < 1e-12) {
        fprintf(stderr, "Error: Denominator is zero or too small in Sigma0 calculation! Check r_max, r_min, and SIGMA_EXP values.\n");
        return 0.0;
    }

    return (long double)init_opts->disk_mass_dust /
           (2.0 * M_PI * init_opts->dust_to_gas_ratio * denominator);
}


static long double calculateGasSurfaceDensityInitTool(double r_au, InitializeDefaultOptions *init_opts, long double current_sigma0) {

    return current_sigma0 * pow(r_au, init_opts->sigma_exponent);
}

static long double calculateDustSurfaceDensityInitTool(double r_au, InitializeDefaultOptions *init_opts, long double current_sigma0) {

    long double sigma_dust = calculateGasSurfaceDensityInitTool(r_au, init_opts, current_sigma0) * init_opts->dust_to_gas_ratio;

    return sigma_dust;
}

static double findMinimumForThreeNumbersInitTool(double value1, double value2, double value3) {

        return fmin(value1, fmin(value2, value3));
}


static int validateInitializationInputs(const InitializeDefaultOptions *opts) {

    if (opts->r_inner <= 0.0) {
        fprintf(stderr,
            "ERROR [runInitialization]: Inner radius (r_inner) must be positive. Current value: %lg\n",
            opts->r_inner);
        return 1;
    }

    if (opts->r_outer <= opts->r_inner) {
        fprintf(stderr,
            "ERROR [runInitialization]: Outer radius (r_outer) must be greater than inner radius (r_inner). "
            "R_inner: %lg, R_outer: %lg\n",
            opts->r_inner, opts->r_outer);
        return 1;
    }

    if (opts->n_grid_points <= 0) {
        fprintf(stderr,
            "ERROR [runInitialization]: Number of gas grid points (n_grid_points) must be positive. "
            "Current value: %d\n",
            opts->n_grid_points);
        return 1;
    }

    if (opts->n_dust_particles <= 0) {
        fprintf(stderr,
            "ERROR [runInitialization]: Number of dust particles (n_dust_particles) must be positive. "
            "Current value: %d\n",
            opts->n_dust_particles);
        return 1;
    }

    if (opts->deadzone_r_inner > 0.0 && opts->deadzone_dr_inner <= 0.0) {
        fprintf(stderr,
            "ERROR [runInitialization]: deadzone_r_inner is set (%lg), "
            "but deadzone_dr_inner is zero or missing. Please provide a non-zero width.\n",
            opts->deadzone_r_inner);
        return 1;
    }

    if (opts->deadzone_r_outer > 0.0 && opts->deadzone_dr_outer <= 0.0) {
        fprintf(stderr,
            "ERROR [runInitialization]: deadzone_r_outer is set (%lg), "
            "but deadzone_dr_outer is zero or missing. Please provide a non-zero width.\n",
            opts->deadzone_r_outer);
        return 1;
    }


    return 0; 
}

int runInitialization(InitializeDefaultOptions *default_options, DiskParameters *disk_params) {

    FILE *dust_ouputput_file = NULL; 
    FILE *disk_parameters_output_file = NULL; 
    FILE *gas_parameters_output_file = NULL; 

    long double current_sigma0_gas;

    if (validateInitializationInputs(default_options) != 0) { 
        return 1; 
    }

    double drdze_inner_calculated =
        pow(default_options->deadzone_r_inner, 1.0 + default_options->flaring_index)
        * default_options->aspect_ratio
        * default_options->deadzone_dr_inner;

    double drdze_outer_calculated =
        pow(default_options->deadzone_r_outer, 1.0 + default_options->flaring_index)
        * default_options->aspect_ratio
        * default_options->deadzone_dr_outer;

    const double DEFAULT_disk_mass_DUST = 0.01;

    if (fabs(default_options->disk_mass_dust - DEFAULT_disk_mass_DUST) > 1e-9) {
        current_sigma0_gas = calculateSigm0FromDiskMass(default_options);
        fprintf(stderr,"Sigma0 calculated from total dust disk mass (Md): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    } else {
        current_sigma0_gas = default_options->sigma0_gas_au;
        fprintf(stderr,"Using explicit Sigma0 (gas surface density at 1 AU): %Lg M_Sun/AU^2\n", current_sigma0_gas);
    }

    const double DEFAULT_ONE_SIZE = 1.0;

    if (fabs(default_options->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9 && default_options->one_size_particle_cm > 0) {
        default_options->two_pop_ratio = 1.0;
    }

    char *full_init_dust_profile_path = NULL;
    char *full_disk_param_path = NULL;
    char *full_init_density_path = NULL;

    asprintf(&full_init_dust_profile_path, "%s/%s%s",default_options->output_base_path,kInitialDustProfileFileName,kFileNamesSuffix);
    asprintf(&full_disk_param_path, "%s/%s%s",default_options->output_base_path,kDiskConfigFile,kFileNamesSuffix);
    asprintf(&full_init_density_path, "%s/%s%s",default_options->output_base_path,kInitialGasProfileFileName,kFileNamesSuffix);

    dust_ouputput_file = fopen(full_init_dust_profile_path, "w");
    if (dust_ouputput_file == NULL) {
        perror("Error opening initial dust profile file in init_tool_module");
        return 1;
    }

    disk_parameters_output_file = fopen(full_disk_param_path, "w");
    if (disk_parameters_output_file == NULL) {
        perror("Error opening disk parameters file in init_tool_module");
        fclose(dust_ouputput_file);
        return 1;
    }

    gas_parameters_output_file = fopen(full_init_density_path, "w");
    if (gas_parameters_output_file == NULL) {
        perror("Error opening initial density profile file in init_tool_module");
        fclose(dust_ouputput_file);
        fclose(disk_parameters_output_file);
        return 1;
    }

    fprintf(stderr,"\n--- Simulation Parameters ---\n");
    fprintf(stderr,"Total dust disk mass (Solar Mass): %lg\n", default_options->disk_mass_dust);
    fprintf(stderr,"Inner disk edge (AU): %lg\n", default_options->r_inner);
    fprintf(stderr,"Outer disk edge (AU): %lg\n", default_options->r_outer);
    fprintf(stderr,"Surface density profile exponent: %lg\n", -default_options->sigma_exponent);
    fprintf(stderr,"Gas surface density at 1 AU (Solar Mass/AU^2): %Lg\n", current_sigma0_gas);
    fprintf(stderr,"Dust to gas ratio: %lg\n", default_options->dust_to_gas_ratio);
    fprintf(stderr,"Number of gas grid points: %d\n", default_options->n_grid_points); 
    fprintf(stderr,"Number of dust particles to generate: %d\n", default_options->n_dust_particles);
    fprintf(stderr,"Dust particle density (g/cm^3): %lg\n", default_options->dust_density_g_cm3);
    fprintf(stderr,"------------------------------\n\n");

    HeaderData initial_header_data;
    initial_header_data.current_time = 0.0;
    initial_header_data.is_initial_data = 1;
    initial_header_data.R_in = default_options->r_inner;
    initial_header_data.R_out = default_options->r_outer;
    initial_header_data.sigma_exponent = default_options->sigma_exponent;
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

    if (dust_ouputput_file != NULL) {
        printFileHeader(dust_ouputput_file, FILE_TYPE_INTIIAL_DUST_PROFILE, &initial_header_data);
    }
    if (gas_parameters_output_file != NULL) {
        printFileHeader(gas_parameters_output_file, FILE_TYPE_GAS_DENSITY, &initial_header_data);
    }
    if (disk_parameters_output_file != NULL) {
        printFileHeader(disk_parameters_output_file, FILE_TYPE_DISK_PARAM, &initial_header_data);
    }

    double fragmentation_velocity_cm_s = 1000.0;
    double fragmentation_velocity_au_yr2pi = fragmentation_velocity_cm_s * CM_PER_SEC_TO_AU_PER_YEAR_2PI;
    double fragmentation_velocity_sq_au_yr2pi_sq = fragmentation_velocity_au_yr2pi * fragmentation_velocity_au_yr2pi;

    disk_params->grid_number = default_options->n_grid_points; // Gas grid resolution
    disk_params->r_min = default_options->r_inner;
    disk_params->r_max = default_options->r_outer;
    disk_params->sigma_0 = current_sigma0_gas;
    disk_params->sigma_power_law_index = default_options->sigma_exponent;
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


    if (!disk_params->radial_grid || !disk_params->gas_surface_density_vector || !disk_params->gas_pressure_vector || !disk_params->gas_pressure_gradient_vector || !disk_params->gas_velocity_vector) {
        fprintf(stderr, "ERROR [runInitialization]: Failed to allocate disk arrays. Exiting.\n");
        if (dust_ouputput_file) fclose(dust_ouputput_file);
        if (disk_parameters_output_file) fclose(disk_parameters_output_file);
        if (gas_parameters_output_file) fclose(gas_parameters_output_file);
        if (disk_params->radial_grid) free(disk_params->radial_grid);
        if (disk_params->gas_surface_density_vector) free(disk_params->gas_surface_density_vector);
        if (disk_params->gas_pressure_vector) free(disk_params->gas_pressure_vector);
        if (disk_params->gas_pressure_gradient_vector) free(disk_params->gas_pressure_gradient_vector);
        if (disk_params->gas_velocity_vector) free(disk_params->gas_velocity_vector);
        return 1;
    }

    readDiskParameters(disk_params);
    createRadialGrid(disk_params);
    createInitialGasSurfaceDensity(disk_params);
    createInitialGasPressure(disk_params);
    createInitialGasPressureGradient(disk_params);
    createInitialGasVelocity(disk_params);

    for (int i_loop = 0; i_loop < default_options->n_grid_points; i_loop++) {
        double r_gas_grid_au = disk_params->radial_grid[i_loop + 1];

        if (r_gas_grid_au <= 0) {
            fprintf(stderr, "ERROR: Calculated gas grid radial position is non-positive at index %d (%lg AU). Skipping this point.\n", i_loop, r_gas_grid_au);
            continue;
        }

        long double sigma_gas_local_val = disk_params->gas_surface_density_vector[i_loop + 1];
        double pressure_local_val = disk_params->gas_pressure_vector[i_loop + 1];
        double dPdr_local_val = disk_params->gas_pressure_gradient_vector[i_loop + 1];

        fprintf(gas_parameters_output_file, "%-5d %-15.6e %-15.6Le %-15.6e %15.6e\n",
            i_loop, r_gas_grid_au,
            sigma_gas_local_val,
            pressure_local_val,
            dPdr_local_val);
    }
    fflush(gas_parameters_output_file); 

    /* allocate vertical grid and dust scaleheight */
    if (default_options->vertical_grid_number > 0) {
        default_options->vertical_grid = (double *)malloc(default_options->vertical_grid_number * sizeof(double));
        default_options->dust_scaleheight = (double *)malloc(default_options->n_grid_points * sizeof(double));
        
        if (!default_options->vertical_grid || !default_options->dust_scaleheight) {
            fprintf(stderr, "ERROR [runInitialization]: Failed to allocate vertical grid arrays.\n");
            if(default_options->vertical_grid) free(default_options->vertical_grid);
            if(default_options->dust_scaleheight) free(default_options->dust_scaleheight);
            return 1;
        }

        for(int i=0;i<default_options->n_grid_points;i++){
            double r = default_options->r_inner + i*(default_options->r_outer - default_options->r_inner)/(default_options->n_grid_points-1);
            double particle_radius = default_options->one_size_particle_cm / AU_IN_CM;
            default_options->dust_scaleheight[i] = calculateDustScaleHeight(r, particle_radius, disk_params);
        }

    }



    for (int i_loop = 0; i_loop < default_options->n_dust_particles; i_loop++) {

        double r_dust_particle_au;

        if (default_options->n_dust_particles > 1) {
            r_dust_particle_au = default_options->r_inner + (default_options->r_outer - default_options->r_inner) * i_loop / ((double)default_options->n_dust_particles - 1.0);
        } else { 
            r_dust_particle_au = default_options->r_inner; 
        }

        if (r_dust_particle_au <= 0) {
            fprintf(stderr, "ERROR: Calculated radial position is non-positive at dust particle index %d (%lg AU). Skipping this point.\n", i_loop, r_dust_particle_au);
            continue;
        }

        double temp_sigma, temp_pressure, temp_dPdr;
        
        linearInterpolation(disk_params->gas_surface_density_vector, disk_params->radial_grid, r_dust_particle_au, &temp_sigma, disk_params->delta_r,disk_params);
        
        long double sigma_gas_local = temp_sigma;
        
        linearInterpolation(disk_params->gas_pressure_vector, disk_params->radial_grid, r_dust_particle_au, &temp_pressure, disk_params->delta_r, disk_params);
        
        double pressure_local = temp_pressure;
        
        linearInterpolation(disk_params->gas_pressure_gradient_vector, disk_params->radial_grid, r_dust_particle_au, &temp_dPdr, disk_params->delta_r, disk_params);
        
        double dPdr_local = temp_dPdr;
        double s_max_cm;

        if (fabs(default_options->one_size_particle_cm - DEFAULT_ONE_SIZE) > 1e-9 && default_options->one_size_particle_cm > 0) {
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
                fprintf(stderr, "Error: Pressure is near zero in dlnPdlnr calculation at r = %lg. Check input parameters. Setting dlnPdlnr = 0.\n", r_dust_particle_au);
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
                fprintf(stderr, "Error: Denominator is near zero in s_df calculation at r = %lg. Check dlnPdlnr value. Setting s_df to a large value.\n", r_dust_particle_au);
                s_df = 1e99;
            } else {
                s_df = fragmentation_velocity_au_yr2pi * calculateKeplerianVelocity_au_yr2pi / dlnPdlnr_abs_cs2_half * 2.0 * sigma_gas_local_cgs / (M_PI * default_options->dust_density_g_cm3);
            }

            s_max_cm = findMinimumForThreeNumbersInitTool(s_drift, s_frag, s_df);
        }

        if (s_max_cm <= 0) {
            fprintf(stderr, "Warning: s_max_cm <= 0 at r = %lg. This might indicate problematic physical parameters. Setting to a small positive value.\n", r_dust_particle_au);
            s_max_cm = 1e-10;
        }

        long double representative_mass_total_in_cell = 2.0 * M_PI * r_dust_particle_au *
                                                        ((default_options->r_outer - default_options->r_inner) / ((double)default_options->n_dust_particles - 1.0)) *
                                                        calculateDustSurfaceDensityInitTool(r_dust_particle_au, default_options, current_sigma0_gas);

        long double repr_mass_pop1 = representative_mass_total_in_cell * default_options->two_pop_ratio;
        long double repr_mass_pop2 = representative_mass_total_in_cell * (1.0 - default_options->two_pop_ratio);

        fprintf(dust_ouputput_file, "%-5d %-15.6e %-20.12Lg %-20.12Lg %-15.6e %-15.6e\n",
                i_loop, r_dust_particle_au, 
                repr_mass_pop1,
                repr_mass_pop2,
                s_max_cm, default_options->micro_size_cm);
    }

    fflush(dust_ouputput_file);
    fclose(dust_ouputput_file);
    dust_ouputput_file = NULL;
    fclose(gas_parameters_output_file);
    gas_parameters_output_file = NULL;

    fprintf(stderr,"Particle data file created (%s). Writing disk parameters file!\n\n", full_init_dust_profile_path);

    fprintf(disk_parameters_output_file, "%-15.6e %-15.6e %-10d %-15.6e %-20.12Lg %-15.6e %-15.6e %-15.6e %-20.12e %-20.12e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
            default_options->r_inner, default_options->r_outer, default_options->n_grid_points, -default_options->sigma_exponent, current_sigma0_gas,
            G_DIMENSIONLESS, default_options->deadzone_r_inner, default_options->deadzone_r_outer,
            drdze_inner_calculated, drdze_outer_calculated,
            default_options->deadzone_alpha_mod, default_options->dust_density_g_cm3, default_options->alpha_viscosity, default_options->star_mass, default_options->flaring_index);

    fflush(disk_parameters_output_file);
    fclose(disk_parameters_output_file);
    disk_parameters_output_file = NULL;

    fprintf(stderr,"Disk parameters file created (%s).\n\n", full_disk_param_path);


/********************************************************/
/********************************************************/
/********************************************************/    

  // allocate structured particle array
    StructuredParticleData structured_data;
    structured_data.n_r = default_options->n_dust_particles;
    structured_data.n_z = 100;  // most 10 z-cell

    // allocate 2D pointer array
    structured_data.particles = malloc(structured_data.n_r * sizeof(DustParticle *));
    if (!structured_data.particles) {
        fprintf(stderr, "ERROR: Failed to allocate structured particles array!\n");
        return 1;
    }

    // allocate each row (n_z z-cells per r)
    for (size_t i = 0; i < structured_data.n_r; i++) {
        structured_data.particles[i] = malloc(structured_data.n_z * sizeof(DustParticle));
        if (!structured_data.particles[i]) {
            fprintf(stderr, "ERROR: Failed to allocate particles row %zu!\n", i);
            for (size_t j = 0; j < i; j++) free(structured_data.particles[j]);
            free(structured_data.particles);
            return 1;
        }
    }

    // open output files
    FILE *mass_file = fopen("mass_matrix.dat", "w");
    FILE *z_file    = fopen("z_matrix.dat", "w");
    if (!mass_file || !z_file) {
        perror("ERROR opening output files");
        for (size_t i = 0; i < structured_data.n_r; i++) free(structured_data.particles[i]);
        free(structured_data.particles);
        return 1;
    }

    // fill particles
    for (size_t i_r = 0; i_r < structured_data.n_r; i_r++) {
        double r_dust_particle_au;
        if (structured_data.n_r > 1) {
            r_dust_particle_au = default_options->r_inner + 
                                 (default_options->r_outer - default_options->r_inner) * i_r / ((double)structured_data.n_r - 1.0);
        } else { 
            r_dust_particle_au = default_options->r_inner; 
        }

        long double sigma_dust_local = calculateDustSurfaceDensityInitTool(r_dust_particle_au, default_options, current_sigma0_gas);
        long double representative_mass_total_in_cell = 2.0 * M_PI * r_dust_particle_au *
                                                        ((default_options->r_outer - default_options->r_inner) / ((double)structured_data.n_r - 1.0)) *
                                                        sigma_dust_local;
        long double repr_mass_pop1 = representative_mass_total_in_cell * default_options->two_pop_ratio;

        // compute scale height for this r
        double H = default_options->dust_scaleheight[i_r]; // most már tömb

        // sum of vertical weights
        long double weight_sum = 0.0;
        double z_array[structured_data.n_z];
        for (size_t i_z = 0; i_z < structured_data.n_z; i_z++) {
            double z = -3.0*H + 6.0*H * i_z / (structured_data.n_z - 1); // [-3H, +3H]
            z_array[i_z] = z;
            weight_sum += exp(-0.5*(z/H)*(z/H));
        }

        // fill each z-slice
        for (size_t i_z = 0; i_z < structured_data.n_z; i_z++) {
            double z = z_array[i_z];
            double vertical_factor = exp(-0.5*(z/H)*(z/H)) / weight_sum; // normált tömeg
            structured_data.particles[i_r][i_z].index  = i_r;
            structured_data.particles[i_r][i_z].r_au   = r_dust_particle_au;
            structured_data.particles[i_r][i_z].z_au   = z;
            structured_data.particles[i_r][i_z].mass_g = (double)(repr_mass_pop1 * vertical_factor);
        }

        // write row for mass matrix
        for (size_t i_z = 0; i_z < structured_data.n_z; i_z++) {
            fprintf(mass_file, "%e ", structured_data.particles[i_r][i_z].mass_g);
        }
        fprintf(mass_file, "\n");

        // write row for z matrix: first column = r, then z-slices
        fprintf(z_file, "%e ", r_dust_particle_au);
        for (size_t i_z = 0; i_z < structured_data.n_z; i_z++) {
            fprintf(z_file, "%e ", structured_data.particles[i_r][i_z].z_au);
        }
        fprintf(z_file, "\n");
    }

    fflush(mass_file);  fclose(mass_file);
    fflush(z_file);     fclose(z_file);

    fprintf(stderr, "DEBUG: 2D structured particle array allocated and written to mass_matrix.dat + z_matrix.dat\n");

    return 0;
}