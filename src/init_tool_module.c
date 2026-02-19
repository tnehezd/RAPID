#include "init_tool_module.h"
#include "config.h" 
#include "disk_model.h" 
#include "dust_physics.h" 
#include "utils.h" 
#include "ascii_output.h" 
#include "particle_data.h"
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

    def->n_grid_points              = 1000; 
    def->n_dust_particles           = 2000;
    def->r_inner                    = 0.1;
    def->r_outer                    = 5.0;
    def->sigma0_gas_au              = 0.01; 
    def->sigma_exponent             = 0.5; 
    def->alpha_viscosity            = 1.0e-2;
    def->star_mass                  = 1.0; 
    def->aspect_ratio               = 5.0e-2;
    def->flaring_index              = 0.0;  

    def->deadzone_r_inner           = 0.0;
    def->deadzone_r_outer           = 0.0;
    def->deadzone_dr_inner          = 0.0;
    def->deadzone_dr_outer          = 0.0;
    def->deadzone_alpha_mod         = 0.01;

    def->dust_to_gas_ratio          = 0.01;
    def->disk_mass_dust             = 0.0100000001; 
    def->one_size_particle_cm       = 1.0;
    def->two_pop_ratio              = 0.85; 
    def->micro_size_cm              = 1e-4; 
    def->drift_factor               = 1.0; 
    def->fragmentation_factor       = 1.0;

    def->output_base_path[0]        = '\0';
    def->dust_density_g_cm3         = 1.6;

    def->vertical_grid_number       = 100;     
    def->vertical_grid_max_height   = 4.0;    
    def->vertical_grid              = NULL;   
    def->dust_scaleheight           = NULL;   

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

int initializeOneDimensions(InitializeDefaultOptions *default_options, DiskParameters *disk_params, long double current_sigma0_gas, StructuredParticleData *structured_data, FILE *dust_output_file, FILE *gas_output_file) {
    
    // 1. Struktúra alapadatok beállítása (n_z = 1 az 1D-hez)
    structured_data->n_r = default_options->n_dust_particles;
    structured_data->n_z = 1; 

    // 2. Memória foglalása
    structured_data->particles = (DustParticle**)malloc(structured_data->n_r * sizeof(DustParticle*));
    if (!structured_data->particles) return 1;
    for (size_t i = 0; i < structured_data->n_r; i++) {
        structured_data->particles[i] = (DustParticle*)malloc(structured_data->n_z * sizeof(DustParticle));
        if (!structured_data->particles[i]) return 1;
    }

    double dr_particles = (default_options->r_outer - default_options->r_inner) / (double)structured_data->n_r;

    // 3. Feltöltés és kiírás
    for (size_t i = 0; i < structured_data->n_r; i++) {
        double r = default_options->r_inner + (i + 0.5) * dr_particles;

        // Tömegszámítás (precíz bin-integrál)
        double r_bin_in  = r - 0.5 * dr_particles;
        double r_bin_out = r + 0.5 * dr_particles;
        if (r_bin_in < default_options->r_inner) r_bin_in = default_options->r_inner;
        if (r_bin_out > default_options->r_outer) r_bin_out = default_options->r_outer;

        long double sigma_dust = calculateDustSurfaceDensityInitTool(r, default_options, current_sigma0_gas);
        long double cell_mass = M_PI * (r_bin_out * r_bin_out - r_bin_in * r_bin_in) * sigma_dust;

        // Két populációs tömegbontás az eredeti fájlformátumhoz
        long double repr_mass_pop1 = cell_mass * default_options->two_pop_ratio;
        long double repr_mass_pop2 = cell_mass * (1.0 - default_options->two_pop_ratio);

        // --- Memória feltöltése (AoS) ---
        structured_data->particles[i][0].index   = i;
        structured_data->particles[i][0].r_au    = r;
        structured_data->particles[i][0].z_au    = 0.0;
        structured_data->particles[i][0].mass_g  = (double)cell_mass;
        structured_data->particles[i][0].radius  = default_options->one_size_particle_cm;

        // --- Eredeti fájlba írás (Dust) ---
        // Megtartva a pontos formázást: index, r, mass_pop1, mass_pop2, s_max, s_mic
        fprintf(dust_output_file, "%-5zu %-15.6e %-20.12Lg %-20.12Lg %-15.6e %-15.6e\n",
                i, r, repr_mass_pop1, repr_mass_pop2, 
                default_options->one_size_particle_cm, default_options->micro_size_cm);
    }

    fprintf(stderr, "1D dust particle distribution written to file and memory.\n");

    // --- Eredeti gáz profil kiírása ---
    for (int i = 0; i < disk_params->grid_number; i++) {
        double r_gas = disk_params->radial_grid[i + 1];
        fprintf(gas_output_file, "%-5d %-15.6e %-15.6e %-15.6e %15.6e\n",
                i, r_gas,
                disk_params->gas_surface_density_vector[i + 1],
                disk_params->gas_pressure_vector[i + 1],
                disk_params->gas_pressure_gradient_vector[i + 1]);
    }

    fprintf(stderr, "1D gas profiles written to file.\n");
    fflush(dust_output_file);
    fflush(gas_output_file);

    return 0;
}

int initializeTwoDimensions(InitializeDefaultOptions *default_options, long double current_sigma0_gas, StructuredParticleData *structured_data) {
    // 1. Memória foglalása a 2D struktúrának
    structured_data->n_r = default_options->n_dust_particles;
    structured_data->n_z = default_options->vertical_grid_number;

    structured_data->particles = malloc(structured_data->n_r * sizeof(DustParticle*));
    if (!structured_data->particles) {
        fprintf(stderr, "ERROR: Failed to allocate structured_data->particles\n");
        return 1;
    }
    for (size_t i = 0; i < structured_data->n_r; i++) {
        structured_data->particles[i] = malloc(structured_data->n_z * sizeof(DustParticle));
        if (!structured_data->particles[i]) {
            fprintf(stderr, "ERROR: Failed to allocate structured_data->particles[%zu]\n", i);
            return 1;
        }
    }

    double *z_array = malloc(structured_data->n_z * sizeof(double));
    if (!z_array) {
        fprintf(stderr, "ERROR: Failed to allocate z_array\n");
        return 1;
    }

    double n_sigma = default_options->vertical_grid_max_height; 
    double dr_particles = (default_options->r_outer - default_options->r_inner) 
                          / (double)(structured_data->n_r);

    // 2. Feltöltés radiális irányban
    for (size_t i_r = 0; i_r < structured_data->n_r; i_r++) {
        // Részecske sugara (szabályos kiosztás, jitter nélkül a tiszta profilért)
        double r = default_options->r_inner + (i_r + 0.5) * dr_particles;

        // --- PRECÍZ BIN-INTEGRÁLT TÖMEGSZÁMÍTÁS ---
        // Kiszámoljuk a részecske "saját" gyűrűjének határait (fél úton a szomszédokig)
        double r_bin_in  = r - 0.5 * dr_particles;
        double r_bin_out = r + 0.5 * dr_particles;

        // Határkorrekciók
        if (r_bin_in < default_options->r_inner) r_bin_in = default_options->r_inner;
        if (r_bin_out > default_options->r_outer) r_bin_out = default_options->r_outer;

        // Lokális skálamagasság és felületi sűrűség
        double H = default_options->dust_scaleheight[i_r];
        long double sigma_dust = calculateDustSurfaceDensityInitTool(r, default_options, current_sigma0_gas);

        // A teljes gyűrű tömege: Area * Sigma_dust
        // Ez a pontos geometriai integrál, ami illeszkedik a CIC-hez
        long double cell_mass = M_PI * (r_bin_out * r_bin_out - r_bin_in * r_bin_in) * sigma_dust;

        // 3. Vertikális eloszlás (Gauss) előkészítése
        long double weight_sum = 0.0;
        double z_min = -n_sigma * H;
        double z_max = +n_sigma * H;

        // Súlyok kiszámítása a normalizáláshoz
        for (size_t iz = 0; iz < structured_data->n_z; iz++) {
            double z = (structured_data->n_z > 1) ? z_min + (z_max - z_min) * iz / (structured_data->n_z - 1) : 0.0;
            z_array[iz] = z;
            weight_sum += exp(-0.5 * (z / H) * (z / H));
        }

        // 4. Részecskék adatainak mentése
        for (size_t iz = 0; iz < structured_data->n_z; iz++) {
            // Gauss-súly (tömegmegmaradással)
            double f = exp(-0.5 * (z_array[iz] / H) * (z_array[iz] / H)) / weight_sum;
            double m = (double)(cell_mass * f);

            structured_data->particles[i_r][iz].mass_g = m;
            structured_data->particles[i_r][iz].radius = default_options->one_size_particle_cm / AU_IN_CM;
            structured_data->particles[i_r][iz].r_au   = r;
            structured_data->particles[i_r][iz].z_au   = z_array[iz];
        }
    }

    free(z_array);
    return 0;
}   

int runInitialization(InitializeDefaultOptions *default_options, DiskParameters *disk_params, StructuredParticleData *structured_data) {

    FILE *dust_output_file = NULL; 
    FILE *disk_parameters_output_file = NULL; 
    FILE *gas_parameters_output_file = NULL; 

    long double current_sigma0_gas;

    if (validateInitializationInputs(default_options) != 0) { 
        return 1; 
    }

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

    asprintf(&full_init_dust_profile_path, "%s/%s%s",
             default_options->output_base_path, kInitialDustProfileFileName, kFileNamesSuffix);
    asprintf(&full_disk_param_path, "%s/%s%s",
             default_options->output_base_path, kDiskConfigFile, kFileNamesSuffix);
    asprintf(&full_init_density_path, "%s/%s%s",
             default_options->output_base_path, kInitialGasProfileFileName, kFileNamesSuffix);

    dust_output_file = fopen(full_init_dust_profile_path, "w");
    disk_parameters_output_file = fopen(full_disk_param_path, "w");
    gas_parameters_output_file = fopen(full_init_density_path, "w");

    if (!dust_output_file || !disk_parameters_output_file || !gas_parameters_output_file) {
        perror("Error opening initialization files");
        if(dust_output_file) fclose(dust_output_file);
        if(disk_parameters_output_file) fclose(disk_parameters_output_file);
        if(gas_parameters_output_file) fclose(gas_parameters_output_file);
        return 1;
    }

    HeaderData initial_header_data = {0};
    initial_header_data.current_time = 0.0;
    initial_header_data.is_initial_data = 1;
    initial_header_data.R_in = default_options->r_inner;
    initial_header_data.R_out = default_options->r_outer;
    initial_header_data.sigma_exponent = default_options->sigma_exponent;
    initial_header_data.sigma0_gas_au = current_sigma0_gas;
    initial_header_data.grav_const = G_DIMENSIONLESS;
    initial_header_data.dz_r_inner = default_options->deadzone_r_inner;
    initial_header_data.dz_r_outer = default_options->deadzone_r_outer;
    initial_header_data.dz_alpha_mod = default_options->deadzone_alpha_mod;
    initial_header_data.dust_density_g_cm3 = default_options->dust_density_g_cm3;
    initial_header_data.alpha_viscosity = default_options->alpha_viscosity;
    initial_header_data.star_mass = default_options->star_mass;
    initial_header_data.flaring_index = default_options->flaring_index;
    initial_header_data.n_grid_points = default_options->n_grid_points;

    printFileHeader(dust_output_file, FILE_TYPE_INTIIAL_DUST_PROFILE, &initial_header_data);
    printFileHeader(gas_parameters_output_file, FILE_TYPE_GAS_DENSITY, &initial_header_data);
    printFileHeader(disk_parameters_output_file, FILE_TYPE_DISK_PARAM, &initial_header_data);

    disk_params->grid_number = default_options->n_grid_points;
    disk_params->r_min = default_options->r_inner;
    disk_params->r_max = default_options->r_outer;
    disk_params->sigma_0 = current_sigma0_gas;
    disk_params->sigma_power_law_index = default_options->sigma_exponent;
    disk_params->r_dze_i = default_options->deadzone_r_inner;
    disk_params->r_dze_o = default_options->deadzone_r_outer;
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

    disk_params->dust_surface_density_euler = (double *)malloc(disk_params->grid_number * sizeof(double));

    if (!disk_params->dust_surface_density_euler) {
        fprintf(stderr, "ERROR: Failed to allocate Eulerian dust surface density array.\n");
        return 1;
    }


    if (!disk_params->radial_grid || !disk_params->gas_surface_density_vector ||
        !disk_params->gas_pressure_vector || !disk_params->gas_pressure_gradient_vector ||
        !disk_params->gas_velocity_vector) 
    {
        fprintf(stderr, "ERROR: Failed to allocate disk arrays.\n");
        return 1;
    }

    readDiskParameters(disk_params);
    createRadialGrid(disk_params);
    createInitialGasSurfaceDensity(disk_params);
    createInitialGasPressure(disk_params);
    createInitialGasPressureGradient(disk_params);
    createInitialGasVelocity(disk_params);

    if (default_options->vertical_grid_number > 0) {
        default_options->vertical_grid = (double *)malloc(default_options->vertical_grid_number * sizeof(double));
        default_options->dust_scaleheight = (double *)malloc(default_options->n_grid_points * sizeof(double));
        for(int i=0;i<default_options->n_grid_points;i++){
            double r = default_options->r_inner + i*(default_options->r_outer - default_options->r_inner)/(default_options->n_grid_points-1);
            double particle_radius = default_options->one_size_particle_cm / AU_IN_CM;
            default_options->dust_scaleheight[i] = calculateDustScaleHeight(r, particle_radius, disk_params);
        }
    }


    fprintf(disk_parameters_output_file, "%-15.6e %-15.6e %-10d %-15.6e %-20.12Lg %-15.6e %-15.6e %-15.6e %-20.12e %-20.12e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n",
            default_options->r_inner, default_options->r_outer, default_options->n_grid_points, 
            -default_options->sigma_exponent, current_sigma0_gas,
            G_DIMENSIONLESS, default_options->deadzone_r_inner, default_options->deadzone_r_outer,
            default_options->deadzone_dr_inner, default_options->deadzone_dr_outer,
            default_options->deadzone_alpha_mod, default_options->dust_density_g_cm3, default_options->alpha_viscosity, 
            default_options->star_mass, default_options->flaring_index);



    if (default_options->dimension == 1) {
        initializeOneDimensions(default_options, disk_params, current_sigma0_gas, structured_data, dust_output_file, gas_parameters_output_file);
    } else if (default_options->dimension == 2) {
        initializeTwoDimensions(default_options, current_sigma0_gas, structured_data);
        writeDustField2D(structured_data, default_options->output_base_path, -1, "initial");

        int N = disk_params->grid_number;

        // 1) nullázás
        for (int i = 0; i < N; i++) {
            disk_params->dust_surface_density_euler[i] = 0.0;
        }

        // 2) részecskék tömegének összegzése radiális cellákba
        // --- ANALYTIKUS KEZDETI PORFELÜLETISŰRŰSÉG (Euler) ---
        for (int i = 0; i < disk_params->grid_number; i++) {
            double r = disk_params->radial_grid[i];

            // analitikus porfelületisűrűség
            disk_params->dust_surface_density_euler[i] = (double) calculateDustSurfaceDensityInitTool(r, default_options, current_sigma0_gas);
        }
        fprintf(stderr, "Initial Euler dust surface density set analytically.\n");
    }


    fclose(dust_output_file);
    fclose(gas_parameters_output_file);
    fclose(disk_parameters_output_file);

    return 0;
}
