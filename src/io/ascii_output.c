#include "utils.h"
#include "print_panels.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>  
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h> 
#include <pwd.h>    
#include "logger.h"


#ifdef _OPENMP
    #include <omp.h>
#endif

#ifdef _WIN32
    #include <direct.h>
    #define MKDIR_CALL(path) _mkdir(path)
    #define access _access
    #define F_OK 0
    #include <io.h>
#else
    #include <unistd.h>
    #define MKDIR_CALL(path) mkdir(path, 0755)
#endif

#include "ascii_output.h"
#include "config.h"         
#include "dust_physics.h"  
#include "utils.h"         
#include "simulation_types.h" 
#include "boundary_conditions.h"

#define INIT_DATA_HEADER_LINES 8
#define HEADER_WIDTH 75 

int calculateNumbersOfParticles(const char *particle_data_file_name) {
    FILE *particle_file = NULL;
    char line_buffer[1024];
    int line_count = 0; 

    particle_file = fopen(particle_data_file_name, "r");
    if (particle_file == NULL) {
        LOG_ERROR("Could not open file '%s'.\n", particle_data_file_name);
        perror("Reason"); 
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < INIT_DATA_HEADER_LINES; i++) {
        if (fgets(line_buffer, sizeof(line_buffer), particle_file) == NULL) {
            LOG_ERROR("Unexpected end of file while skipping %d header lines in '%s'.\n", INIT_DATA_HEADER_LINES, particle_data_file_name);
            fclose(particle_file);
            exit(EXIT_FAILURE);
        }
    }

    while (fgets(line_buffer, sizeof(line_buffer), particle_file) != NULL) {
        if (line_buffer[0] != '#' && line_buffer[0] != '\n' && line_buffer[0] != '\r') {
             line_count++;
        }
    }

    fclose(particle_file); 
    return line_count;
}

void loadDustParticlesFromFile(ParticleData *particle_data, const char *particle_data_file_name) {
    int i, particle_index;
    double distance, particle_radius, micron_particle_radius;
    long double representative_mass;
    long double micron_representative_mass;

    load_dust_particles_file = fopen(particle_data_file_name, "r");

    if (load_dust_particles_file == NULL) {
        LOG_ERROR("Could not open file '%s'.\n", particle_data_file_name);
        perror("Reason");
        exit(EXIT_FAILURE);
    }

    char line_buffer[1024];
    for (int k = 0; k < INIT_DATA_HEADER_LINES; k++) {
        if (fgets(line_buffer, sizeof(line_buffer), load_dust_particles_file) == NULL) {
            LOG_ERROR("Unexpected end of file while skipping headers in '%s'.\n", particle_data_file_name);
            fclose(load_dust_particles_file);
            exit(EXIT_FAILURE);
        }
    }

    for (i = 0; i < particle_number; i++) {
        if (fscanf(load_dust_particles_file, "%d %lg %Lg %Lg %lg %lg", 
                   &particle_index, &distance, &representative_mass, &micron_representative_mass, 
                   &particle_radius, &micron_particle_radius) == 6) {
            
            // Primary population (always valid)
            particle_data->particle_distance_array[i][0] = distance;
            particle_data->particle_distance_array[i][1] = particle_radius / AU_IN_CM; 
            particle_data->dust_particle_mass_grid[i] = representative_mass;

            // SAFE GUARD: Only write to secondary population if arrays are allocated
            if (particle_data->micron_particle_distance_array != NULL) {
                particle_data->micron_particle_distance_array[i][0] = distance;
                particle_data->micron_particle_distance_array[i][1] = micron_particle_radius / AU_IN_CM; 
            }
            if (particle_data->massmicradial_grid != NULL) {
                particle_data->massmicradial_grid[i] = micron_representative_mass;
            }
            
        } else {
            LOG_ERROR("Failed to read line %d from particle data file '%s'!\n", i, particle_data_file_name);
            fclose(load_dust_particles_file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(load_dust_particles_file);
}


void loadGasSurfaceDensityFromFile(DiskParameters *disk_params, const char *disk_file_name) {
    const char *input_disk_file_name = disk_file_name;

    FILE *input_file = fopen(input_disk_file_name, "r");
    if (input_file == NULL) {
        LOG_ERROR("Could not open input file '%s'.\n", input_disk_file_name);
        perror("Reason");
        exit(EXIT_FAILURE);
    }

    char line[512];

    while (fgets(line, sizeof(line), input_file) != NULL) {
        if (line[0] == '#' || strncmp(line, "---", 3) == 0) {
            continue;
        } else {
            fseek(input_file, -strlen(line), SEEK_CUR); 
            break;
        }
    }

    if (feof(input_file) && (line[0] == '#' || strncmp(line, "---", 3) == 0)) {
        LOG_ERROR("File '%s' is empty or only contains comments/headers.\n", input_disk_file_name);
        fclose(input_file);
        exit(EXIT_FAILURE);
    }

    int index;
    double r_value;
    double surfacedensity_gas_value;
    double gas_pressure_value;
    double gas_pressure_gradient_value;

    for (int i = 0; i < disk_params->grid_number; i++) {
        if (fscanf(input_file, "%d %lf %lf %lf %lf",
                         &index, &r_value, &surfacedensity_gas_value, &gas_pressure_value, &gas_pressure_gradient_value) != 5) {
            LOG_ERROR("Failed to read 4 values for row %d from file '%s'. File may be malformed or ended unexpectedly.\n", i, input_disk_file_name);
            fclose(input_file);
            exit(EXIT_FAILURE);
        }

        if ((i + 1) >= 0 && (i + 1) <= disk_params->grid_number + 1) { 
            disk_params->radial_grid[i + 1] = r_value;
            disk_params->gas_surface_density_vector[i + 1] = surfacedensity_gas_value;
            disk_params->gas_pressure_vector[i + 1] = gas_pressure_value;
            disk_params->gas_pressure_gradient_vector[i + 1] = gas_pressure_gradient_value;
        } else {
            LOG_WARN("Attempted to write to out-of-bounds index %d. Max allowed index: %d (grid_number+1).\n", i + 1, disk_params->grid_number + 1);
        }
    }

    fclose(input_file);
}


char *createRunDirectory(const char *dir_path) {

    char *temporary_path = NULL;
    int counter = 0;

    asprintf(&temporary_path, "%s", dir_path);

    while (access(temporary_path, F_OK) == 0) {
        free(temporary_path);
        asprintf(&temporary_path, "%s_%04d", dir_path, ++counter);

        if (counter > 99) {
            LOG_ERROR("Too many directories.\n");
            exit(1);
        }
    }

    if (MKDIR_CALL(temporary_path) != 0) {
        LOG_ERROR("mkdir failed");
        exit(1);
    }

    return temporary_path;
}

void printCentered(FILE *file, const char *text) {

    int length = strlen(text);
    int padding = (HEADER_WIDTH - length) / 2;
    
    fprintf(file, "#");
    for (int i = 0; i < padding; i++) fprintf(file, " ");
    fprintf(file, "%s", text);
    for (int i = 0; i < (HEADER_WIDTH - length - padding); i++) fprintf(file, " ");
    fprintf(file, "#\n"); 
}

void printCurrentInformationAboutRun(const char *directory_name, const DiskParameters *disk_params, const SimulationOptions *sim_opts) {

    char *full_path = NULL;
    char file_name[100];
    
    time_t rawtime;
    struct tm *timeinfo;
    char time_buffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", timeinfo);

    char hostname[1024];
    gethostname(hostname, 1024);
    struct passwd *pw = getpwuid(getuid());
    const char *user = pw ? pw->pw_name : "unknown";

    sprintf(file_name, "%s%s", kCurrentInfoFile, kFileNamesSuffix);
    asprintf(&full_path, "%s/%s", directory_name, file_name);

    LOG_INFO("Writing run info to: '%s'\n", full_path);

    FILE *info_file = fopen(full_path, "w");
    if (info_file == NULL) {
        LOG_ERROR("Could not open info file '%s'.\n", full_path);
        if (full_path) free(full_path);
        return;
    }

    fprintf(info_file, "==========================================================\n");
    fprintf(info_file, "        DUST DRIFT SIMULATION - RUN SNAPSHOT\n");
    fprintf(info_file, "==========================================================\n");
    fprintf(info_file, "  Run Started:       %s\n", time_buffer);
    fprintf(info_file, "  Host Machine:      %s\n", hostname);
    fprintf(info_file, "  User:              %s\n", user);
    
#ifdef _OPENMP
    fprintf(info_file, "  Parallel Threads:  %d\n", omp_get_max_threads());
#else
    fprintf(info_file, "  Parallel Threads:  1 (Serial mode)\n");
#endif

    fprintf(info_file, "  Binary Compiled:   %s %s\n", __DATE__, __TIME__);
    fprintf(info_file, "  Output Directory:  %s\n", directory_name);
    fprintf(info_file, "==========================================================\n\n");
    fprintf(info_file, "\n--- [ Time Parameter ] ---\n");
    fprintf(info_file, "  Total Time:        %.3e yr\n\n", sim_opts->maximum_simulation_time);

    fprintf(info_file, "\n--- [ Central Star ] ---\n");
    fprintf(info_file, "  Stellar Mass:      %.4f M_Sun\n\n", disk_params->stellar_mass);

    fprintf(info_file, "\n--- [ Disk Geometry & Gas ] ---\n");
    fprintf(info_file, "  Radial Range:      %.2f - %.2f AU\n", disk_params->r_min, disk_params->r_max);
    fprintf(info_file, "  Gas Grid Points:   %d\n", disk_params->grid_number);
    fprintf(info_file, "  Sigma_0 (1 AU):    %.4e M_Sun/AU^2\n", disk_params->sigma_0);
    fprintf(info_file, "  Sigma Exponent:    %.4f\n", disk_params->sigma_power_law_index);
    fprintf(info_file, "  Aspect Ratio (H/R): %.4f\n", disk_params->h_aspect_ratio);
    fprintf(info_file, "  Density Floor:     %.3e M_Sun/AU^2\n", disk_params->density_floor);
    fprintf(info_file, "  Flaring Index:     %.4f\n", disk_params->flaring_index);
    fprintf(info_file, "  Alpha Viscosity:   %.4e\n", disk_params->alpha_parameter);
    fprintf(info_file, "  Total Disk Mass:   %.4e M_Sun\n\n", disk_params->total_disk_mass);


    fprintf(info_file, "\n--- [ Gas Cutoff Profile ] ---\n");
    fprintf(info_file, "  Use Cutoff:        %s\n", disk_params->cutoff ? "YES" : "NO");
    if (disk_params->cutoff) {
        fprintf(info_file, "  Cutoff Radius:     %.2f AU\n", disk_params->r_cutoff);
        fprintf(info_file, "  Cutoff Sharpness:  %.2f\n", disk_params->n_cutoff);
    }

    fprintf(info_file, "\n");

    fprintf(info_file, "\n--- [ Dead Zone Configuration ] ---\n");
    if (disk_params->r_dze_i > 0.0 || disk_params->r_dze_o > 0.0) {
        fprintf(info_file, "  Status:            ENABLED\n");
        fprintf(info_file, "  Inner DZE Radius:  %.2f AU (Trans. width: %.2f)\n",
                disk_params->r_dze_i, disk_params->dr_dze_i);
        fprintf(info_file, "  Outer DZE Radius:  %.2f AU (Trans. width: %.2f)\n",
                disk_params->r_dze_o, disk_params->dr_dze_o);
        fprintf(info_file, "  Alpha Mod Factor:  %.4e\n", disk_params->alpha_parameter_modification);
    } else {
        fprintf(info_file, "  Status:            DISABLED (Uniform alpha)\n");
    }

    fprintf(info_file, "\n");

    fprintf(info_file, "\n--- [ Boundary Conditions ] ---\n");
    fprintf(info_file, "  Inner BC:          %s\n", disk_params->inner_bc_string);
    fprintf(info_file, "  Outer BC:          %s\n\n", disk_params->outer_bc_string);

    fprintf(info_file, "\n--- [ Photoevaporation ] ---\n");
    if (strcmp(disk_params->photoevaporation_mode_string, "none") != 0) {
        fprintf(info_file, "  Mode:              DISABLED\n");
    } else {
        fprintf(info_file, "  Mode:              %s\n", disk_params->photoevaporation_mode_string);
        fprintf(info_file, "  X-ray Luminosity:  %.3e erg/s\n", disk_params->xray_luminosity);
    }
    fprintf(info_file, "\n");

    fprintf(info_file, "\n--- [ Dust Properties ] ---\n");
    fprintf(info_file, "  Particle Density:  %.2f g/cm^3\n", disk_params->particle_density);
    fprintf(info_file, "  Fragmentation Vel: %.2f cm/s\n", disk_params->fragmentation_velocity);
    fprintf(info_file, "  Global Dust Count: %d\n", particle_number);
    fprintf(info_file, "  Dust Smoothing:    %s\n", disk_params->dust_smoothing_mode_string);
    if (sim_opts->dust_smoothing_mode == 3) {
        fprintf(info_file, "    Gaussian Sigma:    %.2f grid units\n", disk_params->gaussian_sigma);
        fprintf(info_file, "    Gaussian Cutoff:   %.2f sigma\n", disk_params->gaussian_cutoff);
    }

    fprintf(info_file, "\n");


    fprintf(info_file, "\n==========================================================\n");
    fprintf(info_file, "         End of Configuration Summary\n");
    fprintf(info_file, "==========================================================\n");
    fclose(info_file);
    if (full_path) free(full_path);
}

void printTrapMassEvolution(double output_time, int num_found, const PressureTrap *traps, OutputFiles *output_files) {
    
    if (output_files->mass_file == NULL) return;

    // Idő (Év)

    fprintf(output_files->mass_file, "%-10.0f", output_time);

    // Végigmegyünk az 5 lehetséges sloton
    for (int i = 0; i < 3; i++) {
        if (i < num_found && traps[i].radial_position > 0.0) {
            // Távolság (AU) | Pillanatnyi tömeg (M_sun)
            fprintf(output_files->mass_file, " %-12.4f %-15.6e", 
                    traps[i].radial_position, 
                    traps[i].total_dust_mass);
        } else {
            // Ha nincs csapda, nullákat írunk, hogy az oszlopszám fix maradjon
            fprintf(output_files->mass_file, " %-12.4f %-15.6e", 0.0, 0.0);
        }
    }
    fprintf(output_files->mass_file, "\n");
    fflush(output_files->mass_file);
}


void printGasSurfaceDensityPressurePressureDerivateFile(const DiskParameters *disk_params, OutputFiles *output_files) {

    int i;

    if (output_files->surface_file == NULL) {
        LOG_ERROR("output_files->surface_file is NULL in printGasSurfaceDensityPressurePressureDerivateFile! Cannot write sigma data.\n");
        return;
    }

    for(i = 1; i <= disk_params->grid_number; i++) { 
        fprintf(output_files->surface_file, "%-5d %-15.6e %-15.6le %-15.6e %15.6e\n",i, disk_params->radial_grid[i], disk_params->gas_surface_density_vector[i], disk_params->gas_pressure_vector[i], disk_params->gas_pressure_gradient_vector[i]);
    }

    fflush(output_files->surface_file);
}

void printDustSurfaceDensityPressurePressureDerivateFile(const double *r, const double *rm, const double *dust_surfacedensity, const double *micron_dust_surfacedensity, const DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files, double step) {

    int i;

    if (output_files->dust_file == NULL) {
        LOG_ERROR("output_files->dust_file is NULL in printDustSurfaceDensityPressurePressureDerivateFile! Cannot write main dust surface density.\n");
        return;
    }

    for(i=0;i<particle_number;i++){

        if (r[i] >= disk_params->r_min) { 
            fprintf(output_files->dust_file,"%-16lg %-23.5f %-20.6e\n",(double)step,r[i],dust_surfacedensity[i]);
        }

        if(sim_opts->option_for_dust_secondary_population == 1.0 && output_files->micron_dust_file != NULL) {
            if (rm[i] >= disk_params->r_min) {
                fprintf(output_files->micron_dust_file,"%-16lg %-23.5f %-20.6e\n",(double)step,rm[i],micron_dust_surfacedensity[i]);
            }
        }
    }

    fflush(output_files->dust_file);
    if(sim_opts->option_for_dust_secondary_population == 1.0 && output_files->micron_dust_file != NULL) {
        fflush(output_files->micron_dust_file); 
    }
}

/**
 * @brief Prints dust particle radial positions and sizes to ASCII output files.
 * 
 * Uses the ParticleData structure directly to ensure the most recent particle
 * coordinates are accessed.
 */
void printDustParticleSizeFile(char *size_name, char *size_name2, int step, 
                               ParticleData *particle_data, 
                               const DiskParameters *disk_params, 
                               const SimulationOptions *sim_opts, 
                               OutputFiles *output_files, 
                               SnapshotMode mode) {

    FILE *fout_size = NULL;
    FILE *fout_size2 = NULL;
    int i;
    int is_twopop = isSecondaryPopulationEnabled(mode);

    // Open output files if dust is enabled in the current mode.
    // Note: Removed 'isDustGrowthEnabled' dependency so files are updated during drift.
    if (isDustEnabled(mode)) {
        fout_size = openSnapshotFile(size_name, FILE_TYPE_PARTICLE_SIZE, (double)step / (2.0 * M_PI));        
        if (fout_size == NULL) {
            LOG_ERROR("Could not open size file '%s' in printDustParticleSizeFile!\n", size_name);
            return;
        }

        if (is_twopop && particle_data->micron_particle_distance_array != NULL) {
            fout_size2 = openSnapshotFile(size_name2, FILE_TYPE_MICRON_PARTICLE_SIZE, (double)step / (2.0 * M_PI));
            if (fout_size2 == NULL) {
                LOG_ERROR("Could not open micron size file '%s' in printDustParticleSizeFile!\n", size_name2);
                fclose(fout_size);
                return;
            }
        }
    }

    // Iterate through particles using the direct ParticleData pointers
    for (i = 0; i < particle_number; i++) { 
        
        // Write primary population size data
        if (fout_size != NULL && particle_data->particle_distance_array[i][0] >= disk_params->r_min) {
            fprintf(fout_size, "%-10d %-10d %-20.6e %-20.6e\n", 
                    step, i, 
                    particle_data->particle_distance_array[i][0], 
                    particle_data->particle_distance_array[i][1] * AU_IN_CM);
        }

        // Write secondary population size data if enabled
        if (is_twopop && fout_size2 != NULL && particle_data->micron_particle_distance_array != NULL) {
            if (particle_data->micron_particle_distance_array[i][0] >= disk_params->r_min) {
                fprintf(fout_size2, "%-10d %-10d %-20.6e %-20.6e\n", 
                        step, i, 
                        particle_data->micron_particle_distance_array[i][0], 
                        particle_data->micron_particle_distance_array[i][1] * AU_IN_CM);
            }
        }
    }

    // Close files safely
    if (fout_size != NULL) fclose(fout_size);
    if (fout_size2 != NULL) fclose(fout_size2);
}

void printFileHeader(FILE *file, FileType_e file_type, const HeaderData *header_data) {
    
    if (file == NULL) {
        LOG_ERROR("Attempted to write header to a NULL file pointer!\n");
        return;
    }

    char buffer[128]; 
    
    fprintf(file, "#############################################################################\n");
    printCentered(file, "Generated by Dust Drift Simulation with the RAPID code");
    sprintf(buffer, "Version: %s | Compiled: %s %s", SIM_VERSION, __DATE__, __TIME__);
    printCentered(file, buffer);
    fprintf(file, "#############################################################################\n");

    switch (file_type) {

/*        case FILE_TYPE_MASS_ACCUMULATION:
            fprintf(file, "# This file contains the time evolution of dust mass within specified disk regions.\n");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-5s %-15s %-15s %-15s %-15s\n",
                    "Time", "Inner_DZE_AU", "Accum._mass_inner_MSun", "Outer_DZE_AU", "Accum._mass_outer_MSun");
            fprintf(file, "#--------------------------------------------------------------------------\n");            
            break;
*/


        case FILE_TYPE_MASS_ACCUMULATION:
            fprintf(file, "# Instantaneous dust mass inventory in identified pressure traps (0.5H width)\n");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-13s  %-12s %-15s  %-12s %-15s  %-12s %-15s\n",
                    "Time_yr", "RDzei_AU", "MassDzei_Msun", "RDzeo_AU", "MassDzeo_Msun", "RSecondary_AU", "MassSecondary_Msun");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            break;

        case FILE_TYPE_GAS_DENSITY:
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial gas profile\n");
            } else {
                fprintf(file, "# This file contains the time evolution of the gas component at %e years\n", header_data ? header_data->current_time : 0.0);
            }
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-5s %-15s %-15s %-15s %-15s\n",
                    "Index", "Radius_AU", "GasSurfDensity", "GasPressure", "GasPressureDeriv");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            break;

        case FILE_TYPE_DUST_DENSITY:
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial dust profile\n");
            } else {
                fprintf(file, "# This file contains the time evolution of the dust component at %e years\n", header_data ? header_data->current_time : 0.0);
            }
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-10s %-20s %-20s\n",
                    "Time_step","Radial_distance_AU", "Dust_surfacedensity_M_Sun/AU^2");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            break;

        case FILE_TYPE_DUST_MICRON_DENSITY:
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial micron-sized dust profile\n");
            } else {
                fprintf(file, "# This file contains the time evolution of the micron-sized dust component at %e years\n", header_data ? header_data->current_time : 0.0);
            }
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-10s %-20s %-20s\n",
                    "Time_step","Radial_distance_AU", "Dust_surfacedensity_M_Sun/AU^2");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            break;

        case FILE_TYPE_PARTICLE_SIZE:
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial particle size profile\n");
            } else {
                fprintf(file, "# This file contains the size of each dust particle at %e years\n", header_data ? header_data->current_time : 0.0);
            }
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-10s %-10s  %-20s %-20s\n",
                    "Time_step", "Index" ,"Radial_distance_AU", "Particle_size_cm");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            break;            

        case FILE_TYPE_MICRON_PARTICLE_SIZE:
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial micron particle size profile\n");
            } else {
                fprintf(file, "# This file contains the size of each micron dust particle at %e years\n", header_data ? header_data->current_time : 0.0);
            }
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-10s %-10s  %-20s %-20s\n",
                    "Time_step", "Index" ,"Radial_distance_AU", "Micron_particle_size_cm");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            break;

        case FILE_TYPE_INTIIAL_DUST_PROFILE:
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial particle distribution\n");
            } else {
                fprintf(file, "# Particle distribution, Time: %e years\n", header_data ? header_data->current_time : 0.0);
            }
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-5s %-15s %-20s %-20s %-15s %-15s\n",
                    "Index", "Radius_AU", "RepMass_Pop1_Msun", "RepMass_Pop2_Msun", "MaxPartSize_cm", "MicroSize_cm");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            break;

        case FILE_TYPE_DISK_PARAM:
            fprintf(file, "# Disk Parameters\n");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            fprintf(file, "# %-15s %-15s %-10s %-15s %-20s %-15s %-15s %-15s %-20s %-20s %-15s %-15s %-15s %-15s %-15s\n",
                    "R_Min_AU", "R_Max_AU", "N_Grid", "SigmaExp", "Sigma0_gas_Msun_AU2",
                    "G_GravConst", "DzR_Inner_AU", "DzR_Outer_AU", "DzDr_Inner_Calc_AU", "DzDr_Outer_Calc_AU",
                    "DzAlphaMod", "DustDensity_g_cm3", "AlphaViscosity", "StarMass_Msun", "FlaringIndex");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            break;

        default:
            LOG_WARN("Unknown file type for header generation: %d!\n", file_type);
            break;
    }
    fflush(file);
}



void printFinalSimulationSummary(const char *directory_name, double elapsed_seconds, const SimulationOptions *sim_opts) {

    char *full_path = NULL;
    asprintf(&full_path, "%s/%s%s", directory_name, kCurrentRuntimeInfoFile,kFileNamesSuffix);

    FILE *f = fopen(full_path, "w");
    if (f == NULL) {
        if (full_path) free(full_path);
        return;
    }

    int h = (int)elapsed_seconds / 3600;
    int m = ((int)elapsed_seconds % 3600) / 60;
    int s = (int)elapsed_seconds % 60;
    double speed = sim_opts->maximum_simulation_time / elapsed_seconds;

    fprintf(f, "==========================================================\n");
    fprintf(f, "             SIMULATION PERFORMANCE SUMMARY\n");
    fprintf(f, "==========================================================\n");
    fprintf(f, "  Total Wall-Clock Time:  %02d:%02d:%02d (H:M:S)\n", h, m, s);
    fprintf(f, "  Total Simulation Time:  %.2e years\n", sim_opts->maximum_simulation_time);
    fprintf(f, "  Execution Speed:        %.2f simulated years / wall-sec\n", speed);
    fprintf(f, "==========================================================\n");
    fprintf(f, "Status: Finished Normally\n");
    fclose(f);
    
    if (full_path) free(full_path);
    
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "Simulation finished in %02d:%02d:%02d. Summary written to logs.", h, m, s);
    printHeader(buffer);
}


int setupInitialOutputFiles(OutputFiles *output_files, const SimulationOptions *sim_opts,
                               const DiskParameters *disk_params, HeaderData *header_data_for_files) {

    if(OUTPUT_HDF5){ 
        return 0;
    } else {
        char *mass_output = NULL;

        header_data_for_files->current_time = 0.0;
        header_data_for_files->is_initial_data = 1;
        header_data_for_files->R_in = disk_params->r_min;
        header_data_for_files->R_out = disk_params->r_max;

        asprintf(&mass_output, "%s/%s/%s%s", sim_opts->output_dir_name, kLogFilesDirectory, kDustAccumulationFileName, kFileNamesSuffix);
        LOG_INFO("Writing initial mass output to: %s\n", mass_output);

        output_files->mass_file = fopen(mass_output, "w");
        if (output_files->mass_file == NULL) {
            LOG_ERROR("Could not open %s\n", mass_output);
            goto cleanup_error;
        }

        printFileHeader(output_files->mass_file, FILE_TYPE_MASS_ACCUMULATION, header_data_for_files);
        free(mass_output);

        return 0;

        

    cleanup_error:
        if (output_files->mass_file) fclose(output_files->mass_file);
        
        if (mass_output) free(mass_output);
        
        return 1;

    }
}

void cleanupSimulationResources(ParticleData *particle_data, OutputFiles *output_files) {

    if (particle_number > 0) {
        free(particle_data->particle_distance_array); particle_data->particle_distance_array = NULL;
        free(particle_data->micron_particle_distance_array); particle_data->micron_particle_distance_array = NULL;
        free(particle_data->dust_particle_mass_grid); particle_data->dust_particle_mass_grid = NULL;
        free(particle_data->massmicradial_grid); particle_data->massmicradial_grid = NULL;
        free(particle_data->dust_particle_mass_array); particle_data->dust_particle_mass_array = NULL;
        free(particle_data->micron_dust_particle_mass_array); particle_data->micron_dust_particle_mass_array = NULL;
        free(particle_data->dust_surfacedensity); particle_data->dust_surfacedensity = NULL;
        free(particle_data->micron_dust_surfacedensity); particle_data->micron_dust_surfacedensity = NULL;
        free(particle_data->particle_distance_grid); particle_data->particle_distance_grid = NULL;
        free(particle_data->micron_particle_distance_grid); particle_data->micron_particle_distance_grid = NULL;
        LOG_DEBUG("All dynamically allocated particle arrays freed.\n");
    }

    if (output_files->mass_file != NULL) {
        fclose(output_files->mass_file);
        output_files->mass_file = NULL;
        LOG_DEBUG("Closed %s%s\n", kDustAccumulationFileName,kFileNamesSuffix);
    }
}


FILE *openSnapshotFile(const char *file_name, FileType_e file_type, double current_time_years){

    FILE *file = fopen(file_name, "w");
    if (file == NULL) {
        LOG_ERROR("Could not open %s for writing.\n", file_name);
        return NULL;
    }

    HeaderData header = {
        .current_time = current_time_years,
        .is_initial_data = (current_time_years == 0.0)
    };

    printFileHeader(file, file_type, &header);

    return file;
}


void buildSnapshotFilenames(char *dens_name, char *dust_name, char *dust_name2, char *size_name, char *size_name2, 
                            const SimulationOptions *sim_opts, int snapshot_id, SnapshotMode mode) {

    sprintf(dens_name, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory, kGasDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix);

    sprintf(dust_name, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory, kDustDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix);

    // Használd a segédfüggvényt a sim_opts flag helyett
    if (isSecondaryPopulationEnabled(mode)) {
        sprintf(dust_name2, "%s/%s/%s_%08d%s",
                 sim_opts->output_dir_name, kLogFilesDirectory, kMicronDustDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix); 
    }

    sprintf(size_name, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory, kDustParticleSizeFileName, snapshot_id, kFileNamesSuffix);

    if (isSecondaryPopulationEnabled(mode)) {
        sprintf(size_name2, "%s/%s/%s_%08d%s",
                 sim_opts->output_dir_name, kLogFilesDirectory, kMicronDustParticleSizeFileName, snapshot_id, kFileNamesSuffix);
    }
}

void closeSnapshotFiles(OutputFiles *output_files, SnapshotMode mode) {

    if (output_files->surface_file != NULL) {
        fclose(output_files->surface_file);
        output_files->surface_file = NULL;
    }
    if (output_files->dust_file != NULL) {
        fclose(output_files->dust_file);
        output_files->dust_file = NULL;
    }

    if (isSecondaryPopulationEnabled(mode)) {
        if (output_files->micron_dust_file != NULL) {
            fclose(output_files->micron_dust_file);
            output_files->micron_dust_file = NULL;
        }
        if (output_files->size_micron_file != NULL) {
            fclose(output_files->size_micron_file);
            output_files->size_micron_file = NULL;
        }
    }
    
    if (output_files->size_file != NULL) {
        fclose(output_files->size_file);
        output_files->size_file = NULL;
    }
}