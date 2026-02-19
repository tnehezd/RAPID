#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>  
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h> 
#include <pwd.h>    


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
        fprintf(stderr, "ERROR [calculateNumbersOfParticles]: Could not open file '%s'.\n", particle_data_file_name);
        perror("Reason"); 
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < INIT_DATA_HEADER_LINES; i++) {
        if (fgets(line_buffer, sizeof(line_buffer), particle_file) == NULL) {
            fprintf(stderr, "ERROR [calculateNumbersOfParticles]: Unexpected end of file while skipping %d header lines in '%s'.\n", INIT_DATA_HEADER_LINES, particle_data_file_name);
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

    load_dust_particles_file = fopen(particle_data_file_name,"r");

    if (load_dust_particles_file == NULL) {
        fprintf(stderr, "ERROR [loadDustParticlesFromFile]: Could not open file '%s'.\n", particle_data_file_name);
        perror("Reason");
        exit(EXIT_FAILURE);
    }


    char line_buffer[1024];
    for (int k = 0; k < INIT_DATA_HEADER_LINES; k++) {
        if (fgets(line_buffer, sizeof(line_buffer), load_dust_particles_file) == NULL) {
            fprintf(stderr, "ERROR [loadDustParticlesFromFile]: Unexpected end of file while skipping headers in '%s'.\n", particle_data_file_name);
            fclose(load_dust_particles_file);
            exit(EXIT_FAILURE);
        }
    }


    for (i = 0; i < particle_number; i++) {
        if(fscanf(load_dust_particles_file,"%d %lg %Lg %Lg %lg %lg",&particle_index,&distance,&representative_mass,&micron_representative_mass,&particle_radius,&micron_particle_radius) == 6) {
            particle_data->particle_distance_array[i][0] = distance;
            particle_data->particle_distance_array[i][1] = particle_radius / AU_IN_CM; 
            particle_data->dust_particle_mass_grid[i] = representative_mass;

            particle_data->micron_particle_distance_array[i][0] = distance;
            particle_data->micron_particle_distance_array[i][1] = micron_particle_radius / AU_IN_CM; 
            particle_data->massmicradial_grid[i] = micron_representative_mass;
        } else {
            fprintf(stderr, "\n\n******************* ERROR!      *********************\n\n");
            fprintf(stderr, "   Failed to read line %d from particle data file '%s'!\n", i, particle_data_file_name);
            fprintf(stderr, "   Expected 6 values, but fscanf failed. Program will exit.\n");
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
        fprintf(stderr, "ERROR [loadGasSurfaceDensityFromFile]: Could not open input file '%s'.\n", input_disk_file_name);
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
        fprintf(stderr, "ERROR [loadGasSurfaceDensityFromFile]: File '%s' is empty or only contains comments/headers.\n", input_disk_file_name);
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
            fprintf(stderr, "ERROR [loadGasSurfaceDensityFromFile]: Failed to read 4 values for row %d from file '%s'. File may be malformed or ended unexpectedly.\n", i, input_disk_file_name);
            fclose(input_file);
            exit(EXIT_FAILURE);
        }

        if ((i + 1) >= 0 && (i + 1) <= disk_params->grid_number + 1) { 
            disk_params->radial_grid[i + 1] = r_value;
            disk_params->gas_surface_density_vector[i + 1] = surfacedensity_gas_value;
            disk_params->gas_pressure_vector[i + 1] = gas_pressure_value;
            disk_params->gas_pressure_gradient_vector[i + 1] = gas_pressure_gradient_value;
        } else {
            fprintf(stderr, "WARNING [loadGasSurfaceDensityFromFile]: Attempted to write to out-of-bounds index %d. Max allowed index: %d (grid_number+1).\n", i + 1, disk_params->grid_number + 1);
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
            fprintf(stderr, "ERROR: Too many directories.\n");
            exit(1);
        }
    }

    if (MKDIR_CALL(temporary_path) != 0) {
        perror("mkdir failed");
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

void printCurrentInformationAboutRun(const char *directory_name, const DiskParameters *disk_params) {

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

    fprintf(stderr, "DEBUG [printCurrentInformationAboutRun]: Writing run info to: '%s'\n", full_path);

    FILE *info_file = fopen(full_path, "w");
    if (info_file == NULL) {
        fprintf(stderr, "ERROR: Could not open info file '%s'.\n", full_path);
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
    fprintf(info_file, "--- [ Central Star ] ---\n");
    fprintf(info_file, "  Stellar Mass:      %.4f M_Sun\n\n", disk_params->stellar_mass);
    fprintf(info_file, "--- [ Disk Geometry & Gas ] ---\n");
    fprintf(info_file, "  Radial Range:      %.2f - %.2f AU\n", disk_params->r_min, disk_params->r_max);
    fprintf(info_file, "  Gas Grid Points:   %d\n", disk_params->grid_number);
    fprintf(info_file, "  Sigma_0 (1 AU):    %.4e M_Sun/AU^2\n", disk_params->sigma_0);
    fprintf(info_file, "  Sigma Exponent:    %.4f\n", disk_params->sigma_power_law_index);
    fprintf(info_file, "  Aspect Ratio (H/R): %.4f\n", disk_params->h_aspect_ratio);
    fprintf(info_file, "  Flaring Index:     %.4f\n", disk_params->flaring_index);
    fprintf(info_file, "  Alpha Viscosity:   %.4e\n\n", disk_params->alpha_parameter);

    fprintf(info_file, "--- [ Dead Zone Configuration ] ---\n");
    if (disk_params->r_dze_i > 0.0 || disk_params->r_dze_o > 0.0) {
        fprintf(info_file, "  Status:            ENABLED\n");
        fprintf(info_file, "  Inner DZE Radius:  %.2f AU (Trans. width: %.2f)\n", disk_params->r_dze_i, disk_params->dr_dze_i);
        fprintf(info_file, "  Outer DZE Radius:  %.2f AU (Trans. width: %.2f)\n", disk_params->r_dze_o, disk_params->dr_dze_o);
        fprintf(info_file, "  Alpha Mod Factor:  %.4e\n", disk_params->alpha_parameter_modification);
    } else {
        fprintf(info_file, "  Status:            DISABLED (Uniform alpha)\n");
    }

    fprintf(info_file, "\n--- [ Dust Properties ] ---\n");
    fprintf(info_file, "  Particle Density:  %.2f g/cm^3\n", disk_params->particle_density);
    fprintf(info_file, "  Fragmentation Vel: %.2f cm/s\n", disk_params->fragmentation_velocity);
    fprintf(info_file, "  Global Dust Count: %d\n", particle_number);
    fprintf(info_file, "\n==========================================================\n");
    fprintf(info_file, "         End of Configuration Summary\n");
    fprintf(info_file, "==========================================================\n");
    fclose(info_file);
    if (full_path) free(full_path);
}

void printTrapMassEvolution(double output_time, int num_found, const PressureTrap *traps, OutputFiles *output_files) {
    
    if (output_files->mass_file == NULL) {
        fprintf(stderr, "NO MASS FILE WRITTEN\n");
        return;
    }

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
        fprintf(stderr, "ERROR: output_files->surface_file is NULL in printGasSurfaceDensityPressurePressureDerivateFile! Cannot write sigma data.\n");
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
        fprintf(stderr, "ERROR: output_files->dust_file is NULL in printDustSurfaceDensityPressurePressureDerivateFile! Cannot write main dust surface density.\n");
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

void printDustParticleSizeFile(char *size_name, int step, double (*rad)[3], double (*micron_particle_radius)[3], const DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files) {

    FILE *fout_size = NULL;
    int i;

    if (sim_opts->option_for_dust_growth == 1.0) {
        fout_size = openSnapshotFile(size_name, FILE_TYPE_PARTICLE_SIZE, (double)step / (2.0 * M_PI));        
        if (fout_size == NULL) {
            fprintf(stderr, "ERROR: Could not open size file '%s' in printDustParticleSizeFile!\n", size_name);
            return;
        }
    }

    for (i = 0; i < particle_number; i++) { 
        if (sim_opts->option_for_dust_secondary_population == 1.0 && output_files->micron_motion_file != NULL) {
            if (micron_particle_radius[i][0] >= disk_params->r_min) { 
                fprintf(output_files->micron_motion_file, "%-16lg %-10d %-20.6e %-20.6e\n", (double)step, i, micron_particle_radius[i][0],micron_particle_radius[i][1] * AU_IN_CM);
            }
        }
        if (sim_opts->option_for_dust_growth == 1.0 && fout_size != NULL) {
            if (rad[i][0] >= disk_params->r_min) {
                fprintf(fout_size, "%-10lg %-10d %-20.6e %-20.6e\n", (double)step, i, rad[i][0], rad[i][1] * AU_IN_CM);
            }
        }
    }

    if (sim_opts->option_for_dust_growth == 1.0 && fout_size != NULL) {
        fclose(fout_size);
    }
}

void printFileHeader(FILE *file, FileType_e file_type, const HeaderData *header_data) {
    
    if (file == NULL) {
        fprintf(stderr, "ERROR [printFileHeader]: Attempted to write header to a NULL file pointer!\n");
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
            fprintf(stderr, "WARNING [printFileHeader]: Unknown file type for header generation: %d!\n", file_type);
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
    
    printf("\n>> Simulation finished in %02d:%02d:%02d. Summary written to logs.\n", h, m, s);
}


int setupInitialOutputFiles(OutputFiles *output_files, const SimulationOptions *sim_opts,
                               const DiskParameters *disk_params, HeaderData *header_data_for_files) {

    if (sim_opts->output_format == OUTPUT_HDF5) { 
            return 0; // Ha HDF5-öt kért a júzer, akkor nem kell ASCII fájlokat nyitni
        
    } else {
        char *mass_output = NULL;

        header_data_for_files->current_time = 0.0;
        header_data_for_files->is_initial_data = 1;
        header_data_for_files->R_in = disk_params->r_min;
        header_data_for_files->R_out = disk_params->r_max;

        asprintf(&mass_output, "%s/%s/%s%s", sim_opts->output_dir_name, kLogFilesDirectory, kDustAccumulationFileName, kFileNamesSuffix);
        fprintf(stderr, "DEBUG [setupInitialOutputFiles]: Paths:\n  Mass: %s\n", mass_output);

        output_files->mass_file = fopen(mass_output, "w");
        if (output_files->mass_file == NULL) {
            fprintf(stderr, "ERROR: Could not open %s\n", mass_output);
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
        fprintf(stderr, "DEBUG [cleanupSimulationResources]: All dynamically allocated particle arrays freed.\n");
    }

    if (output_files->mass_file != NULL) {
        fclose(output_files->mass_file);
        output_files->mass_file = NULL;
        fprintf(stderr, "DEBUG [cleanupSimulationResources]: Closed %s%s\n", kDustAccumulationFileName,kFileNamesSuffix);
    }
}


FILE *openSnapshotFile(const char *file_name, FileType_e file_type, double current_time_years){

    FILE *file = fopen(file_name, "w");
    if (file == NULL) {
        fprintf(stderr, "ERROR: Could not open %s for writing.\n", file_name);
        return NULL;
    }

    HeaderData header = {
        .current_time = current_time_years,
        .is_initial_data = (current_time_years == 0.0)
    };

    printFileHeader(file, file_type, &header);

    return file;
}


void buildSnapshotFilenames(char *dens_name, char *dust_name, char *dust_name2, char *size_name, const SimulationOptions *sim_opts, int snapshot_id){

    sprintf(dens_name, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kGasDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix);

    sprintf(dust_name, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kDustDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix);

    if (sim_opts->option_for_dust_secondary_population == 1) {
        sprintf(dust_name2, "%s/%s/%s_%08d%s",
                 sim_opts->output_dir_name, kLogFilesDirectory,
                 kMicronDustDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix); 
    }

    sprintf(size_name, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kDustParticleSizeFileName, snapshot_id, kFileNamesSuffix);
}

void closeSnapshotFiles(OutputFiles *output_files, const SimulationOptions *sim_opts) {

    if (output_files->surface_file != NULL) {
        fclose(output_files->surface_file);
        output_files->surface_file = NULL;
    }
    if (output_files->dust_file != NULL) {
        fclose(output_files->dust_file);
        output_files->dust_file = NULL;
    }

    if (sim_opts->option_for_dust_secondary_population == 1 && output_files->micron_dust_file != NULL) {
        fclose(output_files->micron_dust_file);
        output_files->micron_dust_file = NULL;
    }
}



void writeDustField2D(const StructuredParticleData *sdata, const char *directory, int snapshot_index, const char *label) {
    if (!sdata || !sdata->particles) {
        fprintf(stderr,"ERROR [writeDustField2D]: No structured dust field allocated\n");
        return;
    }

    char *mass_path = NULL;
    char *grid_path = NULL;

    if (snapshot_index >= 0) {
        asprintf(&mass_path,"%s/%s/%s_%08d%s",directory,kLogFilesDirectory,kMassFieldNameFile,snapshot_index,kFileNamesSuffix);
        asprintf(&grid_path,"%s/%s/%s_%08d%s",directory,kLogFilesDirectory,kGridFieldNameFile,snapshot_index,kFileNamesSuffix);
    }
    else {
        if (label)
        {
            asprintf(&mass_path,"%s/%s_%s%s",directory,label,kMassFieldNameFile,kFileNamesSuffix);
            asprintf(&grid_path,"%s/%s_%s%s",directory,label,kGridFieldNameFile,kFileNamesSuffix);
        }
        else
        {
            asprintf(&mass_path,"%s/%s%s",directory,kMassFieldNameFile,kFileNamesSuffix);
            asprintf(&grid_path,"%s/%s%s",directory,kGridFieldNameFile,kFileNamesSuffix);
        }
    }

    FILE *fm = fopen(mass_path,"w");
    FILE *fg = fopen(grid_path,"w");

    if(!fm || !fg){
        fprintf(stderr,"ERROR [writeDustField2D]: cannot open output files\n");
        free(mass_path); free(grid_path);
        return;
    }

    for (size_t i = 0; i < sdata->n_r; i++) {
//        const DustParticle *p0 = &sdata->particles[i][0];

        // fg fájl: első oszlopban r_au, utána a z_au értékek sorban
//        fprintf(fg, "%e", p0->r_au);
        for (size_t j = 0; j < sdata->n_z; j++) {
            const DustParticle *p = &sdata->particles[i][j];
            fprintf(fg, " %e %e", p->r_au, p->z_au);
        }
        fprintf(fg, "\n");

        // fm fájl: első oszlopban r_au, utána a mass_g értékek sorban
//        fprintf(fm, "%e", p0->r_au);
        for (size_t j = 0; j < sdata->n_z; j++) {
            const DustParticle *p = &sdata->particles[i][j];
            fprintf(fm, " %e", p->mass_g);
        }
        fprintf(fm, "\n");
    }

    fclose(fm);
    fclose(fg);


    fprintf(stderr,"DEBUG: 2D dust field written (%s)\n",
            snapshot_index>=0 ? "snapshot" : (label?label:"plain"));

    free(mass_path);
    free(grid_path);
}
