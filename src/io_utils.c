#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

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

#include "config.h"
#include "io_utils.h"
#include "dust_physics.h"
#include "simulation_types.h"
#include "globals.h"
#include "particle_data.h"
#include "dust_particle.h"

#include <sys/stat.h>
#include <sys/types.h>


#define INIT_DATA_HEADER_LINES 5

// Returns the number of data lines in a file, which corresponds to the number of particles to be read.
int get_particle_count(const char *filenev) {
    FILE *fp = NULL;
    char line_buffer[1024];
    int line_count = 0; // Counter for data lines

    fp = fopen(filenev, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR [get_particle_count]: Could not open file '%s'.\n", filenev);
        perror("Reason"); // Prints system error message
        exit(EXIT_FAILURE);
    }

    // Skip header lines
    for (int i = 0; i < INIT_DATA_HEADER_LINES; i++) {
        if (fgets(line_buffer, sizeof(line_buffer), fp) == NULL) {
            // If file ends before all header lines are skipped, it's an error
            fprintf(stderr, "ERROR [get_particle_count]: Unexpected end of file while skipping %d header lines in '%s'.\n", INIT_DATA_HEADER_LINES, filenev);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
    }

    // Count remaining data lines
    while (fgets(line_buffer, sizeof(line_buffer), fp) != NULL) {
        // Check to ensure the line is not empty or a comment
        if (line_buffer[0] != '#' && line_buffer[0] != '\n' && line_buffer[0] != '\r') {
             line_count++;
        }
    }

    fclose(fp); // Close the file after reading
    return line_count;
}

// Reads dust particle data (NEW version with ParticleData_t and PARTICLE_NUMBER).
void load_dust_particles(ParticleData_t *p_data, const char *filename,
                       const disk_t *disk_params, simulation_options_t *sim_opts) {
    FILE *fp = NULL;
    int i;
    int dummy_id;
    double distance;
    double particle_radius_pop1;  // Size read for Pop1 (MaxPartSize_cm)
    long double reprmass_pop1;    // Mass read for Pop1 (RepMass_Pop1)

    // Variables for Pop2, only used if twopop is enabled
    double particle_radius_pop2 = 0.0; // Size read for Pop2 (MicroSize_cm)
    long double reprmass_pop2 = 0.0L;  // Mass read for Pop2 (RepMass_Pop2)

    // 1. Determine number of particles (total lines in the file)
    // 'get_particle_count' is assumed to return the count for ONE population.
    // The global PARTICLE_NUMBER will be set by this.
    PARTICLE_NUMBER = get_particle_count(filename);

    if (PARTICLE_NUMBER <= 0) {
        fprintf(stderr, "ERROR [load_dust_particles]: No particles found or error in counting for file '%s'. PARTICLE_NUMBER = %d.\n", filename, PARTICLE_NUMBER);
        exit(EXIT_FAILURE);
    }

    // Based on PARTICLE_NUMBER, set the particle counts in p_data.
    // If twopop is enabled, both populations get PARTICLE_NUMBER particles.
    p_data->num_particles_pop1 = PARTICLE_NUMBER;
    p_data->num_particles_pop2 = (sim_opts->twopop == 1.0) ? PARTICLE_NUMBER : 0;

    // Update sim_opts->num_dust_particles. This usually refers to the main loop count, which is PARTICLE_NUMBER.
    sim_opts->num_dust_particles = PARTICLE_NUMBER;

    // 2. Allocate memory using the dedicated allocation function.
    // NOTE: The allocation function call is currently commented out and inactive.
/*
    allocate_particle_data(p_data, p_data->num_particles_pop1, p_data->num_particles_pop2, (int)(sim_opts->twopop + 0.5));
*/

    // DEBUG: Memory allocation check
    fprintf(stderr, "DEBUG [load_dust_particles]: Particle data arrays allocated. Pop1 allocated: %s, Pop2 allocated: %s (if twopop enabled).\n",
            (p_data->particles_pop1 != NULL ? "YES" : "NO"), (p_data->particles_pop2 != NULL ? "YES" : "NO"));


    // 3. Open the file for reading actual data
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR [load_dust_particles]: Could not open file '%s' for data reading (second pass).\n", filename);
        perror("Reason"); // Print system error message
        // Free previously allocated memory before exiting on error
        free_particle_data(p_data);
        exit(EXIT_FAILURE);
    }

    // Skip header lines
    char line_buffer[1024];
    for (int k = 0; k < INIT_DATA_HEADER_LINES; k++) {
        if (fgets(line_buffer, sizeof(line_buffer), fp) == NULL) {
            fprintf(stderr, "ERROR [load_dust_particles]: Unexpected end of file while skipping headers in '%s'.\n", filename);
            fclose(fp);
            free_particle_data(p_data);
            exit(EXIT_FAILURE);
        }
    }

// 4. Read particle data and populate ParticleData_t
for (i = 0; i < PARTICLE_NUMBER; i++) {
    // The file always has 6 columns, so we always expect 6 values
    int expected_values = 6;

    // Temporary variables to hold ALL 6 columns read from the file
    int temp_dummy_id;
    double temp_distance_au;
    long double temp_reprmass_pop1;
    long double temp_reprmass_pop2;
    double temp_particle_radius_pop1;
    double temp_particle_radius_pop2;

    // Use a single fscanf call to read all 6 columns into temporary variables
    if (fscanf(fp, "%d %lg %Lg %Lg %lg %lg",
               &temp_dummy_id,
               &temp_distance_au,
               &temp_reprmass_pop1,
               &temp_reprmass_pop2,
               &temp_particle_radius_pop1,
               &temp_particle_radius_pop2) != expected_values) {
        fprintf(stderr, "\n\n******************* ERROR! *********************\n\n");
        fprintf(stderr, "   Failed to read line %d from particle data file '%s'! Expected %d values.\n", i, filename, expected_values);
        fprintf(stderr, "   Ensure the file has exactly 6 columns as expected.\n");
        fclose(fp);
        free_particle_data(p_data);
        exit(EXIT_FAILURE);
    }

    // --- Populate the dust_particle_t array for population 1 ---
    p_data->particles_pop1[i].id = temp_dummy_id;
    p_data->particles_pop1[i].distance_au = temp_distance_au;
    p_data->particles_pop1[i].current_size_au = temp_particle_radius_pop1 / AU_TO_CM; // This now correctly gets MaxPartSize_cm
    p_data->particles_pop1[i].initial_mass_msun = temp_reprmass_pop1;

    // Initialize size_reciprocal for Pop1
    if (p_data->particles_pop1[i].current_size_au > 0.0) {
        p_data->particles_pop1[i].size_reciprocal = 1.0 / p_data->particles_pop1[i].current_size_au;
    } else {
        p_data->particles_pop1[i].size_reciprocal = 0.0;
    }

    // --- Populate the dust_particle_t array for population 2 (only if enabled) ---
    if (sim_opts->twopop == 1.0) {
        p_data->particles_pop2[i].id = temp_dummy_id; // Same ID as Pop1 particle
        p_data->particles_pop2[i].distance_au = temp_distance_au; // Same distance as Pop1 particle
        p_data->particles_pop2[i].current_size_au = temp_particle_radius_pop2 /AU_TO_CM; // This correctly gets MicroSize_cm
        p_data->particles_pop2[i].initial_mass_msun = temp_reprmass_pop2; // This correctly gets RepMass_Pop2_Msun

        // Initialize size_reciprocal for Pop2
        if (p_data->particles_pop2[i].current_size_au > 0.0) {
            p_data->particles_pop2[i].size_reciprocal = 1.0 / p_data->particles_pop2[i].current_size_au;
        } else {
            p_data->particles_pop2[i].size_reciprocal = 0.0;
        }
        }
    }

    // 5. Close file
    fclose(fp);
//    fprintf(stderr, "DEBUG [load_dust_particles]: Successfully read %d particle entries from '%s'.\n", PARTICLE_NUMBER, filename);
}


// Reads the gas profile from a file.
void read_gas_profile(disk_t *disk_params, const char *filename) {
    const char *input_filename = filename;

    FILE *fp = fopen(input_filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR [read_gas_profile]: Could not open input file '%s'.\n", input_filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }

    char line[512];

    // Skip header lines and comments.
    // Read lines until a line that does not start with '#' or '---' is found.
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (line[0] == '#' || strncmp(line, "---", 3) == 0) {
            continue; // Skip to the next line
        } else {
            // This is the first data line. Reset the file pointer to the beginning of the line.
            fseek(fp, -strlen(line), SEEK_CUR);
            break;
        }
    }

    // If the file is empty or only contains comments/headers
    if (feof(fp) && (line[0] == '#' || strncmp(line, "---", 3) == 0)) {
        fprintf(stderr, "ERROR [read_gas_profile]: File '%s' is empty or only contains comments/headers.\n", input_filename);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    double r_val;
    double sigma_gas_val;
    double pressure_gas_val;
    double dpressure_dr_val;

    // The loop runs up to disk_params->NGRID
    for (int i = 0; i < disk_params->NGRID; i++) {
        // Read the 4 expected columns: Radius_AU, GasSurfDensity, GasPressure, GasPressureDeriv
        if (fscanf(fp, "%lf %lf %lf %lf",
                                     &r_val, &sigma_gas_val, &pressure_gas_val, &dpressure_dr_val) != 4) {
            // Handle error if 4 double values cannot be read
            fprintf(stderr, "ERROR [read_gas_profile]: Failed to read 4 values for row %d from file '%s'. File may be malformed or ended unexpectedly.\n", i, input_filename);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        // Assign to disk_params arrays
        // The indexing 'i + 1' is used because of the 0-indexed ghost cell.
        // The arrays have a size of disk_params->NGRID + 2. Valid indices are 0 to NGRID+1.
        // The 'real' data goes from 1 to NGRID.
        if ((i + 1) >= 0 && (i + 1) <= disk_params->NGRID + 1) {
            disk_params->rvec[i + 1] = r_val;
            disk_params->sigmavec[i + 1] = sigma_gas_val;
            disk_params->pressvec[i + 1] = pressure_gas_val;
            disk_params->dpressvec[i + 1] = dpressure_dr_val;
        } else {
            fprintf(stderr, "WARNING [read_gas_profile]: Attempted to write to out-of-bounds index %d. Max allowed index: %d (NGRID+1).\n", i + 1, disk_params->NGRID + 1);
        }

    }

    fclose(fp);

}

// Creates a new output directory, adding a numbered suffix if the directory already exists.
void create_output_directory(char *dir_path) {
    char tmp_path[MAX_PATH_LEN];
    int counter = 0;

    // Copy the original name
    strncpy(tmp_path, dir_path, MAX_PATH_LEN - 1);
    tmp_path[MAX_PATH_LEN - 1] = '\0';

    while (access(tmp_path, F_OK) == 0) {   // File or directory exists
        snprintf(tmp_path, MAX_PATH_LEN, "%s_%04d", dir_path, ++counter);
        if (counter > 999) {
            fprintf(stderr, "ERROR [create_output_directory]: Too many existing directories with similar names: %s.\n",counter);
            exit(1);
        }
    }

    int result = MKDIR_CALL(tmp_path);
    if (result != 0) {
        perror("ERROR [create_output_directory]: mkdir failed");
        fprintf(stderr, "ERROR [create_output_directory]: Could not create directory: '%s'\n", tmp_path);
        exit(1);
    }

    fprintf(stderr, "Directory '%s' created successfully.\n", tmp_path);

    // Copy the created directory name back to the input
    strncpy(dir_path, tmp_path, MAX_PATH_LEN - 1);
    dir_path[MAX_PATH_LEN - 1] = '\0';

    fflush(stderr);
}

// Generates a summary log file for the current simulation run.
void write_summary_log(const char *output_dir_name, const disk_t *disk_params, const simulation_options_t *sim_opts) {

    char full_path[MAX_PATH_LEN];
    // char file_name[100]; // Unused variable


    snprintf(full_path, sizeof(full_path), "%s/%s", output_dir_name, FILE_SUMMARY);


//    fprintf(stderr, "DEBUG [write_summary_log]: Attempting to open file: '%s'\n", full_path);
    fprintf(stderr, "Attempting to open file: '%s'\n", full_path);

    // Open the file using the global info_current_file pointer
    info_current_file = fopen(full_path, "w");

    if (info_current_file == NULL) {
        fprintf(stderr, "ERROR [write_summary_log]: Could not open file '%s'.\n", full_path);
        perror("Reason"); // Prints system error message
        // Don't exit here, as specified, just warn and return
        return;
    }

    // --- Write Simulation Summary to File (English and Formatted) ---
    fprintf(info_current_file, "--- Simulation Configuration Summary ---\n");
    fprintf(info_current_file, "Timestamp: %s %s\n", __DATE__, __TIME__); // Add compilation date/time
    fprintf(info_current_file, "Output files for this run are located in: %s\n\n", output_dir_name);

    fprintf(info_current_file, "--- Disk Parameters ---\n");
    fprintf(info_current_file, "   Inner Radius (RMIN): %.4g AU\n", disk_params->RMIN);
    fprintf(info_current_file, "   Outer Radius (RMAX): %.4g AU\n", disk_params->RMAX);
    fprintf(info_current_file, "   Surface Density at R=1 AU (SIGMA0): %.4e\n", disk_params->SIGMA0);
    fprintf(info_current_file, "   Surface Density Exponent (SIGMA_EXP): %.4g\n", disk_params->SIGMAP_EXP); // Note: SIGMAP_EXP for disk_params
    fprintf(info_current_file, "   Flaring Index: %.4g\n", disk_params->FLIND); // Note: FLIND for disk_params
    fprintf(info_current_file, "   Alpha Viscosity (ALPHA_VISC): %.4g\n", disk_params->alpha_visc);
    fprintf(info_current_file, "   Alpha Modified (ALPHA_MOD): %.4g\n", disk_params->a_mod);
    fprintf(info_current_file, "\n");

    fprintf(info_current_file, "--- Dead Zone Edge (DZE) Parameters ---\n");
    fprintf(info_current_file, "   Inner DZE Radius (R_DZE_I): %.4g AU\n", disk_params->r_dze_i);
    fprintf(info_current_file, "   Outer DZE Radius (R_DZE_O): %.4g AU\n", disk_params->r_dze_o);
    fprintf(info_current_file, "   Inner DZE Width (DR_DZEI): %.4g AU\n", disk_params->Dr_dze_i);
    fprintf(info_current_file, "   Outer DZE Width (DR_DZE_O): %.4g AU\n\n", disk_params->Dr_dze_o);
    fprintf(info_current_file, "   **Note: R_DZE_I/O = 0, or R_DZE_I/O = 0 indicates that the corresponding DZE is not simulated.**\n");
    fprintf(info_current_file, "\n\n");

    fprintf(info_current_file, "--- Central Star Parameters ---\n");
    fprintf(info_current_file, "   Central Star Mass: %.4g Solar Masses\n", disk_params->STAR_MASS);
    fprintf(info_current_file, "\n");

    // Close the file
    fclose(info_current_file);
}


// Writes sigma, pressure, and dpress to a file.
void write_gas_profile_to_file(const disk_t *disk_params, output_files_t *output_files) {

    int i;

    if (output_files->surface_file == NULL) {
        fprintf(stderr, "ERROR: output_files->surface_file is NULL in write_gas_profile_to_file! Cannot write sigma data.\n");
        return;
    }

//%-15.6e %-15.6Lg %-15.6e %-15.6e\n",
    for(i = 1; i <= disk_params->NGRID; i++) { // Using disk_params->NGRID
        fprintf(output_files->surface_file, "%-15.6e %-15.6lg %-15.6e %15.6e\n", disk_params->rvec[i], disk_params->sigmavec[i], disk_params->pressvec[i], disk_params->dpressvec[i]);
    }

    fflush(output_files->surface_file);
}

// Writes dust surface density profile to a file.
void write_dust_profile_to_file(int step,
                       const ParticleData_t *p_data, // The structure has changed
                       const disk_t *disk_params,
                       const simulation_options_t *sim_opts,
                       output_files_t *output_files) {

    int i;
    double interpolated_sigmad_at_r = 0.0;
    double interpolated_sigmadm_at_r = 0.0;

    if (output_files->dust_file == NULL) {
        fprintf(stderr, "ERROR: output_files->dust_file is NULL in write_dust_profile_to_file! Cannot write main dust surface density.\n");
        return;
    }

//    fprintf(stderr, "DEBUG [write_dust_profile_to_file]: Entering write_dust_profile_to_file for step %d.\n", step);

    // Main dust population (particles_pop1)
    if (p_data->particles_pop1 == NULL || p_data->num_particles_pop1 <= 0) {
        fprintf(stderr, "WARNING [write_dust_profile_to_file]: No main dust particles (pop1) to print or p_data->particles_pop1 is NULL.\n");
    } else {
        for(i = 0; i < p_data->num_particles_pop1; i++){
            if (p_data->particles_pop1[i].distance_au >= disk_params->RMIN) {
                // Interpolate from disk_params->sigmadustvec
                // Assumed that disk_params->rvec contains the grid points.
                // If the interpol function expects other parameters, modify this!
                interpol(disk_params->sigmadustvec, disk_params->rvec,
                          p_data->particles_pop1[i].distance_au, &interpolated_sigmad_at_r,
                          disk_params->DD, 0, disk_params); // Check the last 2 parameters!

                fprintf(output_files->dust_file, "%-8d %-15d %-15.8lg %-20.15lg %-20.15lg %-15.15lg \n",
                          step,
                          p_data->particles_pop1[i].id,          // Use the ID
                          p_data->particles_pop1[i].distance_au ,
                          p_data->particles_pop1[i].current_size_au * AU_TO_CM,
                          interpolated_sigmad_at_r,
                          p_data->particles_pop1[i].initial_mass_msun); // initial_mass_msun instead of mass_g
            }
        }
        fflush(output_files->dust_file);
    }

    // Micron dust population (particles_pop2), if enabled
    if(sim_opts->twopop == 1.0) {
        if (output_files->micron_dust_file == NULL) {
            fprintf(stderr, "ERROR: output_files->micron_dust_file is NULL in write_dust_profile_to_file (two-pop enabled)! Cannot write micron dust surface density.\n");
            return;
        }

        if (p_data->particles_pop2 == NULL || p_data->num_particles_pop2 <= 0) {
            fprintf(stderr, "WARNING [write_dust_profile_to_file]: No micron dust particles (pop2) to print or p_data->particles_pop2 is NULL.\n");
        } else {
            for(i = 0; i < p_data->num_particles_pop2; i++){
                if (p_data->particles_pop2[i].distance_au >= disk_params->RMIN) {
                    interpol(disk_params->sigmadustmicrvec, disk_params->rvec,
                              p_data->particles_pop2[i].distance_au, &interpolated_sigmadm_at_r,
                              disk_params->DD, 0, disk_params); // Check the last 2 parameters!

                    fprintf(output_files->micron_dust_file,"%-8d %-15d %-15.8lg %-20.15lg %-20.15lg %-15.15lg \n",
                              step,
                              p_data->particles_pop2[i].id,
                              p_data->particles_pop2[i].distance_au,
                              p_data->particles_pop2[i].current_size_au,
                              interpolated_sigmadm_at_r,
                              p_data->particles_pop2[i].initial_mass_msun);
                }
            }
            fflush(output_files->micron_dust_file);
        }
    }
//    fprintf(stderr, "DEBUG [write_dust_profile_to_file]: Exiting write_dust_profile_to_file.\n");
}



// Writes file headers.
void write_file_header(FILE *file, FileType_e file_type, const HeaderData_t *header_data) {
    if (file == NULL) {
        fprintf(stderr, "ERROR [write_file_header]: Attempted to write header to a NULL file pointer!\n");
        return;
    }

    fprintf(file, "# Generated by Dust Drift Simulation (Date: %s %s)\n", __DATE__, __TIME__);

    switch (file_type) {

        case FILE_TYPE_MASS_ACCUMULATION:
            fprintf(file, "# This file contains the time evolution of dust mass within specified disk regions.\n");
            fprintf(file, "# Columns:\n");
            fprintf(file, "# 1. Time [years]\n");
            fprintf(file, "# 2. Inner Boundary Radius (R_in) of the inner region [AU] (%.2f)\n", header_data ? header_data->R_in : 0.0);
            fprintf(file, "# 3. Total Dust Mass (non-micron + micron) within R < R_in [M_Jupiter]\n");
            fprintf(file, "# 4. Outer Boundary Radius (R_out) of the outer region [AU] (%.2f)\n", header_data ? header_data->R_out : 0.0);
            fprintf(file, "# 5. Total Dust Mass (non-micron + micron) within R > R_out [M_Jupiter]\n");
            fprintf(file, "# All masses include both 'cm-sized' (larger) and 'micron-sized' dust populations.\n");
            break;

        case FILE_TYPE_GAS_DENSITY:
            // Additional info if t=0 (is_initial_data)
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial gas profile\n");
            } else {
                fprintf(file, "# Time: %e years\n", header_data ? header_data->current_time : 0.0);
            }
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-15s %-15s %-15s %-15s\n",
                              "Radius_AU", "GasSurfDensity", "GasPressure", "GasPressureDeriv");
            fprintf(file, "#--------------------------------------------------------------------------\n");

            break;

        case FILE_TYPE_DUST_DENSITY:
            fprintf(file, "# Main Dust Surface Density Profile\n");
            fprintf(file, "# Time: %e years\n", header_data ? header_data->current_time : 0.0);
            fprintf(file, "# Columns: 1. Radius [AU], 2. Sigma_dust [M_Sun/AU^2]\n");
            break;

        case FILE_TYPE_DUST_EVOL:
            fprintf(file, "# Evolution of Dust Particles\n");
            fprintf(file, "# Time: %e years\n", header_data ? header_data->current_time : 0.0);
            fprintf(file, "#--------------------------------------------------------------------------------------------\n");
            fprintf(file, "# %-5s %-15s %-20s %-20s %-20s %-15s\n",
                              "Time", "ID", "Distance_AU", "PartSize_cm", "DustSurfDensity","RepMass");
            fprintf(file, "#--------------------------------------------------------------------------------------------\n");
            break;

        case FILE_TYPE_MICRON_DUST_EVOL:
            fprintf(file, "# Evolution of Micron Sized Dust Particles\n");
            fprintf(file, "# Time: %e years\n", header_data ? header_data->current_time : 0.0);
            fprintf(file, "#--------------------------------------------------------------------------------------------\n");
            fprintf(file, "# %-5s %-15s %-20s %-20s %-20s %-15s\n",
                              "Time", "ID", "Distance_AU", "PartSize_cm", "DustSurfDensity","RepMass");
            fprintf(file, "#--------------------------------------------------------------------------------------------\n");
            break;

        case FILE_TYPE_PARTICLE_SIZE:
            // Additional info if t=0 (is_initial_data)
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
            // The actual parameter values are not written in the header, but in the main data line.
            break;

        case FILE_TYPE_TIMESCALE:
            fprintf(file, "# Dust Depletion Timescale Profile\n");
            fprintf(file, "# Note: Timescale calculated based on a unit system where 2*PI simulation time units = 1 year.\n");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            fprintf(file, "#   %-15s %-15s\n","R_dust [AU]", "Timesc_depletion [years]");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            break;

        default:
            fprintf(stderr, "WARNING [write_file_header]: Unknown file type for header generation: %d!\n", file_type);
            break;
    }
    fflush(file);
}




// This function is currently **inactive** and performs no actions.
// It is intended to initialize the mass accumulation file.
int initialize_mass_accumulation_file(output_files_t *output_files, const simulation_options_t *sim_opts,
                                     const disk_t *disk_params, HeaderData_t *header_data_for_files) {
    char massout[MAX_PATH_LEN] = "";

    // Prepare the header_data_for_files structure with specific data
    header_data_for_files->current_time = 0.0;
    header_data_for_files->is_initial_data = 1;
    header_data_for_files->R_in = disk_params->RMIN;
    header_data_for_files->R_out = disk_params->RMAX;

    // In the write_summary_log function, the full_path was created with a combination of output_dir_name and LOGS_DIR.
    // Here, I assume that the mass_file is also under LOGS_DIR, and massout is a combination of the main output folder and LOGS_DIR.
    // Correction for the massout path to be consistent with write_summary_log.
    // sim_opts->output_dir_name is the root output folder.
    // LOGS_DIR can be a subfolder within it (as used by write_summary_log).
    // If the mass_file is directly under sim_opts->output_dir_name, then
    // snprintf(massout, MAX_PATH_LEN, "%s/%s.dat", sim_opts->output_dir_name, FILE_MASS_ACCUMULATE);
    // If it's under the LOGS_DIR subfolder, then the current snprintf is correct.
/* snprintf(massout, MAX_PATH_LEN, "%s/%s/%s.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_MASS_ACCUMULATE);


//     fprintf(stderr, "DEBUG [initialize_mass_accumulation_file]: Opening output file: %s\n", massout);


    output_files->mass_file = fopen(massout, "w");
    if (output_files->mass_file == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", massout);
        return 1; // Error
    }
    write_file_header(output_files->mass_file, FILE_TYPE_MASS_ACCUMULATION, header_data_for_files);
*/
    return 0; // Success
}



void cleanup_simulation_resources(ParticleData_t *p_data, output_files_t *output_files, const simulation_options_t *sim_opts) {
    // Use the dedicated free function for particle data
    free_particle_data(p_data);

    // Note: The `mass_file` cleanup is currently commented out and inactive.
/* if (output_files->mass_file != NULL) {
        fclose(output_files->mass_file);
        output_files->mass_file = NULL;
        fprintf(stderr, "DEBUG [cleanup_simulation_resources]: Closed %s\n", FILE_MASS_ACCUMULATE);
    }
*/
}

// Helper function to close snapshot files.
void close_snapshot_files(output_files_t *output_files, const char *dens_name, const char *dust_name, const char *dust_name2, const simulation_options_t *sim_opts) {
    if (output_files->surface_file != NULL) {
        fclose(output_files->surface_file);
        output_files->surface_file = NULL;
    }
    if (output_files->dust_file != NULL) {
        fclose(output_files->dust_file);
        output_files->dust_file = NULL;
    }
    if (sim_opts->twopop == 1 && output_files->micron_dust_file != NULL) {
        fclose(output_files->micron_dust_file);
        output_files->micron_dust_file = NULL;
    }
}