#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>     // For errno

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

// Local includes
#include "config.h"          // Biztosítja, hogy az info_current_file látható legyen
#include "io_utils.h"
#include "dust_physics.h"    // If needed for any specific function interactions
#include "simulation_types.h" // For disk_t, simulation_options_t, output_files_t
#include "globals.h"
#include "particle_data.h"   // A ParticleData_t definíciójához és a dust_particle_t típushoz
#include "dust_particle.h"   // Explicit include for dust_particle_t definition

#include <sys/stat.h> // A mkdir-hez (ha ezt a fájl használja)
#include <sys/types.h> // A mkdir-hez (ha ezt a fájl használja)


#define INIT_DATA_HEADER_LINES 5
// --- GLOBAL FILE POINTERS ---
// DO NOT DEFINE THEM HERE IF THEY ARE DEFINED ELSEWHERE (e.g., in main.c or config.c)
// They are simply used here because they are declared as 'extern' in config.h.

// --- FÜGGVÉNY DEFINÍCIÓK ---

/* Visszaadja, hogy hány sora van a beolvasandó fájlnak,
 * ez jelen esetben megadja a beolvasandó részecskék számát. */
int reszecskek_szama(const char *filenev) {
    FILE *fp = NULL;
    char line_buffer[1024];
    int line_count = 0; // Counter for data lines

    fp = fopen(filenev, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR [reszecskek_szama]: Could not open file '%s'.\n", filenev);
        perror("Reason"); // Prints system error message
        exit(EXIT_FAILURE);
    }

    // Skip header lines
    for (int i = 0; i < INIT_DATA_HEADER_LINES; i++) {
        if (fgets(line_buffer, sizeof(line_buffer), fp) == NULL) {
            // If file ends before all header lines are skipped, it's an error
            fprintf(stderr, "ERROR [reszecskek_szama]: Unexpected end of file while skipping %d header lines in '%s'.\n", INIT_DATA_HEADER_LINES, filenev);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
    }

    // Count remaining data lines
    while (fgets(line_buffer, sizeof(line_buffer), fp) != NULL) {
        // You might want to add a check here to ensure the line is not empty or a comment
        // For example, if lines starting with '#' are comments:
        if (line_buffer[0] != '#' && line_buffer[0] != '\n' && line_buffer[0] != '\r') {
             line_count++;
        }
    }

    fclose(fp); // Close the file after reading
    return line_count;
}

/* A porreszecskék adatainak beolvasása (ÚJ verzió, ParticleData_t struktúrával, PARTICLE_NUMBER-rel) */
void ReadDustFile_V2(ParticleData_t *p_data, const char *filename,
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

    // DEBUG: Entering function
//    fprintf(stderr, "DEBUG [ReadDustFile_V2]: Entering function. Reading from '%s'. sim_opts->twopop = %.2f\n", filename, sim_opts->twopop);

    // 1. Determine number of particles (total lines in the file)
    // IMPORTANT: 'reszecskek_szama' is assumed to return the count for ONE population.
    // The global PARTICLE_NUMBER will be set by this.
    PARTICLE_NUMBER = reszecskek_szama(filename);
    
    if (PARTICLE_NUMBER <= 0) {
        fprintf(stderr, "ERROR [ReadDustFile_V2]: No particles found or error in counting for file '%s'. PARTICLE_NUMBER = %d.\n", filename, PARTICLE_NUMBER);
        exit(EXIT_FAILURE);
    }

    // Based on PARTICLE_NUMBER, set the particle counts in p_data.
    // If twopop is enabled, both populations get PARTICLE_NUMBER particles.
    p_data->num_particles_pop1 = PARTICLE_NUMBER;
    p_data->num_particles_pop2 = (sim_opts->twopop == 1.0) ? PARTICLE_NUMBER : 0; 

    // Update sim_opts->num_dust_particles. This usually refers to the main loop count, which is PARTICLE_NUMBER.
    sim_opts->num_dust_particles = PARTICLE_NUMBER;
//    fprintf(stderr, "DEBUG [ReadDustFile_V2]: Global PARTICLE_NUMBER = %d. Setting p_data->num_particles_pop1 = %d, p_data->num_particles_pop2 = %d.\n",
//            PARTICLE_NUMBER, p_data->num_particles_pop1, p_data->num_particles_pop2);
//    fprintf(stderr, "DEBUG [ReadDustFile_V2]: sim_opts->num_dust_particles set to %d.\n", sim_opts->num_dust_particles);


    // 2. Allocate memory using the dedicated allocation function
    // Pass the twopop flag (as int) from sim_opts to allocate_particle_data
//    allocate_particle_data(p_data, p_data->num_particles_pop1, p_data->num_particles_pop2, (int)(sim_opts->twopop + 0.5)); 

    // DEBUG: Memory allocation check
    fprintf(stderr, "DEBUG [ReadDustFile_V2]: Particle data arrays allocated. Pop1 allocated: %s, Pop2 allocated: %s (if twopop enabled).\n",
            (p_data->particles_pop1 != NULL ? "YES" : "NO"), (p_data->particles_pop2 != NULL ? "YES" : "NO"));


    // 3. Open the file for reading actual data
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR [ReadDustFile_V2]: Could not open file '%s' for data reading (second pass).\n", filename);
        perror("Reason"); // Print system error message
        // Free previously allocated memory before exiting on error
        free_particle_data(p_data); 
        exit(EXIT_FAILURE);
    }

    // Skip header lines
    char line_buffer[1024];
    for (int k = 0; k < INIT_DATA_HEADER_LINES; k++) {
        if (fgets(line_buffer, sizeof(line_buffer), fp) == NULL) {
            fprintf(stderr, "ERROR [ReadDustFile_V2]: Unexpected end of file while skipping headers in '%s'.\n", filename);
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
    fprintf(stderr, "DEBUG [ReadDustFile_V2]: Successfully read %d particle entries from '%s'.\n", PARTICLE_NUMBER, filename);
}


void ReadSigmaFile(disk_t *disk_params, const char *filename) {
    const char *input_filename = filename;

    FILE *fp = fopen(input_filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR [ReadSigmaFile]: Could not open input file '%s'.\n", input_filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }

    char line[512];

    // Fejléc sorok átugrása
    // Addig olvassuk a sorokat, amíg komment (#) karakterrel kezdődnek
    // VAGY amíg nem a '---' elválasztó sort találjuk, ami a valódi adatok kezdetét jelzi
    while (fgets(line, sizeof(line), fp) != NULL) {
        // Ellenőrizzük, hogy a sor komment-e vagy a '---' elválasztó
        if (line[0] == '#' || strncmp(line, "---", 3) == 0) { // strncmp "---" karakterekre
            continue; // Ugrás a következő sorra
        } else {
            // Ez az első valós adatsor. Visszaállítjuk a fájlmutatót a sor elejére.
            fseek(fp, -strlen(line), SEEK_CUR);
            break;
        }
    }

    // Ha a fájl üres vagy csak kommenteket tartalmaz
    if (feof(fp) && (line[0] == '#' || strncmp(line, "---", 3) == 0)) {
        fprintf(stderr, "ERROR [ReadSigmaFile]: File '%s' is empty or only contains comments/headers.\n", input_filename);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    double r_val;
    double sigma_gas_val;
    double pressure_gas_val;
    double dpressure_dr_val;

    // A ciklus disk_params->NGRID-ig megy
    for (int i = 0; i < disk_params->NGRID; i++) {
        // MOST MÁR PONTOSAN A FÁJLOD FORMÁTUMÁT OLVASSUK BE:
        // Radius_AU, GasSurfDensity, GasPressure, GasPressureDeriv (4 oszlop, mind double)
        if (fscanf(fp, "%lf %lf %lf %lf",
                            &r_val, &sigma_gas_val, &pressure_gas_val, &dpressure_dr_val) != 4) {
            // Hiba kezelése, ha nem tudunk 4 double értéket beolvasni
            fprintf(stderr, "ERROR [ReadSigmaFile]: Failed to read 4 values for row %d from file '%s'. File may be malformed or ended unexpectedly.\n", i, input_filename);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        // Hozzárendelés a disk_params tömbökhöz
        // Az indexelés 'i + 1' a 0-ás indexű szellemcella miatt (ahogy a disk_t definíciója és a calculate_boundary függvény valószínűsíti).
        // Fontos: ellenőrizzük, hogy az 'i + 1' index a tömb határain belül van-e.
        // A tömbök mérete disk_params->NGRID + 2, tehát az érvényes indexek 0-tól NGRID+1-ig mennek.
        // A "valós" adatok 1-től NGRID-ig kerülnek, a 0 és NGRID+1 pedig a calculate_boundary-hez.
        if ((i + 1) >= 0 && (i + 1) <= disk_params->NGRID + 1) {
            disk_params->rvec[i + 1] = r_val;
            disk_params->sigmavec[i + 1] = sigma_gas_val;
            disk_params->pressvec[i + 1] = pressure_gas_val;
            disk_params->dpressvec[i + 1] = dpressure_dr_val;
        } else {
            fprintf(stderr, "WARNING [ReadSigmaFile]: Attempted to write to out-of-bounds index %d. Max allowed index: %d (NGRID+1).\n", i + 1, disk_params->NGRID + 1);
        }

    }

    fclose(fp);

}


void Mk_Dir(char *dir_path) {
    char tmp_path[MAX_PATH_LEN];
    int counter = 0;

    // Másolat az eredeti névről
    strncpy(tmp_path, dir_path, MAX_PATH_LEN - 1);
    tmp_path[MAX_PATH_LEN - 1] = '\0';

    while (access(tmp_path, F_OK) == 0) {   // Fájl vagy mappa létezik
        snprintf(tmp_path, MAX_PATH_LEN, "%s_%04d", dir_path, ++counter);
        if (counter > 99) {
            fprintf(stderr, "ERROR [Mk_Dir]: Too many existing directories with similar names.\n");
            exit(1);
        }
    }

    int result = MKDIR_CALL(tmp_path);
    if (result != 0) {
        perror("ERROR [Mk_Dir]: mkdir failed");
        fprintf(stderr, "ERROR [Mk_Dir]: Could not create directory: '%s'\n", tmp_path);
        exit(1);
    }

    fprintf(stderr, "DEBUG [Mk_Dir]: Directory '%s' created successfully.\n", tmp_path);

    // Másold vissza a létrehozott mappa nevét a bemenetbe
    strncpy(dir_path, tmp_path, MAX_PATH_LEN - 1);
    dir_path[MAX_PATH_LEN - 1] = '\0';

    fflush(stderr);
}

/**
 * @brief Generates a summary log file for the current simulation run.
 * @details This function creates a file (named `summary.dat`)
 * within the specified output directory. It logs the initial simulation parameters,
 * including disk properties, central star mass, and Dead Zone Edge (DZE) configurations.
 *
 * @param output_dir_name [in] The name of the main output directory for this simulation run.
 * This directory is typically located under the `LOGS_DIR`.
 * @param disk_params [in] A pointer to the `disk_t` structure containing disk-specific parameters.
 * @param sim_opts [in] A pointer to the `simulation_options_t` structure, containing
 * general simulation options.
 *
 * @note This function uses the global `info_current_file` file pointer, which is declared in `config.h`
 * and defined in `config.c`.
 */
void infoCurrent(const char *output_dir_name, const disk_t *disk_params, const simulation_options_t *sim_opts) {

    char full_path[MAX_PATH_LEN];
    // char file_name[100]; // Unused variable, commented out


    snprintf(full_path, sizeof(full_path), "%s/%s", output_dir_name, FILE_SUMMARY);


    fprintf(stderr, "DEBUG [infoCurrent]: Attempting to open file: '%s'\n", full_path);

    // Open the file using the global info_current_file pointer
    info_current_file = fopen(full_path, "w");

    if (info_current_file == NULL) {
        fprintf(stderr, "ERROR [infoCurrent]: Could not open file '%s'.\n", full_path);
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


/* Függvény a sigma, p, dp kiíratására */
void Print_Sigma(const disk_t *disk_params, output_files_t *output_files) {

    int i;

    if (output_files->surface_file == NULL) {
        fprintf(stderr, "ERROR: output_files->surface_file is NULL in Print_Sigma! Cannot write sigma data.\n");
        return;
    }

//%-15.6e %-15.6Lg %-15.6e %-15.6e\n",
    for(i = 1; i <= disk_params->NGRID; i++) { // Using disk_params->NGRID
        fprintf(output_files->surface_file, "%-15.6e %-15.6lg %-15.6e %15.6e\n", disk_params->rvec[i], disk_params->sigmavec[i], disk_params->pressvec[i], disk_params->dpressvec[i]);
    }

    fflush(output_files->surface_file);
}

/* Függvény a por felületisűrűségének kiíratására */
/* Function to write dust surface density */
void Print_Sigmad(int step,
                  const ParticleData_t *p_data, // Itt is változik
                  const disk_t *disk_params,
                  const simulation_options_t *sim_opts,
                  output_files_t *output_files) {

    int i;
    double interpolated_sigmad_at_r = 0.0;
    double interpolated_sigmadm_at_r = 0.0;

    if (output_files->dust_file == NULL) {
        fprintf(stderr, "ERROR: output_files->dust_file is NULL in Print_Sigmad! Cannot write main dust surface density.\n");
        return;
    }

    fprintf(stderr, "DEBUG [Print_Sigmad]: Entering Print_Sigmad for step %d.\n", step);

    // Fő porpopuláció (particles_pop1)
    if (p_data->particles_pop1 == NULL || p_data->num_particles_pop1 <= 0) {
        fprintf(stderr, "WARNING [Print_Sigmad]: No main dust particles (pop1) to print or p_data->particles_pop1 is NULL.\n");
    } else {
        for(i = 0; i < p_data->num_particles_pop1; i++){
            if (p_data->particles_pop1[i].distance_au >= disk_params->RMIN) {
                // Interpolálás a disk_params->sigmadustvec-ből
                // Feltételezve, hogy disk_params->rvec is tartalmazza a rács pontjait.
                // A 'DD' paramétert a disk_t struktúrából veszem, ha az a rács pontok száma.
                // Ha az interpol függvény nem ezt várja, módosítsd!
                interpol(disk_params->sigmadustvec, disk_params->rvec,
                         p_data->particles_pop1[i].distance_au, &interpolated_sigmad_at_r,
                         disk_params->DD, 0, disk_params); // utolsó 2 paramétert ellenőrizd!

                fprintf(output_files->dust_file, "%-8d %-15d %-15.8lg %-20.15lg %-20.15lg %-15.15lg \n",
                        step,
                        p_data->particles_pop1[i].id,               // Használjuk az ID-t
                        p_data->particles_pop1[i].distance_au ,
                        p_data->particles_pop1[i].current_size_au * AU_TO_CM,
                        interpolated_sigmad_at_r,
                        p_data->particles_pop1[i].initial_mass_msun); // mass_g helyett initial_mass_msun
            }
        }
        fflush(output_files->dust_file);
    }

    // Mikron porpopuláció (particles_pop2), ha engedélyezve van
    if(sim_opts->twopop == 1.0) {
        if (output_files->micron_dust_file == NULL) {
            fprintf(stderr, "ERROR: output_files->micron_dust_file is NULL in Print_Sigmad (two-pop enabled)! Cannot write micron dust surface density.\n");
            return;
        }

        if (p_data->particles_pop2 == NULL || p_data->num_particles_pop2 <= 0) {
            fprintf(stderr, "WARNING [Print_Sigmad]: No micron dust particles (pop2) to print or p_data->particles_pop2 is NULL.\n");
        } else {
            for(i = 0; i < p_data->num_particles_pop2; i++){
                if (p_data->particles_pop2[i].distance_au >= disk_params->RMIN) {
// ide majd a sigmamicrdustvec kellene inkáb!!!
                    interpol(disk_params->sigmadustvec, disk_params->rvec,
                             p_data->particles_pop2[i].distance_au, &interpolated_sigmadm_at_r,
                             disk_params->DD, 0, disk_params); // utolsó 2 paramétert ellenőrizd!

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
    fprintf(stderr, "DEBUG [Print_Sigmad]: Exiting Print_Sigmad.\n");
}


/* Az időt tartalmazó fájl paramétereinek beolvasása (vagy beállítása) */
void timePar(double tMax_val, double stepping_val, double current_val, simulation_options_t *sim_opts) {


    sim_opts->TMAX = tMax_val;
    sim_opts->WO = tMax_val / stepping_val;
    sim_opts->TCURR = current_val;
}


// Függvény a fájl fejlécek kiírására
void print_file_header(FILE *file, FileType_e file_type, const HeaderData_t *header_data) {
    if (file == NULL) {
        fprintf(stderr, "ERROR [print_file_header]: Attempted to write header to a NULL file pointer!\n");
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
            // Kiegészítő infók, ha t=0 (is_initial_data)
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial gas profile\n"); // Removed %s, was trying to print current_time here
            } else {
                fprintf(file, "# Time: %e years\n", header_data ? header_data->current_time : 0.0);
            }
            // A fejléc az init_tool_module-ból, módosítva a HeaderData_t használatára
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
            // A fejléc az init_tool_module-ból, módosítva a HeaderData_t használatára
            // Kiegészítő infók, ha t=0 (is_initial_data)
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial particle distribution\n"); // Removed %s
            } else {
                fprintf(file, "# Particle distribution, Time: %e years\n", header_data ? header_data->current_time : 0.0);
            }
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-5s %-15s %-20s %-20s %-15s %-15s\n",
                            "Index", "Radius_AU", "RepMass_Pop1_Msun", "RepMass_Pop2_Msun", "MaxPartSize_cm", "MicroSize_cm");
            fprintf(file, "#--------------------------------------------------------------------------\n");

            break;

        case FILE_TYPE_DISK_PARAM:
            // A fejléc az init_tool_module-ból, módosítva a HeaderData_t használatára
            fprintf(file, "# Disk Parameters\n");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            fprintf(file, "# %-15s %-15s %-10s %-15s %-20s %-15s %-15s %-15s %-20s %-20s %-15s %-15s %-15s %-15s %-15s\n",
                            "R_Min_AU", "R_Max_AU", "N_Grid", "SigmaExp", "Sigma0_gas_Msun_AU2",
                            "G_GravConst", "DzR_Inner_AU", "DzR_Outer_AU", "DzDr_Inner_Calc_AU", "DzDr_Outer_Calc_AU",
                            "DzAlphaMod", "DustDensity_g_cm3", "AlphaViscosity", "StarMass_Msun", "FlaringIndex");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            // A tényleges paraméter értékeket nem a fejlécbe írjuk, hanem a fő adatsorba.
            break;

        case FILE_TYPE_TIMESCALE:
            fprintf(file, "# Dust Depletion Timescale Profile\n");
            fprintf(file, "# Note: Timescale calculated based on a unit system where 2*PI simulation time units = 1 year.\n");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            fprintf(file, "#   %-15s %-15s\n","R_dust [AU]", "Timesc_depletion [years]");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            break;

        default:
            fprintf(stderr, "WARNING [print_file_header]: Unknown file type for header generation: %d!\n", file_type);
            break;
    }
    fflush(file);
}




int setup_initial_output_files(output_files_t *output_files, const simulation_options_t *sim_opts,
                                 const disk_t *disk_params, HeaderData_t *header_data_for_files) {
    char massout[MAX_PATH_LEN] = "";

    // Készítsük elő a header_data_for_files struktúrát a specifikus adatokkal
    header_data_for_files->current_time = 0.0;
    header_data_for_files->is_initial_data = 1;
    header_data_for_files->R_in = disk_params->RMIN;
    header_data_for_files->R_out = disk_params->RMAX;

    // Az infoCurrent függvényben a full_path-ot az output_dir_name és LOGS_DIR kombinációjával hoztad létre.
    // Itt feltételezem, hogy a mass_file is a LOGS_DIR alatt van, és a massout a fő output mappa és a LOGS_DIR kombinációja.
    // Javítás a massout elérési útvonalára, hogy konzisztens legyen a infoCurrent-tel.
    // A sim_opts->output_dir_name a gyökér kimeneti mappa.
    // A LOGS_DIR pedig egy almappa lehet ezen belül (ahogy az infoCurrent is használja).
    // Ha a mass_file közvetlenül a sim_opts->output_dir_name alatt van, akkor
    // snprintf(massout, MAX_PATH_LEN, "%s/%s.dat", sim_opts->output_dir_name, FILE_MASS_ACCUMULATE);
    // Ha a LOGS_DIR almappa alatt van, akkor a mostani snprintf helyes.
/*    snprintf(massout, MAX_PATH_LEN, "%s/%s/%s.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_MASS_ACCUMULATE);


//    fprintf(stderr, "DEBUG [setup_initial_output_files]: Opening output file: %s\n", massout);


    output_files->mass_file = fopen(massout, "w");
    if (output_files->mass_file == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", massout);
        return 1; // Hiba
    }
    print_file_header(output_files->mass_file, FILE_TYPE_MASS_ACCUMULATION, header_data_for_files);
*/
    return 0; // Siker
}




void cleanup_simulation_resources(ParticleData_t *p_data, output_files_t *output_files, const simulation_options_t *sim_opts) {
    // Use the dedicated free function for particle data
    free_particle_data(p_data);

/*    if (output_files->mass_file != NULL) {
        fclose(output_files->mass_file);
        output_files->mass_file = NULL;
        fprintf(stderr, "DEBUG [cleanup_simulation_resources]: Closed %s\n", FILE_MASS_ACCUMULATE);
    }
*/    
}

// Segédfüggvény a pillanatfelvételek fájljainak bezárására
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