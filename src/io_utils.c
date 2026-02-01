#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>    // For errno


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
#include "io_utils.h"
#include "config.h"         
#include "dust_physics.h"   // If needed for any specific function interactions
#include "utils.h"          // For countZeroPoints, findZeroPoint, findRAnnulusAroundDZE
#include "simulation_types.h" // For DiskParameters, SimulationOptions, OutputFiles
#include "boundary_conditions.h"


#define INIT_DATA_HEADER_LINES 5
// --- GLOBAL FILE POINTERS ---
// DO NOT DEFINE THEM HERE IF THEY ARE DEFINED ELSEWHERE (e.g., in main.c or config.c)
// They are simply used here because they are declared as 'extern' in config.h.

// --- FÜGGVÉNY DEFINÍCIÓK ---

/* Visszaadja, hogy hány sora van a beolvasandó fájlnak,
 * ez jelen esetben megadja a beolvasandó részecskék számát. */
int calculateNumbersOfParticles(const char *filenev) {
    FILE *particle_file = NULL;
    char line_buffer[1024];
    int line_count = 0; // Counter for data lines

    particle_file = fopen(filenev, "r");
    if (particle_file == NULL) {
        fprintf(stderr, "ERROR [calculateNumbersOfParticles]: Could not open file '%s'.\n", filenev);
        perror("Reason"); // Prints system error message
        exit(EXIT_FAILURE);
    }

    // Skip header lines
    for (int i = 0; i < INIT_DATA_HEADER_LINES; i++) {
        if (fgets(line_buffer, sizeof(line_buffer), particle_file) == NULL) {
            // If file ends before all header lines are skipped, it's an error
            fprintf(stderr, "ERROR [calculateNumbersOfParticles]: Unexpected end of file while skipping %d header lines in '%s'.\n", INIT_DATA_HEADER_LINES, filenev);
            fclose(particle_file);
            exit(EXIT_FAILURE);
        }
    }

    // Count remaining data lines
    while (fgets(line_buffer, sizeof(line_buffer), particle_file) != NULL) {
        // You might want to add a check here to ensure the line is not empty or a comment
        // For example, if lines starting with '#' are comments:
        if (line_buffer[0] != '#' && line_buffer[0] != '\n' && line_buffer[0] != '\r') {
             line_count++;
        }
    }

    fclose(particle_file); // Close the file after reading
    return line_count;
}

/* A porreszecskek adatainak beolvasasa */
void loadDustParticlesFromFile(ParticleData *particle_data, const char *filename) {

    int i, dummy;
    double distance, particle_radius, radmicr;
    long double reprmass;
    long double reprmassmicr;

    // Use the global 'load_dust_particles_file' which is declared extern in config.h
    load_dust_particles_file = fopen(filename,"r"); // Use the passed filename

    if (load_dust_particles_file == NULL) {
        fprintf(stderr, "ERROR [loadDustParticlesFromFile]: Could not open file '%s'.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }


    char line_buffer[1024];
    for (int k = 0; k < INIT_DATA_HEADER_LINES; k++) {
        if (fgets(line_buffer, sizeof(line_buffer), load_dust_particles_file) == NULL) {
            fprintf(stderr, "ERROR [loadDustParticlesFromFile]: Unexpected end of file while skipping headers in '%s'.\n", filename);
            fclose(load_dust_particles_file);
            exit(EXIT_FAILURE);
        }
    }


    for (i = 0; i < particle_number; i++) {
        if(fscanf(load_dust_particles_file,"%d %lg %Lg %Lg %lg %lg",&dummy,&distance,&reprmass,&reprmassmicr,&particle_radius,&radmicr) == 6) {
            particle_data->particle_distance_array[i][0] = distance;
            particle_data->particle_distance_array[i][1] = particle_radius / AU_IN_CM; // AU_IN_CM from config.h
            particle_data->particle_mass_array[i] = reprmass;

            particle_data->micron_particle_distance_array[i][0] = distance;
            particle_data->micron_particle_distance_array[i][1] = radmicr / AU_IN_CM; // AU_IN_CM from config.h
            particle_data->massmicradial_grid[i] = reprmassmicr;
        } else {
            fprintf(stderr, "\n\n******************* ERROR!      *********************\n\n");
            fprintf(stderr, "   Failed to read line %d from particle data file '%s'!\n", i, filename);
            fprintf(stderr, "   Expected 6 values, but fscanf failed. Program will exit.\n");
            fclose(load_dust_particles_file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(load_dust_particles_file);
}


void loadGasSurfaceDensityFromFile(DiskParameters *disk_params, const char *filename) {
    const char *input_filename = filename;

    FILE *fp = fopen(input_filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR [loadGasSurfaceDensityFromFile]: Could not open input file '%s'.\n", input_filename);
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
        fprintf(stderr, "ERROR [loadGasSurfaceDensityFromFile]: File '%s' is empty or only contains comments/headers.\n", input_filename);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    double r_val;
    double sigma_gas_val;
    double pressure_gas_val;
    double dpressure_dr_val;

    // A ciklus disk_params->grid_number-ig megy
    for (int i = 0; i < disk_params->grid_number; i++) {
        // MOST MÁR PONTOSAN A FÁJLOD FORMÁTUMÁT OLVASSUK BE:
        // Radius_AU, GasSurfDensity, GasPressure, GasPressureDeriv (4 oszlop, mind double)
        if (fscanf(fp, "%lf %lf %lf %lf",
                         &r_val, &sigma_gas_val, &pressure_gas_val, &dpressure_dr_val) != 4) {
            // Hiba kezelése, ha nem tudunk 4 double értéket beolvasni
            fprintf(stderr, "ERROR [loadGasSurfaceDensityFromFile]: Failed to read 4 values for row %d from file '%s'. File may be malformed or ended unexpectedly.\n", i, input_filename);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        // Hozzárendelés a disk_params tömbökhöz
        // Az indexelés 'i + 1' a 0-ás indexű szellemcella miatt (ahogy a DiskParameters definíciója és a applyBoundaryConditions függvény valószínűsíti).
        // Fontos: ellenőrizzük, hogy az 'i + 1' index a tömb határain belül van-e.
        // A tömbök mérete disk_params->grid_number + 2, tehát az érvényes indexek 0-tól grid_number+1-ig mennek.
        // A "valós" adatok 1-től grid_number-ig kerülnek, a 0 és grid_number+1 pedig a applyBoundaryConditions-hez.
        if ((i + 1) >= 0 && (i + 1) <= disk_params->grid_number + 1) { 
            disk_params->radial_grid[i + 1] = r_val;
            disk_params->gas_surface_density_vector[i + 1] = sigma_gas_val;
            disk_params->gas_pressure_vector[i + 1] = pressure_gas_val;
            disk_params->gas_pressure_gradient_vector[i + 1] = dpressure_dr_val;
        } else {
            fprintf(stderr, "WARNING [loadGasSurfaceDensityFromFile]: Attempted to write to out-of-bounds index %d. Max allowed index: %d (grid_number+1).\n", i + 1, disk_params->grid_number + 1);
        }

    }

    fclose(fp);

}


void createRunDirectory(char *dir_path) {
    char tmp_path[MAX_PATH_LEN];
    int counter = 0;

    // Másolat az eredeti névről
    strncpy(tmp_path, dir_path, MAX_PATH_LEN - 1);
    tmp_path[MAX_PATH_LEN - 1] = '\0';

    while (access(tmp_path, F_OK) == 0) {  // Fájl vagy mappa létezik
        snprintf(tmp_path, MAX_PATH_LEN, "%s_%04d", dir_path, ++counter);
        if (counter > 99) {
            fprintf(stderr, "ERROR [createRunDirectory]: Too many existing directories with similar names.\n");
            exit(1);
        }
    }

    int result = MKDIR_CALL(tmp_path);
    if (result != 0) {
        perror("ERROR [createRunDirectory]: mkdir failed");
        fprintf(stderr, "ERROR [createRunDirectory]: Could not create directory: '%s'\n", tmp_path);
        exit(1);
    }

    fprintf(stderr, "DEBUG [createRunDirectory]: Directory '%s' created successfully.\n", tmp_path);

    // Másold vissza a létrehozott mappa nevét a bemenetbe
    strncpy(dir_path, tmp_path, MAX_PATH_LEN - 1);
    dir_path[MAX_PATH_LEN - 1] = '\0';

    fflush(stderr);
}

/* Elkészít egy fájlt, ami tartalmazza a jelenlegi futás paramétereit,
 * és hogy melyik mappában találhatóak a fájlok */
void printCurrentInformationAboutRun(const char *nev, const DiskParameters *disk_params, const SimulationOptions *sim_opts) {

    char full_path[MAX_PATH_LEN]; // Használjuk a MAX_PATH_LEN-t a biztonságos puffereléshez
    char file_name[100]; 

    sprintf(file_name, "%s%s", kCurrentInfoFile,kFileNamesSuffix);
    
    // Építsük fel a teljes elérési utat: <nev>/<file_name>
    snprintf(full_path, sizeof(full_path), "%s/%s", nev, file_name);

    fprintf(stderr, "DEBUG [printCurrentInformationAboutRun]: Attempting to open file: '%s'\n", full_path);

    current_info_file = fopen(full_path, "w"); // Most már a teljes elérési utat használja

    if (current_info_file == NULL) {
        fprintf(stderr, "ERROR [printCurrentInformationAboutRun]: Could not open file '%s'.\n", full_path);
        perror("Reason");
        // Don't exit here, it's not critical, just warn and return
        return;
    }

    fprintf(current_info_file,"The current run is in the %s directory!\n",nev);
    fprintf(current_info_file,"\n\nThe parameters of the disk:\nr_min: %lg, r_max: %lg\nsigma_0: %lg, SIGMA_EXP: %lg, flaring index: %lg\nALPHA_VISC: %lg, ALPHA_MOD: %lg\nR_DZE_I: %lg, R_DZE_O: %lg, DR_DZEI: %lg, DR_DZE_O: %lg   (*** R_DZE_I/O = 0, akkor azt a DZE-t nem szimulálja a futás! ***)\n\n\n",
              disk_params->r_min, disk_params->r_max,
              disk_params->sigma_0, disk_params->sigma_power_law_index, disk_params->flaring_index,
              disk_params->alpha_parameter, disk_params->alpha_parameter_modification,
              disk_params->r_dze_i, disk_params->r_dze_o, disk_params->dr_dze_i, disk_params->dr_dze_o);
    fprintf(current_info_file,"The mass of the central star: %lg M_Sun\n", disk_params->stellar_mass);
    fclose(current_info_file);
}


void printMassGrowthAtDZEFile(double step, double (*partmassind)[5], double (*partmassmicrind)[5], double t, double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, 
                double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, double *tavin, double *tavout, 
                const DiskParameters *disk_params, const SimulationOptions *sim_opts,OutputFiles *output_files) {


    double ind_ii, ind_io, ind_oi, ind_oo, tav, tav2;

    tav = disk_params->r_dze_o;
    tav2 = disk_params->r_dze_i;

    int dim = countZeroPoints(disk_params); 
    double *r_count = (double *)malloc(sizeof(double) * dim); 

    if (dim > 0 && r_count == NULL) {
        fprintf(stderr, "ERROR [printMassGrowthAtDZEFile]: Failed to allocate memory for r_count. Exiting.\n");
        exit(EXIT_FAILURE);
    }

    double temp_new = 0.;
    double temp = 0.;
    double rin = disk_params->r_dze_i;
    double rout = disk_params->r_dze_o;
    double rin_new = 0.0;
    double rout_new = 0.0;




    int j = 0, i;

    if(dim != 0) {
        for(i = 0; i < disk_params->grid_number; i++) {
            // Itt a findZeroPoint valószínűleg disk_params->radial_grid és disk_params->gas_pressure_gradient_vector-et használ
            temp_new = findZeroPoint(i,disk_params->radial_grid,disk_params->gas_pressure_gradient_vector); 
            if(temp != temp_new && i > 3 && temp_new != 0.0) {
                if (j < dim) { 
                    r_count[j] = temp_new;
                    j++;
                } else {
                    fprintf(stderr, "WARNING [printMassGrowthAtDZEFile]: r_count array overflow, skipping data. dim: %d, j: %d\n", dim, j);
                }
            }
            if(sim_opts->flag_for_deadzone == 0.0) { 
                if(temp_new > 0.) {
                    temp = temp_new;
                    rout_new = temp;
                }
            }
        }
    }

    if(sim_opts->flag_for_deadzone == 1.0) {
        if(dim > 0) {
            if (dim == 1) { 
                rin_new = r_count[0]; 
                rout_new = rout; 
            } else {
                rin_new = r_count[0]; 
                rout_new = r_count[1]; 
            } 
        }
        if(dim == 0) { 
            rin_new = rin; 
            rout_new = rout; 
        }
    }

    rin = rin_new;
    if(sim_opts->flag_for_deadzone == 0.0) rin = 0;
//    double rout_current = rout_new;
    rout = rout_new;
    
//    *tavin = rin;  
//    *tavout = rout_current; 

    tav2 = rin;
    tav = rout;

    // findRAnnulusAroundDZE hívása: EZ KISZÁMOLJA AZ INDEX-HATÁROKAT AZ AKTUÁLIS SUGARAK ALAPJÁN
    // Ezt már a disk_params->radial_grid és disk_params->gas_pressure_gradient_vector alapján kellene, nem pedig külön paraméterekből.
    // Ha a findRAnnulusAroundDZE is disk_params-ot kapott, akkor rendben van.
    findRAnnulusAroundDZE(tav2, &ind_ii, &ind_io, tav, &ind_oi, &ind_oo, sim_opts, disk_params);

    double massii = 0, massoi = 0;
    double massiim = 0, massoim = 0;
    double massis = 0, massos = 0;

    calculateParticleMass(particle_number, partmassind, 
            (int)ind_ii, (int)ind_io, 
            (int)ind_oi, (int)ind_oo, 
            &massii, &massoi, sim_opts); 

    if(sim_opts->option_for_dust_secondary_population == 1.0) {
        calculateParticleMass(particle_number, partmassmicrind, 
                (int)ind_ii, (int)ind_io, 
                (int)ind_oi, (int)ind_oo, 
                &massiim, &massoim, sim_opts);
    }

    double massi, massim, masso, massom;

    if(tav2 != disk_params->r_dze_i) {
        massi = massii + massbtempii + massis;
        massim = massiim + massmtempii;
    } else {
        massi = massii + massis;
        massim = massiim;
    }
    if(tav != disk_params->r_dze_o) {
        masso = massoi + massbtempoi + massos;
        massom = massoim + massmtempoi;
    } else {
        masso = massoi + massos;
        massom = massoim;
    }

    *massbtempio = massi;
    *massbtempoo = masso;
    *massmtempio = massim;
    *massmtempoo = massom;

    *tavin = tav2;
    *tavout = tav;

    if (output_files->mass_file != NULL) {
        fprintf(output_files->mass_file, "%lg %lg %lg %lg %lg\n", step, *tavin, massi+massim, *tavout, masso+massom);
        fflush(output_files->mass_file);
    } else {
        fprintf(stderr, "WARNING: output_files->mass_file is NULL in printMassGrowthAtDZEFile. Cannot write mass data or fflush.\n");
    }

    if (dim > 0) {
        free(r_count);
    }
}


/* Függvény a sigma, p, dp kiíratására */
void printGasSurfaceDensityPressurePressureDerivateFile(const DiskParameters *disk_params, OutputFiles *output_files) {

    int i;

    if (output_files->surface_file == NULL) {
        fprintf(stderr, "ERROR: output_files->surface_file is NULL in printGasSurfaceDensityPressurePressureDerivateFile! Cannot write sigma data.\n");
        return;
    }

//%-15.6e %-15.6Lg %-15.6e %-15.6e\n",
    for(i = 1; i <= disk_params->grid_number; i++) { // Using disk_params->grid_number
        fprintf(output_files->surface_file, "%-15.6e %-15.6lg %-15.6e %15.6e\n", disk_params->radial_grid[i], disk_params->gas_surface_density_vector[i], disk_params->gas_pressure_vector[i], disk_params->gas_pressure_gradient_vector[i]);
    }

    fflush(output_files->surface_file);
}

/* Függvény a por felületisűrűségének kiíratására */
void printDustSurfaceDensityPressurePressureDerivateFile(const double *r, const double *rm, const double *dust_surfacedensity, const double *micron_dust_surfacedensity, const DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files, double step) {

    int i;

    if (output_files->dust_file == NULL) {
        fprintf(stderr, "ERROR: output_files->dust_file is NULL in printDustSurfaceDensityPressurePressureDerivateFile! Cannot write main dust surface density.\n");
        return;
    }

    for(i=0;i<particle_number;i++){ // particle_number from config.h
        if (r[i] >= disk_params->r_min) { // Using disk_params->r_min
            fprintf(output_files->dust_file,"%lg %.11lg  %lg \n",(double)step,r[i],dust_surfacedensity[i]);
        }

        if(sim_opts->option_for_dust_secondary_population == 1.0 && output_files->micron_dust_file != NULL) {
            if (rm[i] >= disk_params->r_min) { // Using disk_params->r_min
                fprintf(output_files->micron_dust_file,"%lg %lg  %lg \n",(double)step,rm[i],micron_dust_surfacedensity[i]);
            }
        }
    }

    fflush(output_files->dust_file);
    if(sim_opts->option_for_dust_secondary_population == 1.0 && output_files->micron_dust_file != NULL) {
        fflush(output_files->micron_dust_file); // Corrected this, it was micron_motion_file before
    }
}

/* Függvény a pormozgás és részecskeméret kiíratására */
void printDustParticleSizeFile(char *size_name, int step, double (*rad)[2], double (*radmicr)[2], const DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files) {

    FILE *fout_size = NULL;

    int i;

    if (sim_opts->option_for_dust_growth == 1.0) {
        fout_size = fopen(size_name, "w");
        if (fout_size == NULL) {
            fprintf(stderr, "ERROR: Could not open size file '%s' in printDustParticleSizeFile!\n", size_name);
            return;
        }
    }

    for (i = 0; i < particle_number; i++) { // particle_number from config.h

        if (output_files->dust_motion_file != NULL) {
            if (rad[i][0] >= disk_params->r_min) { // Using disk_params->r_min
                fprintf(output_files->dust_motion_file, "%lg %d %lg\n", (double)step, i, rad[i][0]);
            }
        } else {
            fprintf(stderr, "WARNING: output_files->dust_motion_file is NULL. Cannot write main particle motion.\n");
        }

        if (sim_opts->option_for_dust_secondary_population == 1.0 && output_files->micron_motion_file != NULL) {
            if (radmicr[i][0] >= disk_params->r_min) { // Using disk_params->r_min
                fprintf(output_files->micron_motion_file, "%lg %d %lg\n", (double)step, i, radmicr[i][0]);
            }
        }

        if (sim_opts->option_for_dust_growth == 1.0 && fout_size != NULL) {
            if (rad[i][0] >= disk_params->r_min) { // Using disk_params->r_min
                fprintf(fout_size, "%lg %lg %lg \n", (double)step, rad[i][0], rad[i][1] * AU_IN_CM); // AU_IN_CM from config.h
            }
        }
    }

    if (output_files->dust_motion_file != NULL) {
        fflush(output_files->dust_motion_file);
    } else {
        fprintf(stderr, "WARNING: Cannot fflush dust_motion_file, as it is NULL.\n");
    }

    if (sim_opts->option_for_dust_secondary_population == 1.0 && output_files->micron_motion_file != NULL) {
        fflush(output_files->micron_motion_file);
    }

    if (sim_opts->option_for_dust_growth == 1.0 && fout_size != NULL) {
        fclose(fout_size);
    }
}



// Függvény a fájl fejlécek kiírására
void printFileHeader(FILE *file, FileType_e file_type, const HeaderData *header_data) {
    if (file == NULL) {
        fprintf(stderr, "ERROR [printFileHeader]: Attempted to write header to a NULL file pointer!\n");
        return;
    }

    fprintf(file, "# Generated by Dust Drift Simulation (Date: %s %s)\n", __DATE__, __TIME__);

    switch (file_type) {
        case FILE_TYPE_DUST_MOTION:
            fprintf(file, "# Main Dust Particle Motion and Properties\n");
            fprintf(file, "# Columns: 1. Time Step (Simulation Time), 2. Particle Index, 3. Radius [AU]\n");
            break;

        case FILE_TYPE_MICRON_MOTION:
            fprintf(file, "# Micron Dust Particle Motion and Properties\n");
            fprintf(file, "# Columns: 1. Time Step (Simulation Time), 2. Particle Index, 3. Radius [AU]\n");
            fprintf(file, "# (Two-Population Model)\n");
            break;

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
                fprintf(file, "# Initial gas profile\n", header_data->current_time);
            } else {
                fprintf(file, "# Time: %e years\n", header_data ? header_data->current_time : 0.0);
            }
            // A fejléc az init_tool_module-ból, módosítva a HeaderData használatára
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

        case FILE_TYPE_DUST_MICRON_DENSITY:
            fprintf(file, "# Micron Dust Surface Density Profile\n");
            fprintf(file, "# Time: %e years\n", header_data ? header_data->current_time : 0.0);
            fprintf(file, "# Columns: 1. Radius [AU], 2. Sigma_micron_dust [M_Sun/AU^2]\n");
            break;

        case FILE_TYPE_PARTICLE_SIZE:
            // A fejléc az init_tool_module-ból, módosítva a HeaderData használatára
            // Kiegészítő infók, ha t=0 (is_initial_data)
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial particle distribution\n", header_data->current_time);
            } else {
                fprintf(file, "# Particle distribution, Time: %e years\n", header_data ? header_data->current_time : 0.0);
            }
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-5s %-15s %-20s %-20s %-15s %-15s\n",
                    "Index", "Radius_AU", "RepMass_Pop1_Msun", "RepMass_Pop2_Msun", "MaxPartSize_cm", "MicroSize_cm");
            fprintf(file, "#--------------------------------------------------------------------------\n");

            break;

        case FILE_TYPE_DISK_PARAM:
            // A fejléc az init_tool_module-ból, módosítva a HeaderData használatára
            fprintf(file, "# Disk Parameters\n");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            fprintf(file, "# %-15s %-15s %-10s %-15s %-20s %-15s %-15s %-15s %-20s %-20s %-15s %-15s %-15s %-15s %-15s\n",
                    "R_Min_AU", "R_Max_AU", "N_Grid", "SigmaExp", "Sigma0_gas_Msun_AU2",
                    "G_GravConst", "DzR_Inner_AU", "DzR_Outer_AU", "DzDr_Inner_Calc_AU", "DzDr_Outer_Calc_AU",
                    "DzAlphaMod", "DustDensity_g_cm3", "AlphaViscosity", "StarMass_Msun", "FlaringIndex");
            fprintf(file, "#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
            // A tényleges paraméter értékeket nem a fejlécbe írjuk, hanem a fő adatsorba.
            break;

        default:
            fprintf(stderr, "WARNING [printFileHeader]: Unknown file type for header generation: %d!\n", file_type);
            break;
    }
    fflush(file);
}




int setupInitialOutputFiles(OutputFiles *output_files, const SimulationOptions *sim_opts,
                               const DiskParameters *disk_params, HeaderData *header_data_for_files) {
    char porout[MAX_PATH_LEN] = "";
    char poroutmicr[MAX_PATH_LEN] = "";
    char massout[MAX_PATH_LEN] = "";

    // Készítsük elő a header_data_for_files struktúrát a specifikus adatokkal
    header_data_for_files->current_time = 0.0;
    header_data_for_files->is_initial_data = 1;
    header_data_for_files->R_in = disk_params->r_min;
    header_data_for_files->R_out = disk_params->r_max;

    // Fájlnevek generálása
    snprintf(porout, MAX_PATH_LEN, "%s/%s/%s%s", sim_opts->output_dir_name, kLogFilesDirectory, kDustAccumulationFileName,kFileNamesSuffix);

    if (sim_opts->option_for_dust_secondary_population == 1.0) {
        snprintf(poroutmicr, MAX_PATH_LEN, "%s/%s/%s%st", sim_opts->output_dir_name, kLogFilesDirectory,kDustMicronParticleEvolutionFile,kFileNamesSuffix);
    }
    snprintf(massout, MAX_PATH_LEN, "%s/%s/%s%s", sim_opts->output_dir_name, kLogFilesDirectory, kDustAccumulationFileName,kFileNamesSuffix);

    fprintf(stderr, "DEBUG [setupInitialOutputFiles]: Opening output files: %s, %s (if 2pop), %s\n", porout, poroutmicr, massout);

    // Fájlok megnyitása és fejlécek írása
    output_files->dust_motion_file = fopen(porout, "w");
    if (output_files->dust_motion_file == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", porout);
        return 1; // Hiba
    }
    printFileHeader(output_files->dust_motion_file, FILE_TYPE_DUST_MOTION, header_data_for_files);

    if (sim_opts->option_for_dust_secondary_population == 1.0) {
        output_files->micron_motion_file = fopen(poroutmicr, "w");
        if (output_files->micron_motion_file == NULL) {
            fprintf(stderr, "ERROR: Could not open %s\n", poroutmicr);
            // Itt fontos lehet, hogy felszabadítsuk az eddig megnyitott fájlokat
            fclose(output_files->dust_motion_file);
            output_files->dust_motion_file = NULL;
            return 1; // Hiba
        }
        printFileHeader(output_files->micron_motion_file, FILE_TYPE_MICRON_MOTION, header_data_for_files);
    }

    output_files->mass_file = fopen(massout, "w");
    if (output_files->mass_file == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", massout);
        // Itt is felszabadítjuk a már megnyitott fájlokat
        fclose(output_files->dust_motion_file);
        output_files->dust_motion_file = NULL;
        if (sim_opts->option_for_dust_secondary_population == 1.0 && output_files->micron_motion_file != NULL) {
            fclose(output_files->micron_motion_file);
            output_files->micron_motion_file = NULL;
        }
        return 1; // Hiba
    }
    printFileHeader(output_files->mass_file, FILE_TYPE_MASS_ACCUMULATION, header_data_for_files);

    return 0; // Siker
}


void cleanupSimulationResources(ParticleData *particle_data, OutputFiles *output_files, const SimulationOptions *sim_opts) {
    if (particle_number > 0) {
        free(particle_data->particle_distance_array); particle_data->particle_distance_array = NULL;
        free(particle_data->micron_particle_distance_array); particle_data->micron_particle_distance_array = NULL;
        free(particle_data->particle_mass_array); particle_data->particle_mass_array = NULL;
        free(particle_data->massmicradial_grid); particle_data->massmicradial_grid = NULL;
        free(particle_data->partmassind); particle_data->partmassind = NULL;
        free(particle_data->partmassmicrind); particle_data->partmassmicrind = NULL;
        free(particle_data->dust_surfacedensity); particle_data->dust_surfacedensity = NULL;
        free(particle_data->micron_dust_surfacedensity); particle_data->micron_dust_surfacedensity = NULL;
        free(particle_data->rdvec); particle_data->rdvec = NULL;
        free(particle_data->rmicvec); particle_data->rmicvec = NULL;

        fprintf(stderr, "DEBUG [cleanupSimulationResources]: All dynamically allocated particle arrays freed.\n");
    }

    if (output_files->dust_motion_file != NULL) {
        fclose(output_files->dust_motion_file);
        output_files->dust_motion_file = NULL;
        fprintf(stderr, "DEBUG [cleanupSimulationResources]: Closed %s%s\n", kDustParticleEvolutionFile,kFileNamesSuffix);
    }
    if (output_files->micron_motion_file != NULL) { // Ellenőrzés option_for_dust_secondary_population-ra itt is
        fclose(output_files->micron_motion_file);
        output_files->micron_motion_file = NULL;
        fprintf(stderr, "DEBUG [cleanupSimulationResources]: Closed %s%s\n",kDustMicronParticleEvolutionFile,kFileNamesSuffix);
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
    // Gáz sűrűség fájl
    snprintf(dens_name, MAX_PATH_LEN, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kGasDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix);

    // Por sűrűség fájl
    snprintf(dust_name, MAX_PATH_LEN, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kDustDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix);

    // Mikron por fájl
/*    snprintf(dust_name2, MAX_PATH_LEN, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kDustDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix);
*/
    // Részecskeméret fájl
    snprintf(size_name, MAX_PATH_LEN, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kDustParticleSizeFileName, snapshot_id, kFileNamesSuffix);
}



// Segédfüggvény a pillanatfelvételek fájljainak bezárására
void closeSnapshotFiles(OutputFiles *output_files, const char *dens_name, const char *dust_name, const char *dust_name2, const SimulationOptions *sim_opts) {
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