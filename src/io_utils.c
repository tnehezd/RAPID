#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>    // For errno
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h> // A gethostname-hez (Linux/Mac)
#include <pwd.h>    // A felhasználónévhez


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

// Local includes
#include "io_utils.h"
#include "config.h"         
#include "dust_physics.h"  
#include "utils.h"         
#include "simulation_types.h" 
#include "boundary_conditions.h"


#define INIT_DATA_HEADER_LINES 8
#define HEADER_WIDTH 75 // A '#' keret belső szélessége
// --- GLOBAL FILE POINTERS ---
// DO NOT DEFINE THEM HERE IF THEY ARE DEFINED ELSEWHERE (e.g., in main.c or config.c)
// They are simply used here because they are declared as 'extern' in config.h.

// --- FÜGGVÉNY DEFINÍCIÓK ---

/* Visszaadja, hogy hány sora van a beolvasandó fájlnak,
 * ez jelen esetben megadja a beolvasandó részecskék számát. */
int calculateNumbersOfParticles(const char *particle_data_file_name) {
    FILE *particle_file = NULL;
    char line_buffer[1024];
    int line_count = 0; // Counter for data lines

    particle_file = fopen(particle_data_file_name, "r");
    if (particle_file == NULL) {
        fprintf(stderr, "ERROR [calculateNumbersOfParticles]: Could not open file '%s'.\n", particle_data_file_name);
        perror("Reason"); // Prints system error message
        exit(EXIT_FAILURE);
    }

    // Skip header lines
    for (int i = 0; i < INIT_DATA_HEADER_LINES; i++) {
        if (fgets(line_buffer, sizeof(line_buffer), particle_file) == NULL) {
            // If file ends before all header lines are skipped, it's an error
            fprintf(stderr, "ERROR [calculateNumbersOfParticles]: Unexpected end of file while skipping %d header lines in '%s'.\n", INIT_DATA_HEADER_LINES, particle_data_file_name);
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
void loadDustParticlesFromFile(ParticleData *particle_data, const char *particle_data_file_name) {

    int i, dummy;
    double distance, particle_radius, micron_particle_radius;
    long double representative_mass;
    long double micron_representative_mass;

    // Use the global 'load_dust_particles_file' which is declared extern in config.h
    load_dust_particles_file = fopen(particle_data_file_name,"r"); // Use the passed particle_file_name

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
        if(fscanf(load_dust_particles_file,"%d %lg %Lg %Lg %lg %lg",&dummy,&distance,&representative_mass,&micron_representative_mass,&particle_radius,&micron_particle_radius) == 6) {
            particle_data->particle_distance_array[i][0] = distance;
            particle_data->particle_distance_array[i][1] = particle_radius / AU_IN_CM; // AU_IN_CM from config.h
            particle_data->dust_particle_mass_grid[i] = representative_mass;

            particle_data->micron_particle_distance_array[i][0] = distance;
            particle_data->micron_particle_distance_array[i][1] = micron_particle_radius / AU_IN_CM; // AU_IN_CM from config.h
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

    // Fejléc sorok átugrása
    // Addig olvassuk a sorokat, amíg komment (#) karakterrel kezdődnek
    // VAGY amíg nem a '---' elválasztó sort találjuk, ami a valódi adatok kezdetét jelzi
    while (fgets(line, sizeof(line), input_file) != NULL) {
        // Ellenőrizzük, hogy a sor komment-e vagy a '---' elválasztó
        if (line[0] == '#' || strncmp(line, "---", 3) == 0) { // strncmp "---" karakterekre
            continue; // Ugrás a következő sorra
        } else {
            // Ez az első valós adatsor. Visszaállítjuk a fájlmutatót a sor elejére.
            fseek(input_file, -strlen(line), SEEK_CUR); 
            break;
        }
    }

    // Ha a fájl üres vagy csak kommenteket tartalmaz
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

    // A ciklus disk_params->grid_number-ig megy
    for (int i = 0; i < disk_params->grid_number; i++) {
        // MOST MÁR PONTOSAN A FÁJLOD FORMÁTUMÁT OLVASSUK BE:
        // Radius_AU, GasSurfDensity, GasPressure, GasPressureDeriv (4 oszlop, mind double)
        if (fscanf(input_file, "%d %lf %lf %lf %lf",
                         &index, &r_value, &surfacedensity_gas_value, &gas_pressure_value, &gas_pressure_gradient_value) != 5) {
            // Hiba kezelése, ha nem tudunk 4 double értéket beolvasni
            fprintf(stderr, "ERROR [loadGasSurfaceDensityFromFile]: Failed to read 4 values for row %d from file '%s'. File may be malformed or ended unexpectedly.\n", i, input_disk_file_name);
            fclose(input_file);
            exit(EXIT_FAILURE);
        }

        // Hozzárendelés a disk_params tömbökhöz
        // Az indexelés 'i + 1' a 0-ás indexű szellemcella miatt (ahogy a DiskParameters definíciója és a applyBoundaryConditions függvény valószínűsíti).
        // Fontos: ellenőrizzük, hogy az 'i + 1' index a tömb határain belül van-e.
        // A tömbök mérete disk_params->grid_number + 2, tehát az érvényes indexek 0-tól grid_number+1-ig mennek.
        // A "valós" adatok 1-től grid_number-ig kerülnek, a 0 és grid_number+1 pedig a applyBoundaryConditions-hez.
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

    // első próbálkozás
    asprintf(&temporary_path, "%s", dir_path);

    // ha létezik, generáljunk új neveket
    while (access(temporary_path, F_OK) == 0) {
        free(temporary_path);
        asprintf(&temporary_path, "%s_%04d", dir_path, ++counter);

        if (counter > 99) {
            fprintf(stderr, "ERROR: Too many directories.\n");
            exit(1);
        }
    }

    // létrehozás
    if (MKDIR_CALL(temporary_path) != 0) {
        perror("mkdir failed");
        exit(1);
    }

    return temporary_path;  // <-- VISSZAADJUK AZ ÚJ STRINGET
}

void printCentered(FILE *file, const char *text) {
    int len = strlen(text);
    int padding = (HEADER_WIDTH - len) / 2;
    
    fprintf(file, "#"); // Keret széle
    for (int i = 0; i < padding; i++) fprintf(file, " ");
    fprintf(file, "%s", text);
    for (int i = 0; i < (HEADER_WIDTH - len - padding); i++) fprintf(file, " ");
    fprintf(file, "#\n"); // Keret másik széle
}

/* Elkészít egy részletes információs fájlt a futás paramétereivel és a környezettel */
void printCurrentInformationAboutRun(const char *directory_name, const DiskParameters *disk_params) {
    char *full_path = NULL;
    char file_name[100];
    
    // 1. Időbélyeg generálása
    time_t rawtime;
    struct tm *timeinfo;
    char time_buffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", timeinfo);

    // 2. Gépnév és Felhasználó lekérése
    char hostname[1024];
    gethostname(hostname, 1024);
    struct passwd *pw = getpwuid(getuid());
    const char *user = pw ? pw->pw_name : "unknown";

    // 3. Fájlútvonal összeállítása
    sprintf(file_name, "%s%s", kCurrentInfoFile, kFileNamesSuffix);
    asprintf(&full_path, "%s/%s", directory_name, file_name);

    fprintf(stderr, "DEBUG [printCurrentInformationAboutRun]: Writing run info to: '%s'\n", full_path);

    FILE *info_f = fopen(full_path, "w");
    if (info_f == NULL) {
        fprintf(stderr, "ERROR: Could not open info file '%s'.\n", full_path);
        if (full_path) free(full_path);
        return;
    }

    // --- FEJLÉC ÉS RENDSZERINFÓ ---
    fprintf(info_f, "==========================================================\n");
    fprintf(info_f, "        DUST DRIFT SIMULATION - RUN SNAPSHOT\n");
    fprintf(info_f, "==========================================================\n");
    fprintf(info_f, "  Run Started:       %s\n", time_buffer);
    fprintf(info_f, "  Host Machine:      %s\n", hostname);
    fprintf(info_f, "  User:              %s\n", user);
    
#ifdef _OPENMP
    fprintf(info_f, "  Parallel Threads:  %d\n", omp_get_max_threads());
#else
    fprintf(info_f, "  Parallel Threads:  1 (Serial mode)\n");
#endif
    fprintf(info_f, "  Binary Compiled:   %s %s\n", __DATE__, __TIME__);
    fprintf(info_f, "  Output Directory:  %s\n", directory_name);
    fprintf(info_f, "==========================================================\n\n");

    // --- FIZIKAI PARAMÉTEREK ---
    fprintf(info_f, "--- [ Central Star ] ---\n");
    fprintf(info_f, "  Stellar Mass:      %.4f M_Sun\n\n", disk_params->stellar_mass);

    fprintf(info_f, "--- [ Disk Geometry & Gas ] ---\n");
    fprintf(info_f, "  Radial Range:      %.2f - %.2f AU\n", disk_params->r_min, disk_params->r_max);
    fprintf(info_f, "  Gas Grid Points:   %d\n", disk_params->grid_number);
    fprintf(info_f, "  Sigma_0 (1 AU):    %.4e M_Sun/AU^2\n", disk_params->sigma_0);
    fprintf(info_f, "  Sigma Exponent:    %.4f\n", disk_params->sigma_power_law_index);
    fprintf(info_f, "  Aspect Ratio (H/R): %.4f\n", disk_params->h_aspect_ratio);
    fprintf(info_f, "  Flaring Index:     %.4f\n", disk_params->flaring_index);
    fprintf(info_f, "  Alpha Viscosity:   %.4e\n\n", disk_params->alpha_parameter);

    // --- DEAD ZONE ---
    fprintf(info_f, "--- [ Dead Zone Configuration ] ---\n");
    if (disk_params->r_dze_i > 0.0 || disk_params->r_dze_o > 0.0) {
        fprintf(info_f, "  Status:            ENABLED\n");
        fprintf(info_f, "  Inner DZE Radius:  %.2f AU (Trans. width: %.2f)\n", disk_params->r_dze_i, disk_params->dr_dze_i);
        fprintf(info_f, "  Outer DZE Radius:  %.2f AU (Trans. width: %.2f)\n", disk_params->r_dze_o, disk_params->dr_dze_o);
        fprintf(info_f, "  Alpha Mod Factor:  %.4e\n", disk_params->alpha_parameter_modification);
    } else {
        fprintf(info_f, "  Status:            DISABLED (Uniform alpha)\n");
    }

    // --- POR PARAMÉTEREK ---
    fprintf(info_f, "\n--- [ Dust Properties ] ---\n");
    fprintf(info_f, "  Particle Density:  %.2f g/cm^3\n", disk_params->particle_density);
    fprintf(info_f, "  Fragmentation Vel: %.2f cm/s\n", disk_params->fragmentation_velocity);
    fprintf(info_f, "  Global Dust Count: %d\n", particle_number);

    fprintf(info_f, "\n==========================================================\n");
    fprintf(info_f, "         End of Configuration Summary\n");
    fprintf(info_f, "==========================================================\n");

    fclose(info_f);
    if (full_path) free(full_path);
}


void printMassGrowthAtDZEFile(double step, double (*dust_particle_mass_array)[5], double (*micron_dust_particle_mass_array)[5], double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, 
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

    calculateParticleMass(particle_number, dust_particle_mass_array, 
            (int)ind_ii, (int)ind_io, 
            (int)ind_oi, (int)ind_oo, 
            &massii, &massoi, sim_opts); 

    if(sim_opts->option_for_dust_secondary_population == 1.0) {
        calculateParticleMass(particle_number, micron_dust_particle_mass_array, 
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
        fprintf(output_files->mass_file, "%-10d %-15.6f %-20.6e %-15.6f %-15.6e\n", 
                    (int)step, 
                    (double)*tavin, 
                    (double)(massi + massim), 
                    (double)*tavout, 
                    (double)(masso + massom));
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
        fprintf(output_files->surface_file, "%-5d %-15.6e %-15.6le %-15.6e %15.6e\n",i, disk_params->radial_grid[i], disk_params->gas_surface_density_vector[i], disk_params->gas_pressure_vector[i], disk_params->gas_pressure_gradient_vector[i]);
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
void printDustParticleSizeFile(char *size_name, int step, double (*rad)[2], double (*micron_particle_radius)[2], const DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files) {

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

        if (sim_opts->option_for_dust_secondary_population == 1.0 && output_files->micron_motion_file != NULL) {
            if (micron_particle_radius[i][0] >= disk_params->r_min) { // Using disk_params->r_min
                fprintf(output_files->micron_motion_file, "%lg %d %lg %lg\n", (double)step, i, micron_particle_radius[i][0],micron_particle_radius[i][1] * AU_IN_CM);
            }
        }

        if (sim_opts->option_for_dust_growth == 1.0 && fout_size != NULL) {
            if (rad[i][0] >= disk_params->r_min) { // Using disk_params->r_min
                fprintf(fout_size, "%lg %d %lg %lg \n", (double)step, i, rad[i][0], rad[i][1] * AU_IN_CM); // AU_IN_CM from config.h
            }
        }
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

    char buffer[128]; // Ideiglenes tároló a formázott szövegeknek
    
    // Felső keret
    fprintf(file, "#############################################################################\n");
    // 1. Sor: Program név
    printCentered(file, "Generated by Dust Drift Simulation with the RAPID code");
    // 2. Sor: Verzió
    sprintf(buffer, "Version: %s | Compiled: %s %s", SIM_VERSION, __DATE__, __TIME__);
    printCentered(file, buffer);
    // Alsó keret
    fprintf(file, "#############################################################################\n");

    switch (file_type) {

        case FILE_TYPE_MASS_ACCUMULATION:
            fprintf(file, "# This file contains the time evolution of dust mass within specified disk regions.\n");
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-5s %-15s %-15s %-15s %-15s\n",
                    "Time", "Inner_DZE_AU", "Accum._mass_inner_MSun", "Outer_DZE_AU", "Accum._mass_outer_MSun");
            fprintf(file, "#--------------------------------------------------------------------------\n");            
            break;

        case FILE_TYPE_GAS_DENSITY:
            // Kiegészítő infók, ha t=0 (is_initial_data)
            if (header_data && header_data->is_initial_data) {
                fprintf(file, "# Initial gas profile\n");
            } else {
                fprintf(file, "# Time: %e years\n", header_data ? header_data->current_time : 0.0);
            }
            // A fejléc az init_tool_module-ból, módosítva a HeaderData használatára
            fprintf(file, "#--------------------------------------------------------------------------\n");
            fprintf(file, "# %-5s %-15s %-15s %-15s %-15s\n",
                    "Index", "Radius_AU", "GasSurfDensity", "GasPressure", "GasPressureDeriv");
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



void printFinalSimulationSummary(const char *directory_name, double elapsed_seconds, const SimulationOptions *sim_opts) {
    char *full_path = NULL;
    asprintf(&full_path, "%s/%s%s", directory_name, kCurrentRuntimeInfoFile,kFileNamesSuffix);

    FILE *f = fopen(full_path, "w");
    if (f == NULL) {
        if (full_path) free(full_path);
        return;
    }

    // Idő átváltása olvasható formátumba
    int h = (int)elapsed_seconds / 3600;
    int m = ((int)elapsed_seconds % 3600) / 60;
    int s = (int)elapsed_seconds % 60;

    // Szimulációs sebesség (szimulált év / valódi másodperc)
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
    char *mass_output = NULL;

    // Fejléc adatok előkészítése
    header_data_for_files->current_time = 0.0;
    header_data_for_files->is_initial_data = 1;
    header_data_for_files->R_in = disk_params->r_min;
    header_data_for_files->R_out = disk_params->r_max;

    // --- Fájlnevek generálása ---
    
    
    // Itt marad a kDustAccumulationFileName
    asprintf(&mass_output, "%s/%s/%s%s", sim_opts->output_dir_name, kLogFilesDirectory, kDustAccumulationFileName, kFileNamesSuffix);

    fprintf(stderr, "DEBUG [setupInitialOutputFiles]: Paths:\n  Mass: %s\n", mass_output);

    // --- Fájlok megnyitása ---

    


    // 3. Tömeggyarapodás fájl
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
    // sprintf-et használunk asprintf helyett, mert a memóriát a hívó fél foglalta le!
    sprintf(dens_name, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kGasDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix);

    sprintf(dust_name, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kDustDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix);

    if (sim_opts->option_for_dust_secondary_population == 1) {
        sprintf(dust_name2, "%s/%s/%s_%08d%s",
                 sim_opts->output_dir_name, kLogFilesDirectory,
                 kMicronDustDensityProfileFilePrefix, snapshot_id, kFileNamesSuffix); // Javított prefix
    }

    sprintf(size_name, "%s/%s/%s_%08d%s",
             sim_opts->output_dir_name, kLogFilesDirectory,
             kDustParticleSizeFileName, snapshot_id, kFileNamesSuffix);
}



// Segédfüggvény a pillanatfelvételek fájljainak bezárására
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