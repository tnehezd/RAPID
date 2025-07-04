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
#include "config.h"         // Defines PARTICLE_NUMBER, AU2CM, FILENAME_INIT_PROFILE, and declares extern FILE *fin1, extern FILE *jelfut
#include "dust_physics.h"   // If needed for any specific function interactions
#include "utils.h"          // For find_num_zero, find_zero, find_r_annulus
#include "simulation_types.h" // For disk_t, simulation_options_t, output_files_t


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

/* A porreszecskek adatainak beolvasasa */
void por_be(double radius[][2], double radiusmicr[][2], double *mass, double *massmicr, const char *filename) {

    int i, dummy;
    double distance, particle_radius, radmicr;
    long double reprmass;
    long double reprmassmicr;

    // Use the global 'fin1' which is declared extern in config.h
    fin1 = fopen(filename,"r"); // Use the passed filename

    if (fin1 == NULL) {
        fprintf(stderr, "ERROR [por_be]: Could not open file '%s'.\n", filename);
        perror("Reason");
        exit(EXIT_FAILURE);
    }

    char line_buffer[1024];
    for (int k = 0; k < INIT_DATA_HEADER_LINES; k++) {
        if (fgets(line_buffer, sizeof(line_buffer), fin1) == NULL) {
            fprintf(stderr, "ERROR [por_be]: Unexpected end of file while skipping headers in '%s'.\n", filename);
            fclose(fin1);
            exit(EXIT_FAILURE);
        }
    }

    for (i = 0; i < PARTICLE_NUMBER; i++) {
        if(fscanf(fin1,"%d %lg %Lg %Lg %lg %lg",&dummy,&distance,&reprmass,&reprmassmicr,&particle_radius,&radmicr) == 6) {
            radius[i][0] = distance;
            radius[i][1] = particle_radius / AU2CM; // AU2CM from config.h
            mass[i] = reprmass;

            radiusmicr[i][0] = distance;
            radiusmicr[i][1] = radmicr / AU2CM; // AU2CM from config.h
            massmicr[i] = reprmassmicr;
        } else {
            fprintf(stderr, "\n\n******************* ERROR!      *********************\n\n");
            fprintf(stderr, "   Failed to read line %d from particle data file '%s'!\n", i, filename);
            fprintf(stderr, "   Expected 6 values, but fscanf failed. Program will exit.\n");
            fclose(fin1);
            exit(EXIT_FAILURE);
        }
    }

    fclose(fin1);
    printf("\n\n ******* A file beolvasasa sikerult!     ******* \n ******* Uss egy ENTER-t a folytatashoz! ******* \n\n ");
}


void sigIn(disk_t *disk_params, const char *filename) {
    const char *input_filename = filename;

    FILE *fp = fopen(input_filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR [sigIn]: Could not open input file '%s'.\n", input_filename);
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
        fprintf(stderr, "ERROR [sigIn]: File '%s' is empty or only contains comments/headers.\n", input_filename);
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
            fprintf(stderr, "ERROR [sigIn]: Failed to read 4 values for row %d from file '%s'. File may be malformed or ended unexpectedly.\n", i, input_filename);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        // Hozzárendelés a disk_params tömbökhöz
        // Az indexelés 'i + 1' a 0-ás indexű szellemcella miatt (ahogy a disk_t definíciója és a Perem függvény valószínűsíti).
        // Fontos: ellenőrizzük, hogy az 'i + 1' index a tömb határain belül van-e.
        // A tömbök mérete disk_params->NGRID + 2, tehát az érvényes indexek 0-tól NGRID+1-ig mennek.
        // A "valós" adatok 1-től NGRID-ig kerülnek, a 0 és NGRID+1 pedig a Perem-hez.
        if ((i + 1) >= 0 && (i + 1) <= disk_params->NGRID + 1) { 
            disk_params->rvec[i + 1] = r_val;
            disk_params->sigmavec[i + 1] = sigma_gas_val;
            disk_params->pressvec[i + 1] = pressure_gas_val;
            disk_params->dpressvec[i + 1] = dpressure_dr_val;
        } else {
            fprintf(stderr, "WARNING [sigIn]: Attempted to write to out-of-bounds index %d. Max allowed index: %d (NGRID+1).\n", i + 1, disk_params->NGRID + 1);
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

    while (access(tmp_path, F_OK) == 0) {  // Fájl vagy mappa létezik
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

/* Elkészít egy fájlt, ami tartalmazza a jelenlegi futás paramétereit,
 * és hogy melyik mappában találhatóak a fájlok */
void infoCurrent(const char *nev, const disk_t *disk_params, const simulation_options_t *sim_opts) {

    char full_path[MAX_PATH_LEN]; // Használjuk a MAX_PATH_LEN-t a biztonságos puffereléshez
    char file_name[100]; // Csak a fájlnév, pl. "run_0.dat"

    sprintf(file_name, "run_%i.dat", (int)sim_opts->TCURR);
    
    // Építsük fel a teljes elérési utat: <nev>/<file_name>
    snprintf(full_path, sizeof(full_path), "%s/%s", nev, file_name);

    fprintf(stderr, "DEBUG [infoCurrent]: Attempting to open file: '%s'\n", full_path);

    jelfut = fopen(full_path, "w"); // Most már a teljes elérési utat használja

    if (jelfut == NULL) {
        fprintf(stderr, "ERROR [infoCurrent]: Could not open file '%s'.\n", full_path);
        perror("Reason");
        // Don't exit here, it's not critical, just warn and return
        return;
    }

    fprintf(jelfut,"A jelenlegi futás a %s mappaban taláható!\n",nev);
    fprintf(jelfut,"\n\nA korong paraméterei:\nRMIN: %lg, RMAX: %lg\nSIGMA0: %lg, SIGMA_EXP: %lg, flaring index: %lg\nALPHA_VISC: %lg, ALPHA_MOD: %lg\nR_DZE_I: %lg, R_DZE_O: %lg, DR_DZEI: %lg, DR_DZE_O: %lg   (*** R_DZE_I/O = 0, akkor azt a DZE-t nem szimulálja a futás! ***)\n\n\n",
              disk_params->RMIN, disk_params->RMAX,
              disk_params->SIGMA0, disk_params->SIGMAP_EXP, disk_params->FLIND,
              disk_params->alpha_visc, disk_params->a_mod,
              disk_params->r_dze_i, disk_params->r_dze_o, disk_params->Dr_dze_i, disk_params->Dr_dze_o);
    fprintf(jelfut,"A központi csillag tömege: %lg M_Sun\n", disk_params->STAR_MASS);
    fclose(jelfut);
}


void Print_Mass(double step, double (*partmassind)[5], double (*partmassmicrind)[5], double (*partmasssecind)[5], double t, double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, 
                double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, double *tavin, double *tavout, 
                const disk_t *disk_params, const simulation_options_t *sim_opts,output_files_t *output_files) {


    double ind_ii, ind_io, ind_oi, ind_oo, tav, tav2;

    tav = disk_params->r_dze_o;
    tav2 = disk_params->r_dze_i;

    int dim = find_num_zero(disk_params); 
    double *r_count = (double *)malloc(sizeof(double) * dim); 

    if (dim > 0 && r_count == NULL) {
        fprintf(stderr, "ERROR [Print_Mass]: Failed to allocate memory for r_count. Exiting.\n");
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
        for(i = 0; i < disk_params->NGRID; i++) {
            // Itt a find_zero valószínűleg disk_params->rvec és disk_params->dpressvec-et használ
            temp_new = find_zero(i,disk_params->rvec,disk_params->dpressvec); 
            if(temp != temp_new && i > 3 && temp_new != 0.0) {
                if (j < dim) { 
                    r_count[j] = temp_new;
                    j++;
                } else {
                    fprintf(stderr, "WARNING [Print_Mass]: r_count array overflow, skipping data. dim: %d, j: %d\n", dim, j);
                }
            }
            if(sim_opts->dzone == 0.0) { 
                if(temp_new > 0.) {
                    temp = temp_new;
                    rout_new = temp;
                }
            }
        }
    }

    if(sim_opts->dzone == 1.0) {
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
    if(sim_opts->dzone == 0.0) rin = 0;
//    double rout_current = rout_new;
    rout = rout_new;
    
//    *tavin = rin;  
//    *tavout = rout_current; 

    tav2 = rin;
    tav = rout;

    // find_r_annulus hívása: EZ KISZÁMOLJA AZ INDEX-HATÁROKAT AZ AKTUÁLIS SUGARAK ALAPJÁN
    // Ezt már a disk_params->rvec és disk_params->dpressvec alapján kellene, nem pedig külön paraméterekből.
    // Ha a find_r_annulus is disk_params-ot kapott, akkor rendben van.
    find_r_annulus(tav2, &ind_ii, &ind_io, tav, &ind_oi, &ind_oo, sim_opts, disk_params);

    double massii = 0, massoi = 0;
    double massiim = 0, massoim = 0;
    double massis = 0, massos = 0;

    GetMass(PARTICLE_NUMBER, partmassind, 
            (int)ind_ii, (int)ind_io, 
            (int)ind_oi, (int)ind_oo, 
            &massii, &massoi, sim_opts); 

    if(sim_opts->twopop == 1.0) {
        GetMass(4*PARTICLE_NUMBER, partmasssecind, 
                (int)ind_ii, (int)ind_io, 
                (int)ind_oi, (int)ind_oo, 
                &massis, &massos, sim_opts); 
        GetMass(PARTICLE_NUMBER, partmassmicrind, 
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
        fprintf(stderr, "WARNING: output_files->mass_file is NULL in Print_Mass. Cannot write mass data or fflush.\n");
    }

    if (dim > 0) {
        free(r_count);
    }
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
void Print_Sigmad(const double *r, const double *rm, const double *sigmad, const double *sigmadm, const disk_t *disk_params, const simulation_options_t *sim_opts, output_files_t *output_files) {

    int i;

    if (output_files->dust_file == NULL) {
        fprintf(stderr, "ERROR: output_files->dust_file is NULL in Print_Sigmad! Cannot write main dust surface density.\n");
        return;
    }

    for(i=0;i<PARTICLE_NUMBER;i++){ // PARTICLE_NUMBER from config.h
        if (r[i] >= disk_params->RMIN) { // Using disk_params->RMIN
            fprintf(output_files->dust_file,"%.11lg  %lg \n",r[i],sigmad[i]);
        }

        if(sim_opts->twopop == 1.0 && output_files->micron_dust_file != NULL) {
            if (rm[i] >= disk_params->RMIN) { // Using disk_params->RMIN
                fprintf(output_files->micron_dust_file,"%lg  %lg \n",rm[i],sigmadm[i]);
            }
        }
    }

    fflush(output_files->dust_file);
    if(sim_opts->twopop == 1.0 && output_files->micron_dust_file != NULL) {
        fflush(output_files->micron_dust_file); // Corrected this, it was micron_motion_file before
    }
}

/* Függvény a pormozgás és részecskeméret kiíratására */
void Print_Pormozg_Size(char *size_name, int step, double (*rad)[2], double (*radmicr)[2], const disk_t *disk_params, const simulation_options_t *sim_opts, output_files_t *output_files) {

    FILE *fout_size = NULL;

    int i;

    if (sim_opts->growth == 1.0) {
        fout_size = fopen(size_name, "w");
        if (fout_size == NULL) {
            fprintf(stderr, "ERROR: Could not open size file '%s' in Print_Pormozg_Size!\n", size_name);
            return;
        }
    }

    for (i = 0; i < PARTICLE_NUMBER; i++) { // PARTICLE_NUMBER from config.h

        if (output_files->por_motion_file != NULL) {
            if (rad[i][0] >= disk_params->RMIN) { // Using disk_params->RMIN
                fprintf(output_files->por_motion_file, "%lg %d %lg\n", (double)step, i, rad[i][0]);
            }
        } else {
            fprintf(stderr, "WARNING: output_files->por_motion_file is NULL. Cannot write main particle motion.\n");
        }

        if (sim_opts->twopop == 1.0 && output_files->micron_motion_file != NULL) {
            if (radmicr[i][0] >= disk_params->RMIN) { // Using disk_params->RMIN
                fprintf(output_files->micron_motion_file, "%lg %d %lg\n", (double)step, i, radmicr[i][0]);
            }
        }

        if (sim_opts->growth == 1.0 && fout_size != NULL) {
            if (rad[i][0] >= disk_params->RMIN) { // Using disk_params->RMIN
                fprintf(fout_size, "%lg %lg %lg \n", (double)step, rad[i][0], rad[i][1] * AU2CM); // AU2CM from config.h
            }
        }
    }

    if (output_files->por_motion_file != NULL) {
        fflush(output_files->por_motion_file);
    } else {
        fprintf(stderr, "WARNING: Cannot fflush por_motion_file, as it is NULL.\n");
    }

    if (sim_opts->twopop == 1.0 && output_files->micron_motion_file != NULL) {
        fflush(output_files->micron_motion_file);
    }

    if (sim_opts->growth == 1.0 && fout_size != NULL) {
        fclose(fout_size);
    }
}

/* Az időt tartalmazó fájl paramétereinek beolvasása (vagy beállítása) */
void timePar(double tMax_val, double stepping_val, double current_val, simulation_options_t *sim_opts) {

    printf("\n\n *********** A korong parameterei sikeresen beolvasva!  *********** \n      tmax: %lg, a program %lg evenként írja ki a file-okat\n\n\n", tMax_val, stepping_val);

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

        case FILE_TYPE_DUST_MICRON_DENSITY:
            fprintf(file, "# Micron Dust Surface Density Profile\n");
            fprintf(file, "# Time: %e years\n", header_data ? header_data->current_time : 0.0);
            fprintf(file, "# Columns: 1. Radius [AU], 2. Sigma_micron_dust [M_Sun/AU^2]\n");
            break;

        case FILE_TYPE_PARTICLE_SIZE:
            // A fejléc az init_tool_module-ból, módosítva a HeaderData_t használatára
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

        default:
            fprintf(stderr, "WARNING [print_file_header]: Unknown file type for header generation: %d!\n", file_type);
            break;
    }
    fflush(file);
}




int setup_initial_output_files(output_files_t *output_files, const simulation_options_t *sim_opts,
                               const disk_t *disk_params, HeaderData_t *header_data_for_files) {
    char porout[MAX_PATH_LEN] = "";
    char poroutmicr[MAX_PATH_LEN] = "";
    char massout[MAX_PATH_LEN] = "";

    // Készítsük elő a header_data_for_files struktúrát a specifikus adatokkal
    header_data_for_files->current_time = 0.0;
    header_data_for_files->is_initial_data = 1;
    header_data_for_files->R_in = disk_params->RMIN;
    header_data_for_files->R_out = disk_params->RMAX;

    // Fájlnevek generálása
    snprintf(porout, MAX_PATH_LEN, "%s/%s/%s", sim_opts->output_dir_name, LOGS_DIR, FILE_DUST_EVOLUTION);

    if (sim_opts->twopop == 1.0) {
        snprintf(poroutmicr, MAX_PATH_LEN, "%s/%s/micron_particle_evolution.dat", sim_opts->output_dir_name, LOGS_DIR);
    }
    snprintf(massout, MAX_PATH_LEN, "%s/%s/%s.dat", sim_opts->output_dir_name, LOGS_DIR, FILE_MASS_ACCUMULATE);

    fprintf(stderr, "DEBUG [setup_initial_output_files]: Opening output files: %s, %s (if 2pop), %s\n", porout, poroutmicr, massout);

    // Fájlok megnyitása és fejlécek írása
    output_files->por_motion_file = fopen(porout, "w");
    if (output_files->por_motion_file == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", porout);
        return 1; // Hiba
    }
    print_file_header(output_files->por_motion_file, FILE_TYPE_DUST_MOTION, header_data_for_files);

    if (sim_opts->twopop == 1.0) {
        output_files->micron_motion_file = fopen(poroutmicr, "w");
        if (output_files->micron_motion_file == NULL) {
            fprintf(stderr, "ERROR: Could not open %s\n", poroutmicr);
            // Itt fontos lehet, hogy felszabadítsuk az eddig megnyitott fájlokat
            fclose(output_files->por_motion_file);
            output_files->por_motion_file = NULL;
            return 1; // Hiba
        }
        print_file_header(output_files->micron_motion_file, FILE_TYPE_MICRON_MOTION, header_data_for_files);
    }

    output_files->mass_file = fopen(massout, "w");
    if (output_files->mass_file == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", massout);
        // Itt is felszabadítjuk a már megnyitott fájlokat
        fclose(output_files->por_motion_file);
        output_files->por_motion_file = NULL;
        if (sim_opts->twopop == 1.0 && output_files->micron_motion_file != NULL) {
            fclose(output_files->micron_motion_file);
            output_files->micron_motion_file = NULL;
        }
        return 1; // Hiba
    }
    print_file_header(output_files->mass_file, FILE_TYPE_MASS_ACCUMULATION, header_data_for_files);

    return 0; // Siker
}


void cleanup_simulation_resources(ParticleData_t *p_data, output_files_t *output_files, const simulation_options_t *sim_opts) {
    if (PARTICLE_NUMBER > 0) {
        free(p_data->radius); p_data->radius = NULL;
        free(p_data->radiusmicr); p_data->radiusmicr = NULL;
        free(p_data->radius_rec); p_data->radius_rec = NULL;
        free(p_data->massvec); p_data->massvec = NULL;
        free(p_data->massmicrvec); p_data->massmicrvec = NULL;
        free(p_data->partmassind); p_data->partmassind = NULL;
        free(p_data->partmassmicrind); p_data->partmassmicrind = NULL;
        free(p_data->sigmad); p_data->sigmad = NULL;
        free(p_data->sigmadm); p_data->sigmadm = NULL;
        free(p_data->rdvec); p_data->rdvec = NULL;
        free(p_data->rmicvec); p_data->rmicvec = NULL;

        // Csak akkor szabadítjuk fel, ha twopop módban vannak allokálva
        if (sim_opts->twopop == 1) { // Feltételezve, hogy twopop == 1 esetén vannak ezek allokálva 4*PARTICLE_NUMBER méretben
            free(p_data->radiussec); p_data->radiussec = NULL;
            free(p_data->masssecvec); p_data->masssecvec = NULL;
            free(p_data->partmasssecind); p_data->partmasssecind = NULL;
            free(p_data->sigmads); p_data->sigmads = NULL;
            free(p_data->rsvec); p_data->rsvec = NULL;
        }

        fprintf(stderr, "DEBUG [cleanup_simulation_resources]: All dynamically allocated particle arrays freed.\n");
    }

    if (output_files->por_motion_file != NULL) {
        fclose(output_files->por_motion_file);
        output_files->por_motion_file = NULL;
        fprintf(stderr, "DEBUG [cleanup_simulation_resources]: Closed %s\n", FILE_DUST_EVOLUTION);
    }
    if (output_files->micron_motion_file != NULL) { // Ellenőrzés twopop-ra itt is
        fclose(output_files->micron_motion_file);
        output_files->micron_motion_file = NULL;
        fprintf(stderr, "DEBUG [cleanup_simulation_resources]: Closed micron_particle_evolution.dat\n");
    }
    if (output_files->mass_file != NULL) {
        fclose(output_files->mass_file);
        output_files->mass_file = NULL;
        fprintf(stderr, "DEBUG [cleanup_simulation_resources]: Closed %s\n", FILE_MASS_ACCUMULATE);
    }
}

// Segédfüggvény a pillanatfelvételek fájljainak bezárására
void close_snapshot_files(output_files_t *output_files, const char *dens_name, const char *dust_name, const char *dust_name2, const simulation_options_t *sim_opts) {
    if (output_files->surface_file != NULL) {
        fclose(output_files->surface_file);
        output_files->surface_file = NULL;
        fprintf(stderr, "DEBUG [close_snapshot_files]: Closed %s.\n", dens_name);
    }
    if (output_files->dust_file != NULL) {
        fclose(output_files->dust_file);
        output_files->dust_file = NULL;
        fprintf(stderr, "DEBUG [close_snapshot_files]: Closed %s.\n", dust_name);
    }
    if (sim_opts->twopop == 1 && output_files->micron_dust_file != NULL) {
        fclose(output_files->micron_dust_file);
        output_files->micron_dust_file = NULL;
        fprintf(stderr, "DEBUG [close_snapshot_files]: Closed %s.\n", dust_name2);
    }
}