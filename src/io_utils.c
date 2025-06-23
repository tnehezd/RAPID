#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>    // For errno

#ifdef _WIN32
#include <direct.h> // For _mkdir on Windows
#define MKDIR_CALL(dir) _mkdir(dir)
#else
#include <sys/stat.h> // For mkdir on Linux/macOS
#include <sys/types.h>
#define MKDIR_CALL(dir) mkdir(dir, 0755) // 0755 permissions: owner rwx, group rx, others rx
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
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (line[0] == '#') {
            continue;
        } else {
            fseek(fp, -strlen(line), SEEK_CUR);
            break;
        }
    }

    if (feof(fp) && line[0] == '#') {
        fprintf(stderr, "ERROR [sigIn]: File '%s' is empty or only contains comments.\n", input_filename);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    int i_read;
    double r_read;
    long double repmass_pop1_read;
    long double repmass_pop2_read;
    double max_part_size_read;
    double micro_size_read;

    for (int i = 0; i < disk_params->NGRID; i++) {
        if (fscanf(fp, "%d %lf %Lg %Lg %lf %lf",
                   &i_read, &r_read, &repmass_pop1_read, &repmass_pop2_read,
                   &max_part_size_read, &micro_size_read) != 6) {
            fprintf(stderr, "ERROR [sigIn]: Failed to read data for particle %d from file '%s'. Expected 6 values.\n", i, input_filename);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        // Adjust index for 1-based indexing in array (if NGRID is 0-based for loops)
        // Your code uses 1-based indexing for r_arr and disk_params->sigmavec
        if ((i_read + 1) >= 0 && (i_read + 1) <= disk_params->NGRID) {
             disk_params->rvec[i_read + 1] = r_read;
             disk_params->sigmavec[i_read + 1] = repmass_pop1_read + repmass_pop2_read;
        } else {
            fprintf(stderr, "WARNING [sigIn]: Skipping data for out-of-bounds index %d.\n", i_read);
        }
    }

    fclose(fp);
    fprintf(stderr, "DEBUG [sigIn]: Successfully loaded profile from %s.\n", input_filename);
}


void Mk_Dir(const char *dir_path) {
    fprintf(stderr, "DEBUG [Mk_Dir]: Attempting to create directory: '%s'\n", dir_path);
    int result = MKDIR_CALL(dir_path);

    if (result == 0) {
        fprintf(stderr, "DEBUG [Mk_Dir]: Directory '%s' created successfully.\n", dir_path);
    } else {
        // Hiba esetén ellenőrizzük, hogy már létezik-e a mappa
        #ifdef _WIN32
        // Windows alatt az _mkdir() -1-et ad vissza, ha már létezik, de nincs egyszerű GetLastError-os ellenőrzés itt
        // Egyszerűen csak figyelmeztetést adunk
        fprintf(stderr, "WARNING [Mk_Dir]: Failed to create directory '%s'. It might already exist or there's a permission issue.\n", dir_path);
        #else
        if (errno == EEXIST) {
            fprintf(stderr, "DEBUG [Mk_Dir]: Directory '%s' already exists.\n", dir_path);
        } else {
            perror("ERROR [Mk_Dir]: Failed to create directory"); // Kiírja a rendszerspecifikus hibaüzenetet
            fprintf(stderr, "ERROR [Mk_Dir]: Could not create directory: '%s'. Exiting.\n", dir_path);
            exit(1); // Kilépés súlyos hiba esetén
        }
        #endif
    }
    fflush(stderr); // Biztosítsuk, hogy a debug üzenetek azonnal kiíródjanak
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


/* Függvény a tömegfájl kiíratására */
void Print_Mass(double step, const double *rvec, double (*partmassind)[4], double (*partmassmicrind)[4],
                double (*partmasssecind)[4], const double *dpressvec,
                double massbtempii, double massbtempoi, double massmtempii, double massmtempoi,
                double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo,
                double *tavin, double *tavout,
                const disk_t *disk_params, const simulation_options_t *sim_opts,
                output_files_t *output_files) {

//    fprintf(stderr, "DEBUG [Print_Mass]: Entry. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
//            (void*)disk_params, disk_params->FLIND, disk_params->HASP);

    double ind_ii, ind_io, ind_oi, ind_oo, tav, tav2;

    tav = disk_params->r_dze_o;
    tav2 = disk_params->r_dze_i;

//    fprintf(stderr, "DEBUG [Print_Mass]: Before find_num_zero. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
//            (void*)disk_params, disk_params->FLIND, disk_params->HASP);
    int dim = find_num_zero(disk_params); // Assuming find_num_zero is const-correct
//    fprintf(stderr, "DEBUG [Print_Mass]: After find_num_zero. dim=%d. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
//            dim, (void*)disk_params, disk_params->FLIND, disk_params->HASP);

    // double r_count[dim]; // VLA - C99 standard. If compiling with C11 or later and not using GNU extensions, consider dynamic allocation.
    double *r_count = (double *)malloc(sizeof(double) * dim); // Safer for larger dim
//    fprintf(stderr, "DEBUG [Print_Mass]: After malloc for r_count (size %zu). r_count address=%p. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
//            sizeof(double) * dim, (void*)r_count, (void*)disk_params, disk_params->FLIND, disk_params->HASP);

    if (dim > 0 && r_count == NULL) {
        fprintf(stderr, "ERROR [Print_Mass]: Failed to allocate memory for r_count. Exiting.\n");
        exit(EXIT_FAILURE);
    }

    double temp_new = 0.;
    double temp = 0.;
    double rin_new = 0.0;
    double rout_new = 0.0;

    int j = 0, i;

    if(dim != 0) {
        for(i = 0; i < disk_params->NGRID; i++) { // Using disk_params->NGRID
            temp_new = find_zero(i,rvec,dpressvec); // Assuming find_zero is const-correct

            if(temp != temp_new && i > 3 && temp_new != 0.0) {
                if (j < dim) { // Prevent out-of-bounds write if dim calculation is off
                    r_count[j] = temp_new;
                    j++;
                } else {
                    fprintf(stderr, "WARNING [Print_Mass]: r_count array overflow, skipping data. dim: %d, j: %d\n", dim, j);
                }
            }

            if(sim_opts->dzone == 0.0) { // Using sim_opts->dzone
                if(temp_new > 0.) {
                    temp = temp_new;
                    rout_new = temp;
                }
            }
        }
    }
//    fprintf(stderr, "DEBUG [Print_Mass]: After r_count population loop. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
//            (void*)disk_params, disk_params->FLIND, disk_params->HASP);


    if(sim_opts->dzone == 1.0) { // Using sim_opts->dzone
        if(dim > 0) {
            if (dim == 1) {
                rin_new = r_count[0];
                rout_new = tav; // Using tav (disk_params->r_dze_o)
            } else if (dim >= 2) { // Ensure there are at least two elements
                rin_new = r_count[0];
                rout_new = r_count[1];
            } else { // Should not happen if dim > 0 and not 1
                fprintf(stderr, "WARNING [Print_Mass]: Unexpected dim value %d for sim_opts->dzone == 1.0. Using default r_dze_i/o.\n", dim);
                rin_new = tav2;
                rout_new = tav;
            }
        }
        if(dim == 0) {
            rin_new = tav2; // Using tav2 (disk_params->r_dze_i)
            rout_new = tav; // Using tav (disk_params->r_dze_o)
        }
    }
//    fprintf(stderr, "DEBUG [Print_Mass]: After dzone logic. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
//            (void*)disk_params, disk_params->FLIND, disk_params->HASP);


    // The rin and rout variables are now consistently updated
    double rin_current = rin_new;
    if(sim_opts->dzone == 0.0) rin_current = 0;
    double rout_current = rout_new;
    *tavin = rin_current; // Update output parameter directly
    *tavout = rout_current; // Update output parameter directly

//    fprintf(stderr, "DEBUG [Print_Mass]: Before find_r_annulus call. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
//            (void*)disk_params, disk_params->FLIND, disk_params->HASP);
    // Passing updated tav2 and tav to find_r_annulus
    find_r_annulus(disk_params->rvec,*tavin,&ind_ii,&ind_io,*tavout,&ind_oi,&ind_oo,sim_opts,disk_params);
//    fprintf(stderr, "DEBUG [Print_Mass]: After find_r_annulus call. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
//            (void*)disk_params, disk_params->FLIND, disk_params->HASP);


    double massii = 0, massoi = 0;
    double massiim = 0,massoim = 0;
    double massis = 0, massos = 0;

    // The crash occurred on the *first* call to scale_height which is within find_r_annulus.
    // So the problem is likely before or in find_r_annulus.
    // The following calls are not reached during the crash.

    GetMass(PARTICLE_NUMBER,partmassind,(int)ind_ii,(int)ind_io,*tavin,disk_params->r_dze_i,&massii,(int)ind_oi,(int)ind_oo,*tavout,disk_params->r_dze_o,&massoi, sim_opts);
    if(sim_opts->twopop == 1.0) {
        GetMass(4*PARTICLE_NUMBER,partmasssecind,(int)ind_ii,(int)ind_io,*tavin,disk_params->r_dze_i,&massis,(int)ind_oi,(int)ind_oo,*tavout,disk_params->r_dze_o,&massos,sim_opts);
        GetMass(PARTICLE_NUMBER,partmassmicrind,(int)ind_ii,(int)ind_io,*tavin,disk_params->r_dze_i,&massiim,(int)ind_oi,(int)ind_oo,*tavout,disk_params->r_dze_o,&massoim,sim_opts);
    }

    double massi, massim, masso, massom;

    if(*tavin != disk_params->r_dze_i) {
        massi = massii + massbtempii + massis;
        massim = massiim + massmtempii;
    } else {
        massi = massii + massis;
        massim = massiim;
    }
    if(*tavout != disk_params->r_dze_o) {
        masso = massoi + massbtempoi + massos;
        massom = massoim + massmtempoi;
    } else {
        masso = massoi + massos;
        massom = massoim;
    }

    *massbtempio = massi;
    *massbtempoo = masso;
    *massmtempio = massim; // FIX: Corrected typo from massmtempiout to massmtempio
    *massmtempoo = massom;

    if (output_files->mass_file != NULL) {
        fprintf(output_files->mass_file, "%lg %lg %lg %lg %lg\n", step, *tavin, massi + massim, *tavout, masso + massom);
    } else {
        fprintf(stderr, "WARNING: output_files->mass_file is NULL in Print_Mass. Cannot write mass data.\n");
    }

    if (output_files->mass_file != NULL) {
        fflush(output_files->mass_file);
    } else {
        fprintf(stderr, "WARNING: Cannot fflush output_files->mass_file, as it is NULL.\n");
    }

    // Free dynamically allocated memory for r_count if it was allocated
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

    for(i = 1; i <= disk_params->NGRID; i++) { // Using disk_params->NGRID
        fprintf(output_files->surface_file, "%lg %lg %lg %lg\n", disk_params->rvec[i], disk_params->sigmavec[i], disk_params->pressvec[i], disk_params->dpressvec[i]);
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