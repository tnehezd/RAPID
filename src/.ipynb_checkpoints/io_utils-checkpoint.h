#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h>
#include <stdbool.h>

#include "simulation_types.h"
#include "particle_data.h"

// Globális változók deklarációi, ha az io_utils.c fájlban definiálva vannak.
// Ezeknek EGYEZNIÜK KELL a src/config.h-ban lévő extern deklarációkkal.
extern FILE *fin1, *fin2, *fmo, *fout, *foutmicr, *fout3, *massfil, *jelfut;

// --- FÜGGVÉNY DEKLARÁCIÓK ---

/* reszecskek_szama függvény deklaráció */
int reszecskek_szama(const char *filenev);

/* A porreszecskek adatainak beolvasasa */
// FIX: The original had 'void por_be(double radius[][2], double radiusmicr[][2], double *mass, double *massmicr);'
// You are missing the 'const char *filename' parameter in the .h file.
void por_be(long double radius[][2], long double radiusmicr[][2], long double *mass, long double *massmicr, const char *filename); // MÓDOSÍTVA: double -> long double

/* A sigmat tartalmazo file parametereinek beolvasasa */
void sigIn(disk_t *disk_params, const char *filename);

int compareLongDoubles(const void *a, const void *b); // MÓDOSÍTVA: compareDoubles -> compareLongDoubles

/* Fuggveny az adott futashoz mappa letrehozasara */
void Mk_Dir(char *dir_path);

/* Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok */
// FIX: The original had 'void infoCurrent(const char *nev);'
// You are missing 'const disk_t *disk_params' and 'const simulation_options_t *sim_opts'.
void infoCurrent(const char *nev, const disk_t *disk_params, const simulation_options_t *sim_opts);

/* Fuggveny a tomegfile kiiratasara */
// FIX: The original was missing 'const disk_t *disk_params' and 'const simulation_options_t *sim_opts'.
void Print_Mass(long double step, long double (*partmassind)[5], long double (*partmassmicrind)[5], long double t, // MÓDOSÍTVA: double -> long double
                long double massbtempii, long double massbtempoi, long double massmtempii, long double massmtempoi, // MÓDOSÍTVA: double -> long double
                long double *massbtempio, long double *massbtempoo, long double *massmtempio, long double *massmtempoo, // MÓDOSÍTVA: double -> long double
                long double masstempoibuff, long double *masstempoobuff, long double massmtempoibuff, long double *massmtempoobuff, // MÓDOSÍTVA: double -> long double
                long double *tavin, long double *tavout, long double *tavbuff, // MÓDOSÍTVA: double -> long double
                const disk_t *disk_params, const simulation_options_t *sim_opts, output_files_t *output_files);

/* Fuggveny a sigma, p, dp kiiratasara */
// FIX: The original was missing 'const disk_t *disk_params'.
void Print_Sigma(const disk_t *disk_params, output_files_t *output_files);

/* Fuggveny a por feluletisurusegenek kiiratasara */
// FIX: The original was missing 'const disk_t *disk_params' and 'const simulation_options_t *sim_opts'.
void Print_Sigmad(const long double *r, const long double *rm, const long double *sigmad, const long double *sigmadm, // MÓDOSÍTVA: double -> long double
                  const disk_t *disk_params, const simulation_options_t *sim_opts,
                  output_files_t *output_files);
/* Fuggveny a pormozgas es reszecskemeret kiiratasara */
// FIX: The original was missing 'const disk_t *disk_params' and 'const simulation_options_t *sim_opts'.
void Print_Pormozg_Size(char *size_name, long double step, long double (*rad)[2], long double (*radmicr)[2], // MÓDOSÍTVA: int step -> long double step, double -> long double
                        const disk_t *disk_params, const simulation_options_t *sim_opts,
                        output_files_t *output_files);

/* Az idot tartalmazo file parametereinek beolvasasa (vagy beallitasa) */
// FIX: The original was missing 'simulation_options_t *sim_opts'.
void timePar(long double tMax_val, long double stepping_val, long double current_val, simulation_options_t *sim_opts); // MÓDOSÍTVA: double -> long double

// Enumeráció a fájltípusok azonosítására
typedef enum {
    FILE_TYPE_DUST_MOTION,
    FILE_TYPE_MICRON_MOTION,
    FILE_TYPE_MASS_ACCUMULATION,
    FILE_TYPE_GAS_DENSITY,
    FILE_TYPE_DUST_DENSITY,
    FILE_TYPE_DUST_MICRON_DENSITY,
    FILE_TYPE_PARTICLE_SIZE,
    FILE_TYPE_DISK_PARAM // ÚJ: a disk_config.dat fájlhoz

} FileType_e;

// Struktúra a fejléc-specifikus adatoknak
typedef struct {
    long double current_time;  // Jelenlegi szimulációs idő (pl. években) // MÓDOSÍTVA: double -> long double
    int is_initial_data;  // 1, ha t=0, 0, ha szimulált időpont
    // Ide tehetsz más adatokat is, ami a fejléchez kellhet, pl. R_in, R_out
    long double R_in; // MÓDOSÍTVA: double -> long double
    long double R_out; // MÓDOSÍTVA: double -> long double
    long double sigma_exponent; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double sigma0_gas_au; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double grav_const; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double dz_r_inner; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double dz_r_outer; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double dz_dr_inner_calc; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double dz_dr_outer_calc; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double dz_alpha_mod; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double dust_density_g_cm3; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double alpha_viscosity; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double star_mass; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    long double flaring_index; // Kellhet a disk_param fejlécbe // MÓDOSÍTVA: double -> long double
    int n_grid_points; // Kellhet a disk_param fejlécbe
} HeaderData_t;

// Függvény a fejlécek kiírására
// Az 'header_data' opcionális lehet (NULL is átadható), ha az adott fájltípushoz nem kell
void print_file_header(FILE *file, FileType_e file_type, const HeaderData_t *header_data);

// Függvény a kezdeti kimeneti fájlok beállítására és fejlécek írására
int setup_initial_output_files(output_files_t *output_files, const simulation_options_t *sim_opts,
                               const disk_t *disk_params, HeaderData_t *header_data_for_files);

void cleanup_simulation_resources(ParticleData_t *p_data, output_files_t *output_files, const simulation_options_t *sim_opts);

void close_snapshot_files(output_files_t *output_files, const char *dens_name, const char *dust_name, const char *dust_name2, const simulation_options_t *sim_opts);

#endif // IO_UTILS_H