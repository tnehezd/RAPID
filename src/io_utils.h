#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h>
#include <stdbool.h>

#include "simulation_types.h"
#include "particle_data.h"
#include "config.h"
#include "globals.h"

// Globális változók deklarációi, ha az io_utils.c fájlban definiálva vannak.
// Ezeknek EGYEZNIÜK KELL a src/config.h-ban lévő extern deklarációkkal.
extern FILE *fin1; // Valószínűleg a régi ReadDustFile használja

// --- FÜGGVÉNY DEKLARÁCIÓK ---

/* reszecskek_szama függvény deklaráció */
int get_particle_count(const char *filenev);


/* A porreszecskek adatainak beolvasasa (ÚJ verzió, ParticleData_t struktúrával) */
// Ez a függvény tölti fel a ParticleData_t struktúrában lévő adatokat.
// A fájlnevet kapja paraméterül, és a beolvasott részecskeszámot is felhasználhatja (pl. a reszecskek_szama függvényből).
// A disk_params és sim_opts paraméterek opcionálisak lehetnek, de jó, ha a függvény kontextust kap a szimulációról.
void load_dust_particles(ParticleData_t *p_data, const char *filename,
                     const disk_t *disk_params, simulation_options_t *sim_opts); // REMOVED CONST

/* A sigmat tartalmazo file parametereinek beolvasasa */
void read_gas_profile(disk_t *disk_params, const char *filename);

/* Fuggveny az adott futashoz mappa letrehozasara */
void create_output_directory(char *dir_path); // Megfontolandó: const char* dir_path, ha nem módosítja

/* Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok */
// FIX: The original had 'void write_summary_log(const char *nev);'
// You are missing 'const disk_t *disk_params' and 'const simulation_options_t *sim_opts'.
void write_summary_log(const char *nev, const disk_t *disk_params, const simulation_options_t *sim_opts);


/* Fuggveny a sigma, p, dp kiiratasara */
// FIX: The original was missing 'const disk_t *disk_params'.
void write_gas_profile_to_file(const disk_t *disk_params, output_files_t *output_files);


void write_dust_profile_to_file(int step,
                  const ParticleData_t *p_data, // Itt adjuk át a ParticleData_t struktúrát
                  const disk_t *disk_params,
                  const simulation_options_t *sim_opts,
                  output_files_t *output_files);



// Enumeráció a fájltípusok azonosítására
typedef enum {
    FILE_TYPE_MASS_ACCUMULATION,
    FILE_TYPE_GAS_DENSITY,
    FILE_TYPE_DUST_DENSITY,
    FILE_TYPE_DUST_MICRON_DENSITY,
    FILE_TYPE_PARTICLE_SIZE,
    FILE_TYPE_DISK_PARAM,
    FILE_TYPE_TIMESCALE,
    FILE_TYPE_MICRON_DUST_EVOL,
    FILE_TYPE_DUST_EVOL

} FileType_e;

// Struktúra a fejléc-specifikus adatoknak
typedef struct {
    double current_time;     // Jelenlegi szimulációs idő (pl. években)
    int is_initial_data;    // 1, ha t=0, 0, ha szimulált időpont
    // Ide tehetsz más adatokat is, ami a fejléchez kellhet, pl. R_in, R_out
    double R_in;
    double R_out;
    double sigma_exponent; // Kellhet a disk_param fejlécbe
    long double sigma0_gas_au; // Kellhet a disk_param fejlécbe
    double grav_const; // Kellhet a disk_param fejlécbe
    double dz_r_inner; // Kellhet a disk_param fejlécbe
    double dz_r_outer; // Kellhet a disk_param fejlécbe
    double dz_dr_inner_calc; // Kellhet a disk_param fejlécbe
    double dz_dr_outer_calc; // Kellhet a disk_param fejlécbe
    double dz_alpha_mod; // Kellhet a disk_param fejlécbe
    double dust_density_g_cm3; // Kellhet a disk_param fejlécbe
    double alpha_viscosity; // Kellhet a disk_param fejlécbe
    double star_mass; // Kellhet a disk_param fejlécbe
    double flaring_index; // Kellhet a disk_param fejlécbe
    int n_grid_points; // Kellhet a disk_param fejlécbe
} HeaderData_t;


// Függvény a fejlécek kiírására
// Az 'header_data' opcionális lehet (NULL is átadható), ha az adott fájltípushoz nem kell
void write_file_header(FILE *file, FileType_e file_type, const HeaderData_t *header_data);


// Függvény a kezdeti kimeneti fájlok beállítására és fejlécek írására
int initialize_mass_accumulation_file(output_files_t *output_files, const simulation_options_t *sim_opts,
                                 const disk_t *disk_params, HeaderData_t *header_data_for_files);


void cleanup_simulation_resources(ParticleData_t *p_data, output_files_t *output_files, const simulation_options_t *sim_opts);

void close_snapshot_files(output_files_t *output_files, const char *dens_name, const char *dust_name, const char *dust_name2, const simulation_options_t *sim_opts);


#endif // IO_UTILS_H