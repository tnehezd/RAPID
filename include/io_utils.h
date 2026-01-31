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

/* calculateNumbersOfParticles függvény deklaráció */
int calculateNumbersOfParticles(const char *filenev);

/* A porreszecskek adatainak beolvasasa */
// FIX: The original had 'void loadDustParticlesFromFile(double radius[][2], double radiusmicr[][2], double *mass, double *massmicr);'
// You are missing the 'const char *filename' parameter in the .h file.
void loadDustParticlesFromFile(ParticleData *particle_data, const char *filename);

/* A sigmat tartalmazo file parametereinek beolvasasa */
void loadGasSurfaceDensityFromFile(DiskParameters *disk_params, const char *filename);

/* Fuggveny az adott futashoz mappa letrehozasara */
void createRunDirectory(char *dir_path);

/* Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok */
// FIX: The original had 'void printCurrentInformationAboutRun(const char *nev);'
// You are missing 'const DiskParameters *disk_params' and 'const SimulationOptions *sim_opts'.
void printCurrentInformationAboutRun(const char *nev, const DiskParameters *disk_params, const SimulationOptions *sim_opts);

/* Fuggveny a tomegfile kiiratasara */
// FIX: The original was missing 'const DiskParameters *disk_params' and 'const SimulationOptions *sim_opts'.
void printMassGrowthAtDZEFile(double step, 
                double (*partmassind)[5], double (*partmassmicrind)[5], 
                double t, // Ezt továbbra is meghagyjuk, ha az időre szükség van
                double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, 
                double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, 
                double *tavin, double *tavout, 
                const DiskParameters *disk_params, const SimulationOptions *sim_opts,
                OutputFiles *output_files);

/* Fuggveny a sigma, p, dp kiiratasara */
// FIX: The original was missing 'const DiskParameters *disk_params'.
void printGasSurfaceDensityPressurePressureDerivateFile(const DiskParameters *disk_params, OutputFiles *output_files);

/* Fuggveny a por feluletisurusegenek kiiratasara */
// FIX: The original was missing 'const DiskParameters *disk_params' and 'const SimulationOptions *sim_opts'.
void printDustSurfaceDensityPressurePressureDerivateFile(const double *r, const double *rm, const double *sigmad, const double *sigmadm,
                  const DiskParameters *disk_params, const SimulationOptions *sim_opts,
                  OutputFiles *output_files, double step);
/* Fuggveny a pormozgas es reszecskemeret kiiratasara */
// FIX: The original was missing 'const DiskParameters *disk_params' and 'const SimulationOptions *sim_opts'.
void printDustParticleSizeFile(char *size_name, int step, double (*rad)[2], double (*radmicr)[2],
                        const DiskParameters *disk_params, const SimulationOptions *sim_opts,
                        OutputFiles *output_files);




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
    double current_time;    // Jelenlegi szimulációs idő (pl. években)
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
} HeaderData;


// Függvény a fejlécek kiírására
// Az 'header_data' opcionális lehet (NULL is átadható), ha az adott fájltípushoz nem kell
void printFileHeader(FILE *file, FileType_e file_type, const HeaderData *header_data);


// Függvény a kezdeti kimeneti fájlok beállítására és fejlécek írására
int setupInitialOutputFiles(OutputFiles *output_files, const SimulationOptions *sim_opts,
                               const DiskParameters *disk_params, HeaderData *header_data_for_files);


void cleanupSimulationResources(ParticleData *particle_data, OutputFiles *output_files, const SimulationOptions *sim_opts);

FILE *openSnapshotFile(const char *filename,FileType_e file_type,double current_time_years);

void closeSnapshotFiles(OutputFiles *output_files, const char *dens_name, const char *dust_name, const char *dust_name2, const SimulationOptions *sim_opts);


void buildSnapshotFilenames(char *dens_name, char *dust_name, char *dust_name2, char *size_name, const SimulationOptions *sim_opts, int snapshot_id);


#endif // IO_UTILS_H

