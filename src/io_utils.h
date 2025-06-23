#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h>
#include <stdbool.h>

#include "simulation_types.h"

// Globális változók deklarációi, ha az io_utils.c fájlban definiálva vannak.
// Ezeknek EGYEZNIÜK KELL a src/config.h-ban lévő extern deklarációkkal.
extern FILE *fin1, *fin2, *fmo, *fout, *foutmicr, *fout3, *massfil, *jelfut;

// --- FÜGGVÉNY DEKLARÁCIÓK ---

/* reszecskek_szama függvény deklaráció */
int reszecskek_szama(const char *filenev);

/* A porreszecskek adatainak beolvasasa */
// FIX: The original had 'void por_be(double radius[][2], double radiusmicr[][2], double *mass, double *massmicr);'
// You are missing the 'const char *filename' parameter in the .h file.
void por_be(double radius[][2], double radiusmicr[][2], double *mass, double *massmicr, const char *filename);

/* A sigmat tartalmazo file parametereinek beolvasasa */
void sigIn(disk_t *disk_params, const char *filename);

/* Fuggveny az adott futashoz mappa letrehozasara */
void Mk_Dir(const char *dir_path);

/* Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok */
// FIX: The original had 'void infoCurrent(const char *nev);'
// You are missing 'const disk_t *disk_params' and 'const simulation_options_t *sim_opts'.
void infoCurrent(const char *nev, const disk_t *disk_params, const simulation_options_t *sim_opts);

/* Fuggveny a tomegfile kiiratasara */
// FIX: The original was missing 'const disk_t *disk_params' and 'const simulation_options_t *sim_opts'.
void Print_Mass(double step, const double *rvec, double (*partmassind)[4], double (*partmassmicrind)[4],
                double (*partmasssecind)[4], const double *dpressvec,
                double massbtempii, double massbtempoi, double massmtempii, double massmtempoi,
                double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo,
                double *tavin, double *tavout,
                const disk_t *disk_params, const simulation_options_t *sim_opts,
                output_files_t *output_files);

/* Fuggveny a sigma, p, dp kiiratasara */
// FIX: The original was missing 'const disk_t *disk_params'.
void Print_Sigma(const disk_t *disk_params, output_files_t *output_files);

/* Fuggveny a por feluletisurusegenek kiiratasara */
// FIX: The original was missing 'const disk_t *disk_params' and 'const simulation_options_t *sim_opts'.
void Print_Sigmad(const double *r, const double *rm, const double *sigmad, const double *sigmadm,
                  const disk_t *disk_params, const simulation_options_t *sim_opts,
                  output_files_t *output_files);
/* Fuggveny a pormozgas es reszecskemeret kiiratasara */
// FIX: The original was missing 'const disk_t *disk_params' and 'const simulation_options_t *sim_opts'.
void Print_Pormozg_Size(char *size_name, int step, double (*rad)[2], double (*radmicr)[2],
                        const disk_t *disk_params, const simulation_options_t *sim_opts,
                        output_files_t *output_files);

/* Az idot tartalmazo file parametereinek beolvasasa (vagy beallitasa) */
// FIX: The original was missing 'simulation_options_t *sim_opts'.
void timePar(double tMax_val, double stepping_val, double current_val, simulation_options_t *sim_opts);

#endif // IO_UTILS_H