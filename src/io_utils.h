// include/io_utils.h

#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h>   // For FILE* and printf/fprintf
#include <stdbool.h> // For bool type if you decide to use it (optional)

#include "simulation_types.h" // This is correctly included now

// Globális változók deklarációi, ha az io_utils.c fájlban definiálva vannak.
// Ezeknek EGYEZNIÜK KELL a src/config.h-ban lévő extern deklarációkkal.
// Feltételezem, hogy ezeket a config.h-ból kapod meg, így itt csak a specifikus io_utils fájl stream-ek vannak.
extern FILE *fin1, *fin2, *fmo, *fout, *foutmicr, *fout3, *massfil, *jelfut;

// --- FÜGGVÉNY DEKLARÁCIÓK ---

/* reszecskek_szama függvény deklaráció */
// Az numout-ot már nem adjuk át referenciával, hanem visszatérési érték.
// A const char *filenev javasolt, ha nem módosítod a stringet.
int reszecskek_szama(int numout_dummy, const char *filenev); // numout_dummy, mert a funkció maga adja vissza az értéket, nem paraméterként veszi át

/* A porreszecskek adatainak beolvasasa */
// ha tényleg nem használod. Ha igen, akkor tedd vissza.
void por_be(double radius[][2], double radiusmicr[][2], double *mass, double *massmicr);

/* A sigmat tartalmazo file parametereinek beolvasasa */
// A sigvec és rvec paraméterek maradtak, mert feltételezem, hogy ezekbe töltöd az adatokat.
// --- FIX: Add 'const char *filename' here to match the call in main.c ---
void sigIn(double *sigma_arr, double *r_arr, const disk_t *disk_params, const char *filename);

/* Fuggveny az adott futashoz mappa letrehozasara */
// --- FIX: Change to 'const char *nev' ---
void Mk_Dir(char *nev);

/* Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok */
// --- FIX: Change to 'const char *nev' ---
void infoCurrent(const char *nev);

/* Fuggveny a tomegfile kiiratasara */
// A massbtempio, massbtempoo, massmtempio, massmtempoo, tavin, tavout most már mutatók,
// ezért a függvénydeklarációban is annak kell lenniük (ahogy a .c fájlban is láttam).
void Print_Mass(double step, double *rvec, double partmassind[][4], double partmassmicrind[][4], double partmasssecind[][4], double *dpressvec, double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, double *tavin, double *tavout);
/* Fuggveny a sigma, p, dp kiiratasara */
// --- Consider changing 'char *dens_name' to 'const char *dens_name' if it's not modified ---
void Print_Sigma(char *dens_name, double *rvec, double *sigmavec, double *pressvec, double *dpressvec);

/* Fuggveny a por feluletisurusegenek kiiratasara */
// --- Consider changing 'char *dust_name, char *dust_name2' to 'const char *' if not modified ---
void Print_Sigmad(char *dust_name, char *dust_name2, double *r, double *rm, double *sigmad, double *sigmadm);

/* Fuggveny a pormozgas kiiratasara */
// --- Consider changing 'char *size_name' to 'const char *' if not modified ---
void Print_Pormozg_Size(char *size_name, int step, double rad[][2], double radmicr[][2]);

/* Az idot tartalmazo file parametereinek beolvasasa */
void timePar(double tMax_val, double stepping_val, double current_val); // Ez kell az io_utils.h-ba
#endif // IO_UTILS_H