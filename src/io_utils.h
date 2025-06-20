#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h> // For FILE* and printf/fprintf
#include <stdbool.h> // For bool type if you decide to use it (optional)

// Globális változók deklarációi, ha az io_utils.c fájlban definiálva vannak.
// Ezeknek EGYEZNIÜK KELL a src/config.h-ban lévő extern deklarációkkal.
// Feltételezem, hogy ezeket a config.h-ból kapod meg, így itt csak a specifikus io_utils fájl stream-ek vannak.
extern FILE *fin1, *fin2, *fmo, *fout, *foutmicr, *fout3, *massfil, *jelfut;

// --- FÜGGVÉNY DEKLARÁCIÓK ---

/* A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut */
// Feltehetően a disk_model.h tartalmazza a scale_height deklarációját.
void find_r_annulus(double *rvec, double rin_val, double *ind_ii, double *ind_io, double rout_val, double *ind_oi, double *ind_oo);

/* reszecskek_szama függvény deklaráció */
// Az numout-ot már nem adjuk át referenciával, hanem visszatérési érték.
// A const char *filenev javasolt, ha nem módosítod a stringet.
int reszecskek_szama(int numout_dummy, const char *filenev); // numout_dummy, mert a funkció maga adja vissza az értéket, nem paraméterként veszi át

/* A porreszecskek adatainak beolvasasa */
// Az unused warnings miatt a particle_radius és radmicr paramétereket töröltem,
// ha tényleg nem használod. Ha igen, akkor tedd vissza.
void por_be(double radius[][2], double radiusmicr[][2], double *mass, double *massmicr);

/* A sigmat tartalmazo file parametereinek beolvasasa */
// A sigvec és rvec paraméterek maradtak, mert feltételezem, hogy ezekbe töltöd az adatokat.
void sigIn(double *sigvec, double *rvec);

/* Fuggveny az adott futashoz mappa letrehozasara */
void Mk_Dir(char *nev);

/* Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok */
void infoCurrent(char *nev);

/* Fuggveny a tomegfile kiiratasara */
// A "t" paramétert töröltem, mert unused warning volt, ha kell, tedd vissza.
// A massbtempio, massbtempoo, massmtempio, massmtempoo, tavin, tavout most már mutatók,
// ezért a függvénydeklarációban is annak kell lenniük (ahogy a .c fájlban is láttam).
// In src/io_utils.h:
void Print_Mass(double step, double *rvec, double partmassind[][4], double partmassmicrind[][4], double partmasssecind[][4], double t, double *dpressvec, double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, double *tavin, double *tavout);
/* Fuggveny a sigma, p, dp kiiratasara */
void Print_Sigma(char *dens_name, double *rvec, double *sigmavec, double *pressvec, double *dpressvec);

/* Fuggveny a por feluletisurusegenek kiiratasara */
// A 'min' paramétert eltávolítottam, mert unused warning volt, ha kell, tedd vissza.
// A 'r', 'rm', 'sigmad', 'sigmadm' típusát pontosítottam double * -ra
// a .c fájlban látott tömbös használat miatt (azaz mutatók).
void Print_Sigmad(char *dust_name, char *dust_name2, double min, double *r, double *rm, double *sigmad, double *sigmadm);

/* Fuggveny a pormozgas kiiratasara */
// A 'size_name', 'radmicr', 'rvec', 't' paramétereket eltávolítottam, mert unused warningok voltak.
// Ha szükséged van rájuk, tedd vissza őket!
void Print_Pormozg_Size(int step, double rad[][2]);

/* Az idot tartalmazo file parametereinek beolvasasa */
void timePar(double *tMax, double *step, double *current);

#endif // IO_UTILS_H
