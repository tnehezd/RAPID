// src/disk_model.h

#ifndef DISK_MODEL_H
#define DISK_MODEL_H

// Azért include-olom a config.h-t, mert a disk_param_be globális változókat használ
// (pl. SUN2GR, AU2CM), és a Parabola, load_R, Initial_Profile is NGRID, RMIN, DD, SIGMA0, SIGMAP_EXP-et.
// Bár ez nem szigorúan szükséges, mert a disk_model.c direktben include-olja a config.h-t,
// jó gyakorlat, ha egy header is jelzi a függőségeit, ha a benne lévő deklarációk függenek tőlük.
#include "config.h"

// Funkciódeklarációk a disk_model.c-ből
void disk_param_be(double *sigma0, double *sdexp, double *Rmin, double *Rmax,
                   double *r_dzei, double *r_dzeo, double *dr_dzei, double *dr_dzeo,
                   double *alph_mod, double *rho_p, double *rho_p_dimless,
                   double *alphav, double *mStar, double *gamma);

void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd);

void load_R(double *rvec);

void Initial_Profile(double *sigmavec, double *r);

void Initial_Press(double *pressvec, double *sigmavec, double *rvec);

void Initial_dPress(double *dpressvec, double *pressvec);

void Initial_Ugas(double *sigmavec, double *rvec, double *ug);

// A Perem függvényt itt deklaráljuk, mert a disk_model.c-ben van definiálva
// és más függvények (pl. Initial_Profile) is hívják.
// Fontos: Ha a Perem az io_utils.c-ben van definiálva, akkor az io_utils.h-ba kell tenni!
// A korábbi fordítási kimenet alapján a disk_model.c hívja a Perem-et.
// Feltételezem, hogy a Perem definíciója a disk_model.c-ben van, de nem mutattad meg a teljes fájlt.
// Ha nincs a disk_model.c-ben, akkor egy másik .c fájlban kell lennie, és annak a headerében kell lennie.
// De a "disk_model.c: In function 'Initial_Profile': src/disk_model.c:88:9: warning: implicit declaration of function 'Perem'"
// alapján a Perem-et is includolni kell valahonnan, vagy ide tenni, ha a disk_model.c-ben van.
// Mivel a `Perem` függvény nincs a most megadott `disk_model.c` részletben, de a hibák szerint ott van hívva,
// feltételezem, hogy a `disk_model.c` egy másik részén található.
// Ha a Perem definíciója *mégis* az `io_utils.c`-ben van, akkor a `disk_model.c`-nek include-olnia kell az `io_utils.h`-t!
// Most ide teszem, és ha mégsem itt van, majd módosítjuk.
void Perem(double *vec);


#endif // DISK_MODEL_H
