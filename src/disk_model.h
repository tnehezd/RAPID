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

void load_R(double *rvec);

void Initial_Profile(double *sigmavec, double *r);

void Initial_Press(double *pressvec, double *sigmavec, double *rvec);

void Initial_dPress(double *dpressvec, double *pressvec);

void Initial_Ugas(double *sigmavec, double *rvec, double *ug);

void loadSigDust(double radin[][2], double *massin, double out[][3], double dd, int n);

#endif // DISK_MODEL_H
