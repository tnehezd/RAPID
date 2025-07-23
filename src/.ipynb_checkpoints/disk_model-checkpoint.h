// src/disk_model.h

#ifndef DISK_MODEL_H
#define DISK_MODEL_H


// Azért include-olom a config.h-t, mert a disk_param_be globális változókat használ
// (pl. SUN2GR, AU2CM), és a Parabola, load_R, Initial_Profile is NGRID, RMIN, DD, SIGMA0, SIGMAP_EXP-et.
// Bár ez nem szigorúan szükséges, mert a disk_model.c direktben include-olja a config.h-t,
// jó gyakorlat, ha egy header is jelzi a függőségeit, ha a benne lévő deklarációk függenek tőlük.
#include "config.h"
#include "simulation_types.h" // Feltételezve, hogy itt található a disk_t definíciója

// Funkciódeklarációk a disk_model.c-ből
void disk_param_be(disk_t *disk_params);

void load_R(disk_t *disk_params);

void Initial_Profile(disk_t *disk_params);

void Initial_Press(disk_t *disk_params);

void Initial_dPress(disk_t *disk_params);

void Initial_Ugas(disk_t *disk_params);

// --- MÓDOSÍTVA: A loadSigDust függvény deklarációja long double paraméterekkel ---
void loadSigDust(double radin[][2], long double *massin, long double out[][3], int n, const disk_t *disk_params);
// --- MÓDOSÍTÁS VÉGE ---

#endif // DISK_MODEL_H