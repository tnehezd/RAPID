// src/disk_model.h

#ifndef DISK_MODEL_H
#define DISK_MODEL_H


// Azért include-olom a config.h-t, mert a readDiskParameters globális változókat használ
// (pl. SUN2GR, AU2CM), és a Parabola, createRadialGrid, createInitialGasSurfaceDensity is NGRID, RMIN, DD, SIGMA0, SIGMAP_EXP-et.
// Bár ez nem szigorúan szükséges, mert a disk_model.c direktben include-olja a config.h-t,
// jó gyakorlat, ha egy header is jelzi a függőségeit, ha a benne lévő deklarációk függenek tőlük.
#include "config.h"
#include "simulation_types.h" // Feltételezve, hogy itt található a disk_t definíciója

// Funkciódeklarációk a disk_model.c-ből
void readDiskParameters(disk_t *disk_params);

void createRadialGrid(disk_t *disk_params);

void createInitialGasSurfaceDensity(disk_t *disk_params);

void Initial_Press(disk_t *disk_params);

void Initial_dPress(disk_t *disk_params);

void Initial_Ugas(disk_t *disk_params);

void loadSigDust(double radin[][2], double *massin, double out[][3], int n, const disk_t *disk_params);

#endif // DISK_MODEL_H
