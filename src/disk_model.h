// src/disk_model.h

#ifndef DISK_MODEL_H
#define DISK_MODEL_H


#include "config.h"
#include "simulation_types.h" // Feltételezve, hogy itt található a disk_t definíciója
#include "globals.h"

// Funkciódeklarációk a disk_model.c-ből
void read_disk_parameters(disk_t *disk_params);

void load_R(disk_t *disk_params);

void Initial_Profile(disk_t *disk_params);

void Initial_Press(disk_t *disk_params);

void Initial_dPress(disk_t *disk_params);

void Initial_Ugas(disk_t *disk_params);

void loadSigDust(double radin[][2], double *massin, double out[][3], int n, const disk_t *disk_params);

#endif // DISK_MODEL_H
