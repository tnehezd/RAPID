#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "simulation_types.h" 

// Függvény deklarációk (prototípusok)
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd, const disk_t *disk_params);
void Perem(double *vec,const disk_t *disk_params);

#endif // BOUNDARY_CONDITIONS_H
