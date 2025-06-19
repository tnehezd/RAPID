// src/dust_physics.h
#ifndef DUST_PHYSICS_H
#define DUST_PHYSICS_H

double press(double sigma, double r);
void dpress(double *dpressvec, double *pressvec);
void u_gas(double *sigmavec, double *rvec, double *ug);
// ... (egyéb függvények, pl. GetMass) ...

#endif
