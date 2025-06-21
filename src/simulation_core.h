// src/simulation_core.h

#ifndef SIMULATION_CORE_H
#define SIMULATION_CORE_H

// Include any headers needed *by this header itself* for types used in declarations.
// For tIntegrate, it uses 'char' and 'double', which are built-in, so no specific
// includes are strictly necessary *just for the prototype*. However, if this header
// were to define structs or other types that are part of the public interface,
// you'd include their definitions here.

// 	1D drfit: dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs
void eqrhs(double pradius, double dp, double sigma, double ug, double r, double *drdt);


// for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 3*nu = Coeff_1 
double Coeff_1(double r);

// for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 9*nu /(2*r) = Coeff_2
double Coeff_2(double r);

/*	Runge-Kutta4 integrator	*/
// prad bemenet: AU-ban!
void int_step(double time, double prad, double *pressvec, double *dpressvec, double *sigmavec, double *sigmad, double *rdvec, double *rvec, double *ugvec, double step, double y, double *ynew, double *pradnew);


// a diffúziós egyenlet megoldásához szükséges minimális időlépés
double time_step(double *rvec);

// Function Declarations
void tIntegrate(const char *nev, const disk_t *disk_params);

void secondaryGrowth(double rad[][2], double radmicr[][2], double radsec[][2], double partmicind[][4], double partsecind[][4], double *massmicvec, double *masssecvec);

// If you add other functions to simulation_core.c that are called elsewhere,
// declare them here as well.

#endif // SIMULATION_CORE_H
