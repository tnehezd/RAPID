// src/config.h

#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> // For FILE* and other standard definitions
#include <math.h>  // For tanh, sqrt, pow, and where M_PI is typically defined

// --- Fixed Physical Constants (defined directly for immediate access) ---
#define G_GRAV_CONST 1.0     // Gravitational constant G = 1 in your unit system
// REMOVED: #define M_PI 3.14159265358979323846 // M_PI is already in math.h

// --- General Constants ---
const double AU2CM = 1.495978707e13; // Astronomical Unit in centimeters
const double SUN2GR = 1.989e33; // Solar Mass in grams (Ensure no stray chars here)

// --- Global Variables (declared as extern here, defined and initialized in config.c) ---
// Disk parameters
extern double RMIN;
extern double RMAX;
extern int NGRID;
extern double DD; // DD will be calculated in config.c's init function
extern double SIGMA0;
extern double SIGMAP_EXP;
extern double FLIND;      // Flaring index (gamma from disk_param_be)
extern double HASP;       // Scale height constant (H/R at 1 AU)
extern double STAR;       // Mass of the central star in solar masses
extern double alpha_visc; // alpha viscosity measure
extern double a_mod;      // viscosity reduction measure

// Dust parameters
extern double PDENSITY; // Particle density (Ensure no stray chars here)
extern double PDENSITYDIMLESS; // Dimensionless particle density
extern int PARTICLE_NUMBER;  // Total number of particles

// Dead Zone parameters
extern double r_dze_i; // Inner dead zone radius (Ensure no stray chars here)
extern double r_dze_o; // Outer dead zone radius
extern double Dr_dze_i; // Inner dead zone width
extern double Dr_dze_o; // Outer dead zone width
extern int optdze; // Option for dead zone simulation (Ensure no stray chars here)

// Simulation control options
extern int optev; // Evolution option (Ensure no stray chars here)
extern int optdr; // Drift option
extern int optgr; // Growth option
extern int opttwopop; // Two-population model option
extern int fFrag; // Fragmentation flag
extern int uFrag; // Fragmentation update flag
extern int inputsig; // Input sigma file option

// Time parameters
extern double DT; // Time step (Ensure no stray chars here)
extern double TMAX; // Maximum simulation time
extern double WO; // Output frequency parameter
extern double TCURR; // Current simulation time

// --- Global File Pointers ---
extern FILE *fmo;
extern FILE *fout;
extern FILE *fout2;
extern FILE *fout3;
extern FILE *foutmicr;
extern FILE *massfil;
extern FILE *jelfut;
extern FILE *fin1; // (Ensure no stray chars here)
extern FILE *fin2; // (Ensure no stray chars here)

// --- Global Filenames ---
extern char filenev1[1024];
extern char filenev2[1024];
extern char filenev3[1024];

#endif // CONFIG_H
