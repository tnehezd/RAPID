// src/config.h

#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> // For FILE* and other standard definitions
#include <math.h>  // M_PI definition, if not implicitly defined elsewhere

// --- Fixed Physical Constants (defined directly in header for global access) ---
#define G2 39.47841760435743     // Gravitational constant 4pi^2 AU^3 M_sun^-1 yr^-2
#define M_PI 3.14159265358979323846 // Pi (standard math constant)

// --- General Constants (defined directly in header) ---
const double AU2CM = 1.495978707e13; // Astronomical Unit in centimeters
const double SUN2GR = 1.989e33;     // Solar Mass in grams

// --- Global Variables (declared as extern here, defined and initialized in config.c) ---
// Disk parameters
extern double RMIN;
extern double RMAX;
extern int NGRID;
extern double DD; // Derived from RMIN, RMAX, NGRID, usually calculated in config.c
extern double SIGMA0;
extern double SIGMAP_EXP;
extern double FLIND;      // Flaring index (gamma from disk_param_be)
extern double alpha_visc; // alpha viscosity measure
extern double a_mod;      // viscosity reduction measure
extern double STAR;       // Mass of the central star in solar masses (used in dust_physics)
extern double HASP;       // Scale height constant (used in dust_physics)

// Dust parameters
extern double PDENSITY;      // Particle density (cgs)
extern double PDENSITYDIMLESS; // Dimensionless particle density (calculated in disk_param_be)
extern int PARTICLE_NUMBER;  // Changed to extern, as it might be read from input

// Dead Zone parameters
extern double r_dze_i;      // Inner dead zone radius
extern double r_dze_o;      // Outer dead zone radius
extern double Dr_dze_i;     // Inner dead zone width
extern double Dr_dze_o;     // Outer dead zone width
extern int optdze;          // Option for dead zone simulation

// Simulation control options
extern int optev;           // Evolution option
extern int optdr;           // Drift option
extern int optgr;           // Growth option
extern int opttwopop;       // Two-population model option
extern int fFrag;           // Fragmentation flag
extern int uFrag;           // Fragmentation update flag
extern int inputsig;        // Input sigma file option

// Time parameters
extern double DT;            // Time step
extern double TMAX;          // Maximum simulation time
extern double WO;            // Output frequency parameter
extern double TCURR;         // Current simulation time

// --- Global File Pointers ---
extern FILE *fmo;
extern FILE *fout;
extern FILE *fout2;
extern FILE *fout3;
extern FILE *foutmicr;
extern FILE *massfil;
extern FILE *jelfut;
extern FILE *fin1;
extern FILE *fin2;

// --- Global Filenames ---
extern char filenev1[1024];
extern char filenev2[1024];
extern char filenev3[1024];

#endif // CONFIG_H
