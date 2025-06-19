// src/config.h
#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> // Required for FILE* type

// --- Global Variable Declarations (extern) ---
// Disk parameters
extern double RMIN;
extern double RMAX;
extern int NGRID;
extern double DD;

extern double SIGMA0;
extern double SIGMAP_EXP;
extern double FLIND;
extern double alpha_visc;
extern double a_mod;
extern double STAR;
extern double HASP;

// Dust parameters
extern double PDENSITY;
extern double PDENSITYDIMLESS;
extern int PARTICLE_NUMBER;

// Dead Zone parameters
extern double r_dze_i;
extern double r_dze_o;
extern double Dr_dze_i;
extern double Dr_dze_o;
extern int optdze;

// Simulation control options
extern int optev;
extern int optdr;
extern int optgr;
extern int opttwopop;
extern int fFrag;
extern int uFrag;
extern int inputsig;

// Time parameters
extern double DT;
extern double TMAX;
extern double WO;
extern double TCURR;

// --- Global File Pointer Declarations (extern) ---
extern FILE *fmo;
extern FILE *fout;
extern FILE *fout2;
extern FILE *fout3;
extern FILE *foutmicr;
extern FILE *massfil;
extern FILE *jelfut;
extern FILE *fin1;
extern FILE *fin2;



// --- Physical Constants (Add these if they are not already there) ---
#define SUN2GR 1.989e33  // Solar Mass in grams (M_solar -> g) - ELLENŐRIZD AZ ÉRTÉKET!
#define AU2CM  1.496e13  // Astronomical Unit in centimeters (AU -> cm) - ELLENŐRIZD AZ ÉRTÉKET!


// --- Global Filename Declarations (extern) ---
extern char filenev1[1024];
extern char filenev2[1024];
extern char filenev3[1024];

// Function Declaration
void initialize_derived_config_variables();

#endif // CONFIG_H
