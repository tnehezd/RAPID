#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>  // Required for FILE* type
#include <math.h>   // Required for M_PI (for TWOPI macro)

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
extern FILE *fil;


// --- Physical Constants (Macros) ---
// Renamed G_GRAV_CONST to G for brevity and consistency with G2 (G*G)
#define G           1.0
#define G2          (G * G) // G*G as requested

// Removed duplicate SUN2GR. Keep one, and ensure its value is correct.
#define SUN2GR      1.989e33  // Solar Mass in grams (M_solar -> g) - PLEASE VERIFY THIS VALUE!
#define AU2CM       1.496e13  // Astronomical Unit in centimeters (AU -> cm) - PLEASE VERIFY THIS VALUE!
#define TWOPI       (2.0 * M_PI) // Added parentheses for safety with expressions
#define SDCONV      1.12521e-7
#define CMPSECTOAUPYRP2PI 3.35725e-07


// --- Global Filename Declarations (extern) ---
extern char filenev1[1024];
extern char filenev2[1024];
extern char filenev3[1024];

// Parameter controls
extern double optev, optdr, optgr, opttwopop, optdze, optinp; // variables to control different parts of the code
extern char *inputsig; // Or const char *inputsig; - typically const char* for string literals

// Part. physics
extern double fFrag, uFrag;


// Function Declaration (if this function is implemented in a .c file and used elsewhere)
void initialize_derived_config_variables();

#endif // CONFIG_H
