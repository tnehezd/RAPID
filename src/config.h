#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> // For FILE* and other standard definitions

// --- Global Variables (declared as extern here, defined and initialized in config.c) ---
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
extern double STAR; // Mass of the central star in solar masses
extern const double AU2CM;         // Declare it as a const double


// Dust parameters
extern double PDENSITY;      // Particle density
extern double PDENSITYDIMLESS; // Dimensionless particle density
extern const int PARTICLE_NUMBER; // Declare it as a const int

// Dead Zone parameters
extern double r_dze_i;      // Inner dead zone radius
extern double r_dze_o;      // Outer dead zone radius
extern double Dr_dze_i;     // Inner dead zone width
extern double Dr_dze_o;     // Outer dead zone width
extern int optdze;          // Option for dead zone simulation (1 if inner DZE exists, 0 otherwise)

// Simulation control options
extern int optev;           // Evolution option (1 for gas evolution, 0 otherwise)
extern int optdr;           // Drift option (1 for particle drift, 0 otherwise)
extern int optgr;           // Growth option (1 for particle growth, 0 otherwise)
extern int opttwopop;       // Two-population model option (1 for two populations, 0 for one)
extern int fFrag;           // Fragmentation flag
extern int uFrag;           // Fragmentation update flag
extern int inputsig;        // Input sigma file option (0 if input file is used, 1 otherwise)

// Time parameters
extern double DT;            // Time step (from input)
extern double TMAX;          // Maximum simulation time
extern double WO;            // Output frequency parameter (period related)
extern double TCURR;         // Current simulation time (start time)

// --- Global File Pointers (if they are truly global and needed across multiple files) ---
// Note: It's generally better to pass file pointers as arguments to functions,
// but if they are few and used widely, global declaration might be acceptable.
extern FILE *fmo;
extern FILE *fout;
extern FILE *fout2;
extern FILE *fout3;
extern FILE *foutmicr;
extern FILE *massfil;
extern FILE *jelfut;

// --- Global Filenames (if they are used across multiple files) ---
extern char filenev1[1024]; // General purpose filename (e.g., for disk_param_be input)
extern char filenev2[1024]; // General purpose filename (e.g., for initial disk profile)
extern char filenev3[1024]; // General purpose filename (e.g., for time parameters)

#endif // CONFIG_H
