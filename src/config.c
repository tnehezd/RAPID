// src/config.c
#include "config.h"

// --- Global Variable Definitions and Initializations ---
// Disk parameters
double RMIN;
double RMAX;
int NGRID;
double DD; // Declare DD, but DO NOT initialize here

double SIGMA0;
double SIGMAP_EXP;
double FLIND;
double alpha_visc;
double a_mod;
double STAR;
double HASP;

// Dust parameters
double PDENSITY;
double PDENSITYDIMLESS;
int PARTICLE_NUMBER;

// Dead Zone parameters
double r_dze_i;
double r_dze_o;
double Dr_dze_i;
double Dr_dze_o;


// Time parameters
double DT;
double TMAX;
double WO;
double TCURR;

// --- Global File Pointer Definitions ---
FILE *fmo = NULL;
FILE *fout = NULL;
FILE *fout2 = NULL;
FILE *fout3 = NULL;
FILE *foutmicr = NULL;
FILE *massfil = NULL;
FILE *jelfut = NULL;
FILE *fin1 = NULL;
FILE *fin2 = NULL;
FILE *fil = NULL;
const char *inputsig = NULL;

// --- Global Filename Definitions ---
char filenev1[1024] = "init_data.dat";
char filenev2[1024] = "disk_param.dat";
char filenev3[1024] = "time.dat";


// opt parsing things
double optev, optdr, optgr, opttwopop, optdze, optinp; // Inicializáld őket 0.0-ra vagy 1.0-ra ha szükséges
double fFrag, uFrag;


// Function to initialize derived global variables
void initialize_derived_config_variables() {
    DD = (RMAX - RMIN) / NGRID;
}


