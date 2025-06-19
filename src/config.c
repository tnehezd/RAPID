// src/config.c
#include "config.h"

// --- Global Variable Definitions and Initializations ---
// Disk parameters
double RMIN = 0.1;
double RMAX = 100.0;
int NGRID = 100;
double DD; // Declare DD, but DO NOT initialize here

double SIGMA0 = 1.0;
double SIGMAP_EXP = -1.5;
double FLIND = -0.5;
double alpha_visc = 1e-4;
double a_mod = 0.1;
double STAR = 1.0;
double HASP = 0.05;

// Dust parameters
double PDENSITY = 0.0;
double PDENSITYDIMLESS = 0.0;
int PARTICLE_NUMBER = 1000;

// Dead Zone parameters
double r_dze_i = 0.0;
double r_dze_o = 0.0;
double Dr_dze_i = 0.0;
double Dr_dze_o = 0.0;
int optdze = 0;

// Simulation control options
int optev = 0;
int optdr = 0;
int optgr = 0;
int opttwopop = 0;
int fFrag = 0;
int uFrag = 0;
int inputsig = 0;

// Time parameters
double DT = 0.0;
double TMAX = 0.0;
double WO = 0.0;
double TCURR = 0.0;

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

// --- Global Filename Definitions ---
char filenev1[1024] = "param.dat";
char filenev2[1024] = "sigma.dat";
char filenev3[1024] = "time.dat";

// Function to initialize derived global variables
void initialize_derived_config_variables() {
    DD = (RMAX - RMIN) / NGRID;
}


// Function to initialize derived global variables
void initialize_derived_config_variables() {
    DD = (RMAX - RMIN) / NGRID;
}
