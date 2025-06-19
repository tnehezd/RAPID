#include "config.h"

// --- Global Variable Definitions and Initializations ---
// Disk parameters
double RMIN = 0.0;
double RMAX = 0.0;
int NGRID = 0;
double DD = 0.0;
double SIGMA0 = 0.0;
double SIGMAP_EXP = 0.0;
double FLIND = 0.0;
double alpha_visc = 0.0;
double a_mod = 0.0;
double STAR = 0.0;


// Dust parameters
double PDENSITY = 0.0;
double PDENSITYDIMLESS = 0.0;


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

// --- Global Filename Definitions ---
char filenev1[1024] = "param.dat"; // Example default, adjust as needed
char filenev2[1024] = "sigma.dat"; // Example default, adjust as needed
char filenev3[1024] = "time.dat";  // Example default, adjust as needed
