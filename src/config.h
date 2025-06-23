#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h> // Required for FILE* type
#include <math.h>  // Required for M_PI (for TWOPI macro)

// --- Global Variable Declarations (extern) ---
// Disk parameters



// Dust parameters
extern int PARTICLE_NUMBER;


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
// These constants have been moved here from init_tool.c for global access,
// to prevent duplication or static limitations.

#define SDCONV            1.12521e-7          // Surface density conversion factor
#define ICEFACTOR         3.0                 // Factor for dust density beyond snowline
#define SNOWLINE          2.7                 // Snowline radius in AU
#define G_GRAV_CONST 1.0 // Gravitációs Konstans (dimenziótlan G=1 rendszerben)
// Ha az init_tool_module.c is ezt a G-t akarja használni, akkor itt lehet egy alias
#define G_GRAV_CONST2 (G_GRAV_CONST * G_GRAV_CONST) // Vagy (G*G) ha az aliast használod
#define AUPDAY2CMPSEC     1.7314568e8         // AU/Day to cm/sec conversion
#define CMPSECTOAUPYRP2PI 3.35725e-07         // cm/sec to AU/(yr/2pi) conversion
#define GRPCM32MSUNAU3    1.68329e6           // gr/cm^3 to M_sun/AU^3 conversion

#define SUN2GR            1.989e33            // Solar Mass in grams (M_solar -> g) - PLEASE VERIFY THIS VALUE!
#define AU2CM             1.496e13            // Astronomical Unit in centimeters (AU -> cm) - PLEASE VERIFY THIS VALUE!
#define TWOPI             (2.0 * M_PI)        // Added parentheses for safety with expressions
#define KEREK             1.0                 // Make sure this value is correct for your physics model!


// --- Global Filename Declarations (extern) ---
extern const char * const FILENAME_INIT_PROFILE;
extern const char * const FILENAME_DISK_PARAM;


// Parameter controls
extern const char *inputsig;


// Function Declaration (if this function is implemented in a .c file and used elsewhere)
void initialize_derived_config_variables();

#endif // CONFIG_H