#ifndef SIMULATION_TYPES_H
#define SIMULATION_TYPES_H

#include <stdio.h>

#define MAX_PATH_LEN 8192 // Definiálj egy maximális hosszt fájlnevekhez

// --- Particle Structure ---
// Represents a single dust particle
/*typedef struct {
    double r;         // Radial distance from star (AU)
    double prad;      // Particle radius (cm)
    double reppmass;  // Representative particle mass (mass of the material this particle represents)
} DustParticle;

// --- Dust Population Structure ---
// Represents a collection of particles of a certain type (e.g., cm, micron, secondary)
typedef struct {
    DustParticle *particles; // Pointer to an array of particle_t structures
    int num_particles;     // Number of particles in this specific population
} DustPopulation;
*/
// --- Disk Structure ---
// Encapsulates all global disk parameters and dynamic arrays
typedef struct {
    // Geometriai és rács paraméterek
    double r_min;
    double r_max;
    int grid_number;
    double delta_r; // r_min, r_max, grid_number alapján számolva

    // Gáz korong paraméterek
    double sigma_0; // Kezdeti referencia gáz sűrűség
    double sigma_power_law_index; // Sűrűség profil kitevő
    double alpha_parameter; // Alfa viszkozitás
    double stellar_mass; // Központi csillag tömege
    double disk_mass;
    double h_aspect_ratio; // H/R arány (diszk magasság)
    double flaring_index; // Fáklyázási index

    // Holt zóna paraméterek (dead zone)
    double r_dze_i;
    double r_dze_o;
    double dr_dze_i;
    double dr_dze_o;
    double alpha_parameter_modification; 

    // Por paraméterek (ezeket számolja a disk_param_be)
    double particle_density; // Por sűrűség (fizikai egységekben, pl. g/cm^3)
    double particle_density_dimensionless; // Dimenziómentes por sűrűség

    // Dinamikusan allokált tömbök (ezek a main-ben vannak beállítva, vagy itt lesznek allokálva)
    double *radial_grid;          // Rádiusz vektor
    double *sigmavec;      // Gáz felületi sűrűség
    double *pressvec;      // Gáz nyomás
    double *dpressvec;     // Gáz nyomás gradiens (rádiusz szerinti derivált)
    double *ugvec;         // Gáz belső energia vagy hőmérséklet
    double *sigmadustvec;  // Por felületi sűrűség (Új: sigmad)
    double *rhopvec;       // Por sűrűsége a részecskében (Új: rho_p)

    // Por koagulációs/fragmentációs paraméterek
    double fFrag;          // Fragmentációs hatékonysági faktor
    double uFrag;          // Fragmentációs sebesség küszöb
    double fDrift;         // fd in Birnstiel 2012, set to 0.55

    // További általános paraméterek vagy tömbök
    // Ha a 'dp' és 'rho_p' nem tömbök, hanem skalárok, akkor azok is ide jöhetnek,
    // de a kontextusból valószínűbb, hogy rácspontonként változnak.
    // Mivel az 'dp' már dpressvec néven létezik, csak a 'rho_p'-t kell hozzáadni.
    // Látom, a korábbi hibáknál 'dp' néven hivatkoztál rá, az most 'dpressvec'.
    // A 'rho_p' viszont egy új bejegyzés.

} DiskParameters;

// --- Simulation Options/Control Structure ---
// Groups all boolean/flag-like options
typedef struct {
    double evol;    // Gas evolution (replaces optev)
    double drift;   // Particle drift (replaces optdr)
    double growth;  // Particle growth (replaces optgr)
    double twopop;  // Two-population simulation (replaces opttwopop)
    double dzone;   // Dead zone flag (replaces optdze)
    double DT;      // User-defined fixed time step
    double TMAX;    // Maximum simulation time
    double WO;      // Write-out interval (TMAX/WO)
    double TCURR;

    int num_dust_particles;

    char input_filename[MAX_PATH_LEN];  // Input fájl neve (pl. init_data.dat)
    char output_dir_name[MAX_PATH_LEN]; // Kimeneti könyvtár neve
    char dust_input_filename[MAX_PATH_LEN]; // NEW: Input file for dust particles (e.g., initial_dust_profile.dat for por_be)


} SimulationOptions;


// --- Output File Pointers Structure ---
typedef struct {
    FILE *por_motion_file;      
    FILE *micron_motion_file;   
    FILE *mass_file;            
    FILE *surface_file;         
    FILE *dust_file;            
    FILE *micron_dust_file;     
    FILE *size_file; 
    // Add other file pointers here if you need more output files
} OutputFiles;

#endif // SIMULATION_TYPES_H