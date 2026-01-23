#ifndef SIMULATION_TYPES_H
#define SIMULATION_TYPES_H

#include <stdio.h>

#define MAX_PATH_LEN 8192 // Definiálj egy maximális hosszt fájlnevekhez

// Define fundamental physical constants if they're not in config.h
// #define AU2CM 1.496e13 // Example: Astronomical Unit to Centimeters

// --- Particle Structure ---
// Represents a single dust particle
typedef struct {
    double r;         // Radial distance from star (AU)
    double prad;      // Particle radius (cm)
    double reppmass;  // Representative particle mass (mass of the material this particle represents)
    // Add any other per-particle properties here, e.g.:
    // double velocity_r; // Radial velocity if you track it explicitly
    // int    index_in_grid; // Or calculated on the fly
} particle_t;

// --- Dust Population Structure ---
// Represents a collection of particles of a certain type (e.g., cm, micron, secondary)
typedef struct {
    particle_t *particles; // Pointer to an array of particle_t structures
    int num_particles;     // Number of particles in this specific population

    // If needed, specific properties for the population itself:
    // double total_mass; // Total mass represented by this population
    // double initial_mass; // Initial total mass for comparison
} dust_population_t;

// --- Disk Structure ---
// Encapsulates all global disk parameters and dynamic arrays
typedef struct disk_t {
    // Geometriai és rács paraméterek
    double RMIN;
    double RMAX;
    int NGRID;
    double DD; // RMIN, RMAX, NGRID alapján számolva

    // Gáz korong paraméterek
    double SIGMA0; // Kezdeti referencia gáz sűrűség
    double SIGMAP_EXP; // Sűrűség profil kitevő
    double alpha_visc; // Alfa viszkozitás
    double STAR_MASS; // Központi csillag tömege
    double DISK_MASS;
    double HASP; // H/R arány (diszk magasság)
    double FLIND; // Fáklyázási index

    // Holt zóna paraméterek (dead zone)
    double r_dze_i;
    double r_dze_o;
    double Dr_dze_i;
    double Dr_dze_o;
    double a_mod; // Alfa modifikációs faktor

    // Por paraméterek (ezeket számolja a disk_param_be)
    double PDENSITY; // Por sűrűség (fizikai egységekben, pl. g/cm^3)
    double PDENSITYDIMLESS; // Dimenziómentes por sűrűség

    // Dinamikusan allokált tömbök (ezek a main-ben vannak beállítva, vagy itt lesznek allokálva)
    double *rvec;          // Rádiusz vektor
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

} disk_t;

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


} simulation_options_t;


// --- Output File Pointers Structure ---
typedef struct {
    FILE *por_motion_file;      // pormozgas.dat (fout a régi kódokban)
    FILE *micron_motion_file;   // pormozgasmic.dat (foutmicr a régi kódokban)
    FILE *mass_file;            // mass.dat (massfil a régi kódokban)
    FILE *surface_file;         // surface.dat (a Print_Sigma számára)
    FILE *dust_file;            // dust.dat (a Print_Sigmad fő por kimenetéhez)
    FILE *micron_dust_file;     // dustmic.dat (a Print_Sigmad mikronos por kimenetéhez)
    FILE *size_file; 
    // Add other file pointers here if you need more output files
} output_files_t;

#endif // SIMULATION_TYPES_H