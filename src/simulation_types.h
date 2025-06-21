#ifndef SIMULATION_TYPES_H
#define SIMULATION_TYPES_H

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

    // Dinamikusan allokált tömbök (ezek a main-ben vannak beállítva)
    double *rvec;
    double *sigmavec;
    double *pressvec;
    double *dpressvec;
    double *ugvec;

    // Ha vannak további konstansok vagy segédértékek, amiket globálisan használnál,
    // azok is ide kerülhetnek, ha a diskkel kapcsolatosak.

} disk_t;

// --- Simulation Options/Control Structure ---
// Groups all boolean/flag-like options
typedef struct {
    double evol;      // optev in your code (gas evolution)
    double drift;     // optdr (particle drift)
    double growth;    // optgr (particle growth)
    double twopop;    // opttwopop (two-population simulation)
    double fFrag;
    double uFrag;
    double DT;        // User-defined fixed time step
    double TMAX;      // Maximum simulation time
    double WO;        // Write-out interval (TMAX/WO)
    double TCURR; // <-- Add this line

} simulation_options_t;

// --- Output File Pointers Structure (optional but good for organization) ---
typedef struct {
    FILE *por_motion_file;      // fout
    FILE *micron_motion_file;   // foutmicr
    FILE *mass_file;            // massfil
    // Add other file pointers here (e.g., surface density, dust profiles)
} output_files_t;


#endif // SIMULATION_TYPES_H