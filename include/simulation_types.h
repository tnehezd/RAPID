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
    double *gas_surface_density_vector;      // Gáz felületi sűrűség
    double *gas_pressure_vector;      // Gáz nyomás
    double *gas_pressure_gradient_vector;     // Gáz nyomás gradiens (rádiusz szerinti derivált)
    double *gas_velocity_vector;         // Gáz belső energia vagy hőmérséklet
//    double *dust_surface_density_vector;  // Por felületi sűrűség (Új: sigmad)

    // Por koagulációs/fragmentációs paraméterek
    double fragmentation_factor;          // Fragmentációs hatékonysági faktor
    double fragmentation_velocity;          // Fragmentációs sebesség küszöb
    double drift_factor;         // fd in Birnstiel 2012, set to 0.55

    // További általános paraméterek vagy tömbök
    // Ha a 'dp' és 'rho_p' nem tömbök, hanem skalárok, akkor azok is ide jöhetnek,
    // de a kontextusból valószínűbb, hogy rácspontonként változnak.
    // Mivel az 'dp' már gas_pressure_gradient_vector néven létezik, csak a 'rho_p'-t kell hozzáadni.
    // Látom, a korábbi hibáknál 'dp' néven hivatkoztál rá, az most 'gas_pressure_gradient_vector'.
    // A 'rho_p' viszont egy új bejegyzés.

} DiskParameters;

// --- Simulation Options/Control Structure ---
// Groups all boolean/flag-like options
typedef struct {
    double option_for_evolution;    // Gas evolution (replaces optev)
    double option_for_dust_drift;   // Particle drift (replaces optdr)
    double option_for_dust_growth;  // Particle growth (replaces optgr)
    double option_for_dust_secondary_population;  // Two-population simulation (replaces optoption_for_dust_secondary_population)
    double flag_for_deadzone;   // Dead zone flag (replaces optdze)
    double user_defined_time_step;      // User-defined fixed time step
    double maximum_simulation_time;    // Maximum simulation time
    double output_frequency;      // Write-out interval (maximum_simulation_time/output_frequency)

    int number_of_dust_particles;

    char input_filename[MAX_PATH_LEN];  // Input fájl neve (pl. init_data.dat)
    char output_dir_name[MAX_PATH_LEN]; // Kimeneti könyvtár neve
    char dust_input_filename[MAX_PATH_LEN]; // NEW: Input file for dust particles (e.g., initial_dust_profile.dat for por_be)


} SimulationOptions;


// --- Output File Pointers Structure ---
typedef struct {
    FILE *dust_motion_file;      
    FILE *micron_motion_file;   
    FILE *mass_file;            
    FILE *surface_file;         
    FILE *dust_file;            
    FILE *micron_dust_file;     
    FILE *size_file; 
    // Add other file pointers here if you need more output files
} OutputFiles;


typedef enum {
    SnapshotNonevolving = 0,
    SnapshotGas = 1,
    SnapshotDrift = 2,
    SnapshotGrowth = 3,
    SnapshotDriftTwoPop = 4,
    SnapshotGrowthTwoPop = 5
} SnapshotMode;


SnapshotMode determineSnapshotMode(const SimulationOptions *sim_opts);
const char* snapshotModeToString(SnapshotMode mode);



#endif // SIMULATION_TYPES_H