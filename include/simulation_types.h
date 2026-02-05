#ifndef SIMULATION_TYPES_H
#define SIMULATION_TYPES_H

#include <stdio.h>

#define MAX_PATH_LEN 16384 

typedef struct {
    double r_min;
    double r_max;
    int grid_number;
    double delta_r;
    double sigma_0;
    double sigma_power_law_index;
    double alpha_parameter;
    double stellar_mass;
    double disk_mass;
    double h_aspect_ratio; 
    double flaring_index; 
    double r_dze_i;
    double r_dze_o;
    double dr_dze_i;
    double dr_dze_o;
    double alpha_parameter_modification; 
    double particle_density; 
    double particle_density_dimensionless;
    double *radial_grid;          
    double *gas_surface_density_vector;
    double *gas_pressure_vector;     
    double *gas_pressure_gradient_vector; 
    double *gas_velocity_vector;         
    double fragmentation_factor;   
    double fragmentation_velocity; 
    double drift_factor;         
} DiskParameters;

typedef struct {
    double option_for_evolution;  
    double option_for_dust_drift;  
    double option_for_dust_growth; 
    double option_for_dust_secondary_population; 
    double flag_for_deadzone;   
    double user_defined_time_step;
    double maximum_simulation_time; 
    double output_frequency;      
    int number_of_dust_particles;
    char input_filename[MAX_PATH_LEN]; 
    char output_dir_name[MAX_PATH_LEN];
    char dust_input_filename[MAX_PATH_LEN]; 
} SimulationOptions;


typedef struct {
    FILE *dust_motion_file;      
    FILE *micron_motion_file;   
    FILE *mass_file;            
    FILE *surface_file;         
    FILE *dust_file;            
    FILE *micron_dust_file;     
    FILE *size_file; 
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


/**
 * @struct PressureTrap
 * @brief Data structure characterizing gas pressure maxima and the localized dust accumulation.
 *
 * In protoplanetary disks, pressure maxima occur where the radial gas pressure gradient 
 * vanishes (transitioning from positive to negative values). These regions act as 
 * "dust traps" by halting the inward radial drift of solid particles, making them 
 * prime locations for planetesimal formation.
 * * This structure stores the precise location of the trap, its defined measurement 
 * boundaries, and the mass of various dust populations collected within.
 */
typedef struct {
    /** @brief Radial position of the pressure maximum [AU].
     * Determined by finding the root (zero-crossing) of the gas pressure gradient. */
    double radial_position;             

    /** @brief Inner boundary of the measurement annulus [AU].
     * Defines the starting radius for integrating the dust mass around the trap. */
    double inner_boundary;              

    /** @brief Outer boundary of the measurement annulus [AU].
     * Defines the ending radius for integrating the dust mass around the trap. */
    double outer_boundary;              
    
    /** @brief Total mass of the primary (cm-sized) dust population [M_sun]. */
    double primary_dust_mass;           
    
    /** @brief Total mass of the secondary (micron-sized) dust population [M_sun]. */
    double secondary_dust_mass;         
    
    /** @brief Combined mass of all dust populations within the trap [M_sun].
     * Calculated as primary_dust_mass + secondary_dust_mass. */
    double total_dust_mass;             
    
    /** @brief Unique identifier for the trap within the current timestep. 
     * Essential for tracking multiple features (e.g., DZE edges vs. planet-carved gaps). */
    int trap_id;                        
} PressureTrap;

#endif // SIMULATION_TYPES_H