#ifndef INIT_TOOL_MODULE_H
#define INIT_TOOL_MODULE_H

#include <stdbool.h>
#include "simulation_types.h"

typedef struct {
    int     n_grid_points;
    int     n_dust_particles; 
    double  r_inner;
    double  r_outer;
    double  sigma0_gas_au;    
    double  sigma_exponent;   
    double  alpha_viscosity;  
    double  star_mass;        
    double  aspect_ratio;     
    double  flaring_index;    
    double  deadzone_r_inner;
    double  deadzone_r_outer;
    double  deadzone_dr_inner; 
    double  deadzone_dr_outer; 
    double  deadzone_alpha_mod; 
    double  dust_to_gas_ratio;  
    double  disk_mass_dust;     
    double  one_size_particle_cm;  
    double  two_pop_ratio;         
    double  micro_size_cm;         
    double  drift_factor;          
    double  fragmentation_factor;  
    char output_base_path[MAX_PATH_LEN];
    double dust_density_g_cm3; 
} InitializeDefaultOptions;

void initializeDefaultOptions(InitializeDefaultOptions *opt);
int runInitialization(InitializeDefaultOptions *opts, DiskParameters *output_disk_params);

#endif // INIT_TOOL_MODULE_H