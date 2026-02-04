#ifndef PARSER_H
#define PARSER_H

#include <stdbool.h>  
#include "config.h"   
#include "simulation_types.h" 

#define MAX_OUTPUT_DIR_LEN MAX_PATH_LEN

typedef struct {
    double option_for_dust_drift;
    double option_for_dust_growth;
    double option_for_evolution;
    double option_for_dust_secondary_population;
    double fragmenatation_velocity;
    double fragmenatation_factor;
    int    number_of_grid_points;  
    int    number_of_dust_particles;
    double rmin_val;        
    double rmax_val;        
    double sigma0_val;      
    double sigmap_exp_val;  
    double alpha_visc_val;  
    double star_val;        
    double hasp_val;        
    double flind_val;       
    double r_dze_i_val;     
    double r_dze_o_val;     
    double dr_dze_i_val;    
    double dr_dze_o_val;    
    double a_mod_val;       
    const char *input_file;    
    char output_dir_name[MAX_OUTPUT_DIR_LEN];
    double user_defined_time_step; 
    double maximum_simulation_time; 
    double output_frequency;
    double eps_val;         
    double ratio_val;       
    double mic_val;         
    double onesize_val;     
    double pdensity_val;    
} ParserOptions;

void createDefaultOptions(ParserOptions *opt);
int parseCLIOptions(int argc, const char **argv, ParserOptions *opt);
void printUsageToTerminal();

#endif // PARSER_H