#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h>
#include <stdbool.h>
#include "simulation_types.h"
#include "particle_data.h"

int calculateNumbersOfParticles(const char *particle_data_file_name);
void loadDustParticlesFromFile(ParticleData *particle_data, const char *particle_data_file_name);
void loadGasSurfaceDensityFromFile(DiskParameters *disk_params, const char *disk_file_name);
char *createRunDirectory(const char *base_path);
void printCurrentInformationAboutRun(const char *directory_name, const DiskParameters *disk_params);
void printMassGrowthAtDZEFile(double step, 
                double (*dust_particle_mass_array)[5], double (*micron_dust_particle_mass_array)[5], 
                double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, 
                double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, 
                double *tavin, double *tavout, 
                const DiskParameters *disk_params, const SimulationOptions *sim_opts,
                OutputFiles *output_files);
void printGasSurfaceDensityPressurePressureDerivateFile(const DiskParameters *disk_params, OutputFiles *output_files);
void printDustSurfaceDensityPressurePressureDerivateFile(const double *r, const double *rm, const double *dust_surfacedensity, const double *micron_dust_surfacedensity,
                  const DiskParameters *disk_params, const SimulationOptions *sim_opts,
                  OutputFiles *output_files, double step);
void printDustParticleSizeFile(char *size_name, int step, double (*rad)[2], double (*radmicr)[2],
                        const DiskParameters *disk_params, const SimulationOptions *sim_opts,
                        OutputFiles *output_files);

typedef enum {
    FILE_TYPE_MASS_ACCUMULATION,
    FILE_TYPE_GAS_DENSITY,
    FILE_TYPE_DUST_DENSITY,
    FILE_TYPE_DUST_MICRON_DENSITY,
    FILE_TYPE_INTIIAL_DUST_PROFILE,
    FILE_TYPE_PARTICLE_SIZE,
    FILE_TYPE_DISK_PARAM 
} FileType_e;

typedef struct {
    double current_time;  
    int is_initial_data;  
    double R_in;
    double R_out;
    double sigma_exponent; 
    long double sigma0_gas_au;
    double grav_const; 
    double dz_r_inner; 
    double dz_r_outer; 
    double dz_dr_inner_calc;
    double dz_dr_outer_calc;
    double dz_alpha_mod; 
    double dust_density_g_cm3;
    double alpha_viscosity; 
    double star_mass; 
    double flaring_index;
    int n_grid_points; 
} HeaderData;

void printFileHeader(FILE *file, FileType_e file_type, const HeaderData *header_data);
int setupInitialOutputFiles(OutputFiles *output_files, const SimulationOptions *sim_opts,
                               const DiskParameters *disk_params, HeaderData *header_data_for_files);
void cleanupSimulationResources(ParticleData *particle_data, OutputFiles *output_files);
FILE *openSnapshotFile(const char *filename,FileType_e file_type,double current_time_years);
void closeSnapshotFiles(OutputFiles *output_files, const SimulationOptions *sim_opts);
void printFinalSimulationSummary(const char *directory_name, double elapsed_seconds, const SimulationOptions *sim_opts);
void buildSnapshotFilenames(char *dens_name, char *dust_name, char *dust_name2, char *size_name, const SimulationOptions *sim_opts, int snapshot_id);

#endif // IO_UTILS_H

