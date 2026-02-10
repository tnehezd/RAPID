#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <stdio.h>
#include <stdbool.h>
#include "simulation_types.h"
#include "particle_data.h"

/**
 * @brief Calculates the number of particles stored in a particle data file.
 *
 * @param particle_data_file_name Path to the particle data file.
 * @return The number of particles found in the file.
 */
int calculateNumbersOfParticles(const char *particle_data_file_name);

/**
 * @brief Loads dust particle information from a file into a ParticleData structure.
 *
 * @param particle_data Pointer to the ParticleData structure to fill.
 * @param particle_data_file_name Path to the particle data file.
 */
void loadDustParticlesFromFile(ParticleData *particle_data, const char *particle_data_file_name);

/**
 * @brief Loads gas surface density values from a file into DiskParameters.
 *
 * @param disk_params Pointer to the DiskParameters structure to fill.
 * @param disk_file_name Path to the disk parameter file.
 */
void loadGasSurfaceDensityFromFile(DiskParameters *disk_params, const char *disk_file_name);

/**
 * @brief Creates a new run directory based on a base path.
 *
 * @param base_path Base directory path.
 * @return Newly allocated string containing the created directory path.
 */
char *createRunDirectory(const char *base_path);

/**
 * @brief Prints current simulation information to the console or log.
 *
 * @param directory_name Name of the run directory.
 * @param disk_params Pointer to the DiskParameters structure.
 */
void printCurrentInformationAboutRun(const char *directory_name, const DiskParameters *disk_params);

/**
 * @brief Writes gas surface density, pressure, and pressure derivative data to file.
 *
 * @param disk_params Pointer to the DiskParameters structure.
 * @param output_files Pointer to the OutputFiles structure containing file handles.
 */
void printGasSurfaceDensityPressurePressureDerivateFile(const DiskParameters *disk_params, OutputFiles *output_files);

/**
 * @brief Writes dust surface density and pressure derivative data to file.
 *
 * @param r Radial grid array.
 * @param rm Midpoint radial grid array.
 * @param dust_surfacedensity Dust surface density array.
 * @param micron_dust_surfacedensity Micron-sized dust surface density array.
 * @param disk_params Pointer to DiskParameters.
 * @param sim_opts Pointer to SimulationOptions.
 * @param output_files Pointer to OutputFiles.
 * @param step Current simulation step.
 */
void printDustSurfaceDensityPressurePressureDerivateFile(const double *r, const double *rm, const double *dust_surfacedensity, const double *micron_dust_surfacedensity,
                  const DiskParameters *disk_params, const SimulationOptions *sim_opts,
                  OutputFiles *output_files, double step);

/**
 * @brief Writes dust particle size distribution to file.
 *
 * @param size_name Output filename.
 * @param step Current simulation step.
 * @param rad 2D array of particle radii.
 * @param radmicr 2D array of micron-sized particle radii.
 * @param disk_params Pointer to DiskParameters.
 * @param sim_opts Pointer to SimulationOptions.
 * @param output_files Pointer to OutputFiles.
 */
void printDustParticleSizeFile(char *size_name, int step, double (*rad)[2], double (*radmicr)[2],
                        const DiskParameters *disk_params, const SimulationOptions *sim_opts,
                        OutputFiles *output_files);

/**
 * @brief Writes the evolution of trap mass over time to file.
 *
 * @param current_time_years Current simulation time in years.
 * @param num_found Number of traps found.
 * @param traps Pointer to array of PressureTrap structures.
 * @param output_files Pointer to OutputFiles.
 */
void printTrapMassEvolution(double current_time_years, int num_found, const PressureTrap *traps, OutputFiles *output_files);

/**
 * @brief Enumeration of supported file header types.
 */
typedef enum {
    FILE_TYPE_MASS_ACCUMULATION,
    FILE_TYPE_GAS_DENSITY,
    FILE_TYPE_DUST_DENSITY,
    FILE_TYPE_DUST_MICRON_DENSITY,
    FILE_TYPE_INTIIAL_DUST_PROFILE,
    FILE_TYPE_PARTICLE_SIZE,
    FILE_TYPE_DISK_PARAM 
} FileType_e;

/**
 * @brief Structure containing metadata for file headers.
 */
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

/**
 * @brief Prints a formatted header to a file based on file type and header data.
 *
 * @param file File pointer to write to.
 * @param file_type Type of file header to print.
 * @param header_data Pointer to HeaderData structure.
 */
void printFileHeader(FILE *file, FileType_e file_type, const HeaderData *header_data);

/**
 * @brief Initializes output files and writes initial headers.
 *
 * @param output_files Pointer to OutputFiles structure to initialize.
 * @param sim_opts Pointer to SimulationOptions.
 * @param disk_params Pointer to DiskParameters.
 * @param header_data_for_files Pointer to HeaderData structure to populate.
 * @return 0 on success, non-zero on failure.
 */
int setupInitialOutputFiles(OutputFiles *output_files, const SimulationOptions *sim_opts,
                               const DiskParameters *disk_params, HeaderData *header_data_for_files);

/**
 * @brief Frees memory and closes files associated with the simulation.
 *
 * @param particle_data Pointer to ParticleData to free.
 * @param output_files Pointer to OutputFiles to close.
 */
void cleanupSimulationResources(ParticleData *particle_data, OutputFiles *output_files);

/**
 * @brief Opens a snapshot file for writing simulation data.
 *
 * @param filename Name of the snapshot file.
 * @param file_type Type of snapshot file.
 * @param current_time_years Current simulation time in years.
 * @return File pointer to the opened file.
 */
FILE *openSnapshotFile(const char *filename,FileType_e file_type,double current_time_years);

/**
 * @brief Closes all snapshot files associated with the simulation.
 *
 * @param output_files Pointer to OutputFiles.
 * @param sim_opts Pointer to SimulationOptions.
 */
void closeSnapshotFiles(OutputFiles *output_files, const SimulationOptions *sim_opts);

/**
 * @brief Prints a summary of the simulation after completion.
 *
 * @param directory_name Name of the simulation directory.
 * @param elapsed_seconds Total runtime in seconds.
 * @param sim_opts Pointer to SimulationOptions.
 */
void printFinalSimulationSummary(const char *directory_name, double elapsed_seconds, const SimulationOptions *sim_opts);

/**
 * @brief Builds filenames for snapshot output files.
 *
 * @param dens_name Output buffer for gas density filename.
 * @param dust_name Output buffer for dust density filename.
 * @param dust_name2 Output buffer for micron dust density filename.
 * @param size_name Output buffer for particle size filename.
 * @param sim_opts Pointer to SimulationOptions.
 * @param snapshot_id Snapshot index.
 */
void buildSnapshotFilenames(char *dens_name, char *dust_name, char *dust_name2, char *size_name, const SimulationOptions *sim_opts, int snapshot_id);

#endif // IO_UTILS_H

