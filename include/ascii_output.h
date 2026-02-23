/**
 * @file ascii_output.h
 * @brief Utility functions for reading and writing simulation data.
 *
 * Provides routines for loading input profiles, writing snapshot files,
 * generating output filenames, printing standardized headers, and managing
 * simulation output directories and file resources.
 */

#ifndef ASCII_OUTPUT_H
#define ASCII_OUTPUT_H

#include <stdio.h>
#include <stdbool.h>
#include "simulation_types.h"
#include "particle_data.h"

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
 *
 * This enum identifies the different kinds of output files produced
 * during the simulation. Each value corresponds to a specific file
 * format or data category, allowing the I/O system to write the
 * appropriate header information.
 */
typedef enum {
    FILE_TYPE_MASS_ACCUMULATION,      /**< Mass accumulation output file. */
    FILE_TYPE_GAS_DENSITY,            /**< Gas surface density file. */
    FILE_TYPE_DUST_DENSITY,           /**< Dust surface density file. */
    FILE_TYPE_DUST_MICRON_DENSITY,    /**< Micron-sized dust density file. */
    FILE_TYPE_INTIIAL_DUST_PROFILE,   /**< Initial dust distribution profile. */
    FILE_TYPE_PARTICLE_SIZE,          /**< Dust particle size distribution file. */
    FILE_TYPE_DISK_PARAM              /**< Disk parameter summary file. */
} FileType_e;

/**
 * @brief Structure containing metadata for file headers.
 *
 * This structure stores all metadata written into the header section
 * of simulation output files. These values describe the physical
 * parameters of the disk model, numerical setup, and the simulation
 * state at the time of writing.
 */
typedef struct {
    double current_time;          /**< Current simulation time (in years or code units). */
    int    is_initial_data;       /**< Flag indicating whether this is initial-condition output. */
    double R_in;                  /**< Inner radius of the simulation domain. */
    double R_out;                 /**< Outer radius of the simulation domain. */
    double sigma_exponent;        /**< Power-law exponent of the gas surface density profile. */
    long double sigma0_gas_au;    /**< Gas surface density normalization at 1 AU. */
    double grav_const;            /**< Gravitational constant used in the simulation. */
    double dz_r_inner;            /**< Inner dead-zone edge radius. */
    double dz_r_outer;            /**< Outer dead-zone edge radius. */
    double dz_dr_inner_calc;      /**< Calculated width of the inner dead-zone transition. */
    double dz_dr_outer_calc;      /**< Calculated width of the outer dead-zone transition. */
    double dz_alpha_mod;          /**< Alpha-viscosity modification factor inside the dead zone. */
    double dust_density_g_cm3;    /**< Dust material density (g/cm³). */
    double alpha_viscosity;       /**< Disk viscosity parameter α. */
    double star_mass;             /**< Stellar mass (in solar masses). */
    double flaring_index;         /**< Disk flaring index. */
    int    n_grid_points;         /**< Number of radial grid points in the simulation. */
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
int setupInitialOutputFiles(OutputFiles *output_files, const SimulationOptions *sim_opts, const DiskParameters *disk_params, HeaderData *header_data_for_files);

/**
 * @brief Frees memory and closes files associated with the simulation.
 *
 * @param particle_data Pointer to ParticleData to free.
 * @param output_files Pointer to OutputFiles to close.
 */
void cleanupSimulationResources(StructuredParticleData *data, OutputFiles *output_files);

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


void writeGasField2D(const DiskParameters *disk_params, const char *directory, int snapshot_index, const char *label);

void writeDustField2D(const StructuredParticleData *sdata, const char *directory, int snapshot_index, const char *label);


void printDustParticleSizeFileStructured(char *size_name, int step, const StructuredParticleData *spd, const DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files);

#endif // ASCII_OUTPUT_H

