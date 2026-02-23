#ifndef HDF5_OUTPUT_H
#define HDF5_OUTPUT_H

#include "simulation_types.h"
#include "particle_data.h"
#include <hdf5.h>


/**
 * @file hdf5_output.h
 * @brief HDF5 output module for simulation snapshots and time-series data.
 *
 * This module provides:
 *  - Creation and writing of per-snapshot HDF5 files
 *  - Creation and incremental writing of a time-series HDF5 file
 *  - Safe cleanup of all HDF5 resources
 *
 * The snapshot writer stores full gas/dust grid data and particle states.
 * The time-series writer stores the evolution of pressure traps and their
 * associated dust masses.
 */

#define MAX_TRAPS 5   /**< Maximum number of pressure traps tracked in time-series output */


/* ============================================================================
 *  SNAPSHOT HDF5 OUTPUT
 * ==========================================================================*/

/**
 * @brief Create a new HDF5 file for writing a simulation snapshot.
 *
 * This function creates a fresh HDF5 file and stores the file handle inside
 * the OutputFiles structure. The file is overwritten if it already exists.
 *
 * @param filename Path to the HDF5 file to create.
 * @param output_files Pointer to the OutputFiles structure that will store the file handle.
 * @return 0 on success, non-zero on failure.
 */
int initHDF5File(const char *filename, OutputFiles *output_files);

/**
 * @brief Write a complete simulation snapshot into an already opened HDF5 file.
 *
 * This function writes:
 *  - Gas grid (radius, surface density, pressure, gradients, velocity)
 *  - Dust Eulerian grid (surface density)
 *  - Lagrangian particle positions, sizes, and indices
 *  - Frame metadata (time, code version, compile date/time)
 *
 * @param time Simulation time (in internal units).
 * @param file_id HDF5 file identifier returned by initHDF5File().
 * @param sim_opts Simulation options.
 * @param disk_params Disk model parameters.
 * @param particle_data Particle data arrays.
 */
void writeHDF5SnapshotToFile(double time, hid_t file_id, const SimulationOptions *sim_opts,
                             DiskParameters *disk_params, StructuredParticleData *data);
/**
 * @brief Close the currently open snapshot HDF5 file.
 *
 * This function closes the file handle stored in OutputFiles and resets it.
 *
 * @param output_files Pointer to the OutputFiles structure.
 */
void closeHDF5File(OutputFiles *output_files);


/* ============================================================================
 *  TIME-SERIES (TRAP EVOLUTION) HDF5 OUTPUT
 * ==========================================================================*/

/**
 * @brief Initialize the time-series HDF5 file for trap evolution tracking.
 *
 * Creates the following datasets, each with unlimited time dimension:
 *  - time[ts_row][MAX_TRAPS]
 *  - primary_mass[ts_row][MAX_TRAPS]
 *  - secondary_mass[ts_row][MAX_TRAPS]
 *  - total_mass[ts_row][MAX_TRAPS]
 *  - trap_position[ts_row][MAX_TRAPS]
 *
 * @param filename Path to the time-series HDF5 file to create.
 */
void initializeMassTimeSeries(const char *filename);

/**
 * @brief Append one row of trap evolution data to the time-series HDF5 file.
 *
 * Each call appends a new time index (ts_row) containing:
 *  - snapshot time
 *  - per-trap primary dust mass
 *  - per-trap secondary dust mass
 *  - per-trap total dust mass
 *  - per-trap radial positions
 *
 * @param snapshot Snapshot time (in years).
 * @param m1 Array of primary dust masses for each trap.
 * @param m2 Array of secondary dust masses for each trap.
 * @param mtot Array of total dust masses for each trap.
 * @param trap_pos Array of trap radial positions.
 */
void appendMassTimeSeries(double snapshot,
                          double m1[MAX_TRAPS],
                          double m2[MAX_TRAPS],
                          double mtot[MAX_TRAPS],
                          double trap_pos[MAX_TRAPS]);

/**
 * @brief Close the time-series HDF5 file and all associated datasets.
 */
void closeMassTimeSeries(void);




#endif // HDF5_OUTPUT_H
