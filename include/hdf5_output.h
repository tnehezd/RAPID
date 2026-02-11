#ifndef HDF5_OUTPUT_H
#define HDF5_OUTPUT_H

#include "simulation_types.h"
#include "particle_data.h"
#include <hdf5.h>


/**
 * @brief Initialize HDF5 file for simulation output.
 *
 * @param filename Path to the HDF5 file.
 * @param output_files Pointer to OutputFiles struct.
 * @return 0 on success, non-zero on failure.
 */
int initHDF5File(const char *filename, OutputFiles *output_files);

/**
 * @brief Write a snapshot of gas/dust data to HDF5.
 *
 * @param snapshot_index Index of the snapshot.
 * @param output_files Pointer to OutputFiles struct.
 * @return 0 on success, non-zero on failure.
 */
void writeHDF5Snapshot(const SimulationOptions *sim_opts, OutputFiles *output_files,
                       DiskParameters *disk_params, ParticleData *particle_data, int snapshot_index);



void writeHDF5SnapshotToFile(hid_t file_id, const SimulationOptions *sim_opts,
                             DiskParameters *disk_params, ParticleData *particle_data);

/**
 * @brief Close HDF5 file and cleanup resources.
 *
 * @param output_files Pointer to OutputFiles struct.
 */
void closeHDF5File(OutputFiles *output_files);

#endif // HDF5_OUTPUT_H
