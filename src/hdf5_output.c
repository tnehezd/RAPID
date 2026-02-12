#include "hdf5_output.h"
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>

#include "particle_data.h"
#include "config.h"

int initHDF5File(const char *filename, OutputFiles *output_files) {
    if (!filename || !output_files) return 1;

    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "ERROR: Could not create HDF5 file %s\n", filename);
        return 1;
    }

    output_files->hdf5_file = (void *)(intptr_t)file_id;
    output_files->snapshot_index = 0;

    return 0;
}

void writeHDF5SnapshotToFile(hid_t file_id, const SimulationOptions *sim_opts,
                             DiskParameters *disk_params, ParticleData *particle_data) {

    (void)sim_opts; // unused

    // -------------------
    // GAS CSOPORT
    // -------------------
    hid_t group_gas = H5Gcreate2(file_id, "/gas", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_gas);

    hsize_t dims_gas[1] = { disk_params->grid_number };
    hid_t space_gas = H5Screate_simple(1, dims_gas, NULL);

    hid_t dset_gas_grid = H5Dcreate2(file_id, "/gas/radial_grid", H5T_NATIVE_DOUBLE, space_gas,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_gas_surface = H5Dcreate2(file_id, "/gas/surface_density", H5T_NATIVE_DOUBLE, space_gas,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_gas_pressure = H5Dcreate2(file_id, "/gas/pressure", H5T_NATIVE_DOUBLE, space_gas,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_gas_grad = H5Dcreate2(file_id, "/gas/pressure_gradient", H5T_NATIVE_DOUBLE, space_gas,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_gas_vel = H5Dcreate2(file_id, "/gas/radial_velocity", H5T_NATIVE_DOUBLE, space_gas,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dset_gas_grid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk_params->radial_grid);
    H5Dwrite(dset_gas_surface, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk_params->gas_surface_density_vector);
    H5Dwrite(dset_gas_pressure, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk_params->gas_pressure_vector);
    H5Dwrite(dset_gas_grad, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk_params->gas_pressure_gradient_vector);
    H5Dwrite(dset_gas_vel, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk_params->gas_velocity_vector);

    H5Dclose(dset_gas_grid);
    H5Dclose(dset_gas_surface);
    H5Dclose(dset_gas_pressure);
    H5Dclose(dset_gas_grad);
    H5Dclose(dset_gas_vel);
    H5Sclose(space_gas);

    // -------------------
    // DUST CSOPORT
    // -------------------
    hid_t group_dust = H5Gcreate2(file_id, "/dust", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_dust);

    hsize_t dims_particles[1] = { particle_data->allocated_particle_number };
    hid_t space_particles = H5Screate_simple(1, dims_particles, NULL);

    // -- positions (AU -> cm)
    hid_t dset_positions = H5Dcreate2(file_id, "/dust/position", H5T_NATIVE_DOUBLE, space_particles,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double *positions = malloc(dims_particles[0] * sizeof(double));
    for (size_t i = 0; i < dims_particles[0]; i++) {
        positions[i] = particle_data->particle_distance_array[i][1] * AU_IN_CM;
    }
    H5Dwrite(dset_positions, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, positions);
    H5Dclose(dset_positions);
    free(positions);

    // -- surface density
    hid_t dset_surface = H5Dcreate2(file_id, "/dust/surface_density", H5T_NATIVE_DOUBLE, space_particles,
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double *surface_densities = malloc(dims_particles[0] * sizeof(double));
    for (size_t i = 0; i < dims_particles[0]; i++) {
        surface_densities[i] = particle_data->dust_surfacedensity[i];
    }
    H5Dwrite(dset_surface, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, surface_densities);
    H5Dclose(dset_surface);
    free(surface_densities);

    // -- size (1. oszlop a dust_particle_mass_array-ból)
    hid_t dset_size = H5Dcreate2(file_id, "/dust/size", H5T_NATIVE_DOUBLE, space_particles,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double *sizes = malloc(dims_particles[0] * sizeof(double));
    for (size_t i = 0; i < dims_particles[0]; i++) {
        sizes[i] = particle_data->dust_particle_mass_array[i][1];
    }
    H5Dwrite(dset_size, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sizes);
    H5Dclose(dset_size);
    free(sizes);

    // -- original index 0..N-1
    hid_t dset_index = H5Dcreate2(file_id, "/dust/index", H5T_NATIVE_INT, space_particles,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int *indices = malloc(dims_particles[0] * sizeof(int));
    for (size_t i = 0; i < dims_particles[0]; i++) indices[i] = (int)i;
    H5Dwrite(dset_index, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
    H5Dclose(dset_index);
    free(indices);

    H5Sclose(space_particles);
}



void closeHDF5File(OutputFiles *output_files) {
    if (output_files->hdf5_file != NULL) {
        H5Fclose((hid_t)output_files->hdf5_file);
        output_files->hdf5_file = NULL;
    }
}

