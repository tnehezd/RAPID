#include "hdf5_output.h"
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>

#include "particle_data.h"

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


    hid_t group_id;
    hid_t dataspace_gas_id, dataset_gas_surfacedensity, dataset_gas_pressure, dataset_gas_pressure_gradient, dataset_gas_radial_velocity, dataset_gas_grid;
    hid_t dataspace_dust_id, dataset_dust_surfacedensity;

    // gas csoport létrehozása
    group_id = H5Gcreate2(file_id, "/gas", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    // dust csoport létrehozása
    group_id = H5Gcreate2(file_id, "/dust", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    // Gas surface density
    hsize_t dims_gas[1] = {disk_params->grid_number};
    dataspace_gas_id = H5Screate_simple(1, dims_gas, NULL);
    dataset_gas_radial_velocity = H5Dcreate2(file_id, "/gas/radial_velocity", H5T_NATIVE_DOUBLE,
                                  dataspace_gas_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dataset_gas_grid = H5Dcreate2(file_id, "/gas/radial_grid", H5T_NATIVE_DOUBLE,
                                  dataspace_gas_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dataset_gas_surfacedensity = H5Dcreate2(file_id, "/gas/surface_density", H5T_NATIVE_DOUBLE,
                                  dataspace_gas_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dataset_gas_pressure = H5Dcreate2(file_id, "/gas/pressure", H5T_NATIVE_DOUBLE,
                                      dataspace_gas_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dataset_gas_pressure_gradient = H5Dcreate2(file_id, "/gas/pressure_gradient", H5T_NATIVE_DOUBLE,
                                      dataspace_gas_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset_gas_grid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk_params->radial_grid);
    H5Dclose(dataset_gas_grid);
    H5Dwrite(dataset_gas_radial_velocity, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk_params->gas_velocity_vector);
    H5Dclose(dataset_gas_radial_velocity);
    H5Dwrite(dataset_gas_surfacedensity, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, disk_params->gas_surface_density_vector);
    H5Dclose(dataset_gas_surfacedensity);
    H5Dwrite(dataset_gas_pressure, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,disk_params->gas_pressure_vector);
    H5Dclose(dataset_gas_pressure);
    H5Dwrite(dataset_gas_pressure_gradient, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,disk_params->gas_pressure_gradient_vector);
    H5Dclose(dataset_gas_pressure_gradient);

//    fprintf(stderr, "DEBUG: gas_surface_density_vector=%p, first_val=%.3e\n", disk_params->gas_surface_density_vector, disk_params->gas_surface_density_vector[0]);
//    fprintf(stderr, "DEBUG: dust_particle_mass_array=%p, first_val=%.3e\n", particle_data->dust_particle_mass_array, particle_data->dust_particle_mass_array[0][0]);

    H5Sclose(dataspace_gas_id);

    // Dust surface density
    hsize_t dims_dust[1] = {particle_data->allocated_particle_number};
    dataspace_dust_id = H5Screate_simple(1, dims_dust, NULL);
    dataset_dust_surfacedensity = H5Dcreate2(file_id, "/dust/surface_density", H5T_NATIVE_DOUBLE,
                             dataspace_dust_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_dust_surfacedensity, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, particle_data->dust_particle_mass_array);
    H5Dclose(dataset_dust_surfacedensity);
    H5Sclose(dataspace_dust_id);
}


void closeHDF5File(OutputFiles *output_files) {
    if (output_files->hdf5_file != NULL) {
        H5Fclose((hid_t)output_files->hdf5_file);
        output_files->hdf5_file = NULL;
    }
}

