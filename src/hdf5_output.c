#include "hdf5_output.h"
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>


#include "particle_data.h"
#include "config.h"
#include "utils.h"

#include "hdf5_output.h"

#define MAX_TRAPS 5

static hid_t ts_file;
static hid_t ts_group;
static hid_t dset_time, dset_m1, dset_m2, dset_mtot, dset_trap_pos;

void initializeMassTimeSeries(const char *filename)
{
    ts_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    ts_group = H5Gcreate(ts_file, "/trap_evolution",
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dims[2]    = {0, MAX_TRAPS};
    hsize_t maxdims[2] = {H5S_UNLIMITED, MAX_TRAPS};

    hid_t dataspace = H5Screate_simple(2, dims, maxdims);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);

    hsize_t chunk_dims[2] = {128, MAX_TRAPS};
    H5Pset_chunk(plist, 2, chunk_dims);

    dset_time = H5Dcreate(ts_group, "time", H5T_NATIVE_DOUBLE,
                          dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);

    dset_m1 = H5Dcreate(ts_group, "primary_mass", H5T_NATIVE_DOUBLE,
                        dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);

    dset_m2 = H5Dcreate(ts_group, "secondary_mass", H5T_NATIVE_DOUBLE,
                        dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);

    dset_mtot = H5Dcreate(ts_group, "total_mass", H5T_NATIVE_DOUBLE,
                          dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);

    dset_trap_pos = H5Dcreate(ts_group, "trap_position", H5T_NATIVE_DOUBLE,
                              dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);

    H5Sclose(dataspace);
    H5Pclose(plist);
}

void appendMassTimeSeries(double snapshot,
                          double m1[MAX_TRAPS],
                          double m2[MAX_TRAPS],
                          double mtot[MAX_TRAPS],
                          double trap_pos[MAX_TRAPS])
{
    static size_t ts_row = 0;

    hsize_t new_size[2] = { ts_row + 1, MAX_TRAPS };

    H5Dset_extent(dset_time, new_size);
    H5Dset_extent(dset_m1, new_size);
    H5Dset_extent(dset_m2, new_size);
    H5Dset_extent(dset_mtot, new_size);
    H5Dset_extent(dset_trap_pos, new_size);

    hid_t filespace = H5Dget_space(dset_time);

    hsize_t start[2] = { ts_row, 0 };
    hsize_t count[2] = { 1, MAX_TRAPS };

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);

    hid_t memspace = H5Screate_simple(2, count, NULL);

    double time_row[MAX_TRAPS];
    for (int i = 0; i < MAX_TRAPS; i++)
        time_row[i] = snapshot;

    H5Dwrite(dset_time, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, time_row);
    H5Dwrite(dset_m1,   H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, m1);
    H5Dwrite(dset_m2,   H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, m2);
    H5Dwrite(dset_mtot, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, mtot);

    hid_t filespace_trap = H5Dget_space(dset_trap_pos);
    H5Sselect_hyperslab(filespace_trap, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dwrite(dset_trap_pos, H5T_NATIVE_DOUBLE, memspace, filespace_trap, H5P_DEFAULT, trap_pos);

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Sclose(filespace_trap);

    H5Fflush(ts_file, H5F_SCOPE_GLOBAL);

    ts_row++;
}

void closeMassTimeSeries()
{
    H5Dclose(dset_time);
    H5Dclose(dset_m1);
    H5Dclose(dset_m2);
    H5Dclose(dset_mtot);
    H5Dclose(dset_trap_pos);

    H5Gclose(ts_group);
    H5Fclose(ts_file);
}


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

void writeHDF5SnapshotToFile(double time, hid_t file_id, const SimulationOptions *sim_opts,
                             DiskParameters *disk_params, ParticleData *particle_data) {

    (void)sim_opts; // unused

    // -------------------
    // GAS GROUP
    // -------------------
    hid_t group_gas = H5Gcreate2(file_id, "/gas_grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_gas);

    hsize_t dims_gas[1] = { disk_params->grid_number };
    hid_t space_gas = H5Screate_simple(1, dims_gas, NULL);

    hid_t dset_gas_grid = H5Dcreate2(file_id, "/gas_grid/radial_grid", H5T_NATIVE_DOUBLE, space_gas,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_gas_surface = H5Dcreate2(file_id, "/gas_grid/surface_density", H5T_NATIVE_DOUBLE, space_gas,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_gas_pressure = H5Dcreate2(file_id, "/gas_grid/pressure", H5T_NATIVE_DOUBLE, space_gas,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_gas_grad = H5Dcreate2(file_id, "/gas_grid/pressure_gradient", H5T_NATIVE_DOUBLE, space_gas,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_gas_vel = H5Dcreate2(file_id, "/gas_grid/radial_velocity", H5T_NATIVE_DOUBLE, space_gas,
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
    // DUST EULER
    // -------------------
    hid_t group_dust = H5Gcreate2(file_id, "/dust_grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_dust);

    hsize_t dims_grid[1] = { disk_params->grid_number };
    hid_t space_grid = H5Screate_simple(1, dims_grid, NULL);

    hid_t dset_surface = H5Dcreate2(file_id, "/dust_grid/surface_density", H5T_NATIVE_DOUBLE,
                                    space_grid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dset_dust_grid = H5Dcreate2(file_id, "/dust_grid/radial_grid", H5T_NATIVE_DOUBLE,
                                      space_grid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dset_surface, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             particle_data->dust_surfacedensity);
    H5Dwrite(dset_dust_grid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             disk_params->radial_grid);

    H5Dclose(dset_surface);
    H5Dclose(dset_dust_grid);
    H5Sclose(space_grid);


    // --------------------
    // PARTICLE LAGRANGIAN
    // --------------------

    hid_t group_particles = H5Gcreate2(file_id, "/particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_particles);

    hsize_t dims_particles[1] = { particle_data->allocated_particle_number };
    hid_t space_particles = H5Screate_simple(1, dims_particles, NULL);

    hid_t dset_positions = H5Dcreate2(file_id, "/particles/position", H5T_NATIVE_DOUBLE, space_particles,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double *positions = malloc(dims_particles[0] * sizeof(double));
    for (size_t i = 0; i < dims_particles[0]; i++) {
        positions[i] = particle_data->particle_distance_array[i][1] * AU_IN_CM;
    }
    H5Dwrite(dset_positions, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, positions);
    H5Dclose(dset_positions);
    free(positions);

    hid_t dset_size = H5Dcreate2(file_id, "/particles/size", H5T_NATIVE_DOUBLE, space_particles,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double *sizes = malloc(dims_particles[0] * sizeof(double));
    for (size_t i = 0; i < dims_particles[0]; i++) {
        sizes[i] = particle_data->dust_particle_mass_array[i][1];
    }
    H5Dwrite(dset_size, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sizes);
    H5Dclose(dset_size);
    free(sizes);

    // -- original index 0..N-1
    hid_t dset_index = H5Dcreate2(file_id, "/particles/index", H5T_NATIVE_INT, space_particles,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int *indices = malloc(dims_particles[0] * sizeof(int));
    for (size_t i = 0; i < dims_particles[0]; i++) indices[i] = (int)i;
    H5Dwrite(dset_index, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
    H5Dclose(dset_index);
    free(indices);

    H5Sclose(space_particles);


    // -------------------
    // FRAME / TIME
    // -------------------

    hid_t group_frame = H5Gcreate2(file_id, "/frame", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_frame < 0) {
        fprintf(stderr, "ERROR: Could not create /frame group\n");
        return;
    }

    hsize_t dims_time[1] = {1};
    hid_t space_time = H5Screate_simple(1, dims_time, NULL);
    if (space_time < 0) {
        fprintf(stderr, "ERROR: Could not create dataspace for time\n");
        H5Gclose(group_frame);
        return;
    }

    hid_t dset_time = H5Dcreate2(group_frame, "time", H5T_NATIVE_DOUBLE, space_time,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset_time < 0) {
        fprintf(stderr, "ERROR: Could not create /frame/time dataset\n");
        H5Sclose(space_time);
        H5Gclose(group_frame);
        return;
    }

    double t_val = time;
    H5Dwrite(dset_time, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t_val);

    H5Dclose(dset_time);
    H5Sclose(space_time);


    // -------------------
    // FRAME ATTRIBUTE: code_version
    // -------------------

    hid_t attr_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(attr_type, strlen(SIM_VERSION));
    H5Tset_strpad(attr_type, H5T_STR_NULLTERM);

    hid_t attr_space = H5Screate(H5S_SCALAR);

    hid_t attr = H5Acreate2(group_frame,
                            "code_version",
                            attr_type,
                            attr_space,
                            H5P_DEFAULT,
                            H5P_DEFAULT);

    H5Awrite(attr, attr_type, SIM_VERSION);

    H5Aclose(attr);
    H5Sclose(attr_space);
    H5Tclose(attr_type);



    // -------------------
    // FRAME ATTRIBUTE: compile_date
    // -------------------

    hid_t attr_type_date = H5Tcopy(H5T_C_S1);
    H5Tset_size(attr_type_date, strlen(__DATE__));
    H5Tset_strpad(attr_type_date, H5T_STR_NULLTERM);

    hid_t attr_space_date = H5Screate(H5S_SCALAR);

    hid_t attr_date = H5Acreate2(group_frame,
                                 "compile_date",
                                 attr_type_date,
                                 attr_space_date,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);

    H5Awrite(attr_date, attr_type_date, __DATE__);

    H5Aclose(attr_date);
    H5Sclose(attr_space_date);
    H5Tclose(attr_type_date);


    // -------------------
    // FRAME ATTRIBUTE: compile_time
    // -------------------

    hid_t attr_type_time = H5Tcopy(H5T_C_S1);
    H5Tset_size(attr_type_time, strlen(__TIME__));
    H5Tset_strpad(attr_type_time, H5T_STR_NULLTERM);

    hid_t attr_space_time = H5Screate(H5S_SCALAR);

    hid_t attr_time = H5Acreate2(group_frame,
                                 "compile_time",
                                 attr_type_time,
                                 attr_space_time,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);

    H5Awrite(attr_time, attr_type_time, __TIME__);

    H5Aclose(attr_time);
    H5Sclose(attr_space_time);
    H5Tclose(attr_type_time);
    H5Gclose(group_frame);

}

void closeHDF5File(OutputFiles *output_files) {
    if (output_files->hdf5_file != NULL) {
        H5Fclose((hid_t)output_files->hdf5_file);
        output_files->hdf5_file = NULL;
    }
}

