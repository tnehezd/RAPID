// src/disk_model.h

#ifndef DISK_MODEL_H
#define DISK_MODEL_H


#include "config.h"
#include "simulation_types.h" // Feltételezve, hogy itt található a disk_t definíciója
#include "globals.h"

// Funkciódeklarációk a disk_model.c-ből
void read_disk_parameters(disk_t *disk_params);

void initialize_grid_cells(disk_t *disk_params);

void initial_gas_surface_density_profile(disk_t *disk_params);

void initial_gas_pressure_profile(disk_t *disk_params);

void initial_gas_pressure_gradient_profile(disk_t *disk_params);

void initial_gas_velocity_profile(disk_t *disk_params);

/*	Lokalis viszkozitas erteke	*/
double calculate_gas_viscosity(double r, const disk_t *disk_params);

/*	local scale height	*/
double calculate_scale_height(double r, const disk_t *disk_params);

/*	lokális kepleri sebesség	*/
double calculate_keplerian_velocity(double r, const disk_t *disk_params);

/*	lokalis kepleri korfrekvencia	*/
double calculate_keplerian_angular_velocity(double r, const disk_t *disk_params);

/*	local sound speed		*/
double calculate_local_sound_speed(double r, const disk_t *disk_params);

/*	Suruseg a midplane-ben	*/
double calculate_midplane_gas_density(double sigma, double r, const disk_t *disk_params);

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double calculate_gas_pressure(double sigma, double r, const disk_t *disk_params);

/*	a nyomas derivaltja	*/
void calculate_gas_pressure_gradient(disk_t *disk_params);


#endif // DISK_MODEL_H
