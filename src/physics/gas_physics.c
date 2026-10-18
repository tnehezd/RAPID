// src/dust_physics.c
#include "gas_physics.h" 
#include "config.h"       
#include "simulation_types.h" 
#include "boundary_conditions.h"
#include "simulation_core.h" 
#include "utils.h"           
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>       

#include "photoevaporation.h"
#include "disk_thermal.h"
#include "self_gravity.h"
#include "backreaction.h"
#include "logger.h"

#define WIND_CRIT 1e-20  // Critical value for wind profile (below this value, wind is considered zero) [M_sol/AU^2/day]

double calculateTurbulentAlpha(double radial_distance, const DiskParameters *disk_params) {
    
    if (isRadiativeThermodynamics(disk_params->thermodynamics_model)) {
        double inner = disk_params->dr_dze_i > 0.0 ?
            tanh((radial_distance - disk_params->r_dze_i) / disk_params->dr_dze_i) :
            (radial_distance >= disk_params->r_dze_i ? 1.0 : -1.0);
        double outer = disk_params->dr_dze_o > 0.0 ?
            tanh((disk_params->r_dze_o - radial_distance) / disk_params->dr_dze_o) :
            (radial_distance < disk_params->r_dze_o ? 1.0 : -1.0);
        return disk_params->alpha_parameter *
            (1.0 - 0.5 * (1.0 - disk_params->alpha_parameter_modification) * (inner + outer));
    }
    double alpha_r;
    alpha_r = 1.0 - 0.5 * (1.0 - disk_params->alpha_parameter_modification) * (tanh((radial_distance - disk_params->r_dze_i) / disk_params->dr_dze_i) + tanh((disk_params->r_dze_o - radial_distance) / disk_params->dr_dze_o));
    return alpha_r * disk_params->alpha_parameter;
}


double calculateKinematicViscosity(double radial_distance, const DiskParameters *disk_params) {
    
    double gas_viscosity, local_soundspeed, local_pressure_scaleheight;
    local_pressure_scaleheight = calculatePressureScaleHeight(radial_distance,disk_params);
    local_soundspeed = calculateLocalSoundSpeed(radial_distance,disk_params);
    gas_viscosity = calculateTurbulentAlpha(radial_distance, disk_params) * local_soundspeed * local_pressure_scaleheight;
    return gas_viscosity;
}

double calculateKeplerianVelocity(double radial_distance, const DiskParameters *disk_params) {

    return calculateKeplerianFrequency(radial_distance, disk_params) * radial_distance;
}

double calculateKeplerianFrequency(double radial_distance, const DiskParameters *disk_params) {

    double omega_squared = G_DIMENSIONLESS * disk_params->stellar_mass / radial_distance / radial_distance / radial_distance;
    if (disk_params->enable_selfgravity && disk_params->selfgravity_omega2_extra != NULL) {
        omega_squared += interpolateSelfGravityOmega2(radial_distance, disk_params);
        if (omega_squared <= 0.0) omega_squared = 1e-30;
    }
    return sqrt(omega_squared);
}

void calculateGasPressureGradient(DiskParameters *disk_params) {

    int i;
    double temporary_gas_pressure, calculated_pressure_gradient[disk_params->grid_number + 2];

    for (i = 1; i <= disk_params->grid_number; i++) {
        double left_dr = disk_params->radial_grid[i] - disk_params->radial_grid[i - 1];
        double right_dr = disk_params->radial_grid[i + 1] - disk_params->radial_grid[i];
        temporary_gas_pressure = (disk_params->gas_pressure_vector[i + 1] -
                                  disk_params->gas_pressure_vector[i - 1]) /
                                 (left_dr + right_dr);
        calculated_pressure_gradient[i] = temporary_gas_pressure;

    }

    for (i = 1; i <= disk_params->grid_number; i++) {
        disk_params->gas_pressure_gradient_vector[i] = calculated_pressure_gradient[i];
    }

}

double coefficientForGasRadialVelocity(double gas_surface_density, double radial_distance) {

    return -1.0 * (3.0 / (gas_surface_density * sqrt(radial_distance)));
}

void calculateGasRadialVelocity(DiskParameters *disk_params) {

    double temporary_gas_velocity;
    double gas_velocity_vector[disk_params->grid_number + 2];
    double gas_velocity_vectortemp[disk_params->grid_number + 1];

    int i;

    #pragma omp parallel for private(i)
    for (i = 0; i <= disk_params->grid_number + 1; i++) { 
        gas_velocity_vector[i] = disk_params->gas_surface_density_vector[i] * calculateKinematicViscosity(disk_params->radial_grid[i], disk_params) * sqrt(disk_params->radial_grid[i]);
    }

    #pragma omp parallel for private(i, temporary_gas_velocity)
    for (i = 1; i <= disk_params->grid_number; i++) { 
        double left_dr = disk_params->radial_grid[i] - disk_params->radial_grid[i - 1];
        double right_dr = disk_params->radial_grid[i + 1] - disk_params->radial_grid[i];
        temporary_gas_velocity = (gas_velocity_vector[i + 1] -
                                  gas_velocity_vector[i - 1]) /
                                 (left_dr + right_dr);
        gas_velocity_vectortemp[i] = coefficientForGasRadialVelocity(disk_params->gas_surface_density_vector[i], disk_params->radial_grid[i]) * temporary_gas_velocity;
    }

    for (i = 1; i <= disk_params->grid_number; i++) { 
        disk_params->gas_velocity_vector[i] = gas_velocity_vectortemp[i];
    }
}

void refreshGasSurfaceDensityPressurePressureGradient(SimulationOptions *sim_opts, DiskParameters *disk_params) { 
    if (disk_params->thermodynamics_model == THERMODYNAMICS_RADIATIVE_ENERGY) {
        evolveRadiativeEnergy(disk_params, sim_opts);
        return;
    }

    double gas_sigma_dot_viscosity, gas_sigma_dot_viscosity_backwards, gas_sigma_dot_viscosity_forward;
    double gas_surface_density_temp[disk_params->grid_number + 2];
    double gas_velocity_array[disk_params->grid_number + 2]; 
    double dt = sim_opts->user_defined_time_step; // Internal time units.

    int i;

    // Boundary ghost cells setup for the viscous flux array (\Sigma * \nu)
    gas_velocity_array[0] = disk_params->gas_surface_density_vector[0] * calculateKinematicViscosity(disk_params->radial_grid[0], disk_params); 
    gas_velocity_array[disk_params->grid_number + 1] = disk_params->gas_surface_density_vector[disk_params->grid_number + 1] * calculateKinematicViscosity(disk_params->radial_grid[disk_params->grid_number + 1], disk_params); 
    #pragma omp parallel for
    for(i = 1; i <= disk_params->grid_number; i++) { 
        gas_velocity_array[i] = disk_params->gas_surface_density_vector[i] * calculateKinematicViscosity(disk_params->radial_grid[i], disk_params); 
    }

    // =========================================================================
    // 1. VISCOUS TIME EVOLUTION (FTCS scheme for nu * Sigma)
    // =========================================================================
    for (i = 1; i <= disk_params->grid_number; i++) {
        gas_sigma_dot_viscosity           = gas_velocity_array[i];
        gas_sigma_dot_viscosity_backwards = gas_velocity_array[i - 1];
        gas_sigma_dot_viscosity_forward   = gas_velocity_array[i + 1];

        // The FTCS operator acts on viscous_sigma = nu * Sigma.
        double left_dr = disk_params->radial_grid[i] - disk_params->radial_grid[i - 1];
        double right_dr = disk_params->radial_grid[i + 1] - disk_params->radial_grid[i];
        double span_dr = left_dr + right_dr;
        double first_derivative = (gas_sigma_dot_viscosity_forward -
                                   gas_sigma_dot_viscosity_backwards) / span_dr;
        double second_derivative = 2.0 *    
            (gas_sigma_dot_viscosity_backwards / (left_dr * span_dr) -
             gas_sigma_dot_viscosity / (left_dr * right_dr) +
             gas_sigma_dot_viscosity_forward / (right_dr * span_dr));

        double viscosity = calculateKinematicViscosity(disk_params->radial_grid[i], disk_params);
        double viscous_sigma_dot =
            ftcsSecondDerivativeCoefficient(disk_params->radial_grid[i], disk_params) * second_derivative +
            ftcsFirstDerivativeCoefficient(disk_params->radial_grid[i], disk_params) * first_derivative;

        // Convert the evolved nu * Sigma back to physical Sigma before applying sinks.
        if (viscosity == 0.0) {
            gas_surface_density_temp[i] = disk_params->gas_surface_density_vector[i] +
                dt * (3.0 * second_derivative +
                      4.5 / disk_params->radial_grid[i] * first_derivative);
        } else {
            gas_surface_density_temp[i] = (gas_velocity_array[i] + dt * viscous_sigma_dot) /
                                          viscosity;
        }
    }

    // =========================================================================
    // 2. COMPUTE PHOTOEVAPORATION SINK TERM (Anna's module updates disk_params->sigma_dot_photoevap)
    // =========================================================================
    computePhotoevaporationSink(disk_params);



    // =========================================================================
    // 3. APPLY BOTH VISCOUS EVOLUTION AND PHOTOEVAPORATION LOSS COHERENTLY
    // =========================================================================
    for (i = 1; i <= disk_params->grid_number; i++) { 
        gas_surface_density_temp[i] -= dt * disk_params->sigma_dot_photoevap[i];
        if (disk_params->enable_backreaction) {
            if (disk_params->backreaction == NULL) {
                LOG_ERROR("FTCS backreaction correction requested before drag initialization.\n");
                exit(EXIT_FAILURE);
            }
            /* Cached flux contains only u_g - u_nu, on the same old state as FTCS. */
            double correction = disk_params->backreaction->flux[i] -
                                disk_params->backreaction->flux[i - 1];
            if (correction != 0.0) {
                gas_surface_density_temp[i] -= dt * correction / backreactionCellArea(i, disk_params);
            }
            if (!isfinite(gas_surface_density_temp[i]) ||
                gas_surface_density_temp[i] < disk_params->density_floor) {
                LOG_ERROR("Invalid FTCS plus backreaction density at cell %d; reduce timestep.\n", i);
                exit(EXIT_FAILURE);
            }
        }
        if (gas_surface_density_temp[i] < disk_params->density_floor) {
            gas_surface_density_temp[i] = disk_params->density_floor;
        }
    }

    // Boundary conditions are defined for the physical surface density Sigma.
    sim_opts->current_bc_target = 0;
    applySurfaceDensityBoundaryConditions(gas_surface_density_temp, disk_params, sim_opts);

    if (disk_params->thermodynamics_model == THERMODYNAMICS_RADIATIVE) {
        for (i = 0; i <= disk_params->grid_number + 1; ++i) {
            disk_params->gas_surface_density_vector[i] = gas_surface_density_temp[i];
        }
        refreshRadiativeDisk(disk_params, sim_opts);
        return;
    }

    // Commit the physical surface-density state and update gas pressure.
    for (i = 1; i <= disk_params->grid_number; i++) {
        disk_params->gas_surface_density_vector[i] = gas_surface_density_temp[i];

        // Calculate the updated pressure vector using the synchronized density
        disk_params->gas_pressure_vector[i] = calculateEOS(disk_params->gas_surface_density_vector[i], disk_params->radial_grid[i], disk_params);
    }

    // =========================================================================
    // 4. UPDATE DERIVATIVES AND BOUNDARY CONDITIONS
    // =========================================================================
    calculateGasPressureGradient(disk_params);
    sim_opts->current_bc_target = 1;
    applyBoundaryConditions(disk_params->gas_pressure_vector, disk_params, sim_opts);
    sim_opts->current_bc_target = 2;
    applyBoundaryConditions(disk_params->gas_pressure_gradient_vector, disk_params, sim_opts);

    updateSelfGravity(disk_params);
}
