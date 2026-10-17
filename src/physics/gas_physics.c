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

#define WIND_CRIT 1e-20  // Critical value for wind profile (below this value, wind is considered zero) [M_sol/AU^2/day]

double calculateTurbulentAlpha(double radial_distance, const DiskParameters *disk_params) {
    
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

double calculatePressureScaleHeight(double radial_distance, const DiskParameters *disk_params) {

    double local_pressure_scaleheight = pow(radial_distance, 1. + disk_params->flaring_index) * disk_params->h_aspect_ratio;
    return local_pressure_scaleheight;
}

double calculateKeplerianVelocity(double radial_distance, const DiskParameters *disk_params) {

    return calculateKeplerianFrequency(radial_distance, disk_params) * radial_distance;
}

// Linear interpolation of the cached self-gravity term of Omega^2 onto an arbitrary radius.
static double interpolateSelfGravityOmega2(double radial_distance, const DiskParameters *disk_params) {

    const double *grid = disk_params->radial_grid;
    const double *extra = disk_params->selfgravity_omega2_extra;
    int low = 0;
    int high = disk_params->grid_number + 1;

    if (radial_distance <= grid[low]) return extra[low];
    if (radial_distance >= grid[high]) return extra[high];

    while (high - low > 1) {
        int middle = (low + high) / 2;
        if (grid[middle] <= radial_distance) low = middle; else high = middle;
    }
    double fraction = (radial_distance - grid[low]) / (grid[high] - grid[low]);
    return extra[low] + fraction * (extra[high] - extra[low]);
}

double calculateKeplerianFrequency(double radial_distance, const DiskParameters *disk_params) {

    double omega_squared = G_DIMENSIONLESS * disk_params->stellar_mass / radial_distance / radial_distance / radial_distance;
    if (disk_params->enable_selfgravity && disk_params->selfgravity_omega2_extra != NULL) {
        omega_squared += interpolateSelfGravityOmega2(radial_distance, disk_params);
        if (omega_squared <= 0.0) omega_squared = 1e-30;
    }
    return sqrt(omega_squared);
}

// c_s = Omega * H with the prescribed H(r); with self-gravity Omega includes the disk contribution,
// so viscosity, pressure and the midplane density feel the modified rotation curve.
double calculateLocalSoundSpeed(double radial_distance, const DiskParameters *disk_params) {
    
    return calculateKeplerianFrequency(radial_distance,disk_params) * calculatePressureScaleHeight(radial_distance,disk_params);
}

double calcualteMidplaneGasDensity(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
   
    return 1. / sqrt(2.0 * M_PI) * gas_surface_density / calculatePressureScaleHeight(radial_distance,disk_params);
}

double calculateGasPressure(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
  
    return calcualteMidplaneGasDensity(gas_surface_density, radial_distance, disk_params) * calculateLocalSoundSpeed(radial_distance,disk_params) * calculateLocalSoundSpeed(radial_distance, disk_params);
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
        gas_surface_density_temp[i] = (gas_velocity_array[i] + dt * viscous_sigma_dot) /
                                      viscosity;
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
        if (gas_surface_density_temp[i] < disk_params->density_floor) {
            gas_surface_density_temp[i] = disk_params->density_floor;
        }
    }

    // Boundary conditions are defined for the physical surface density Sigma.
    sim_opts->current_bc_target = 0;
    applyBoundaryConditions(gas_surface_density_temp, disk_params, sim_opts);

    // Commit the physical surface-density state and update gas pressure.
    for (i = 1; i <= disk_params->grid_number; i++) {
        disk_params->gas_surface_density_vector[i] = gas_surface_density_temp[i];

        // Calculate the updated pressure vector using the synchronized density
        disk_params->gas_pressure_vector[i] = calculateGasPressure(disk_params->gas_surface_density_vector[i], disk_params->radial_grid[i], disk_params);
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

// Softening length of the ring kernel in units of the local pressure scale height (Mueller et al. 2012).
#define SELFGRAVITY_SOFTENING_SCALEHEIGHT_FACTOR 0.6

// Complete elliptic integrals K(k), E(k) by the arithmetic-geometric mean.
static void completeEllipticIntegrals(double k_squared, double *K, double *E) {

    double a = 1.0;
    double b = sqrt(1.0 - k_squared);
    double c_squared_sum = 0.5 * k_squared;
    double power_of_two = 0.5;

    for (int n = 0; n < 40 && fabs(a - b) > 1e-15 * a; n++) {
        double a_next = 0.5 * (a + b);
        double c = 0.5 * (a - b);
        b = sqrt(a * b);
        a = a_next;
        power_of_two *= 2.0;
        c_squared_sum += power_of_two * c * c;
    }
    *K = M_PI / (2.0 * a);
    *E = *K * (1.0 - c_squared_sum);
}

// Radial gravitational acceleration at radius R of a softened ring of mass m and radius a (negative = inward).
static double ringRadialAcceleration(double R, double a, double m, double softening) {

    double s_squared = (R + a) * (R + a) + softening * softening;
    double s = sqrt(s_squared);
    double k_squared = 4.0 * R * a / s_squared;
    double K, E;
    completeEllipticIntegrals(k_squared, &K, &E);

    double k = sqrt(k_squared);
    double dK_dk = E / (k * (1.0 - k_squared)) - K / k;
    double dk2_dR = 4.0 * a / s_squared - 8.0 * R * a * (R + a) / (s_squared * s_squared);
    double dk_dR = dk2_dR / (2.0 * k);
    double d_KoverS_dR = dK_dk * dk_dR / s - K * (R + a) / (s_squared * s);

    return (2.0 * G_DIMENSIONLESS * m / M_PI) * d_KoverS_dR;
}

void updateSelfGravity(DiskParameters *disk_params) {

    if (!disk_params->enable_selfgravity) return;

    int n = disk_params->grid_number;
    if (disk_params->selfgravity_omega2_extra == NULL) {
        disk_params->selfgravity_omega2_extra = (double *)calloc(n + 2, sizeof(double));
        if (disk_params->selfgravity_omega2_extra == NULL) {
            fprintf(stderr, "Failed to allocate memory for the self-gravity array.\n");
            exit(1);
        }
    }

    const double *r = disk_params->radial_grid;
    const double *sigma = disk_params->gas_surface_density_vector;
    double *extra = disk_params->selfgravity_omega2_extra;

    double ring_mass[n + 2];
    for (int j = 1; j <= n; j++) {
        ring_mass[j] = 2.0 * M_PI * r[j] * sigma[j] * 0.5 * (r[j + 1] - r[j - 1]);
    }

    #pragma omp parallel for
    for (int i = 0; i <= n + 1; i++) {
        double softening = SELFGRAVITY_SOFTENING_SCALEHEIGHT_FACTOR * calculatePressureScaleHeight(r[i], disk_params);
        double radial_acceleration = 0.0;
        for (int j = 1; j <= n; j++) {
            radial_acceleration += ringRadialAcceleration(r[i], r[j], ring_mass[j], softening);
        }
        extra[i] = -radial_acceleration / r[i];
    }
}
