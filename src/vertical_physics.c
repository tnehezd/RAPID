#include "vertical_physics.h"
#include "gas_physics.h"
#include "dust_physics.h"
#include <math.h>

double calculateDustScaleHeight(double radial_distance,
                                double particle_radius,
                                const DiskParameters *disk_params)
{
    double H_gas = calculatePressureScaleHeight(radial_distance, disk_params);
    double alpha_turb = calculateTurbulentAlpha(radial_distance, disk_params);

    double gas_surface_density =
        disk_params->gas_surface_density_vector[
        (int)((radial_distance - disk_params->r_min)/disk_params->delta_r)];

    double St = calculateStokesNumber(particle_radius,
                                      gas_surface_density,
                                      disk_params);

    return H_gas * sqrt(alpha_turb/(alpha_turb + St));
}
