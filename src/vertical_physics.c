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

    double St = calculateStokesNumber(particle_radius, gas_surface_density, disk_params);

    return H_gas * sqrt(alpha_turb/(alpha_turb + St));
}


double calculateLocalStokesNumber(double r, double z, double particle_radius, const DiskParameters *disk_params)
{
    double Hg = calculatePressureScaleHeight(r, disk_params);

    int i_r = (int)((r - disk_params->r_min)/disk_params->delta_r);
    double Sigma_g = disk_params->gas_surface_density_vector[i_r];

    /* midplane gas density */
    double rho_mid = Sigma_g / (sqrt(2.0*M_PI) * Hg);

    /* vertical stratification */
    double rho_local = rho_mid * exp(-(z*z)/(2.0*Hg*Hg));

    /* Epstein stopping time → Stokes number */
    return (M_PI/2.0)
           * disk_params->particle_density_dimensionless
           * particle_radius
           / (rho_local * Hg);
}
