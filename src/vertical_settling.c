#include "vertical_settling.h"
#include <math.h>
#include <stdio.h>
#include "dust_physics.h" 


double calculateDustScaleHeight(double radial_distance, double particle_radius, const DiskParameters *disk_params) {

    double H_gas = calculatePressureScaleHeight(radial_distance, disk_params);
    double alpha_turb = calculateTurbulentAlpha(radial_distance, disk_params);
    double gas_surface_density = disk_params->gas_surface_density_vector[(int)((radial_distance - disk_params->r_min)/disk_params->delta_r)]; 

    double St = calculateStokesNumber(particle_radius, gas_surface_density, disk_params);

    double Hd = H_gas * sqrt(alpha_turb / (alpha_turb + St));

    fprintf(stderr, "Hg %lg  alpha %lg  sd %lg  St %lg  Hd%lg\n", H_gas, alpha_turb, gas_surface_density, St, Hd);
    return Hd;
}

void calculateVerticalDustProfile(double Hd, double dust_surface_density, double *rho_z, int z_points, double z_max) {
    int i;
    double dz = 2.0 * z_max / (z_points - 1);
    double norm = dust_surface_density / (sqrt(2.0 * M_PI) * Hd);

    for(i = 0; i < z_points; i++){
        double z = -z_max + i * dz;
        rho_z[i] = norm * exp(-0.5 * (z/Hd) * (z/Hd));
    }
}

void calculateVerticalDistribution(double radial_distance, double particle_radius, double dust_surface_density, const DiskParameters *disk_params, double *Hd, double *rho_z_array, int z_points, double z_max) {
    *Hd = calculateDustScaleHeight(radial_distance, particle_radius, disk_params);
    calculateVerticalDustProfile(*Hd, dust_surface_density, rho_z_array, z_points, z_max);
}
