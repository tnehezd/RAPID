#include "vertical_settling.h"
#include "vertical_physics.h"
#include "vertical_profile.h"

void calculateVerticalDistribution(double radial_distance,
                                   double particle_radius,
                                   double dust_surface_density,
                                   const DiskParameters *disk_params,
                                   double *Hd,
                                   double *rho_z_array,
                                   int z_points,
                                   double N_Hd)
{
    *Hd = calculateDustScaleHeight(radial_distance,
                                   particle_radius,
                                   disk_params);

    double z_max = N_Hd * (*Hd);

    calculateVerticalDustProfile(*Hd,
                                 dust_surface_density,
                                 rho_z_array,
                                 z_points,
                                 z_max);
}
