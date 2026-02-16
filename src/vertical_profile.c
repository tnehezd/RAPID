#include "vertical_profile.h"
#include <math.h>

void calculateVerticalDustProfile(double Hd,
                                  double dust_surface_density,
                                  double *rho_z,
                                  int z_points,
                                  double z_max)
{
    double dz = 2.0*z_max/(z_points-1);
    double norm = dust_surface_density/(sqrt(2.0*M_PI)*Hd);

    for(int i=0;i<z_points;i++){
        double z = -z_max + i*dz;
        rho_z[i] = norm*exp(-0.5*(z/Hd)*(z/Hd));
    }
}
