#ifndef VERTICAL_PROFILE_H
#define VERTICAL_PROFILE_H

void calculateVerticalDustProfile(double Hd,
                                  double dust_surface_density,
                                  double *rho_z,
                                  int z_points,
                                  double z_max);

#endif
