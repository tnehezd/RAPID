#ifndef EOS_H
#define EOS_H

#include "simulation_types.h"

/* Thermal closure in simulation code units. Lociso prescribes H(r) and uses
 * cs = Omega*H; radiative mode reads H and cs from the equilibrium cache.
 * Internal-energy evolution, when enabled, is handled by energy_equation.c. */
double calculatePressureScaleHeight(double radial_distance, const DiskParameters *disk_params);
double calculateLocalSoundSpeed(double radial_distance, const DiskParameters *disk_params);

/* Gaussian vertical closure: rho_mid = Sigma/(sqrt(2*pi)*H).
 * The existing public spelling is retained for compatibility. */
double calcualteMidplaneGasDensity(double gas_surface_density, double radial_distance,
                                  const DiskParameters *disk_params);

/* Selects the lociso or radiative closure and returns midplane pressure,
 * not vertically integrated pressure: P_mid = rho_mid*cs^2. */
double calculateEOS(double gas_surface_density, double radial_distance,
                    const DiskParameters *disk_params);

#endif
