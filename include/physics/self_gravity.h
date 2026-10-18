#ifndef SELF_GRAVITY_H
#define SELF_GRAVITY_H

#include "simulation_types.h"

/**
 * @brief Updates the disk self-gravity contribution to the rotation curve.
 *
 * Sums the softened radial gravitational pull of all gas rings (complete elliptic
 * integral kernel, softening 0.6 H) and stores -a_r/r in
 * DiskParameters::selfgravity_omega2_extra, so that
 * Omega^2 = Omega_K^2 + selfgravity_omega2_extra.
 * Does nothing if self-gravity is disabled.
 */
void updateSelfGravity(DiskParameters *disk_params);

/* Linear interpolation of the initialized SG contribution to Omega^2 in code
 * units, with endpoint values used outside the radial grid. */
double interpolateSelfGravityOmega2(double radial_distance, const DiskParameters *disk_params);

#endif
