#include <math.h>
#include <stdlib.h>


#include "vertical_settling.h"
#include "vertical_physics.h"
#include "vertical_profile.h"
#include "particle_data.h"
#include "simulation_types.h"
#include "gas_physics.h"

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


/**
 * @brief Apply vertical settling + turbulent diffusion to dust particles.
 *
 * Updates the z-coordinate of each particle in StructuredParticleData according to
 * local settling velocity, Keplerian frequency, and turbulent diffusion.
 *
 * @param sdata Pointer to the 2D structured particle data.
 * @param disk_params Pointer to disk parameters.
 * @param dt Timestep in internal units.
 */
void applyVerticalSettling(StructuredParticleData *sdata, const DiskParameters *disk_params, double dt) {
    if (!sdata || !disk_params) return;

    for (size_t i_r = 0; i_r < sdata->n_r; i_r++) {

        // radial location of this grid column
        double r = sdata->particles[i_r][0].r_au;

        // local turbulent alpha at this radius
        double alpha_local = calculateTurbulentAlpha(r, disk_params);

        for (size_t i_z = 0; i_z < sdata->n_z; i_z++) {

            DustParticle *p = &sdata->particles[i_r][i_z];
            double z = p->z_au;
            double a = p->radius; // particle radius

            // local Stokes number
            double St = calculateLocalStokesNumber(r, z, a, disk_params);

            // vertical settling velocity: v_settle = - St * Omega_k * z
            double omega_k = calculateKeplerianFrequency(r, disk_params);
            double v_settle = - St * omega_k * z;

            // Update z-position with simple explicit Euler
            double z_new = z + v_settle * dt;

            // Turbulent diffusion: random step with stddev ~ sqrt(2*alpha*H^2*dt)
            double Hd = calculateDustScaleHeight(r, a, disk_params);
            double dz_std = sqrt(2.0 * alpha_local * Hd * Hd * dt);
            double rand_diff = ((double)rand() / RAND_MAX - 0.5) * dz_std;
            z_new += rand_diff;

            // store back
            p->z_au = z_new;

            // optional vertical clipping, e.g., ±5 scale heights
            if (p->z_au > 5.0*Hd) p->z_au = 5.0*Hd;
            if (p->z_au < -5.0*Hd) p->z_au = -5.0*Hd;
        }
    }
}