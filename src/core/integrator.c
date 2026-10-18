#include "integrator.h"

#include "dust_physics.h"
#include "logger.h"
#include "simulation_core.h"
#include "utils.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

static double driftAt(double radius, double size, const DiskParameters *disk)
{
    double sigma, gradient, gas_velocity, drift;
    linearInterpolation(disk->gas_surface_density_vector, disk->radial_grid,
                        radius, &sigma, disk->delta_r, disk);
    linearInterpolation(disk->gas_pressure_gradient_vector, disk->radial_grid,
                        radius, &gradient, disk->delta_r, disk);
    linearInterpolation(disk->gas_velocity_vector, disk->radial_grid,
                        radius, &gas_velocity, disk->delta_r, disk);
    calculate1DDustDrift(size, gradient, sigma, gas_velocity, radius, &drift, disk);
    if (!isfinite(drift)) {
        LOG_ERROR("Non-finite particle drift at r=%.9g AU, size=%.9g AU.\n", radius, size);
        exit(EXIT_FAILURE);
    }
    return drift;
}

static int validStage(double radius, double start, double dr, const DiskParameters *disk)
{
    return isfinite(radius) && radius > disk->r_min && radius < disk->r_max &&
           fabs(radius - start) <= .5 * dr;
}

static int reversesDrift(double initial, double trial)
{
    /* An autonomous 1D trajectory cannot cross a stationary drift root. */
    return (initial < 0 && trial > 0) || (initial > 0 && trial < 0);
}

void integrateParticleRungeKutta4(double particle_radius, double actual_timestep,
                                 double particle_distance, double *particle_distance_new,
                                 double *particle_radius_new, const DiskParameters *disk)
{
    if (!isfinite(actual_timestep) || actual_timestep < 0 ||
        !isfinite(particle_distance) || !isfinite(particle_radius) || particle_radius <= 0) {
        LOG_ERROR("Invalid particle RK4 input: r=%.9g, size=%.9g, dt=%.9g.\n",
                  particle_distance, particle_radius, actual_timestep);
        exit(EXIT_FAILURE);
    }
    double radius = particle_distance, remaining = actual_timestep;
    *particle_radius_new = particle_radius;
    *particle_distance_new = radius;
    if (radius <= disk->r_min || radius >= disk->r_max) {
        *particle_distance_new = *particle_radius_new = 0;
        return;
    }
    for (int step = 0; remaining > 0; ++step) {
        if (step >= 1000000) {
            LOG_ERROR("Particle RK4 exceeded substep limit at r=%.9g AU, remaining dt=%.9g.\n",
                      radius, remaining);
            exit(EXIT_FAILURE);
        }
        double fraction;
        int left = findDustLeftCell(radius, disk, &fraction) + 1;
        double dr = disk->radial_grid[left + 1] - disk->radial_grid[left];
        double k1 = driftAt(radius, particle_radius, disk);
        double distance = k1 < 0 ? radius - disk->r_min : disk->r_max - radius;
        double boundary = k1 < 0 ? disk->r_min : disk->r_max;
        if (k1 != 0 && distance <= 32 * DBL_EPSILON * fmax(1, fabs(boundary))) {
            /* Absorbing crossing localized to floating-point radial precision. */
            *particle_distance_new = *particle_radius_new = 0;
            return;
        }
        double h = remaining;
        if (k1 != 0) {
            h = fmin(h, .2 * dr / fabs(k1));
            h = fmin(h, .8 * distance / fabs(k1));
        }
        double next;
        for (;;) {
            if (!isfinite(h) || h <= 0) {
                LOG_ERROR("Particle RK4 cannot advance at r=%.9g AU.\n", radius);
                exit(EXIT_FAILURE);
            }
            double r2 = radius + .5 * h * k1;
            if (!validStage(r2, radius, dr, disk)) { h *= .5; continue; }
            double k2 = driftAt(r2, particle_radius, disk);
            if (reversesDrift(k1, k2)) { h *= .5; continue; }
            double r3 = radius + .5 * h * k2;
            if (!validStage(r3, radius, dr, disk)) { h *= .5; continue; }
            double k3 = driftAt(r3, particle_radius, disk);
            if (reversesDrift(k1, k3)) { h *= .5; continue; }
            double r4 = radius + h * k3;
            if (!validStage(r4, radius, dr, disk)) { h *= .5; continue; }
            double k4 = driftAt(r4, particle_radius, disk);
            if (reversesDrift(k1, k4)) { h *= .5; continue; }
            next = radius + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            if (!validStage(next, radius, dr, disk)) { h *= .5; continue; }
            if (reversesDrift(k1, driftAt(next, particle_radius, disk))) { h *= .5; continue; }
            break;
        }
        radius = next;
        remaining -= h;
    }
    *particle_distance_new = radius;
}
