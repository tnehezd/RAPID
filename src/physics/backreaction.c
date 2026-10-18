#include "backreaction.h"

#include "boundary_conditions.h"
#include "config.h"
#include "dust_physics.h"
#include "gas_physics.h"
#include "logger.h"
#include "utils.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

void solveStationaryDrag(double a, double b, double viscous_velocity,
                         double pressure_velocity, double *radial, double *azimuthal)
{
    if (!isfinite(a) || a < 0 || !isfinite(b) || b < 0 ||
        !isfinite(viscous_velocity) || !isfinite(pressure_velocity)) {
        LOG_ERROR("Invalid stationary drag coefficients.\n");
        exit(EXIT_FAILURE);
    }
    double denominator = (1 + a) * (1 + a) + b * b;
    *radial = ((1 + a) * viscous_velocity - b * pressure_velocity) / denominator;
    *azimuthal = (b * viscous_velocity + (1 + a) * pressure_velocity) / (2 * denominator);
    if (!isfinite(*radial) || !isfinite(*azimuthal)) {
        LOG_ERROR("Non-finite stationary drag velocity.\n");
        exit(EXIT_FAILURE);
    }
}

double backreactionCellArea(int i, const DiskParameters *disk)
{
    double left = .5 * (disk->radial_grid[i - 1] + disk->radial_grid[i]);
    double right = .5 * (disk->radial_grid[i] + disk->radial_grid[i + 1]);
    return M_PI * (right * right - left * left);
}

void refreshBackreaction(const ParticleData *particles, int count,
                         DiskParameters *disk, SimulationOptions *options)
{
    int n = disk->grid_number, m = n + 2;
    if (disk->backreaction == NULL) {
        disk->backreaction = calloc(1, sizeof(BackreactionState) + 6 * (size_t)m * sizeof(double));
        if (disk->backreaction == NULL) {
            LOG_ERROR("Could not allocate stationary drag state.\n");
            exit(EXIT_FAILURE);
        }
        BackreactionState *state = disk->backreaction;
        state->cells = n;
        state->loading_a = state->values;
        state->loading_b = state->values + m;
        state->radial = state->values + 2 * m;
        state->azimuthal = state->values + 3 * m;
        state->flux = state->values + 4 * m;
        state->viscous = state->values + 5 * m;
    }
    BackreactionState *state = disk->backreaction;
    if (state->cells != n) {
        LOG_ERROR("Backreaction grid changed after allocation.\n");
        exit(EXIT_FAILURE);
    }
    options->current_bc_target = 0;
    applySurfaceDensityBoundaryConditions(disk->gas_surface_density_vector, disk, options);
    for (int i = 0; i < m; ++i) {
        double r = disk->radial_grid[i], sigma = disk->gas_surface_density_vector[i];
        if (!isfinite(r) || r <= 0 || !isfinite(sigma) || sigma <= 0) {
            LOG_ERROR("Backreaction requires positive radii and gas density at cell %d.\n", i);
            exit(EXIT_FAILURE);
        }
        disk->gas_pressure_vector[i] = calculateEOS(sigma, r, disk);
        state->loading_a[i] = state->loading_b[i] = 0;
    }
    calculateGasPressureGradient(disk);
    options->current_bc_target = 2;
    applyBoundaryConditions(disk->gas_pressure_gradient_vector, disk, options);
    calculateGasRadialVelocity(disk);

    for (int k = 0; k < count; ++k) {
        double r = particles->particle_distance_array[k][0];
        double size = particles->particle_distance_array[k][1];
        double mass = particles->dust_particle_mass_grid[k];
        if (!isfinite(r) || !isfinite(size) || !isfinite(mass) || mass < 0) {
            LOG_ERROR("Invalid backreaction particle %d.\n", k);
            exit(EXIT_FAILURE);
        }
        if (r < disk->r_min || r > disk->r_max || mass == 0) continue;
        if (size <= 0) {
            LOG_ERROR("Nonpositive active backreaction grain size at particle %d.\n", k);
            exit(EXIT_FAILURE);
        }
        double fraction;
        int left = findDustLeftCell(r, disk, &fraction) + 1;
        for (int j = 0; j < 2; ++j) {
            int i = left + j;
            double sigma = disk->gas_surface_density_vector[i];
            double stokes = calculateStokesNumber(size, sigma, disk);
            double epsilon = mass * (j == 0 ? 1 - fraction : fraction) /
                             (backreactionCellArea(i, disk) * sigma);
            if (!isfinite(stokes) || stokes < 0 || !isfinite(epsilon)) {
                LOG_ERROR("Invalid deposited drag loading at cell %d.\n", i);
                exit(EXIT_FAILURE);
            }
            state->loading_a[i] += epsilon / (1 + stokes * stokes);
            state->loading_b[i] += epsilon * stokes / (1 + stokes * stokes);
        }
    }
    state->loading_a[0] = state->loading_a[1];
    state->loading_a[n + 1] = state->loading_a[n];
    state->loading_b[0] = state->loading_b[1];
    state->loading_b[n + 1] = state->loading_b[n];
    for (int i = 0; i < m; ++i) {
        double r = disk->radial_grid[i];
        double cs = calculateLocalSoundSpeed(r, disk);
        double pressure_velocity = cs * calculatePressureScaleHeight(r, disk) *
            disk->gas_pressure_gradient_vector[i] / disk->gas_pressure_vector[i];
        double viscous_velocity;
        if (i == 0 || i == n + 1) {
            int j = i == 0 ? 1 : n;
            viscous_velocity = disk->gas_velocity_vector[j];
        } else {
            viscous_velocity = disk->gas_velocity_vector[i];
        }
        state->viscous[i] = viscous_velocity;
        solveStationaryDrag(state->loading_a[i], state->loading_b[i], viscous_velocity,
                            pressure_velocity, &state->radial[i], &state->azimuthal[i]);
    }
    for (int i = 0; i < m; ++i) disk->gas_velocity_vector[i] = state->radial[i];
    for (int j = 0; j <= n; ++j) {
        double face = .5 * (disk->radial_grid[j] + disk->radial_grid[j + 1]);
        double velocity = .5 * (state->radial[j] + state->radial[j + 1]);
        int donor = velocity >= 0 ? j : j + 1;
        state->flux[j] = 2 * M_PI * face * velocity * disk->gas_surface_density_vector[donor];
        if ((j == 0 && options->inner_boundary_condition_type == 4) ||
            (j == n && options->outer_boundary_condition_type == 4)) state->flux[j] = 0;
        if (!isfinite(state->flux[j])) {
            LOG_ERROR("Non-finite backreaction gas flux at face %d.\n", j);
            exit(EXIT_FAILURE);
        }
    }
}

double backreactionDustVelocity(double radius, double size, double sigma,
                                const DiskParameters *disk)
{
    if (disk->backreaction == NULL) {
        LOG_ERROR("Dust drift requested before stationary drag initialization.\n");
        exit(EXIT_FAILURE);
    }
    double radial, azimuthal, a, b, viscous, gradient;
    linearInterpolation(disk->backreaction->loading_a, disk->radial_grid, radius,
                        &a, disk->delta_r, disk);
    linearInterpolation(disk->backreaction->loading_b, disk->radial_grid, radius,
                        &b, disk->delta_r, disk);
    linearInterpolation(disk->backreaction->viscous, disk->radial_grid, radius,
                        &viscous, disk->delta_r, disk);
    linearInterpolation(disk->gas_pressure_gradient_vector, disk->radial_grid, radius,
                        &gradient, disk->delta_r, disk);
    double pressure_velocity = calculateLocalSoundSpeed(radius, disk) *
        calculatePressureScaleHeight(radius, disk) * gradient / calculateEOS(sigma, radius, disk);
    solveStationaryDrag(a, b, viscous, pressure_velocity, &radial, &azimuthal);
    double stokes = calculateStokesNumber(size, sigma, disk);
    double velocity = (radial + 2 * stokes * azimuthal) / (1 + stokes * stokes);
    if (!isfinite(velocity)) {
        LOG_ERROR("Non-finite backreaction dust velocity at %.9g AU.\n", radius);
        exit(EXIT_FAILURE);
    }
    return velocity;
}

double backreactionTimeStep(const ParticleData *particles, int count,
                            const DiskParameters *disk, const SimulationOptions *options)
{
    const BackreactionState *state = disk->backreaction;
    double dt = DBL_MAX;
    for (int i = 1; i <= disk->grid_number; ++i) {
        double r = disk->radial_grid[i], a = state->loading_a[i], b = state->loading_b[i];
        double denominator = (1 + a) * (1 + a) + b * b;
        /* Both viscous drift and pressure feedback contain density gradients. */
        double diffusivity = (3 * calculateKinematicViscosity(r, disk) * (1 + a) +
            calculateLocalSoundSpeed(r, disk) * calculatePressureScaleHeight(r, disk) * b) /
            denominator;
        if (!isfinite(diffusivity) || diffusivity < 0) {
            LOG_ERROR("Invalid backreaction diffusivity at cell %d.\n", i);
            exit(EXIT_FAILURE);
        }
        double dr = fmin(r - disk->radial_grid[i - 1], disk->radial_grid[i + 1] - r);
        if (diffusivity > 0) dt = fmin(dt, .2 * dr * dr / (2 * diffusivity));
        double outgoing = fmax(state->flux[i], 0) + fmax(-state->flux[i - 1], 0);
        if (outgoing > 0) {
            double available = disk->gas_surface_density_vector[i] - disk->density_floor;
            if (available <= 0) {
                LOG_ERROR("Backreaction gas transport reached the density floor at cell %d; floor mass injection is not supported.\n", i);
                exit(EXIT_FAILURE);
            }
            dt = fmin(dt, .2 * backreactionCellArea(i, disk) * available / outgoing);
        }
    }
    for (int k = 0; k < count; ++k) {
        double r = particles->particle_distance_array[k][0];
        if (r < disk->r_min || r > disk->r_max) continue;
        double fraction, sigma;
        int left = findDustLeftCell(r, disk, &fraction) + 1;
        linearInterpolation(disk->gas_surface_density_vector, disk->radial_grid, r,
                            &sigma, disk->delta_r, disk);
        double velocity = backreactionDustVelocity(r, particles->particle_distance_array[k][1], sigma, disk);
        double dr = disk->radial_grid[left + 1] - disk->radial_grid[left];
        if (velocity != 0) dt = fmin(dt, .2 * dr / fabs(velocity));
        if (options->option_for_dust_growth == 1) {
            double dust = (1 - fraction) * particles->dust_surfacedensity[left - 1] +
                          fraction * particles->dust_surfacedensity[left];
            if (dust > 0) dt = fmin(dt, .1 * calculateGrowthTimescale(r, dust / sigma, disk));
        }
    }
    return dt;
}

void advanceBackreactionGas(double dt, DiskParameters *disk, SimulationOptions *options)
{
    double next[disk->grid_number + 2];
    for (int i = 1; i <= disk->grid_number; ++i) {
        next[i] = disk->gas_surface_density_vector[i] -
            dt * (disk->backreaction->flux[i] - disk->backreaction->flux[i - 1]) /
            backreactionCellArea(i, disk);
        if (!isfinite(next[i]) || next[i] < disk->density_floor) {
            LOG_ERROR("Invalid backreaction density update at cell %d; reduce timestep.\n", i);
            exit(EXIT_FAILURE);
        }
    }
    for (int i = 1; i <= disk->grid_number; ++i) disk->gas_surface_density_vector[i] = next[i];
    options->current_bc_target = 0;
    applySurfaceDensityBoundaryConditions(disk->gas_surface_density_vector, disk, options);
}
