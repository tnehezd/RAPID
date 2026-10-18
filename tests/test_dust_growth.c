#include "config.h"
#include "dust_physics.h"
#include "integrator.h"
#include "simulation_core.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static void testGrowth(int count, int logarithmic)
{
    enum { N = 8, P = 24 };
    double radii[N + 2], gas[N + 2], pressure[N + 2], gradient[N + 2], velocity[N + 2];
    double dust[N], positions[P][2], expected[P];
    DiskParameters disk = {
        .grid_number = N, .stellar_mass = 1, .r_min = 1, .r_max = 4,
        .h_aspect_ratio = 0.05, .alpha_parameter = 0.003,
        .alpha_parameter_modification = 1, .dr_dze_i = 1, .dr_dze_o = 1,
        .particle_density = 1.6, .fragmentation_velocity = 1000,
        .fragmentation_factor = 0.37, .drift_factor = 0.55,
        .radial_grid = radii, .gas_surface_density_vector = gas,
        .gas_pressure_vector = pressure, .gas_pressure_gradient_vector = gradient,
        .gas_velocity_vector = velocity, .delta_r = 3.0 / (N - 1)
    };
    disk.particle_density_dimensionless = 1.6 * pow(AU_IN_CM, 3) / SOLAR_MASS_IN_GRAMS;
    ParticleData particles = {.particle_distance_array = positions, .dust_surfacedensity = dust};
    for (int i = 0; i <= N + 1; ++i) {
        radii[i] = logarithmic ? pow(4, (i - 1.0) / (N - 1)) : 1 + (i - 1) * disk.delta_r;
        gas[i] = 1e-5 * (1 + radii[i]);
        pressure[i] = 2e-8 * (1 + radii[i]);
        gradient[i] = -1e-8 * (1 + radii[i]);
        velocity[i] = -1e-6;
        if (i >= 1 && i <= N) dust[i - 1] = 1e-7 * (1 + 0.1 * radii[i]);
    }
    double dt = 0.01, size = 1e-4 / AU_IN_CM;
    for (int i = 0; i < count; ++i) {
        double r = 1.01 + 2.98 * i / (count - 1.0);
        positions[i][0] = r;
        positions[i][1] = size;
        expected[i] = calculateDustParticleSize(size, 1.6, 1e-5 * (1 + r),
            1e-7 * (1 + 0.1 * r), r, 2e-8 * (1 + r), -1e-8 * (1 + r), dt, &disk);
        assert(expected[i] > size);
    }
    updateParticleSizes(&particles, count, dt, &disk);
    for (int i = 0; i < count; ++i) {
        assert(fabs(positions[i][1] / expected[i] - 1) < 1e-13);
        double next_r, next_size;
        integrateParticleRungeKutta4(positions[i][1], dt, positions[i][0],
                                    &next_r, &next_size, &disk);
        assert(isfinite(next_r));
        assert(next_size == positions[i][1]);
    }
    SimulationOptions options = {.option_for_dust_growth = 1};
    calculateDustDistance("", &particles, dt, 0, count, &options, &disk);
    for (int i = 0; i < count; ++i) {
        assert(positions[i][1] == expected[i]);
    }
}

static void testRadiativeBoundaryCrossing(int logarithmic, int outward)
{
    enum { N = 8 };
    double radii[N + 2], gas[N + 2], gradient[N + 2], velocity[N + 2];
    ThermalState thermal[N + 2] = {0};
    DiskParameters disk = {
        .grid_number = N, .r_min = 1, .r_max = 4, .delta_r = 3.0 / (N - 1),
        .thermodynamics_model = THERMODYNAMICS_RADIATIVE_ENERGY,
        .thermal_states = thermal, .radial_grid = radii,
        .gas_surface_density_vector = gas, .gas_pressure_gradient_vector = gradient,
        .gas_velocity_vector = velocity
    };
    for (int i = 0; i < N + 2; ++i) {
        radii[i] = logarithmic ? pow(4, (i - 1.0) / (N - 1)) : 1 + (i - 1) * disk.delta_r;
        gas[i] = 1e-20;
        gradient[i] = 0;
        velocity[i] = outward ? 1e8 : -1e8;
        thermal[i].sound_speed = 1e5;
        thermal[i].scale_height = .05 * AU_IN_CM;
    }
    double radius, size;
    integrateParticleRungeKutta4(1e-10 / AU_IN_CM, 1, outward ? 3.99 : 1.01,
                                &radius, &size, &disk);
    assert(radius == 0 && size == 0);
    for (int i = 0; i < N + 2; ++i) velocity[i] = outward ? .01 : -.01;
    double initial = outward ? 3.9 : 1.1;
    integrateParticleRungeKutta4(1e-10 / AU_IN_CM, .1, initial, &radius, &size, &disk);
    assert(fabs(radius - (initial + .1 * velocity[1])) < 1e-14);
    assert(size == 1e-10 / AU_IN_CM);

    /* A steep spatial velocity change must be sampled, not jumped over. */
    for (int i = 0; i < N + 2; ++i) velocity[i] = -100 * (radii[i] - 2);
    double coarse, fine, exact = 2 + .8 * exp(-.1);
    integrateParticleRungeKutta4(1e-10 / AU_IN_CM, .001, 2.8, &coarse, &size, &disk);
    integrateParticleRungeKutta4(1e-10 / AU_IN_CM, .0005, 2.8, &fine, &size, &disk);
    integrateParticleRungeKutta4(size, .0005, fine, &fine, &size, &disk);
    double error_ratio = fabs(coarse - exact) / fabs(fine - exact);
    assert(error_ratio > 12 && error_ratio < 20);
    integrateParticleRungeKutta4(1e-10 / AU_IN_CM, .2, 2.8, &radius, &size, &disk);
    assert(fabs(radius - 2) < .01);
    assert(size > 0);
}

static void testDriftTimeStepSampling(void)
{
    enum { N = 4 };
    double r[] = {.5, 1, 2, 3, 4, 5};
    double sigma[] = {1, 1, 1, 1, 1, 1}, gradient[N + 2] = {0};
    double velocity[] = {1e8, 1, 3, 5, 7, 1e8};
    double positions[][2] = {{1.5, 1e-10}}, micron[][2] = {{3.5, 1e-10}};
    DiskParameters disk = {.grid_number = N, .r_min = 1, .r_max = 4,
        .stellar_mass = 1, .h_aspect_ratio = .05, .radial_grid = r,
        .gas_surface_density_vector = sigma, .gas_pressure_gradient_vector = gradient,
        .gas_velocity_vector = velocity};
    ParticleData particles = {.particle_distance_array = positions,
                              .micron_particle_distance_array = micron};
    assert(getMaximumDriftVelocity(&particles, 1, &disk, SnapshotDrift) == 2);
    assert(getMaximumDriftVelocity(&particles, 1, &disk, SnapshotDriftTwoPop) == 6);
}

int main(void)
{
    testDriftTimeStepSampling();
    for (int grid = 0; grid < 2; ++grid) {
        testGrowth(4, grid);
        testGrowth(8, grid);
        testGrowth(24, grid);
        testRadiativeBoundaryCrossing(grid, 0);
        testRadiativeBoundaryCrossing(grid, 1);
    }
    puts("Dust: growth, interpolated timestep velocities, bounded RK4 and radiative absorbing exits passed");
    return 0;
}
