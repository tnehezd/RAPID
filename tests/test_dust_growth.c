#include "config.h"
#include "dust_physics.h"
#include "eos.h"
#include "gas_physics.h"
#include "integrator.h"
#include "simulation_core.h"
#include "utils.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static void testDriftFragmentationFactor(void)
{
    DiskParameters disk = {
        .stellar_mass = 1, .h_aspect_ratio = 0.05,
        .fragmentation_velocity = 1000, .fragmentation_factor = 1
    };
    const double sigma = 1e-5, r = 2, pressure = 2e-8, gradient = -1e-8;
    const double density = 1.6;
    double baseline = calculateDriftInducedFragmentationBarrier(
        sigma, r, pressure, gradient, density, &disk);
    double sound_speed = calculateLocalSoundSpeed(r, &disk);
    double expected = (disk.fragmentation_velocity * CM_PER_SEC_TO_AU_PER_YEAR_2PI)
        * calculateKeplerianVelocity(r, &disk)
        / fabs((r * gradient / pressure) * sound_speed * sound_speed * 0.5)
        * 2.0 * (sigma / SURFACE_DENSITY_CONVERSION_FACTOR) / (M_PI * density);
    assert(isfinite(baseline) && baseline > 0);
    assert(fabs(baseline / expected - 1) < 1e-13);
    disk.fragmentation_factor = 0.37;
    double scaled = calculateDriftInducedFragmentationBarrier(
        sigma, r, pressure, gradient, density, &disk);
    assert(fabs(scaled / baseline - 0.37) < 1e-13);
    disk.fragmentation_factor = 1;
    assert(calculateDriftInducedFragmentationBarrier(
        sigma, r, pressure, gradient, density, &disk) == baseline);
}

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

static void testFixedStokes(int logarithmic)
{
    enum { N = 8 };
    double r[N + 2], sigma[N + 2], gradient[N + 2], velocity[N + 2];
    double positions[][2] = {{1.5, 1e-10}, {3, 1e-10}, {0, 0}};
    DiskParameters disk = {
        .grid_number = N, .r_min = 1, .r_max = 4, .stellar_mass = 1,
        .h_aspect_ratio = .05, .fixed_stokes = .5,
        .particle_density_dimensionless = 2,
        .radial_grid = r, .gas_surface_density_vector = sigma,
        .gas_pressure_gradient_vector = gradient, .gas_velocity_vector = velocity,
        .delta_r = 3.0 / (N - 1)
    };
    ParticleData particles = {.particle_distance_array = positions};
    for (int i = 0; i < N + 2; ++i) {
        r[i] = logarithmic ? pow(4, (i - 1.0) / (N - 1)) : 1 + (i - 1) * disk.delta_r;
        sigma[i] = 1 + r[i];
        gradient[i] = 0;
        velocity[i] = r[i];
    }
    updateFixedStokesSizes(&particles, 3, &disk);
    for (int i = 0; i < 2; ++i) {
        assert(fabs(calculateStokesNumber(positions[i][1], 1 + positions[i][0], &disk) - .5) < 1e-14);
    }
    assert(positions[2][0] == 0 && positions[2][1] == 0);
    double next, size;
    integrateParticleRungeKutta4(positions[0][1], .1, positions[0][0], &next, &size, &disk);
    assert(fabs(next - 1.5 * exp(.1 / 1.25)) < 1e-6);
    positions[0][0] = next;
    for (int i = 0; i < N + 2; ++i) sigma[i] *= 2;
    updateFixedStokesSizes(&particles, 3, &disk);
    assert(fabs(calculateStokesNumber(positions[0][1], 2 * (1 + next), &disk) - .5) < 1e-14);
    assert(positions[0][1] > size * 2);
}

static double integratedDustMass(const double *surface, const double *r, int n)
{
    double mass = 0;
    for (int i = 0; i < n; ++i)
        mass += surface[i] * 2 * M_PI * r[i + 1] * (r[i + 2] - r[i + 1]);
    return mass;
}

static void testDustBoundaryMass(int outward, int logarithmic)
{
    enum { N = 8, P = 4 };
    double r[N + 2], gas[N + 2], gradient[N + 2] = {0}, velocity[N + 2];
    DiskParameters disk = {
        .grid_number = N, .r_min = 1, .r_max = 4, .stellar_mass = 1,
        .h_aspect_ratio = .05, .delta_r = 3.0 / (N - 1),
        .radial_grid = r, .gas_surface_density_vector = gas,
        .gas_pressure_gradient_vector = gradient, .gas_velocity_vector = velocity
    };
    for (int i = 0; i < N + 2; ++i) {
        r[i] = logarithmic ? pow(4, (i - 1.0) / (N - 1)) : 1 + (i - 1) * disk.delta_r;
        gas[i] = 1;
        velocity[i] = outward ? 1 : -1;
    }
    ParticleData particles;
    assert(allocateParticleData(&particles, P, N, 1) == 0);
    SimulationOptions opts = {.option_for_dust_secondary_population = 1,
                              .gaussian_sigma = 1, .gaussian_cutoff = 3};
    double boundary = outward ? 4 : 1;
    double near = outward ? 3.99 : 1.01;
    double primary[] = {near, 2, boundary, 0};
    double micron[] = {2, near, 2.5, 3};
    for (int i = 0; i < P; ++i) {
        particles.particle_distance_array[i][0] = primary[i];
        particles.particle_distance_array[i][1] = 1e-10;
        particles.micron_particle_distance_array[i][0] = micron[i];
        particles.micron_particle_distance_array[i][1] = 1e-10;
        particles.dust_particle_mass_grid[i] = 2;
        particles.massmicradial_grid[i] = 3;
        particles.dust_particle_mass_array[i][0] = 2;
        particles.micron_dust_particle_mass_array[i][0] = 3;
    }
    particle_number = P;
    calculateDustDistance("", &particles, .1, 0, P, &opts, &disk);
    assert(particles.particle_distance_array[0][0] == 0);
    assert(particles.dust_particle_mass_grid[0] == 0);
    assert(particles.dust_particle_mass_array[0][0] == 0);
    assert(particles.micron_particle_distance_array[0][0] > 1);
    assert(particles.massmicradial_grid[0] == 3);
    assert(particles.micron_particle_distance_array[1][0] == 0);
    assert(particles.massmicradial_grid[1] == 0);
    assert(particles.micron_dust_particle_mass_array[1][0] == 0);
    assert(particles.particle_distance_array[1][0] > 1);
    assert(particles.dust_particle_mass_grid[1] == 2);
    assert(particles.dust_particle_mass_grid[2] == 0);
    assert(particles.massmicradial_grid[2] == 3);
    assert(particles.dust_particle_mass_grid[3] == 0);
    assert(particles.massmicradial_grid[3] == 3);
    calculateDustDistance("", &particles, .01, .1, P, &opts, &disk);
    assert(particles.micron_particle_distance_array[0][0] > 1);
    assert(particles.particle_distance_array[1][0] > 1);
    updateParticleGridIndices(&particles, .11, P, &disk, 1);
    for (int i = 1; i < P; ++i) particles.massmicradial_grid[i] = 0;
    applyCoagulationMassTransfer(&particles, &disk, &opts, .01);
    assert(particles.massmicradial_grid[0] == 3);
    assert(particles.dust_particle_mass_grid[0] == 0);
    particles.massmicradial_grid[2] = 3;
    particles.massmicradial_grid[3] = 3;
    for (int mode = 0; mode < 4; ++mode) {
        opts.dust_smoothing_mode = mode;
        calculateDustSurfaceDensity(&particles, &opts, &disk);
        double before[N], before_mic[N];
        for (int j = 0; j < N; ++j) {
            before[j] = particles.dust_surfacedensity[j];
            before_mic[j] = particles.micron_dust_surfacedensity[j];
        }
        if (mode == 0 || mode == 1) {
            assert(fabs(integratedDustMass(before, r, N) - 2) < 1e-12);
            assert(fabs(integratedDustMass(before_mic, r, N) - 9) < 1e-12);
        }
        particles.dust_particle_mass_grid[0] = 17;
        particles.massmicradial_grid[1] = 19;
        calculateDustSurfaceDensity(&particles, &opts, &disk);
        for (int j = 0; j < N; ++j) {
            assert(particles.dust_surfacedensity[j] == before[j]);
            assert(particles.micron_dust_surfacedensity[j] == before_mic[j]);
        }
        particles.dust_particle_mass_grid[0] = 0;
        particles.massmicradial_grid[1] = 0;
    }
    freeParticleData(&particles);
    particle_number = 0;
}

static void testDustBoundaryMapping(void)
{
    enum { N = 4, P = 6 };
    double r[] = {0, 1, 2, 3, 4, 5};
    DiskParameters disk = {.grid_number = N, .r_min = 1, .r_max = 4, .radial_grid = r};
    ParticleData particles;
    assert(allocateParticleData(&particles, P, N, 1) == 0);
    SimulationOptions opts = {.option_for_dust_secondary_population = 1};
    double positions[] = {0, .5, 1, 2, 4, 4.5};
    for (int i = 0; i < P; ++i) {
        particles.particle_distance_array[i][0] = positions[i];
        particles.micron_particle_distance_array[i][0] = positions[P - 1 - i];
        particles.dust_particle_mass_grid[i] = 1;
        particles.massmicradial_grid[i] = 2;
    }
    particle_number = P;
    for (int mode = 0; mode <= 1; ++mode) {
        opts.dust_smoothing_mode = mode;
        calculateDustSurfaceDensity(&particles, &opts, &disk);
        assert(fabs(integratedDustMass(particles.dust_surfacedensity, r, N) - 3) < 1e-12);
        assert(fabs(integratedDustMass(particles.micron_dust_surfacedensity, r, N) - 6) < 1e-12);
        assert(fabs(particles.dust_surfacedensity[0] * 2 * M_PI - 1) < 1e-12);
        assert(fabs(particles.dust_surfacedensity[3] * 8 * M_PI - 1) < 1e-12);
    }
    freeParticleData(&particles);
    particle_number = 0;
}

static void testCoagulationGasIndex(int logarithmic, int cell)
{
    enum { N = 8 };
    double r[N+2], gas[N+2];
    DiskParameters disk = {
        .grid_number = N, .r_min = 1, .r_max = 4, .stellar_mass = 1,
        .radial_grid = r, .gas_surface_density_vector = gas
    };
    for (int i = 0; i < N+2; ++i) {
        r[i] = logarithmic ? pow(4., (i-1.)/(N-1)) : 1+(i-1.)*3/(N-1);
        gas[i] = 1e-3*(i+1);
    }
    ParticleData particles;
    assert(allocateParticleData(&particles, 1, N, 1) == 0);
    SimulationOptions opts = {.option_for_dust_secondary_population = 1,
                              .dust_smoothing_mode = 1};
    particles.particle_distance_array[0][0] = 2;
    particles.micron_particle_distance_array[0][0] = r[cell+1];
    particles.dust_particle_mass_grid[0] = 1e-6;
    particles.massmicradial_grid[0] = 2e-6;
    particles.micron_dust_particle_mass_array[0][1] = cell;
    particle_number = 1;
    calculateDustSurfaceDensity(&particles, &opts, &disk);
    double ratio = (particles.dust_surfacedensity[cell] +
                    particles.micron_dust_surfacedensity[cell]) / gas[cell+1];
    double dt = .01;
    double tau = calculateGrowthTimescale(r[cell+1], ratio, &disk);
    double expected = 2e-6*dt/tau;
    applyCoagulationMassTransfer(&particles, &disk, &opts, dt);
    assert(fabs((particles.dust_particle_mass_grid[0]-1e-6)-expected) < expected*1e-8);
    assert(fabs(particles.dust_particle_mass_grid[0]+particles.massmicradial_grid[0]-3e-6) < 1e-20);
    freeParticleData(&particles);
}

int main(void)
{
    for (int grid = 0; grid < 2; ++grid) {
        testCoagulationGasIndex(grid, 0);
        testCoagulationGasIndex(grid, 3);
        testCoagulationGasIndex(grid, 7);
    }
    testDriftFragmentationFactor();
    testDriftTimeStepSampling();
    testDustBoundaryMapping();
    for (int grid = 0; grid < 2; ++grid) {
        testGrowth(4, grid);
        testGrowth(8, grid);
        testGrowth(24, grid);
        testFixedStokes(grid);
        testDustBoundaryMass(0, grid);
        testDustBoundaryMass(1, grid);
        testRadiativeBoundaryCrossing(grid, 0);
        testRadiativeBoundaryCrossing(grid, 1);
    }
    puts("Dust: growth, interpolated timestep velocities, bounded RK4 and radiative absorbing exits passed");
    return 0;
}
