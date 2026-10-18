#include "config.h"
#include "dust_physics.h"
#include "integrator.h"

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
}

int main(void)
{
    for (int grid = 0; grid < 2; ++grid) {
        testGrowth(4, grid);
        testGrowth(8, grid);
        testGrowth(24, grid);
    }
    puts("Dust growth: gas-grid interpolation and single size update checks passed");
    return 0;
}
