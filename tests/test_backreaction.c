#include "backreaction.h"
#include "config.h"
#include "dust_physics.h"
#include "gas_physics.h"
#include "simulation_core.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void testDragBalance(void)
{
    const double loadings[] = {0, 1e-8, .01, 1, 100};
    const double stokes_numbers[] = {1e-6, .01, 1, 100};
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 4; ++j) {
            double epsilon = loadings[i], st = stokes_numbers[j];
            double viscous = -1e-5, pressure = -1e-3, ug, vg;
            solveStationaryDrag(epsilon / (1 + st * st),
                                epsilon * st / (1 + st * st), viscous, pressure, &ug, &vg);
            double ud = (ug + 2 * st * vg) / (1 + st * st);
            double denominator = (1 + epsilon) * (1 + epsilon) + st * st;
            double expected = ((1 + epsilon) * viscous + st * pressure) / denominator;
            assert(fabs(ud - expected) < 1e-15);
            assert(fabs(ug + epsilon * ud - viscous) < 1e-15);
            double vd = (vg - .5 * st * ug) / (1 + st * st);
            assert(fabs(2 * vd - (ud - ug) / st) < 1e-12);
            assert(fabs(-.5 * ud - (vd - vg) / st) < 1e-12);
        }
    }
    /* Different grain sizes contribute separate moments, not a mean Stokes number. */
    double eps[] = {.3, .7}, st[] = {.01, 2}, a = 0, b = 0, ug, vg;
    for (int i = 0; i < 2; ++i) {
        a += eps[i] / (1 + st[i] * st[i]);
        b += eps[i] * st[i] / (1 + st[i] * st[i]);
    }
    solveStationaryDrag(a, b, -1e-5, -1e-3, &ug, &vg);
    double mass_velocity = ug;
    double radial_drag = 0, azimuthal_drag = 0;
    for (int i = 0; i < 2; ++i) {
        double ud = (ug + 2 * st[i] * vg) / (1 + st[i] * st[i]);
        double vd = (vg - .5 * st[i] * ug) / (1 + st[i] * st[i]);
        mass_velocity += eps[i] * ud;
        radial_drag += eps[i] * (ud - ug) / st[i];
        azimuthal_drag += eps[i] * (vd - vg) / st[i];
    }
    assert(fabs(mass_velocity + 1e-5) < 1e-15);
    assert(fabs(2 * vg + 1e-3 + radial_drag) < 1e-15);
    assert(fabs(-.5 * ug + .5 * -1e-5 + azimuthal_drag) < 1e-15);
}

static double area(int i, const DiskParameters *disk)
{
    double left = .5 * (disk->radial_grid[i - 1] + disk->radial_grid[i]);
    double right = .5 * (disk->radial_grid[i] + disk->radial_grid[i + 1]);
    return M_PI * (right * right - left * left);
}

static void testTransport(int logarithmic, int count, int reflecting)
{
    enum { N = 8, P = 26 };
    double r[N + 2], sigma[N + 2], pressure[N + 2], grad[N + 2], velocity[N + 2];
    double positions[P][2], masses[P], density[N], sink[N + 2] = {0};
    DiskParameters disk = {
        .grid_number = N, .stellar_mass = 1, .r_min = 1, .r_max = 4, .delta_r = 3.0 / (N - 1),
        .h_aspect_ratio = .05, .alpha_parameter = .003, .alpha_parameter_modification = 1,
        .r_dze_i = 2, .r_dze_o = 3, .dr_dze_i = .1, .dr_dze_o = .1,
        .radial_grid = r, .gas_surface_density_vector = sigma, .gas_pressure_vector = pressure,
        .gas_pressure_gradient_vector = grad, .gas_velocity_vector = velocity,
        .particle_density_dimensionless = 1.6 * pow(AU_IN_CM, 3) / SOLAR_MASS_IN_GRAMS,
        .density_floor = 1e-20, .enable_backreaction = true, .sigma_dot_photoevap = sink
    };
    ParticleData particles = {.particle_distance_array = positions, .dust_particle_mass_grid = masses,
                               .dust_surfacedensity = density};
    SimulationOptions options = {.inner_boundary_condition_type = reflecting ? 4 : 0,
                                 .outer_boundary_condition_type = reflecting ? 4 : 0};
    for (int i = 0; i < N + 2; ++i) {
        r[i] = logarithmic ? pow(4, (i - 1.0) / (N - 1)) : 1 + (i - 1) * disk.delta_r;
        sigma[i] = 1e-5 / sqrt(r[i]);
    }
    for (int i = 0; i < count; ++i) {
        positions[i][0] = 1.1 + 2.8 * i / (count > 1 ? count - 1.0 : 1);
        positions[i][1] = (.1 + .01 * i) / AU_IN_CM;
        masses[i] = 1e-5 / count;
    }
    positions[count][0] = .5;
    positions[count + 1][0] = 5;
    positions[count][1] = positions[count + 1][1] = .1 / AU_IN_CM;
    masses[count] = masses[count + 1] = 1e-3;
    int total_count = count + 2;
    refreshBackreaction(&particles, total_count, &disk, &options);
    double expected_a[N + 2] = {0}, expected_b[N + 2] = {0}, total_mass = 0;
    for (int i = 1; i <= N; ++i) {
        assert(isfinite(disk.gas_velocity_vector[i]));
        total_mass += area(i, &disk) * sigma[i];
    }
    /* Independently locate receiving cells and reconstruct the size moments. */
    for (int k = 0; k < count; ++k) {
        int left = 1;
        while (left < N - 1 && positions[k][0] > r[left + 1]) ++left;
        double fraction = (positions[k][0] - r[left]) / (r[left + 1] - r[left]);
        for (int j = 0; j < 2; ++j) {
            int i = left + j;
            double st = calculateStokesNumber(positions[k][1], sigma[i], &disk);
            double epsilon = masses[k] * (j ? fraction : 1 - fraction) / (area(i, &disk) * sigma[i]);
            expected_a[i] += epsilon / (1 + st * st);
            expected_b[i] += epsilon * st / (1 + st * st);
        }
        double local_sigma = sigma[left];
        double dust_velocity = backreactionDustVelocity(r[left], positions[k][1], local_sigma, &disk);
        double st = calculateStokesNumber(positions[k][1], local_sigma, &disk);
        assert(fabs(dust_velocity - (disk.backreaction->radial[left] +
                    2 * st * disk.backreaction->azimuthal[left]) / (1 + st * st)) < 1e-15);
    }
    for (int i = 1; i <= N; ++i) {
        assert(fabs(expected_a[i] - disk.backreaction->loading_a[i]) < 1e-15);
        assert(fabs(expected_b[i] - disk.backreaction->loading_b[i]) < 1e-15);
    }
    for (int j = 0; j <= N; ++j) {
        double correction_velocity = .5 * (
            disk.backreaction->radial[j] - disk.backreaction->viscous[j] +
            disk.backreaction->radial[j + 1] - disk.backreaction->viscous[j + 1]);
        int donor = correction_velocity >= 0 ? j : j + 1;
        double expected_flux = 2 * M_PI * .5 * (r[j] + r[j + 1]) *
                               correction_velocity * sigma[donor];
        if (reflecting && (j == 0 || j == N)) expected_flux = 0;
        assert(fabs(expected_flux - disk.backreaction->flux[j]) < 1e-20);
    }
    if (count == 0) {
        double expected[N + 2];
        calculateGasRadialVelocity(&disk);
        for (int i = 1; i <= N; ++i) expected[i] = velocity[i];
        refreshBackreaction(&particles, total_count, &disk, &options);
        for (int i = 1; i <= N; ++i) assert(velocity[i] == expected[i]);
        DiskParameters passive = disk;
        passive.enable_backreaction = false;
        double radius = .5 * (r[3] + r[4]), gas_sigma = .5 * (sigma[3] + sigma[4]);
        double gas_grad = .5 * (grad[3] + grad[4]);
        double viscous = .5 * (expected[3] + expected[4]), original;
        calculate1DDustDrift(.1 / AU_IN_CM, gas_grad, gas_sigma, viscous, radius,
                             &original, &passive);
        assert(fabs(backreactionDustVelocity(radius, .1 / AU_IN_CM, gas_sigma, &disk) -
                    original) < 1e-15);
    }
    particle_number = total_count;
    calculateDustSurfaceDensity(&particles, &options, &disk);
    double deposited_mass = 0;
    for (int i = 1; i <= N; ++i) deposited_mass += area(i, &disk) * density[i - 1];
    assert(fabs(deposited_mass - (count ? 1e-5 : 0)) < 1e-18);
    double dt = fmin(1, backreactionTimeStep(&particles, total_count, &disk, &options));
    for (int i = 1; i <= N; ++i) {
        double a = disk.backreaction->loading_a[i], b = disk.backreaction->loading_b[i];
        double denominator = (1 + a) * (1 + a) + b * b;
        double diffusion = 3 * calculateKinematicViscosity(r[i], &disk) *
            (1 + fabs((1 + a) / denominator - 1)) +
            calculateLocalSoundSpeed(r[i], &disk) * calculatePressureScaleHeight(r[i], &disk) * b /
            denominator;
        double dr = fmin(r[i] - r[i - 1], r[i + 1] - r[i]);
        assert(dt <= .2 * dr * dr / (2 * diffusion));
    }
    double boundary = dt * (disk.backreaction->flux[0] - disk.backreaction->flux[N]);
    if (reflecting) assert(boundary == 0);
    double reference_sigma[N + 2], reference_pressure[N + 2], reference_gradient[N + 2];
    DiskParameters reference = disk;
    reference.enable_backreaction = false;
    reference.gas_surface_density_vector = reference_sigma;
    reference.gas_pressure_vector = reference_pressure;
    reference.gas_pressure_gradient_vector = reference_gradient;
    memcpy(reference_sigma, sigma, sizeof(sigma));
    memcpy(reference_pressure, pressure, sizeof(pressure));
    memcpy(reference_gradient, grad, sizeof(grad));
    options.user_defined_time_step = dt;
    refreshGasSurfaceDensityPressurePressureGradient(&options, &reference);
    double reference_mass = 0;
    for (int i = 1; i <= N; ++i) reference_mass += area(i, &disk) * reference_sigma[i];
    advanceBackreactionGas(dt, &disk, &options);
    double final_mass = 0;
    for (int i = 1; i <= N; ++i) {
        final_mass += area(i, &disk) * sigma[i];
        assert(sigma[i] >= disk.density_floor);
        double expected = reference_sigma[i] -
            dt * (disk.backreaction->flux[i] - disk.backreaction->flux[i - 1]) / area(i, &disk);
        assert(sigma[i] == expected);
    }
    assert(fabs(final_mass - reference_mass - boundary) < 1e-15 * total_mass);
    if (count == 0) {
        for (int step = 0; step < 10; ++step) {
            refreshBackreaction(&particles, total_count, &disk, &options);
            for (int i = 0; i <= N; ++i) assert(disk.backreaction->flux[i] == 0);
            refreshGasSurfaceDensityPressurePressureGradient(&options, &reference);
            advanceBackreactionGas(dt, &disk, &options);
            for (int i = 0; i < N + 2; ++i) assert(sigma[i] == reference_sigma[i]);
        }
    }
    free(disk.backreaction);
}

static double viscousContinuumError(int n, int logarithmic)
{
    double r[n + 2], sigma[n + 2], pressure[n + 2], gradient[n + 2], velocity[n + 2];
    double original[n + 2], sink[n + 2];
    memset(sink, 0, sizeof(sink));
    DiskParameters disk = {
        .grid_number = n, .stellar_mass = 1, .r_min = 1, .r_max = 4, .delta_r = 3.0 / (n - 1),
        .h_aspect_ratio = .05, .alpha_parameter = .003, .alpha_parameter_modification = 1,
        .r_dze_i = 2, .r_dze_o = 3, .dr_dze_i = .1, .dr_dze_o = .1,
        .radial_grid = r, .gas_surface_density_vector = sigma, .gas_pressure_vector = pressure,
        .gas_pressure_gradient_vector = gradient, .gas_velocity_vector = velocity,
        .density_floor = 1e-20, .enable_backreaction = true, .sigma_dot_photoevap = sink
    };
    ParticleData particles = {0};
    SimulationOptions options = {0};
    for (int i = 0; i < n + 2; ++i) {
        r[i] = logarithmic ? pow(4, (i - 1.0) / (n - 1)) : 1 + (i - 1) * disk.delta_r;
        sigma[i] = 1e-5 * pow(r[i], -.25);
    }
    refreshBackreaction(&particles, 0, &disk, &options);
    memcpy(original, sigma, sizeof(sigma));
    const double dt = .01;
    advanceBackreactionGas(dt, &disk, &options);
    double error = 0, reference = 0;
    for (int i = 2; i < n; ++i) {
        if (r[i] < 1.5 || r[i] > 3.5) continue;
        double exact = 3 * calculateKinematicViscosity(r[i], &disk) * original[i] *
                       .25 * .75 / (r[i] * r[i]);
        double numerical = (sigma[i] - original[i]) / dt;
        error += fabs(numerical - exact) * area(i, &disk);
        reference += exact * area(i, &disk);
    }
    free(disk.backreaction);
    return error / reference;
}

int main(void)
{
    testDragBalance();
    for (int grid = 0; grid < 2; ++grid) {
        double coarse = viscousContinuumError(32, grid);
        double medium = viscousContinuumError(64, grid);
        double fine = viscousContinuumError(128, grid);
        assert(medium < .7 * coarse);
        assert(fine < .7 * medium);
        assert(fine < .02);
        for (int boundary = 0; boundary < 2; ++boundary) {
            testTransport(grid, 0, boundary);
            testTransport(grid, 4, boundary);
            testTransport(grid, 8, boundary);
            testTransport(grid, 24, boundary);
        }
    }
    puts("Backreaction: drag balance, CIC, exact FTCS limit, conservative drag correction and viscous convergence passed");
    return 0;
}
