#include "disk_thermal.h"
#include "config.h"
#include "gas_physics.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 8

typedef struct {
    DiskParameters disk;
    SimulationOptions options;
    double radius[N + 2], sigma[N + 2], pressure[N + 2], gradient[N + 2];
    double velocity[N + 2], sink[N + 2];
    ThermalState thermal[N + 2];
} Fixture;

static void initialize(Fixture *f, double alpha, double gamma, int logarithmic)
{
    *f = (Fixture){0};
    f->disk = (DiskParameters){
        .grid_number = N, .stellar_mass = 1.0,
        .thermodynamics_model = THERMODYNAMICS_RADIATIVE_ENERGY,
        .thermal_mean_molecular_weight = 2.34, .thermal_gamma = gamma,
        .thermal_background_temperature = 10.0,
        .thermal_solver_options = {1.0, 1500.0, 1e-8, 128},
        .alpha_parameter = alpha, .alpha_parameter_modification = 1.0,
        .density_floor = 1e-20,
        .radial_grid = f->radius, .gas_surface_density_vector = f->sigma,
        .gas_pressure_vector = f->pressure, .gas_pressure_gradient_vector = f->gradient,
        .gas_velocity_vector = f->velocity, .sigma_dot_photoevap = f->sink,
        .thermal_states = f->thermal
    };
    for (int i = 0; i <= N + 1; ++i) {
        f->radius[i] = logarithmic ? pow(1.2, i) : 1.0 + 0.25 * i;
        f->sigma[i] = 1e-5 * pow(f->radius[i], -0.7);
        f->thermal[i].temperature = 40.0;
        f->thermal[i].scale_height = 0.03 * f->radius[i] * AU_IN_CM;
    }
    refreshRadiativeDisk(&f->disk, &f->options);
}

static double heatCapacity(const Fixture *f)
{
    double speed_unit = AU_IN_CM / thermalTimeUnitSeconds();
    return BOLTZMANN_CONSTANT_CGS /
        (2.34 * HYDROGEN_MASS_IN_GRAMS * (f->disk.thermal_gamma - 1.0) *
         speed_unit * speed_unit);
}

static void closeEnough(double actual, double expected)
{
    assert(fabs(actual - expected) <= 2e-12 * fmax(fabs(expected), 1e-30));
}

static void testTransport(int logarithmic, int boundary)
{
    Fixture f;
    initialize(&f, 0.003, 1.4, logarithmic);
    for (int i = 1; i <= N; ++i) f.thermal[i].temperature = 40.0 + 15.0 * sin(i);
    f.thermal[0].temperature = f.thermal[1].temperature;
    f.thermal[N + 1].temperature = f.thermal[N].temperature;
    refreshRadiativeDisk(&f.disk, &f.options);
    f.options.inner_boundary_condition_type = boundary;
    f.options.outer_boundary_condition_type = boundary;
    double cv = heatCapacity(&f);
    double dt = 1e-5;
    double flux[N + 1], eflux[N + 1], vflux[N + 1], face[N + 1];
    double expected_sigma[N + 2], expected_energy[N + 2];
    double conversion = pow(thermalTimeUnitSeconds(), 3) / SOLAR_MASS_IN_GRAMS;
    for (int j = 0; j <= N; ++j) {
        double l = f.radius[j], r = f.radius[j + 1];
        face[j] = 0.5 * (l + r);
        /* Independent reconstruction of the viscous face flux. */
        double a = f.thermal[j].viscosity * pow(thermalTimeUnitSeconds() / AU_IN_CM, 2) /
                   thermalTimeUnitSeconds() * f.sigma[j] * sqrt(l);
        double b = f.thermal[j + 1].viscosity * pow(thermalTimeUnitSeconds() / AU_IN_CM, 2) /
                   thermalTimeUnitSeconds() * f.sigma[j + 1] * sqrt(r);
        flux[j] = -3 * sqrt(face[j]) * (b - a) / (r - l);
        if (boundary == 4 && (j == 0 || j == N)) flux[j] = 0;
        if (boundary == 5 && j == 0) flux[j] = fmin(flux[j], 0);
        if (boundary == 5 && j == N) flux[j] = fmax(flux[j], 0);
        int donor = flux[j] >= 0 ? j : j + 1;
        eflux[j] = flux[j] * cv * f.thermal[donor].temperature;
        vflux[j] = flux[j] / (0.5 * (f.sigma[j] + f.sigma[j + 1]));
    }
    double old_mass = 0, new_mass = 0, old_energy = 0, new_energy = 0, sources = 0;
    for (int i = 1; i <= N; ++i) {
        double area = 0.5 * (face[i] * face[i] - face[i - 1] * face[i - 1]);
        double e = f.sigma[i] * cv * f.thermal[i].temperature;
        double source = (f.thermal[i].heating - f.thermal[i].cooling) * conversion -
            (f.disk.thermal_gamma - 1) * e * (vflux[i] - vflux[i - 1]) / area;
        expected_sigma[i] = f.sigma[i] - dt * (flux[i] - flux[i - 1]) / area;
        expected_energy[i] = e - dt * (eflux[i] - eflux[i - 1]) / area + dt * source;
        old_mass += area * f.sigma[i];
        old_energy += area * e;
        sources += area * source;
    }
    f.options.user_defined_time_step = dt;
    evolveRadiativeEnergy(&f.disk, &f.options);
    for (int i = 1; i <= N; ++i) {
        double area = 0.5 * (face[i] * face[i] - face[i - 1] * face[i - 1]);
        closeEnough(f.sigma[i], expected_sigma[i]);
        double e = f.sigma[i] * cv * f.thermal[i].temperature;
        closeEnough(e, expected_energy[i]);
        closeEnough(f.pressure[i], calculateEOS(f.sigma[i], f.radius[i], &f.disk));
        new_mass += area * f.sigma[i];
        new_energy += area * e;
    }
    closeEnough(new_mass, old_mass - dt * (flux[N] - flux[0]));
    closeEnough(new_energy, old_energy - dt * (eflux[N] - eflux[0]) + dt * sources);
    double budget_change = f.disk.energy_budget.heating - f.disk.energy_budget.cooling +
        f.disk.energy_budget.work + f.disk.energy_budget.boundary -
        f.disk.energy_budget.wind + f.disk.energy_budget.floor;
    closeEnough(totalDiskInternalEnergy(&f.disk), f.disk.energy_budget.initial + budget_change);
    closeEnough(f.disk.energy_budget.boundary, -2.0 * M_PI * dt * (eflux[N] - eflux[0]));
    if (boundary == 5) assert(new_mass <= old_mass);
}

static void testOutflowDirections(void)
{
    Fixture f;
    initialize(&f, 0.003, 1.4, 0);
    for (int side = 0; side < 2; ++side) {
        int face = side ? N : 0;
        int ghost = side ? N + 1 : 0;
        for (int inflow = 0; inflow < 2; ++inflow) {
            /* Large ghost density drives inflow; zero ghost density drives outflow. */
            f.sigma[ghost] = inflow ? 1.0 : 0.0;
            f.options.inner_boundary_condition_type = 0;
            f.options.outer_boundary_condition_type = 0;
            double raw = radiativeMassFlux(face, &f.disk, &f.options);
            assert((side ? -raw : raw) * (inflow ? 1.0 : -1.0) > 0);
            f.options.inner_boundary_condition_type = 5;
            f.options.outer_boundary_condition_type = 5;
            double limited = radiativeMassFlux(face, &f.disk, &f.options);
            if (inflow) assert(limited == 0);
            else closeEnough(limited, raw);
        }
    }
}

static void testCooling(double gamma)
{
    Fixture f;
    initialize(&f, 0.0, gamma, 0);
    double cv_cgs = BOLTZMANN_CONSTANT_CGS / (2.34 * HYDROGEN_MASS_IN_GRAMS * (gamma - 1));
    double sigma_unit = SOLAR_MASS_IN_GRAMS / (AU_IN_CM * AU_IN_CM);
    double rate = f.thermal[4].cooling / (f.sigma[4] * sigma_unit * cv_cgs) *
                  thermalTimeUnitSeconds();
    double dt = 0.01 * 40.0 / rate;
    double expected[N + 2];
    for (int i = 1; i <= N; ++i) {
        expected[i] = 40 - dt * thermalTimeUnitSeconds() * f.thermal[i].cooling /
                      (f.sigma[i] * sigma_unit * cv_cgs);
    }
    f.options.user_defined_time_step = dt;
    evolveRadiativeEnergy(&f.disk, &f.options);
    for (int i = 1; i <= N; ++i) {
        closeEnough(f.thermal[i].temperature, expected[i]);
        assert(f.thermal[i].temperature > 10 && f.thermal[i].temperature < 40);
    }
    /* Long source integration must relax, not jump instantaneously or oscillate. */
    f.options.user_defined_time_step = 1000 * dt;
    evolveRadiativeEnergy(&f.disk, &f.options);
    for (int i = 1; i <= N; ++i) {
        assert(f.thermal[i].temperature >= 10 && f.thermal[i].temperature < expected[i]);
    }
}

static double coolingAtTime(double duration, int steps)
{
    Fixture f;
    initialize(&f, 0.0, 1.4, 0);
    f.options.user_defined_time_step = duration / steps;
    for (int i = 0; i < steps; ++i) evolveRadiativeEnergy(&f.disk, &f.options);
    return f.thermal[4].temperature;
}

static void testConvergenceAndGamma(void)
{
    Fixture a, b;
    initialize(&a, 0.0, 1.4, 0);
    initialize(&b, 0.0, 5.0 / 3.0, 0);
    double rate = a.thermal[4].cooling * pow(thermalTimeUnitSeconds(), 3) /
                  (SOLAR_MASS_IN_GRAMS * a.sigma[4] * heatCapacity(&a));
    double duration = 0.1 * 40 / rate;
    double coarse = coolingAtTime(duration, 10);
    double fine = coolingAtTime(duration, 20);
    double reference = coolingAtTime(duration, 200);
    assert(fabs(fine - reference) < 0.6 * fabs(coarse - reference));
    a.options.user_defined_time_step = b.options.user_defined_time_step = duration / 100;
    evolveRadiativeEnergy(&a.disk, &a.options);
    evolveRadiativeEnergy(&b.disk, &b.options);
    closeEnough((40 - b.thermal[4].temperature) / (40 - a.thermal[4].temperature),
                (b.disk.thermal_gamma - 1) / (a.disk.thermal_gamma - 1));
}

static void testHeating(void)
{
    Fixture f;
    initialize(&f, 0.0, 1.4, 0);
    f.disk.irradiation_luminosity_lsun = 1.0;
    f.disk.irradiation_factor = 0.05;
    refreshRadiativeDisk(&f.disk, &f.options);
    assert(f.thermal[4].cooling < 0);
    double rate = -f.thermal[4].cooling * pow(thermalTimeUnitSeconds(), 3) /
                  (SOLAR_MASS_IN_GRAMS * f.sigma[4] * heatCapacity(&f));
    f.options.user_defined_time_step = 0.001 * 40 / rate;
    double expected = 40 + f.options.user_defined_time_step * rate;
    evolveRadiativeEnergy(&f.disk, &f.options);
    closeEnough(f.thermal[4].temperature, expected);
}

static void geometricAtTime(double duration, int steps, double *temperatures)
{
    Fixture f;
    initialize(&f, 0.0, 1.4, 1);
    f.disk.irradiation_model = IRRADIATION_GEOMETRY;
    f.disk.irradiation_surface_height = 4.0;
    f.disk.irradiation_stellar_radius_rsun = 2.0;
    f.disk.irradiation_luminosity_lsun = 1.0;
    refreshRadiativeDisk(&f.disk, &f.options);
    f.options.user_defined_time_step = duration / steps;
    for (int step = 0; step < steps; ++step) evolveRadiativeEnergy(&f.disk, &f.options);
    for (int i = 1; i <= N; ++i) temperatures[i] = f.thermal[i].temperature;
    double expected = f.disk.energy_budget.initial + f.disk.energy_budget.heating -
        f.disk.energy_budget.cooling + f.disk.energy_budget.work + f.disk.energy_budget.boundary -
        f.disk.energy_budget.wind + f.disk.energy_budget.floor;
    closeEnough(totalDiskInternalEnergy(&f.disk), expected);
}

static void testGeometricConvergence(void)
{
    double coarse[N + 2], fine[N + 2], reference[N + 2];
    geometricAtTime(0.1, 10, coarse);
    geometricAtTime(0.1, 20, fine);
    geometricAtTime(0.1, 200, reference);
    double coarse_error = 0, fine_error = 0;
    for (int i = 1; i <= N; ++i) {
        coarse_error += fabs(coarse[i] - reference[i]);
        fine_error += fabs(fine[i] - reference[i]);
    }
    assert(coarse_error > 1e-10);
    assert(fine_error < 0.6 * coarse_error);
}

static void testRadiationFloorStability(void)
{
    Fixture f;
    initialize(&f, 0.0, 1.4, 0);
    for (int i = 0; i <= N + 1; ++i) f.thermal[i].temperature = 10.01;
    refreshRadiativeDisk(&f.disk, &f.options);
    double emission = 2 * STEFAN_BOLTZMANN_CONSTANT_CGS * pow(10.01, 4) /
                      f.thermal[4].effective_optical_depth;
    double conversion = pow(thermalTimeUnitSeconds(), 3) / SOLAR_MASS_IN_GRAMS;
    f.options.user_defined_time_step = f.sigma[4] * heatCapacity(&f) * 10.01 /
                                      (emission * conversion);
    evolveRadiativeEnergy(&f.disk, &f.options);
    for (int i = 1; i <= N; ++i) {
        assert(f.thermal[i].temperature >= 10.0);
        assert(f.thermal[i].temperature < 10.001);
    }
}

int main(void)
{
    for (int grid = 0; grid <= 1; ++grid)
        for (int boundary = 0; boundary <= 5; ++boundary)
            if (boundary == 0 || boundary == 4 || boundary == 5) testTransport(grid, boundary);
    testOutflowDirections();
    testCooling(1.4);
    testCooling(5.0 / 3.0);
    testConvergenceAndGamma();
    testHeating();
    testRadiationFloorStability();
    testGeometricConvergence();
    puts("Energy equation: transport, compression, budgets and finite cooling checks passed");
    return 0;
}
