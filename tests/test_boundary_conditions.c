#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "boundary_conditions.h"
#include "thermal_model.h"

static void checkParabolicDensity(int logarithmic, int inner)
{
    double radii[5];
    double density[5] = {0};
    for (int i = 0; i < 5; ++i) {
        radii[i] = logarithmic ? pow(1.1, i) : 1.0 + i;
    }
    double edge = inner ? radii[0] : radii[4];
    double gap = inner ? radii[1] - edge : edge - radii[3];
    for (int i = 1; i <= 3; ++i) {
        density[i] = 1e-8 * ((inner ? radii[i] - edge : edge - radii[i]) - 0.1 * gap);
        assert(density[i] > 0.0);
    }
    DiskParameters disk = {0};
    disk.grid_number = 3;
    disk.radial_grid = radii;
    disk.density_floor = 1e-20;
    SimulationOptions options = {0};
    options.inner_boundary_condition_type = inner ? 1 : 0;
    options.outer_boundary_condition_type = inner ? 0 : 1;
    int ghost = inner ? 0 : 4;
    double raw[5];
    memcpy(raw, density, sizeof(raw));
    applyBoundaryConditions(raw, &disk, &options);
    assert(raw[ghost] < 0.0);

    ThermalParameters parameters = {
        .surface_density = raw[ghost], .orbital_frequency = 1e-8,
        .alpha = 0.003, .mean_molecular_weight = 2.34,
        .background_temperature = 10.0, .shear_parameter = 1.5
    };
    ThermalSolverOptions solver = defaultThermalSolverOptions();
    ThermalState state;
    assert(solveMidplaneTemperature(&parameters, &solver, &state) == THERMAL_INVALID_INPUT);

    options.current_bc_target = 2;
    applySurfaceDensityBoundaryConditions(density, &disk, &options);
    assert(options.current_bc_target == 2);
    assert(density[ghost] == density[inner ? 1 : 3]);
    for (int i = 1; i <= 3; ++i) assert(density[i] == raw[i]);
    parameters.surface_density = density[ghost] * 1.989e33 / (1.496e13 * 1.496e13);
    ThermalStatus status = solveMidplaneTemperature(&parameters, &solver, &state);
    if (status != THERMAL_SUCCESS) fprintf(stderr, "Ghost equilibrium: %s\n", thermalStatusMessage(status));
    assert(status == THERMAL_SUCCESS);
    assert(isfinite(state.temperature) && state.temperature >= 10.0);

    for (int field = 0; field <= 2; ++field) {
        options.current_bc_target = field;
        applyBoundaryConditions(raw, &disk, &options);
        assert(raw[ghost] < 0.0);
    }

    for (int i = 1; i <= 3; ++i) density[i] = 1.0;
    applySurfaceDensityBoundaryConditions(density, &disk, &options);
    assert(fabs(density[ghost] - 1.0) < 1e-12);

    for (int type = 0; type <= 5; ++type) {
        if (type == 1) continue;
        options.inner_boundary_condition_type = type;
        options.outer_boundary_condition_type = type;
        options.current_bc_target = 0;
        memcpy(raw, density, sizeof(raw));
        applyBoundaryConditions(raw, &disk, &options);
        applySurfaceDensityBoundaryConditions(density, &disk, &options);
        assert(memcmp(raw, density, sizeof(raw)) == 0);
        if (type == 3) assert(density[0] == 0.0 && density[4] == 0.0);
        if (type == 5) assert(density[0] == density[1] && density[4] == density[3]);
    }
}

int main(void)
{
    for (int logarithmic = 0; logarithmic <= 1; ++logarithmic) {
        for (int inner = 0; inner <= 1; ++inner) {
            checkParabolicDensity(logarithmic, inner);
        }
    }
    puts("Density boundaries: all checks passed");
    return 0;
}
