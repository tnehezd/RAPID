#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "simulation_types.h"
#include "ring_test.h"
#include "logger.h"

void setupViscousRingTestCustom(
    SimulationOptions *sim,
    DiskParameters *disk,
    double ring_center,   // r_peak
    double ring_width,    // Gaussian width
    double ring_mass      // M0 (total mass)
) {
    int N = disk->grid_number;
    double *r = disk->radial_grid;
    double *Sigma = disk->gas_surface_density_vector;

    // 1) LB&P-compatible peak value:
    //    sigma_peak = M0 / (2π r_peak)
    double sigma_peak = ring_mass / (2.0 * M_PI * ring_center);

    // 2) Initialize density with LB&P Gaussian
    for (int i = 1; i <= N; i++) {
        double dist = r[i] - ring_center;
        Sigma[i] = sigma_peak * exp(-(dist * dist) / (2.0 * ring_width * ring_width));
    }

    // 3) Disable dust physics for ring test
    sim->option_for_evolution = 1.0;
    sim->option_for_dust_drift = 0.0;
    sim->option_for_dust_growth = 0.0;
    sim->option_for_dust_secondary_population = 0.0;

    LOG_INFO("LB&P-compatible ring test initialized: center=%.3f AU, width=%.3f AU, mass=%.3e",
             ring_center, ring_width, ring_mass);
}