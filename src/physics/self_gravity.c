#include "self_gravity.h"

#include "config.h"
#include "eos.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Softening length of the ring kernel in units of the local pressure scale height (Mueller et al. 2012).
#define SELFGRAVITY_SOFTENING_SCALEHEIGHT_FACTOR 0.6

double interpolateSelfGravityOmega2(double radial_distance, const DiskParameters *disk_params) {

    const double *grid = disk_params->radial_grid;
    const double *extra = disk_params->selfgravity_omega2_extra;
    int low = 0;
    int high = disk_params->grid_number + 1;

    if (radial_distance <= grid[low]) return extra[low];
    if (radial_distance >= grid[high]) return extra[high];

    while (high - low > 1) {
        int middle = (low + high) / 2;
        if (grid[middle] <= radial_distance) low = middle; else high = middle;
    }
    double fraction = (radial_distance - grid[low]) / (grid[high] - grid[low]);
    return extra[low] + fraction * (extra[high] - extra[low]);
}

// Complete elliptic integrals K(k), E(k) by the arithmetic-geometric mean.
static void completeEllipticIntegrals(double k_squared, double *K, double *E) {

    double a = 1.0;
    double b = sqrt(1.0 - k_squared);
    double c_squared_sum = 0.5 * k_squared;
    double power_of_two = 0.5;

    for (int n = 0; n < 40 && fabs(a - b) > 1e-15 * a; n++) {
        double a_next = 0.5 * (a + b);
        double c = 0.5 * (a - b);
        b = sqrt(a * b);
        a = a_next;
        power_of_two *= 2.0;
        c_squared_sum += power_of_two * c * c;
    }
    *K = M_PI / (2.0 * a);
    *E = *K * (1.0 - c_squared_sum);
}

// Radial gravitational acceleration at radius R of a softened ring of mass m and radius a (negative = inward).
static double ringRadialAcceleration(double R, double a, double m, double softening) {

    double s_squared = (R + a) * (R + a) + softening * softening;
    double s = sqrt(s_squared);
    double k_squared = 4.0 * R * a / s_squared;
    double K, E;
    completeEllipticIntegrals(k_squared, &K, &E);

    double k = sqrt(k_squared);
    double dK_dk = E / (k * (1.0 - k_squared)) - K / k;
    double dk2_dR = 4.0 * a / s_squared - 8.0 * R * a * (R + a) / (s_squared * s_squared);
    double dk_dR = dk2_dR / (2.0 * k);
    double d_KoverS_dR = dK_dk * dk_dR / s - K * (R + a) / (s_squared * s);

    return (2.0 * G_DIMENSIONLESS * m / M_PI) * d_KoverS_dR;
}

void updateSelfGravity(DiskParameters *disk_params) {

    if (!disk_params->enable_selfgravity) return;

    int n = disk_params->grid_number;
    if (disk_params->selfgravity_omega2_extra == NULL) {
        disk_params->selfgravity_omega2_extra = (double *)calloc(n + 2, sizeof(double));
        if (disk_params->selfgravity_omega2_extra == NULL) {
            fprintf(stderr, "Failed to allocate memory for the self-gravity array.\n");
            exit(1);
        }
    }

    const double *r = disk_params->radial_grid;
    const double *sigma = disk_params->gas_surface_density_vector;
    double *extra = disk_params->selfgravity_omega2_extra;

    double ring_mass[n + 2];
    for (int j = 1; j <= n; j++) {
        ring_mass[j] = 2.0 * M_PI * r[j] * sigma[j] * 0.5 * (r[j + 1] - r[j - 1]);
    }

    #pragma omp parallel for
    for (int i = 0; i <= n + 1; i++) {
        double softening = SELFGRAVITY_SOFTENING_SCALEHEIGHT_FACTOR * calculatePressureScaleHeight(r[i], disk_params);
        double radial_acceleration = 0.0;
        for (int j = 1; j <= n; j++) {
            radial_acceleration += ringRadialAcceleration(r[i], r[j], ring_mass[j], softening);
        }
        extra[i] = -radial_acceleration / r[i];
    }
}
