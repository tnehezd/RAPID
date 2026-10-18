#include "irradiation.h"

#include <math.h>
#include <stddef.h>

static double visibleFraction(double cutoff)
{
    if (cutoff <= -1.0) return 1.0;
    if (cutoff >= 1.0) return 0.0;
    return (acos(cutoff) - cutoff * sqrt(1.0 - cutoff * cutoff)) / acos(-1.0);
}

int calculateIrradiationGeometry(const IrradiationGeometry *geometry, double height,
                                double *factor, double *visibility)
{
    if (geometry == NULL || geometry->radii == NULL || geometry->heights == NULL ||
        geometry->index < 1 || factor == NULL || visibility == NULL ||
        !isfinite(height) || height <= 0.0 ||
        !isfinite(geometry->surface_height_factor) || geometry->surface_height_factor <= 0.0 ||
        !isfinite(geometry->stellar_radius) || geometry->stellar_radius <= 0.0) return 0;
    int i = geometry->index;
    double radius = geometry->radii[i];
    if (!isfinite(radius) || radius <= geometry->stellar_radius) return 0;
    double surface = geometry->surface_height_factor * height;
    double slope = surface / radius;
    double cutoff = -1.0;
    for (int j = 1; j < i; ++j) {
        double inner_radius = geometry->radii[j];
        double inner_surface = geometry->surface_height_factor * geometry->heights[j];
        if (!isfinite(inner_radius) || inner_radius <= 0.0 || inner_radius >= radius ||
            !isfinite(inner_surface) || inner_surface <= 0.0) return 0;
        /* A stellar point at height y clears this annulus only if y exceeds
         * the intercept of the ray through the two surface points. */
        double intercept = (inner_surface - surface * inner_radius / radius) /
                           (1.0 - inner_radius / radius);
        cutoff = fmax(cutoff, intercept / geometry->stellar_radius);
        if (j == i - 1) slope = (surface - inner_surface) / (radius - inner_radius);
    }
    double distance = hypot(radius, surface);
    double normal = hypot(1.0, slope);
    double grazing = (radius * slope - surface) / (distance * normal);
    double angular_size = geometry->stellar_radius / (distance * normal);
    double lower = fmax(cutoff, -grazing / angular_size);
    double incident = 0.0;
    if (lower < 1.0) {
        lower = fmax(lower, -1.0);
        double segment = sqrt(fmax(0.0, 1.0 - lower * lower));
        /* Integrate positive incidence over the unocculted uniform stellar
         * disk. Half of the intercepted radiation is directed toward midplane. */
        incident = grazing * visibleFraction(lower) +
                   angular_size * 2.0 / (3.0 * acos(-1.0)) * segment * segment * segment;
    }
    if (!isfinite(surface) || !isfinite(incident) || !isfinite(cutoff)) return 0;
    *factor = 0.5 * fmin(1.0, fmax(0.0, incident));
    *visibility = visibleFraction(cutoff);
    return 1;
}
