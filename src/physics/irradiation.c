#include "irradiation.h"

#include <math.h>
#include <stddef.h>

static double visibleFraction(double cutoff)
{
    if (cutoff <= -1.0) return 1.0;
    if (cutoff >= 1.0) return 0.0;
    return (acos(cutoff) - cutoff * sqrt(1.0 - cutoff * cutoff)) / acos(-1.0);
}

int appendIrradiationShadowHull(const IrradiationGeometry *geometry, int *hull, int *count)
{
    if (geometry == NULL || geometry->radii == NULL || geometry->heights == NULL ||
        hull == NULL || count == NULL || geometry->index < 1 ||
        *count < 0 || *count >= geometry->index) return 0;
    int i = geometry->index;
    double radius = geometry->radii[i], height = geometry->heights[i];
    if (!isfinite(radius) || radius <= 0.0 || !isfinite(height) || height <= 0.0 ||
        (*count > 0 && (hull[*count - 1] != i - 1 ||
                       geometry->radii[hull[*count - 1]] >= radius))) return 0;
    while (*count >= 2) {
        int a = hull[*count - 2], b = hull[*count - 1];
        long double ab = (long double)geometry->radii[b] - geometry->radii[a];
        long double bc = (long double)radius - geometry->radii[b];
        long double rise_ab = (long double)geometry->heights[b] - geometry->heights[a];
        long double rise_bc = (long double)height - geometry->heights[b];
        /* Points below an upper chord cannot set the stellar horizon for
         * any exterior query, independent of the trial surface height. */
        if (ab * rise_bc < rise_ab * bc) break;
        --*count;
    }
    hull[(*count)++] = i;
    return 1;
}

static double interceptAt(const IrradiationGeometry *geometry, int j,
                          double radius, double surface)
{
    double inner_radius = geometry->radii[j];
    double inner_surface = geometry->surface_height_factor * geometry->heights[j];
    return (inner_surface - surface * inner_radius / radius) /
           (1.0 - inner_radius / radius);
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
    if (geometry->shadow_hull != NULL) {
        int count = geometry->shadow_hull_count;
        if (count < 0 || count >= i || (i > 1 &&
            (count == 0 || geometry->shadow_hull[count - 1] != i - 1))) return 0;
        if (count > 0) {
            int low = 0, high = count - 1;
            /* Secants to a concave upper hull have a single maximum. */
            while (low < high) {
                int middle = low + (high - low) / 2;
                int a = geometry->shadow_hull[middle], b = geometry->shadow_hull[middle + 1];
                long double za = (long double)geometry->surface_height_factor * geometry->heights[a];
                long double zb = (long double)geometry->surface_height_factor * geometry->heights[b];
                long double left = (za - surface) * ((long double)radius - geometry->radii[b]);
                long double right = (zb - surface) * ((long double)radius - geometry->radii[a]);
                if (left <= right) low = middle + 1;
                else high = middle;
            }
            cutoff = fmax(cutoff, interceptAt(geometry, geometry->shadow_hull[low],
                                             radius, surface) / geometry->stellar_radius);
            int previous = i - 1;
            slope = (surface - geometry->surface_height_factor * geometry->heights[previous]) /
                    (radius - geometry->radii[previous]);
        }
    } else for (int j = 1; j < i; ++j) {
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
