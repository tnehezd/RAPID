#ifndef IRRADIATION_H
#define IRRADIATION_H

typedef enum {
    IRRADIATION_PRESCRIBED = 0,
    IRRADIATION_GEOMETRY
} IrradiationModel;

/* CGS lengths; heights are gas H, not the irradiated surface chi*H.
 * Only physical cells 1..index-1 are used as occulting screens. */
typedef struct {
    const double *radii;
    const double *heights;
    int index;
    double surface_height_factor;
    double stellar_radius;
} IrradiationGeometry;

/* Uniform stellar disk, small stellar angular size, opaque disk surface.
 * Returns zero on invalid geometry; otherwise sets both outputs. */
int calculateIrradiationGeometry(const IrradiationGeometry *geometry, double height,
                                double *factor, double *visibility);

#endif
