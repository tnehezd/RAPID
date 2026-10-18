#include "eos.h"

#include <math.h>
#include <stdlib.h>

#include "disk_thermal.h"
#include "gas_physics.h"
#include "logger.h"

typedef struct {
    double scale_height;
    double sound_speed;
} EOSState;

/* SG initialization needs only its initial H guess, before cs is available. */
static EOSState calculateEOSState(double radial_distance, const DiskParameters *disk_params,
                                 int need_sound_speed) {
    EOSState state = {0};

    switch (disk_params->thermodynamics_model) {
        case THERMODYNAMICS_LOCALLY_ISOTHERMAL:
            state.scale_height = pow(radial_distance, 1. + disk_params->flaring_index) * disk_params->h_aspect_ratio;
            if (need_sound_speed) {
                state.sound_speed = calculateKeplerianFrequency(radial_distance, disk_params) * state.scale_height;
            }
            break;

        case THERMODYNAMICS_RADIATIVE:
        case THERMODYNAMICS_RADIATIVE_ENERGY:
            state.scale_height = interpolateThermalScaleHeight(radial_distance, disk_params);
            if (need_sound_speed) {
                state.sound_speed = interpolateThermalSoundSpeed(radial_distance, disk_params);
            }
            break;

        default:
            LOG_ERROR("Unsupported thermodynamics model %d in calculateEOS.\n",
                      disk_params->thermodynamics_model);
            exit(EXIT_FAILURE);
    }

    return state;
}

double calculatePressureScaleHeight(double radial_distance, const DiskParameters *disk_params) {
    return calculateEOSState(radial_distance, disk_params, 0).scale_height;
}

double calculateLocalSoundSpeed(double radial_distance, const DiskParameters *disk_params) {
    return calculateEOSState(radial_distance, disk_params, 1).sound_speed;
}

double calcualteMidplaneGasDensity(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
    EOSState state = calculateEOSState(radial_distance, disk_params, 0);
    return 1. / sqrt(2.0 * M_PI) * gas_surface_density / state.scale_height;
}

double calculateEOS(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
    EOSState state = calculateEOSState(radial_distance, disk_params, 1);
    double midplane_density = 1. / sqrt(2.0 * M_PI) * gas_surface_density / state.scale_height;
    return midplane_density * state.sound_speed * state.sound_speed;
}
