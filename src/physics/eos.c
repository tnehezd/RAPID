#include "eos.h"

#include <math.h>
#include <stdlib.h>

#include "disk_thermal.h"
#include "gas_physics.h"
#include "logger.h"

double calculatePressureScaleHeight(double radial_distance, const DiskParameters *disk_params) {
    if (disk_params->thermodynamics_model == THERMODYNAMICS_RADIATIVE) {
        return interpolateThermalScaleHeight(radial_distance, disk_params);
    }
    double local_pressure_scaleheight = pow(radial_distance, 1. + disk_params->flaring_index) * disk_params->h_aspect_ratio;
    return local_pressure_scaleheight;
}

double calculateLocalSoundSpeed(double radial_distance, const DiskParameters *disk_params) {
    if (disk_params->thermodynamics_model == THERMODYNAMICS_RADIATIVE) {
        return interpolateThermalSoundSpeed(radial_distance, disk_params);
    }
    return calculateKeplerianFrequency(radial_distance,disk_params) * calculatePressureScaleHeight(radial_distance,disk_params);
}

double calcualteMidplaneGasDensity(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
    return 1. / sqrt(2.0 * M_PI) * gas_surface_density / calculatePressureScaleHeight(radial_distance,disk_params);
}

double calculateEOS(double gas_surface_density, double radial_distance, const DiskParameters *disk_params) {
    double scale_height;
    double sound_speed;

    switch (disk_params->thermodynamics_model) {
        case THERMODYNAMICS_LOCALLY_ISOTHERMAL:
            scale_height = pow(radial_distance, 1. + disk_params->flaring_index) * disk_params->h_aspect_ratio;
            sound_speed = calculateKeplerianFrequency(radial_distance, disk_params) * scale_height;
            break;

        case THERMODYNAMICS_RADIATIVE:
            scale_height = interpolateThermalScaleHeight(radial_distance, disk_params);
            sound_speed = interpolateThermalSoundSpeed(radial_distance, disk_params);
            break;

        default:
            LOG_ERROR("Unsupported thermodynamics model %d in calculateEOS.\n",
                      disk_params->thermodynamics_model);
            exit(EXIT_FAILURE);
    }

    double midplane_density = 1. / sqrt(2.0 * M_PI) * gas_surface_density / scale_height;
    return midplane_density * sound_speed * sound_speed;
}
