#include <stdio.h>
#include "print_panels.h"


void printWarningPanelForPhotoevaporation(const DiskParameters *disk_params)
{
    fprintf(stderr, "\n******************************************************************************************\n");
    fprintf(stderr, "[WARNING:] Photoevaporation is globally set to FALSE, but a specific mode ('%s')\n", disk_params->photoevaporation_mode_string);
    fprintf(stderr, "           and options were provided in the configuration!\n");
    fprintf(stderr, "           The simulation will IGNORE these parameters and run WITHOUT photoevaporation.\n");
    fprintf(stderr, "******************************************************************************************\n\n");
}

void printRunConfigurationHeader(const ParserOptions *def, const DiskParameters *disk_params) {
    fprintf(stderr, "\n==========================================================\n");
    fprintf(stderr, " INITIAL GAS INITIAL CONDITION PROFILE:\n");
    if (def->use_cutoff) {
        fprintf(stderr, "   STYLE:          [Exponential Cutoff / Self-Similar Taper]\n");
        fprintf(stderr, "   R_CUTOFF:       [%.2f AU]\n", def->r_cutoff);
        fprintf(stderr, "   SHARPNESS N:    [%.2f]\n", def->n_for_cutoff);
    } else {
        fprintf(stderr, "   STYLE:          [Standard Pure Power-Law]\n");
    }
    fprintf(stderr, "----------------------------------------------------------\n");
    if (disk_params->enable_photoevaporation) {
        fprintf(stderr, " PHOTOEVAPORATION: [ONLINE]\n");
        fprintf(stderr, " MODEL PROFILE:    [%s]\n", disk_params->photoevaporation_mode_string);
        fprintf(stderr, " X-RAY LUM_LX:     [%.2e erg/s]\n", disk_params->xray_luminosity);
    } else {
        fprintf(stderr, " PHOTOEVAPORATION: [OFFLINE] (Pure viscous hydrodynamics)\n");
    }
    fprintf(stderr, "==========================================================\n\n");
}

void printBoundaryConditionsStatus(const DiskParameters *disk_params) {
    fprintf(stderr, "==========================================================\n");
    fprintf(stderr, " BOUNDARY CONDITIONS:\n");
    fprintf(stderr, "   INNER BC:       [%s]\n", disk_params->inner_bc_string);
    fprintf(stderr, "   OUTER BC:       [%s]\n", disk_params->outer_bc_string);
    fprintf(stderr, "----------------------------------------------------------\n");
}

void printDustSmoothingStatus(const SimulationOptions *sim_opts, const DiskParameters *disk_params) {
    if (sim_opts->option_for_dust_drift) {
        fprintf(stderr, " DUST SMOOTHING:\n");
        fprintf(stderr, "   METHOD:         [%s]\n", disk_params->dust_smoothing_mode_string);
        fprintf(stderr, "   GAUSS SIGMA:    [%.2f grid units]\n", disk_params->gaussian_sigma);
        fprintf(stderr, "   GAUSS CUTOFF:   [%.2f sig]\n", disk_params->gaussian_cutoff);
    } else {
        fprintf(stderr, " DUST SMOOTHING:   [NO DUST PRESENT]\n");
    }
    fprintf(stderr, "==========================================================\n\n");
}