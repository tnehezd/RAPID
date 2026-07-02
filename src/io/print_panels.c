#include <stdio.h>
#include "print_panels.h"
#include "logger.h"
#include <string.h>
#include "init_tool_module.h"


void printCenteredTerminal(const char *text) {
    int len = strlen(text);
    int padding = (kTerminalWidth - len) / 2;
    
    fprintf(stderr, "||%*s%s%*s||\n", padding, "", text, padding + (kTerminalWidth - len) % 2, "");
}

void printCenteredText(const char *text) {
    int len = strlen(text);
    int padding = (kTerminalWidth - len) / 2;
    
    fprintf(stderr, "%*s%s%*s\n", padding, "", text, padding + (kTerminalWidth - len) % 2, "");
}

static void printLine(void) {
    fprintf(stderr, "==========================================================================\n");
}

void printHeader(const char *title) {
    fprintf(stderr, "\n");
    printLine();
    fprintf(stderr, "|| %-68s ||\n", title);
    printLine();
}


void printSimulationWelcomeBanner(void) {
    
    printLine();
    fprintf(stderr, "||                                                                      ||\n");
    printCenteredTerminal("RAPID SIMULATION CORE");
    printCenteredTerminal("Version: " SIM_VERSION);
    fprintf(stderr, "||                                                                      ||\n");
    printCenteredTerminal("Lead Developer: [D. Tarczay-Nehez]");
    printCenteredTerminal("Co-Developer: [A. Csaszar]");
    fprintf(stderr, "||                                                                      ||\n");
    printLine();
    fprintf(stderr, "\n");
    printCenteredText("Initializing simulation core...\n\n");
}

void printWarningPanelForPhotoevaporation(const DiskParameters *disk_params)
{
    printLine();
    printCenteredTerminal("WARNING!!!\n");
    fprintf(stderr, "Photoevaporation is globally set to FALSE, but a specific mode ('%s')\n", disk_params->photoevaporation_mode_string);
    fprintf(stderr, "and options were provided in the configuration!\n");
    fprintf(stderr, "The simulation will IGNORE these parameters and run WITHOUT photoevaporation.\n");
    printLine();
}

void printRunConfigurationHeader(const ParserOptions *def, const DiskParameters *disk_params) {
    printHeader("INITIAL GAS INITIAL CONDITION PROFILE:");
    if (def->use_cutoff) {
        fprintf(stderr, "   STYLE:          [Exponential Cutoff / Self-Similar Taper]\n");
        fprintf(stderr, "   R_CUTOFF:       [%.2f AU]\n", def->r_cutoff);
        fprintf(stderr, "   SHARPNESS N:    [%.2f]\n", def->n_for_cutoff);
    } else {
        fprintf(stderr, "   STYLE:          [Standard Pure Power-Law]\n");
    }
    printLine();
    if (disk_params->enable_photoevaporation) {
        fprintf(stderr, " PHOTOEVAPORATION: [ONLINE]\n");
        fprintf(stderr, " MODEL PROFILE:    [%s]\n", disk_params->photoevaporation_mode_string);
        fprintf(stderr, " X-RAY LUM_LX:     [%.2e erg/s]\n", disk_params->xray_luminosity);
    } else {
        fprintf(stderr, " PHOTOEVAPORATION: [OFFLINE] (Pure viscous hydrodynamics)\n");
    }
    printLine();
}

void printBoundaryConditionsStatus(const DiskParameters *disk_params) {
    printHeader("BOUNDARY CONDITIONS");
    fprintf(stderr, "   INNER BC:       [%s]\n", disk_params->inner_bc_string);
    fprintf(stderr, "   OUTER BC:       [%s]\n", disk_params->outer_bc_string);
    printLine();
}

void printDustSmoothingStatus(const SimulationOptions *sim_opts, const DiskParameters *disk_params) {
    if (sim_opts->option_for_dust_drift) {
        printHeader("DUST SMOOTHING METHOD");
        fprintf(stderr, "   DUST SMOOTHING:  [ONLINE]\n");
        fprintf(stderr, "   METHOD:         [%s]\n", disk_params->dust_smoothing_mode_string);
        fprintf(stderr, "   GAUSS SIGMA:    [%.2f grid units]\n", disk_params->gaussian_sigma);
        fprintf(stderr, "   GAUSS CUTOFF:   [%.2f sig]\n", disk_params->gaussian_cutoff);
    } else {
        printHeader("DUST SMOOTHING METHOD");
        fprintf(stderr, "   DUST SMOOTHING:   [NO DUST PRESENT]\n");
    }
    printLine();
}



void printInitializationParameters(const InitializeDefaultOptions *default_options, long double current_sigma0_gas) {
    printHeader("INITIALIZATION PARAMETERS");
    fprintf(stderr, "  Stellar mass (Solar Mass):          %lg\n", default_options->star_mass);
    fprintf(stderr, "  Disk aspect ratio (H/R):            %lg\n", default_options->aspect_ratio);
    fprintf(stderr, "  Disk flaring index:                 %lg\n", default_options->flaring_index);
    fprintf(stderr, "  Alpha viscosity:                    %lg\n", default_options->alpha_viscosity);
    fprintf(stderr, "  Total disk mass (Solar Mass):       %lg\n", default_options->total_disk_mass); // FIXED: Renamed to total_disk_mass
    fprintf(stderr, "  Inner disk edge (AU):               %lg\n", default_options->r_inner);
    fprintf(stderr, "  Outer disk edge (AU):               %lg\n", default_options->r_outer);
    if (default_options->use_cutoff) {
        fprintf(stderr, "  Profile Style:                      [Exponential Cutoff Taper]\n");
        fprintf(stderr, "  Characteristic Cutoff Radius (AU):  %lg\n", default_options->r_cutoff);
        fprintf(stderr, "  Cutoff Sharpness Exponent (n):      %lg\n", default_options->n_for_cutoff);
    } else {
        fprintf(stderr, "  Profile Style:                      [Standard Pure Power-Law]\n");
    }
    fprintf(stderr, "  Surface density profile exponent:   %lg\n", fabs(default_options->sigma_exponent));
    fprintf(stderr, "  Gas surface density at 1 AU:        %Lg M_Sun/AU^2\n", current_sigma0_gas);
    fprintf(stderr, "  Dust to gas ratio:                  %lg\n", default_options->dust_to_gas_ratio);
    fprintf(stderr, "  Number of gas grid points:          %d\n", default_options->n_grid_points); 
    fprintf(stderr, "  Number of dust particles:           %d\n", default_options->n_dust_particles);
    printLine();
}


void printFatalErrorMessageForDiskMass(void) {
    printHeader("FATAL RUNTIME ERROR");
    fprintf(stderr, "Overdetermined initial conditions detected.\n");
    fprintf(stderr, "You provided BOTH '-disk_mass' and '-sigma0_init'.\n");
    fprintf(stderr, "Please specify only one method to define the disk's initial scale.\n");
    printLine();

}