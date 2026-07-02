#ifndef PRINT_PANELS_H
#define PRINT_PANELS_H

#include "parser.h"
#include "disk_model.h"
#include "simulation_types.h"
#include "init_tool_module.h"

void printCenteredTerminal(const char *message);
void printCenteredText(const char *message);
void printSimulationWelcomeBanner(void);
void printHeader(const char *title);
void printWarningPanelForPhotoevaporation(const DiskParameters *disk_params);
void printRunConfigurationHeader(const ParserOptions *def, const DiskParameters *disk_params);
void printBoundaryConditionsStatus(const DiskParameters *disk_params);
void printDustSmoothingStatus(const SimulationOptions *sim_opts, const DiskParameters *disk_params);
void printInitializationParameters(const InitializeDefaultOptions *default_options, long double current_sigma0_gas);
void printFatalErrorMessageForDiskMass();

#endif // PRINT_PANELS_H