#ifndef PRINT_PANELS_H
#define PRINT_PANELS_H

#include "parser.h"
#include "disk_model.h"
#include "simulation_types.h"

void printWarningPanelForPhotoevaporation(const DiskParameters *disk_params);
void printRunConfigurationHeader(const ParserOptions *def, const DiskParameters *disk_params);
void printBoundaryConditionsStatus(const DiskParameters *disk_params);
void printDustSmoothingStatus(const SimulationOptions *sim_opts, const DiskParameters *disk_params);

#endif // PRINT_PANELS_H