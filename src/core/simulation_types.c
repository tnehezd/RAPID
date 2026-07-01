#include "simulation_types.h"

SnapshotMode determineSnapshotMode(const SimulationOptions *sim_opts) {

    if (sim_opts->option_for_evolution != 1)
        return SnapshotNonevolving;

    if (sim_opts->option_for_dust_drift != 1)
        return SnapshotGas;

    if (sim_opts->option_for_dust_growth != 1) {
        if (sim_opts->option_for_dust_secondary_population == 1)
            return SnapshotDriftTwoPop;
        else
            return SnapshotDrift;
    }

    if (sim_opts->option_for_dust_secondary_population == 1)
        return SnapshotGrowthTwoPop;

    return SnapshotGrowth;
}


const char* snapshotModeToString(SnapshotMode mode) {
    switch (mode) {
        case SnapshotNonevolving:       return "NON-EVOLVING DISK";
        case SnapshotGas:               return "GAS-ONLY EVOLUTION";
        case SnapshotDrift:             return "DUST DRIFT";
        case SnapshotGrowth:            return "DUST SIZE EVOLUTION";
        case SnapshotDriftTwoPop:       return "DUST DRIFT FOR TWOPOP";
        case SnapshotGrowthTwoPop:      return "DUST SIZE EVOLUTION FOR TWOPOP";
        default:                        return "UNKNOWN";
    }
}
