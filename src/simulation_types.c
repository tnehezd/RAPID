#include "simulation_types.h"

SnapshotMode determineSnapshotMode(const SimulationOptions *sim_opts)
{
    // 0) Ha nincs semmilyen fejlődés
    if (sim_opts->option_for_evolution != 1)
        return SnapshotNonevolving;

    // 1) Ha nincs por drift → csak gáz
    if (sim_opts->option_for_dust_drift != 1)
        return SnapshotGas;

    // 2) Drift ON, growth OFF
    if (sim_opts->option_for_dust_growth != 1) {
        if (sim_opts->option_for_dust_secondary_population == 1)
            return SnapshotDriftTwoPop;
        else
            return SnapshotDrift;
    }

    // 3) Drift ON, growth ON
    if (sim_opts->option_for_dust_secondary_population == 1)
        return SnapshotGrowthTwoPop;

    return SnapshotGrowth;
}


const char* snapshotModeToString(SnapshotMode mode)
{
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
