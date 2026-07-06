#include "simulation_types.h"

SnapshotMode determineSnapshotMode(const SimulationOptions *sim_opts) {
    // 1. Legelső ellenőrzés: Egyáltalán van-e evolúció?
    int is_evol = (sim_opts->option_for_evolution == 1);
    int is_drift = (sim_opts->option_for_dust_drift == 1);
    int is_growth = (sim_opts->option_for_dust_growth == 1);
    int is_2pop = (sim_opts->option_for_dust_secondary_population == 1);

    if (!is_evol && !is_drift && !is_growth) return SnapshotNonevolving;
    if (!is_drift && !is_growth) return SnapshotGas;

    // 2. Por-folyamatok (szétválasztva a statikus gázt)
    if (is_drift && !is_growth) {
        if (!is_evol) return is_2pop ? SnapshotDriftStaticGasTwoPop : SnapshotDriftStaticGas;
        return is_2pop ? SnapshotDriftTwoPop : SnapshotDrift;
    }

    if (is_growth) {
        if (!is_evol) return is_2pop ? SnapshotGrowthStaticGasTwoPop : SnapshotGrowthStaticGas;
        return is_2pop ? SnapshotGrowthTwoPop : SnapshotGrowth;
    }

    return SnapshotNonevolving; // Biztonsági visszatérés
}


const char* snapshotModeToString(SnapshotMode mode) {
    switch (mode) {
        case SnapshotNonevolving:           return "NON-EVOLVING DISK";
        case SnapshotGas:                   return "GAS-ONLY EVOLUTION";
        case SnapshotDrift:                 return "DUST DRIFT";
        case SnapshotGrowth:                return "DUST SIZE EVOLUTION";
        case SnapshotDriftTwoPop:           return "DUST DRIFT FOR TWOPOP";
        case SnapshotGrowthTwoPop:          return "DUST SIZE EVOLUTION FOR TWOPOP";
        case SnapshotDriftStaticGas:        return "DUST DRIFT IN STATIC GAS";
        case SnapshotGrowthStaticGas:       return "DUST SIZE EVOLUTION IN STATIC GAS";
        case SnapshotDriftStaticGasTwoPop:  return "DUST DRIFT IN STATIC GAS FOR TWOPOP";
        case SnapshotGrowthStaticGasTwoPop: return "DUST SIZE EVOLUTION IN STATIC GAS FOR TWOPOP";
        default:                            return "UNKNOWN";
    }
}
