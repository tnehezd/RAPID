#ifndef SIMULATION_CORE_H
#define SIMULATION_CORE_H

#include "disk_model.h"         // DiskParameters struktúra miatt
#include "simulation_types.h"   // SimulationOptions és OutputFiles struktúrák miatt

// FÜGGVÉNY PROTOTÍPUSOK

/**
 * @brief Kiszámítja az 1D-s részecske drift sebességét (dr/dt).
 * @param particle_radius A részecske sugara [AU].
 * @param pressure_gradient Nyomásgradiens (dP/dr) [dimenziótlan].
 * @param gas_surface_density A gáz felszíni sűrűsége [CGS].
 * @param gas_velocity A gáz sebessége (valószínűleg a viszkózus ráta) [cm/s].
 * @param radial_distance Aktuális sugárpozíció [AU].
 * @param drift_velocity A kimeneti drift sebesség [AU/s].
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @return void (az eredményt a drift_velocity pointerbe írja).
 */
void calculate1DDustDrift(double particle_radius, double pressure_gradient, double gas_surface_density, double gas_velocity, double radial_distance, double *drift_velocity, const DiskParameters *disk_params);


/**
 * @brief Kiszámítja a minimális időlépést a szimulációhoz.
 * @param radial_grid A diszk sugárkoordinátáinak vektora.
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @return A számított időlépés.
 */
double calculateTimeStep(const DiskParameters *disk_params);


/**
 * @brief A fő időintegrációs ciklus.
 * Koordinálja a részecskék mozgását, növekedését,
 * és a gázfelszíni sűrűség option_for_evolutionúcióját, valamint a kimeneti adatok írását.
 * @param output_dir_name A kimeneti könyvtár neve.
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @param sim_opts A szimulációs opciókat tartalmazó struktúra.
 * @param output_files A kimeneti fájl mutatókat tartalmazó struktúra.
 */
void timeIntegrationForTheSystem(DiskParameters *disk_params, const SimulationOptions *sim_opts, OutputFiles *output_files);



#endif // SIMULATION_CORE_H