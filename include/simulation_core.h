#ifndef SIMULATION_CORE_H
#define SIMULATION_CORE_H

#include "disk_model.h"         // DiskParameters struktúra miatt
#include "simulation_types.h"   // SimulationOptions és output_files_t struktúrák miatt

// FÜGGVÉNY PROTOTÍPUSOK

/**
 * @brief Kiszámítja az 1D-s részecske drift sebességét (dr/dt).
 * @param pradius A részecske sugara [AU].
 * @param dp Nyomásgradiens (dP/dr) [dimenziótlan].
 * @param sigma A gáz felszíni sűrűsége [CGS].
 * @param ug A gáz sebessége (valószínűleg a viszkózus ráta) [cm/s].
 * @param r Aktuális sugárpozíció [AU].
 * @param drdt A kimeneti drift sebesség [AU/s].
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @return void (az eredményt a drdt pointerbe írja).
 */
void calculate1DDustDrift(double pradius, double dp, double sigma, double ug, double r, double *drdt, const DiskParameters *disk_params);

/**
 * @brief Kiszámítja a viszkozitással kapcsolatos ftcsSecondDerivativeCoefficient együtthatót.
 * A diffúziós egyenlet 3 * nu tagja.
 * @param r Aktuális sugárpozíció [AU].
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @return Az ftcsSecondDerivativeCoefficient értéke.
 */
double ftcsSecondDerivativeCoefficient(double r, const DiskParameters *disk_params);

/**
 * @brief Kiszámítja a viszkozitással kapcsolatos ftcsFirstDerivativeCoefficient együtthatót.
 * A diffúziós egyenlet 9 * nu / (2 * r) tagja.
 * @param r Aktuális sugárpozíció [AU].
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @return Az ftcsFirstDerivativeCoefficient értéke.
 */
double ftcsFirstDerivativeCoefficient(double r, const DiskParameters *disk_params);

/**
 * @brief Kiszámítja a minimális időlépést a szimulációhoz.
 * @param rvec A diszk sugárkoordinátáinak vektora.
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @return A számított időlépés.
 */
double calculateTimeStep(const DiskParameters *disk_params);


/**
 * @brief A fő időintegrációs ciklus.
 * Koordinálja a részecskék mozgását, növekedését,
 * és a gázfelszíni sűrűség evolúcióját, valamint a kimeneti adatok írását.
 * @param output_dir_name A kimeneti könyvtár neve.
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @param sim_opts A szimulációs opciókat tartalmazó struktúra.
 * @param output_files A kimeneti fájl mutatókat tartalmazó struktúra.
 */
void timeIntegrationForTheSystem(DiskParameters *disk_params, const SimulationOptions *sim_opts, output_files_t *output_files);



#endif // SIMULATION_CORE_H