#ifndef SIMULATION_CORE_H
#define SIMULATION_CORE_H

#include "disk_model.h"         // disk_t struktúra miatt
#include "simulation_types.h"   // simulation_options_t és output_files_t struktúrák miatt

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
void eqrhs(double pradius, double dp, double sigma, double ug, double r, double *drdt, const disk_t *disk_params);

/**
 * @brief Kiszámítja a viszkozitással kapcsolatos Coeff_1 együtthatót.
 * A diffúziós egyenlet 3 * nu tagja.
 * @param r Aktuális sugárpozíció [AU].
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @return Az Coeff_1 értéke.
 */
double Coeff_1(double r, const disk_t *disk_params);

/**
 * @brief Kiszámítja a viszkozitással kapcsolatos Coeff_2 együtthatót.
 * A diffúziós egyenlet 9 * nu / (2 * r) tagja.
 * @param r Aktuális sugárpozíció [AU].
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @return Az Coeff_2 értéke.
 */
double Coeff_2(double r, const disk_t *disk_params);

/**
 * @brief Kiszámítja a minimális időlépést a szimulációhoz.
 * @param rvec A diszk sugárkoordinátáinak vektora.
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @return A számított időlépés.
 */
double time_step(const disk_t *disk_params);

/**
 * @brief Egyetlen Runge-Kutta 4 lépést hajt végre a részecske mozgására.
 * @param time Aktuális szimulációs idő.
 * @param prad Aktuális részecskesugár.
 * @param pressvec Nyomás vektor a diszk rácspontjain.
 * @param dpressvec Nyomásgradiens vektor a diszk rácspontjain.
 * @param sigmavec Gáz felszíni sűrűség vektor a diszk rácspontjain.
 * @param sigmad Por felszíni sűrűség (adott részecske rácspontján).
 * @param rdvec A por részecskék radiális koordinátái.
 * @param rvec A diszk rácspontjainak radiális koordinátái.
 * @param ugvec Gáz sebesség vektor a diszk rácspontjain.
 * @param step Az időlépés mérete.
 * @param y A részecske aktuális radiális pozíciója.
 * @param ynew A részecske új radiális pozíciója (kimenet).
 * @param pradnew A részecske új sugara (kimenet, növekedés után).
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @param sim_opts A szimulációs opciókat tartalmazó struktúra.
 */
void int_step(double time, double prad, const double *sigmad, const double *rdvec, double step, double y, double *ynew, double *pradnew, const disk_t *disk_params, const simulation_options_t *sim_opts);

/**
 * @brief A fő időintegrációs ciklus.
 * Koordinálja a részecskék mozgását, növekedését,
 * és a gázfelszíni sűrűség evolúcióját, valamint a kimeneti adatok írását.
 * @param output_dir_name A kimeneti könyvtár neve.
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @param sim_opts A szimulációs opciókat tartalmazó struktúra.
 * @param output_files A kimeneti fájl mutatókat tartalmazó struktúra.
 */
void tIntegrate(disk_t *disk_params, const simulation_options_t *sim_opts, output_files_t *output_files);


void secondaryGrowth(double rad[][2], double radmicr[][2], double radsec[][2], double partmicind[][4], double partsecind[][4], double *massmicvec, double *masssecvec, const disk_t *disk_params);

#endif // SIMULATION_CORE_H