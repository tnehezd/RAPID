#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "disk_model.h"         // DiskParameters struktúra miatt
#include "simulation_types.h"   // SimulationOptions és output_files_t struktúrák miatt


/**
 * @brief Egyetlen Runge-Kutta 4 lépést hajt végre a részecske mozgására.
 * @param time Aktuális szimulációs idő.
 * @param prad Aktuális részecskesugár.
 * @param pressvec Nyomás vektor a diszk rácspontjain.
 * @param dpressvec Nyomásgradiens vektor a diszk rácspontjain.
 * @param sigmavec Gáz felszíni sűrűség vektor a diszk rácspontjain.
 * @param sigmad Por felszíni sűrűség (adott részecske rácspontján).
 * @param rdvec A por részecskék radiális koordinátái.
 * @param radial_grid A diszk rácspontjainak radiális koordinátái.
 * @param ugvec Gáz sebesség vektor a diszk rácspontjain.
 * @param step Az időlépés mérete.
 * @param y A részecske aktuális radiális pozíciója.
 * @param ynew A részecske új radiális pozíciója (kimenet).
 * @param pradnew A részecske új sugara (kimenet, növekedés után).
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @param sim_opts A szimulációs opciókat tartalmazó struktúra.
 */
void integrateParticleRungeKutta4(double time, double prad, const double *sigmad, const double *rdvec, double step, double y, double *ynew, double *pradnew, const DiskParameters *disk_params, const SimulationOptions *sim_opts);



#endif // INTEGRATOR_H