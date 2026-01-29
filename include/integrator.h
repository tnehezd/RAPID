#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "disk_model.h"         // DiskParameters struktúra miatt
#include "simulation_types.h"   // SimulationOptions és output_files_t struktúrák miatt


/**
 * @brief Egyetlen Runge-Kutta 4 lépést hajt végre a részecske mozgására.
 * @param time Aktuális szimulációs idő.
 * @param prad Aktuális részecskesugár.
 * @param gas_pressure_vector Nyomás vektor a diszk rácspontjain.
 * @param gas_pressure_gradient_vector Nyomásgradiens vektor a diszk rácspontjain.
 * @param gas_surface_density_vector Gáz felszíni sűrűség vektor a diszk rácspontjain.
 * @param sigmad Por felszíni sűrűség (adott részecske rácspontján).
 * @param rdvec A por részecskék radiális koordinátái.
 * @param radial_grid A diszk rácspontjainak radiális koordinátái.
 * @param gas_velocity_vector Gáz sebesség vektor a diszk rácspontjain.
 * @param step Az időlépés mérete.
 * @param y A részecske aktuális radiális pozíciója.
 * @param ynew A részecske új radiális pozíciója (kimenet).
 * @param pradnew A részecske új sugara (kimenet, növekedés után).
 * @param disk_params A diszk paramétereit tartalmazó struktúra.
 * @param sim_opts A szimulációs opciókat tartalmazó struktúra.
 */
void integrateParticleRungeKutta4(double time, double prad, const double *sigmad, const double *rdvec, double step, double y, double *ynew, double *pradnew, const DiskParameters *disk_params, const SimulationOptions *sim_opts);



#endif // INTEGRATOR_H