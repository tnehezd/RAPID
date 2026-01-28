#include "utils.h"    // Ezt kell includolni, mert ebben lesz a parabolicExtrapolationToGhostCells deklarációja
#include "config.h"   // Szükséges a r_min és DD makrók miatt, amiket a parabolicExtrapolationToGhostCells használ
#include <math.h>     // Bár a parabolicExtrapolationToGhostCells most nem használ math.h függvényt,
#include <stdlib.h>                      // más utility függvényeknek szüksége lehet rá.
                      // Jó gyakorlat ide tenni.

#include "simulation_types.h" 
#include "dust_physics.h"
#include "gas_physics.h"
#include "boundary_conditions.h"


void parabolicExtrapolationToGhostCells(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd, const disk_t *disk_params) {


/**
 * @brief Parabolic ghost-cell extrapolation for smooth outflow boundaries.
 *
 * This routine computes a second-order polynomial (parabolic) extrapolation
 * based on three interior grid points and evaluates the resulting polynomial
 * at the ghost-cell location. The extrapolated value is then used to fill
 * the ghost cell.
 *
 * This is NOT a physical boundary condition by itself. It does not impose
 * Dirichlet or Neumann constraints. Instead, it provides a smooth numerical
 * continuation of the interior solution beyond the computational domain.
 *
 * In practice, this acts as a non-reflecting, smooth outflow boundary:
 * the field is allowed to "leave" the domain by following its interior trend,
 * without enforcing a fixed value (Dirichlet) or a strict zero gradient
 * (Neumann). This reduces spurious reflections and maintains second-order
 * accuracy at the boundary.
 *
 * Use this method whenever a numerically stable, trend-following outflow
 * closure is desired. For true Dirichlet or Neumann boundaries, use the
 * corresponding explicit ghost-cell formulas instead.
 */


	double x1, x2, x3;	/*	meghatározott x pontok, ahol illesztünk					*/
	double y1, y2, y3;	/*	amit illesztünk a meghatározott pontokban				*/
	double av, bv, cv;	/*	illesztéshez szükséges együtthatók --> ezt adja vissza a függvény	*/

	x1 = disk_params->r_min + (i1-1) * dd;
	x2 = disk_params->r_min + (i2-1) * dd;
	x3 = disk_params->r_min + (i3-1) * dd;
 
	y1 = vec[i1];
	y2 = vec[i2];
	y3 = vec[i3];

	av = ((y1 - y3) / (x1 - x3) - (y1 - y2) / (x1 - x2)) / (x3 - x2);
	bv = (y1 - y2) / (x1 - x2) - av * (x1 + x2);
	cv = y1 - av * x1 * x1 - bv * x1;

	*a = av;
	*b = bv;
	*c = cv;

}


void applyBoundaryConditions(double *vec, const disk_t *disk_params) {					/*	boundary condition for sigma, p, dp...	*/


// OPEN BOUNDARY: mind a sebességre, mind a többi fizikai mennyiségre parabola illesztést használunk
	

//	parabolicExtrapolationToGhostCells(vec, 1, 2, 3, &a, &b, &c, disk_params->DD,disk_params);
//	vec[0] =  a * (disk_params->r_min - disk_params->DD) * (disk_params->r_min - disk_params->DD) + b * (disk_params->r_min - disk_params->DD) + c;
	vec[0] = vec[1];
//	parabolicExtrapolationToGhostCells(vec, disk_params->grid_number - 2, disk_params->grid_number - 1, disk_params->grid_number, &a, &b, &c, disk_params->DD,disk_params);
//	vec[disk_params->grid_number+1] = a * (disk_params->r_max + disk_params->DD) * (disk_params->r_max + disk_params->DD) + b * (disk_params->r_max + disk_params->DD) + c;
	vec[disk_params->grid_number+1] = vec[disk_params->grid_number];
}
