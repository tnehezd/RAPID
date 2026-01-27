#include "utils.h"    // Ezt kell includolni, mert ebben lesz a Parabola deklarációja
#include "config.h"   // Szükséges a RMIN és DD makrók miatt, amiket a Parabola használ
#include <math.h>     // Bár a Parabola most nem használ math.h függvényt,
#include <stdlib.h>                      // más utility függvényeknek szüksége lehet rá.
                      // Jó gyakorlat ide tenni.

#include "simulation_types.h" 
#include "dust_physics.h"
#include "gas_physics.h"
#include "boundary_conditions.h"


/*	Parabola illesztés a peremen	*/
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd, const disk_t *disk_params) {

	double x1, x2, x3;	/*	meghatározott x pontok, ahol illesztünk					*/
	double y1, y2, y3;	/*	amit illesztünk a meghatározott pontokban				*/
	double av, bv, cv;	/*	illesztéshez szükséges együtthatók --> ezt adja vissza a függvény	*/

	x1 = disk_params->RMIN + (i1-1) * dd;
	x2 = disk_params->RMIN + (i2-1) * dd;
	x3 = disk_params->RMIN + (i3-1) * dd;
 
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


/*	A peremen parabolat illeszt	*/
void Perem(double *vec, const disk_t *disk_params) {					/*	boundary condition for sigma, p, dp...	*/

//	Parabola(vec, 1, 2, 3, &a, &b, &c, disk_params->DD,disk_params);
//	vec[0] =  a * (disk_params->RMIN - disk_params->DD) * (disk_params->RMIN - disk_params->DD) + b * (disk_params->RMIN - disk_params->DD) + c;
	vec[0] = vec[1];
//	Parabola(vec, disk_params->NGRID - 2, disk_params->NGRID - 1, disk_params->NGRID, &a, &b, &c, disk_params->DD,disk_params);
//	vec[disk_params->NGRID+1] = a * (disk_params->RMAX + disk_params->DD) * (disk_params->RMAX + disk_params->DD) + b * (disk_params->RMAX + disk_params->DD) + c;
	vec[disk_params->NGRID+1] = vec[disk_params->NGRID];
}
