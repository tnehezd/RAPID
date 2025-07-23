// src/disk_model.c

#include "disk_model.h"   // Saját header
#include "config.h"       // Globális változók és konstansok
#include "simulation_types.h" // ide kell majd a long double-ös disk_t struktúra deklarációja

#include "dust_physics.h" // press, dpress, u_gas függvények deklarációi
#include "io_utils.h"     // sigIn és egyéb I/O függvények deklarációi (ha használja)
#include "utils.h" // Hogy a disk_model.c lássa a Parabola prototípusát

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// *****************************************************************************
// FONTOS: Ahhoz, hogy ez a fájl helyesen működjön, a 'disk_t' struktúrát
//         és az összes kapcsolódó függvénydeklarációt is long double-re kell
//         módosítani a 'disk_model.h' (vagy 'simulation_types.h') fájlban!
// *****************************************************************************


/*	A korong parametereinek beolvasasa	*/
void disk_param_be(disk_t *disk_params) {
    // Ellenőrzés, ha a pointer NULL (jó gyakorlat)
    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [disk_param_be]: Received NULL disk_params pointer.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "DEBUG [disk_param_be]: Calculating derived disk parameters and writing to output file.\n");

    // A PDENSITYDIMLESS számítása a PDENSITY, csillagtömeg és más konstansok alapján
    // Minden long double-re konvertálva
    disk_params->PDENSITYDIMLESS = disk_params->PDENSITY / (long double)SUN2GR * (long double)AU2CM * (long double)AU2CM * (long double)AU2CM;

    fprintf(stderr, "DEBUG [disk_param_be]: Calculated PDENSITY = %.2Le, PDENSITYDIMLESS = %.2Le.\n", // %Le for long double
           disk_params->PDENSITY, disk_params->PDENSITYDIMLESS);
}


/*	r vektor (gridcellák) inicializálása	*/
void load_R(disk_t *disk_params) {
	
	int i;
 	for(i = 0; i <= disk_params->NGRID + 1; i++) {						/*	load an array of radii	*/
 		disk_params->rvec[i] = disk_params->RMIN + (long double)(i-1) * disk_params->DD; // Casting i-1 to long double
	}
}

/*	a sigmara kezdeti profil betoltese	*/
void Initial_Profile(disk_t *disk_params){		/*	initial profile of sigma		*/

 	int i;
 	
 	for(i = 1; i <= disk_params->NGRID; i++) {
 		// pow függvény long double verziója: powl()
 		disk_params->sigmavec[i] = disk_params->SIGMA0 * powl(disk_params->rvec[i], disk_params->SIGMAP_EXP);		/*	sigma0*r^x (x could be eg. -1/2)	*/
 	}
 	
 	// Feltételezve, hogy a Perem függvény is long double tömböt vár
 	Perem(disk_params->sigmavec, disk_params);

}

void Initial_Press(disk_t *disk_params){		/*	initial profile of pressure		*/

 	int i;
 	
 	for(i = 1; i <= disk_params->NGRID; i++) {
 		// press függvény is long double paramétereket fogad és long double-t ad vissza
 		disk_params->pressvec[i] = press(disk_params->sigmavec[i], disk_params->rvec[i], disk_params);
 	}
 	// Feltételezve, hogy a Perem függvény is long double tömböt vár
 	Perem(disk_params->pressvec, disk_params);
}

void Initial_dPress(disk_t *disk_params){		/*	initial profile of pressure		*/

	// dpress függvény is long double-t kezel
	dpress(disk_params);
 	// Feltételezve, hogy a Perem függvény is long double tömböt vár
 	Perem(disk_params->dpressvec, disk_params);
}

/*	ug vektor feltoltese az u_gas ertekevel	*/
void Initial_Ugas(disk_t *disk_params){		/*	initial profile of pressure		*/
	
	// u_gas függvény is long double-t kezel
	u_gas(disk_params);
 	// Feltételezve, hogy a Perem függvény is long double tömböt vár
 	Perem(disk_params->ugvec, disk_params);
}



// loadSigDust függvény, már a korábbiakban módosítottuk,
// hogy long double paramétereket fogadjon (massin és out)
void loadSigDust(double radin[][2], long double *massin, long double out[][3], int n, const disk_t *disk_params) {
    int i;

    for(i = 0; i < n; i++) {
        if((radin[i][0] >= disk_params->RMIN)) {
            // Minden literálhoz hozzáadtam az 'L' suffix-et, ami `long double` literált jelöl.
            // Ez biztosítja, hogy a számítások long double pontossággal történjenek.
            // M_PI is általában double, így érdemes lehet M_PIl-t használni, ha elérhető (vagy saját konstansot definiálni).
            // A disk_params->DD is feltehetően long double lesz a disk_t struktúra módosítása után.
            out[i][0] = massin[i] / (2.0L * (radin[i][0] - disk_params->DD / 2.0L) * (long double)M_PI * disk_params->DD);
            out[i][1] = radin[i][0]; // radin[i][0] továbbra is double, de ide (out[i][1]) long double-ként mentődik
            
            // rmid számítása long double-ben, majd lefelé kerekítés
            long double rmid = ((long double)radin[i][0] - disk_params->RMIN) / disk_params->DD;
            int rindex = (int) floorl(rmid); // floorl a long double-re
            out[i][2] = (long double)rindex; // tárolás long double-ként
        } else {
            out[i][0] = 0.0L;
            out[i][1] = 0.0L;
            out[i][2] = 0.0L;
        }
    }
}