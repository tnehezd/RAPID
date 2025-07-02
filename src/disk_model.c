// src/disk_model.c

#include "disk_model.h"   // Saját header
#include "config.h"       // Globális változók és konstansok
#include "simulation_types.h"

#include "dust_physics.h" // press, dpress, u_gas függvények deklarációi
#include "io_utils.h"     // sigIn és egyéb I/O függvények deklarációi (ha használja)
#include "utils.h" // Hogy a disk_model.c lássa a Parabola prototípusát

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



/*	A korong parametereinek beolvasasa	*/
void disk_param_be(disk_t *disk_params) {
    // Ellenőrzés, ha a pointer NULL (jó gyakorlat)
    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [disk_param_be]: Received NULL disk_params pointer.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "DEBUG [disk_param_be]: Calculating derived disk parameters and writing to output file.\n");


    // A PDENSITYDIMLESS számítása a PDENSITY, csillagtömeg és más konstansok alapján
    disk_params->PDENSITYDIMLESS = disk_params->PDENSITY / (G_GRAV_CONST * disk_params->STAR_MASS * SDCONV / (AU2CM * disk_params->RMIN * disk_params->RMIN));
    // Az eredeti kódodban volt egy másik képlet is, ami a SUN2GR-t használta:
    // disk_params->PDENSITYDIMLESS = disk_params->PDENSITY / SUN2GR * AU2CM * AU2CM * AU2CM;
    // Kérlek, ellenőrizd, melyik a helyes dimenziómentesítés a te modellben!
    // A fenti verziót használtam, mert az tűnik konzisztensebbnek azzal, ahogy a G_GRAV_CONST-t is használtad.

    fprintf(stderr, "DEBUG [disk_param_be]: Calculated PDENSITY = %.2e, PDENSITYDIMLESS = %.2e.\n",
           disk_params->PDENSITY, disk_params->PDENSITYDIMLESS);
}




/*	r vektor (gridcellák) inicializálása	*/
void load_R(disk_t *disk_params) {
	
	int i;
 	for(i = 0; i <= disk_params->NGRID+1; i++) {						/*	load an array of radii	*/
 		disk_params->rvec[i] = disk_params->RMIN + (i-1) * disk_params->DD;
//        fprintf(stderr, "DEBUG [load_R]: r: %lg\n", disk_params->rvec[i]);
	}
}

/*	a sigmara kezdeti profil betoltese	*/
void Initial_Profile(disk_t *disk_params){		/*	initial profile of sigma		*/

  	int i;
  
  	for(i = 1; i <= disk_params->NGRID; i++) {
    		disk_params->sigmavec[i] = disk_params->SIGMA0 * pow(disk_params->rvec[i],disk_params->SIGMAP_EXP);		/*	sigma0*r^x (x could be eg. -1/2)	*/
    }
  

  	Perem(disk_params->sigmavec,disk_params);

}

void Initial_Press(disk_t *disk_params){		/*	initial profile of pressure		*/
    fprintf(stderr, "DEBUG [Initial_Press]: Entry. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
            (void*)disk_params, disk_params->FLIND, disk_params->HASP);

  	int i;
  
  	for(i = 1; i <= disk_params->NGRID; i++) {
    		disk_params->pressvec[i] = press(disk_params->sigmavec[i],disk_params->rvec[i],disk_params);
  	}
  	Perem(disk_params->pressvec,disk_params);


}

void Initial_dPress(disk_t *disk_params){		/*	initial profile of pressure		*/
    fprintf(stderr, "DEBUG [Initial_dPress]: Entry. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
            (void*)disk_params, disk_params->FLIND, disk_params->HASP);

	dpress(disk_params);

   	Perem(disk_params->dpressvec,disk_params);


}

/*	ug vektor feltoltese az u_gas ertekevel	*/
void Initial_Ugas(disk_t *disk_params){		/*	initial profile of pressure		*/
      fprintf(stderr, "DEBUG [Initial_Ugas]: Entry. disk_params address=%p, FLIND=%.2f, HASP=%.2f\n",
            (void*)disk_params, disk_params->FLIND, disk_params->HASP);
 	
	u_gas(disk_params);
  	Perem(disk_params->ugvec,disk_params);
}



void loadSigDust(double radin[][2], double *massin, double out[][3], int n, const disk_t *disk_params) {

	int i;

	for(i=0;i<n;i++){

/*	cm-es por feluletisurusegenek kiszamolasa	*/
/*	ha a reszecske tavolsaga nagyobb, mint disk_params->RMIN, azaz a szamolas tartomanyan belul van, a feluletisuruseget az altala kepviselt tomegbol szamolja vissza a program	*/
		if((radin[i][0] >= disk_params->RMIN)) {
			out[i][0] = massin[i] / (2. * (radin[i][0]-disk_params->DD/2.) * M_PI * disk_params->DD);	// sigma = m /(2 * r * pi * dr) --> itt a dr az a tavolsag, ami porreszecske "generalo" programban az eredeti gridfelbontas
			out[i][1] = radin[i][0];					// elmenti a reszecske tavolsagat is

  			double rmid = (radin[i][0] - disk_params->RMIN) / disk_params->DD;     						/* 	The integer part of this gives at which index is the body			*/
			int rindex = (int) floor(rmid);							/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
			out[i][2] = (double) rindex;

/*	ha a reszecske disk_params->RMIN-en belul van, akkor az o "tavolsagaban" a feluletisuruseg 0		*/	
		} else {
			out[i][0] = 0;
			out[i][1] = 0;					// r = 0, mert ha disk_params->RMIN-en belulre kerult a reszecske, akkor a program automatikusan kinullazza a reszecske tavolsagat. Itt tehat a sigdtemp[i][1] = 0 lesz!
			out[i][2] = 0;
		}
	}
}
