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


// Fizikai konstansok (ideális esetben ezek egy külön "constants.h" vagy "config.h" fájlban lennének,
// vagy egy dedikált struktúrában, amit átadsz a függvényeknek)
// Most itt definiáljuk őket a példa kedvéért, ha még nincsenek máshol elérhetővé téve.
#ifndef G_CGS
    #define G_CGS 6.674e-8 // Gravitációs állandó cgs-ben
#endif
#ifndef AU2CM
    #define AU2CM 1.496e13 // AU centiméterben
#endif
#ifndef SDCONV
    #define SDCONV 1.0e-3 // SigmaDensityConversion (ha CGS-be konvertálunk) - ellenőrizd az értékét!
#endif
#ifndef SUN2GR
    #define SUN2GR 1.989e33 // Nap tömege grammban
#endif

// A korábbi PDENSITY_DEF_CGS alapjául szolgáló érték.
// Ha ez egy fix érték, akkor maradhat define. Ha számítás, akkor azt a számítást kell ide beírni.
// A korábbi logod alapján ez 1.6e+00-ra jött ki. Ellenőrizd a pontos forrását!
#ifndef PDENSITY_DEFAULT_CGS
    #define PDENSITY_DEFAULT_CGS 1.0 // **FONTOS: Ezt az értéket a VALÓDI alapértelmezett PDENSITY-re cseréld!**
#endif



/*	A korong parametereinek beolvasasa	*/
void disk_param_be(disk_t *disk_params) {
    // Ellenőrzés, ha a pointer NULL (jó gyakorlat)
    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [disk_param_be]: Received NULL disk_params pointer.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "DEBUG [disk_param_be]: Calculating derived disk parameters and writing to output file.\n");

    // <<< FONTOS VÁLTOZTATÁS: Nincs fájlbeolvasás (fopen, fscanf) itt! >>>
    // A függvény most már feltételezi, hogy a `disk_params` struktúra már tartalmazza az összes
    // bemeneti adatot (pl. disk_params->RMIN, STAR_MASS stb.), amit a `main` függvény már feltöltött.

    // A PDENSITY és PDENSITYDIMLESS értékek kiszámítása,
    // a disk_params struktúrában tárolt paraméterek alapján.
    disk_params->PDENSITY = PDENSITY_DEFAULT_CGS; // Ezt a sort ellenőrizd a programod logikája szerint!
                                                // Ha a PDENSITY valójában egy számítás eredménye, akkor azt ide kell tenni.

    // A PDENSITYDIMLESS számítása a PDENSITY, csillagtömeg és más konstansok alapján
    disk_params->PDENSITYDIMLESS = disk_params->PDENSITY / (G_CGS * disk_params->STAR_MASS * SDCONV / (AU2CM * disk_params->RMIN * disk_params->RMIN));
    // Az eredeti kódodban volt egy másik képlet is, ami a SUN2GR-t használta:
    // disk_params->PDENSITYDIMLESS = disk_params->PDENSITY / SUN2GR * AU2CM * AU2CM * AU2CM;
    // Kérlek, ellenőrizd, melyik a helyes dimenziómentesítés a te modellben!
    // A fenti verziót használtam, mert az tűnik konzisztensebbnek azzal, ahogy a G_CGS-t is használtad.

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
