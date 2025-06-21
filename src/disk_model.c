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

    printf("DEBUG [disk_param_be]: Calculating derived disk parameters and writing to output file.\n");

    // <<< FONTOS VÁLTOZTATÁS: Nincs fájlbeolvasás (fopen, fscanf) itt! >>>
    // A függvény most már feltételezi, hogy a `disk_params` struktúra már tartalmazza az összes
    // bemeneti adatot (pl. RMIN, STAR_MASS stb.), amit a `main` függvény már feltöltött.

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

    printf("DEBUG [disk_param_be]: Calculated PDENSITY = %.2e, PDENSITYDIMLESS = %.2e.\n",
           disk_params->PDENSITY, disk_params->PDENSITYDIMLESS);

    // Fájlba írás: Most már az aktuális, disk_params-ban lévő értékeket írjuk ki.
    // A `filenev2` globális változó helyett közvetlenül megadjuk a fájlnevet.
    // IDEÁLIS ESETBEN EZT IS PARAMÉTERKÉNT KELLENE ÁTADNI, PL. EGY output_files_t STRUKTÚRÁBAN.
    // Most egy fix fájlnevet használunk, de gondold végig az output fájlok kezelését is.
    const char *output_filename = "disk_param.dat"; // Ez lehetne a sim_opts-ból, ha az output fájlokat is ott tárolod
    FILE *fout = fopen(output_filename, "w");
    if (fout == NULL) {
        fprintf(stderr, "ERROR [disk_param_be]: Could not open output file '%s'.\n", output_filename);
        perror("Reason"); // Kiírja a pontos rendszerhiba üzenetet
        exit(EXIT_FAILURE);
    }

    // A disk_params struktúra tagjait írjuk ki
    fprintf(fout, "RMIN %lf\n", disk_params->RMIN);
    fprintf(fout, "RMAX %lf\n", disk_params->RMAX);
    fprintf(fout, "NGRID %d\n", disk_params->NGRID);
    fprintf(fout, "SIGMA0 %e\n", disk_params->SIGMA0);
    fprintf(fout, "SIGMAP_EXP %lf\n", disk_params->SIGMAP_EXP);
    fprintf(fout, "ALPHA %e\n", disk_params->alpha_visc);
    fprintf(fout, "STAR_MASS %lf\n", disk_params->STAR_MASS);
    fprintf(fout, "HASP %lf\n", disk_params->HASP);
    fprintf(fout, "FLIND %lf\n", disk_params->FLIND);
    fprintf(fout, "R_DZE_I %lf\n", disk_params->r_dze_i);
    fprintf(fout, "R_DZE_O %lf\n", disk_params->r_dze_o);
    fprintf(fout, "DR_DZE_I %lf\n", disk_params->Dr_dze_i);
    fprintf(fout, "DR_DZE_O %lf\n", disk_params->Dr_dze_o);
    fprintf(fout, "A_MOD %lf\n", disk_params->a_mod);
    fprintf(fout, "PDENSITY %e\n", disk_params->PDENSITY);
    fprintf(fout, "PDENSITYDIMLESS %e\n", disk_params->PDENSITYDIMLESS);

    fclose(fout);
    printf("DEBUG [disk_param_be]: Disk parameters written to %s.\n", output_filename);
}




/*	r vektor (gridcellák) inicializálása	*/
void load_R(double *rvec) {
	
	int i;
 	for(i = 0; i <= NGRID+1; i++) {						/*	load an array of radii	*/
 		rvec[i] = RMIN + (i-1) * DD;
	}

}

/*	a sigmara kezdeti profil betoltese	*/
void Initial_Profile(double *sigmavec, double *r){		/*	initial profile of sigma		*/

  	int i;
  
  	for(i = 1; i <= NGRID; i++) {
    		sigmavec[i] = SIGMA0 * pow(r[i],SIGMAP_EXP);		/*	sigma0*r^x (x could be eg. -1/2)	*/
  	}
  
  	Perem(sigmavec);

}

void Initial_Press(double *pressvec, double *sigmavec, double *rvec){		/*	initial profile of pressure		*/

  	int i;
  
  	for(i = 1; i <= NGRID; i++) {
    		pressvec[i] = press(sigmavec[i],rvec[i]);
  	}
  
  	Perem(pressvec);

}

void Initial_dPress(double *dpressvec, double *pressvec){		/*	initial profile of pressure		*/

	dpress(dpressvec,pressvec);
   	Perem(dpressvec);

}

/*	ug vektor feltoltese az u_gas ertekevel	*/
void Initial_Ugas(double *sigmavec, double *rvec, double *ug){		/*	initial profile of pressure		*/

	u_gas(sigmavec,rvec,ug);
  	Perem(ug);
}



void loadSigDust(double radin[][2], double *massin, double out[][3], double dd, int n) {

	int i;

	for(i=0;i<n;i++){

/*	cm-es por feluletisurusegenek kiszamolasa	*/
/*	ha a reszecske tavolsaga nagyobb, mint RMIN, azaz a szamolas tartomanyan belul van, a feluletisuruseget az altala kepviselt tomegbol szamolja vissza a program	*/
		if((radin[i][0] >= RMIN)) {
			out[i][0] = massin[i] / (2. * (radin[i][0]-dd/2.) * M_PI * dd);	// sigma = m /(2 * r * pi * dr) --> itt a dr az a tavolsag, ami porreszecske "generalo" programban az eredeti gridfelbontas
			out[i][1] = radin[i][0];					// elmenti a reszecske tavolsagat is

  			double rmid = (radin[i][0] - RMIN) / dd;     						/* 	The integer part of this gives at which index is the body			*/
			int rindex = (int) floor(rmid);							/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
			out[i][2] = (double) rindex;

/*	ha a reszecske RMIN-en belul van, akkor az o "tavolsagaban" a feluletisuruseg 0		*/	
		} else {
			out[i][0] = 0;
			out[i][1] = 0;					// r = 0, mert ha RMIN-en belulre kerult a reszecske, akkor a program automatikusan kinullazza a reszecske tavolsagat. Itt tehat a sigdtemp[i][1] = 0 lesz!
			out[i][2] = 0;
		}
	}
}
