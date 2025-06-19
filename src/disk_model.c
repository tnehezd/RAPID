// src/disk_model.c

#include "disk_model.h"   // Saját header
#include "config.h"       // Globális változók és konstansok
#include "dust_physics.h" // press, dpress, u_gas függvények deklarációi
#include "io_utils.h"     // sigIn és egyéb I/O függvények deklarációi (ha használja)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/*	A korong parametereinek beolvasasa	*/
void disk_param_be(double *sigma0, double *sdexp, double *Rmin, double *Rmax, double *r_dzei, double *r_dzeo, double *dr_dzei, double *dr_dzeo, double *alph_mod, double *rho_p, double *rho_p_dimless, double *alphav, double *mStar, double *gamma) {

	double dummy, rmin, rmax, drdzei, drdzeo, rdzei, rdzeo, amod, rhop, alpha, mstar, flind;
	int dummy2;
	double sig0, exp;

	fin2 = fopen(filenev2,"r");

           	if(fscanf(fin2,"%lg  %lg %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",&rmin,&rmax,&dummy2,&exp,&sig0,&dummy,&rdzei,&rdzeo,&drdzei,&drdzeo,&amod,&rhop,&alpha,&mstar,&flind) == 15) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 15, mert 15 oszlopot szeretnénk beolvasni.	*/
 			printf("\n\n *********** A korong parameterei sikeresen beolvasva!  *********** \n                  sigma0: %lg, sdexp: %lg\n\n\n", sig0,exp);  
			*sigma0 = sig0;						/*	sigma az r=1AU tavolsagon		*/
			*sdexp = exp; 						/*	a feluleti suruseg profiljanak kitevoje	*/
			*Rmin = rmin;						/*	a korong belso hatara			*/
			*Rmax = rmax;       					/*	a korong kulso hatara			*/
			*r_dzei = rdzei;					/*	a DZE belso hatara			*/
			*r_dzeo = rdzeo;					/*	a DZE kulso harara			*/
			*dr_dzei = drdzei;					/*	a DZE belso hatar atmenet vastagsaga	*/
			*dr_dzeo = drdzeo;					/*	a DZE kulso hatar atmenet vastagsaga	*/
			*alph_mod = amod;					/*	a viszkozitas redukcio merteke		*/
			*rho_p = rhop;						/*	a reszecskek atlagos surusege cgs-ben	*/
			*rho_p_dimless = rhop / SUN2GR * AU2CM * AU2CM * AU2CM;	/*	a reszecskek atlagos surusege (dimtlan)	*/
			*alphav = alpha;					/*	alpha viszkozitas merteke		*/
			*mStar = mstar;						/*	a kozponti csillag tomege		*/
			*gamma = flind;						/*	flaring index				*/

	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
//			getchar();
			exit(EXIT_FAILURE);
   	        }

	fclose(fin2);
	
}


/*	Parabola illesztés a peremen	*/
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd) {

	double x1, x2, x3;	/*	meghatározott x pontok, ahol illesztünk					*/
	double y1, y2, y3;	/*	amit illesztünk a meghatározott pontokban				*/
	double av, bv, cv;	/*	illesztéshez szükséges együtthatók --> ezt adja vissza a függvény	*/

	x1 = RMIN + (i1-1) * dd;
	x2 = RMIN + (i2-1) * dd;
	x3 = RMIN + (i3-1) * dd;
 
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
