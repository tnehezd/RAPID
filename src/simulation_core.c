// src/simulation_core.c

// Standard C Library Includes
#include <stdio.h>    // For printf, fopen, fclose, fscanf, snprintf, sprintf
#include <stdlib.h>   // For exit, EXIT_FAILURE, EXIT_SUCCESS, system
#include <math.h>     // For M_PI, fmod, HUGE_VAL (and pow if used by other functions)
#include <string.h>   // For snprintf, sprintf

// Your Project Header Includes
#include "config.h"       // For PARTICLE_NUMBER, TMAX, WO, RMIN, DT, optdr, opttwopop, optgr, optev, r_dze_i, r_dze_o
#include "io_utils.h"     // For timePar (though not called in tIntegrate, it's io-related), reszecskek_szama, por_be, Print_Sigma, Print_Pormozg_Size, Print_Mass, Print_Sigmad. Also for globals: filenev1, filenev3, fout, foutmicr, massfil
#include "disk_model.h"   // If any disk_model functions are called (e.g., Perem indirectly if sigma/press depend on it) - Though not directly visible in tIntegrate, often needed for global disk parameters. Add if you hit implicit declaration for disk_model functions.
#include "dust_physics.h" // For Count_Mass, secondaryGrowth, find_max, find_min, Get_Sigmad, Get_Radius
#include "utils.h"        // For time_step, Get_Sigma_P_dP, and potentially other utility functions



/*	Kiszamolja az 1D-s driftet	*/
/*    	dr/dt = St/(1+St*St)*H(r)/r*dlnP/dlnr*cs = St/(1+St*St) * (H/r) * (r/P) * (dP/dr) * cs		*/
void eqrhs(double pradius, double dp, double sigma, double ug, double r, double *drdt) {
     
     	double P, H, dPdr, St, csound;
    
	St = Stokes_Number(pradius,sigma);
	H = scale_height(r);  
	P = press(sigma,r);
	dPdr = dp;
	csound = c_sound(r); 

     	*drdt = ug / (1. + St * St) + St / (1. + St * St) * H / P * dPdr * csound;	/* 	bearamlas sebessege: Birnstiel PHD	*/
     
}


/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 3*nu = Coeff_1 	*/
double Coeff_1(double r){					
  
  	double A;
  	A = 3.0 * visc(r);
  	return A;
  
}

/* 	for solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*hu/(2*r)*dsigma/dr	--> 9*nu /(2*r) = Coeff_2 	*/
double Coeff_2(double r){							
  
  	double B;
  	B = 9.0 * visc(r) / (2.0 * r);
  	return B;
    
}

double time_step(double *rvec) {					/*	a diffúziós egyenlet megoldásához szükséges minimális időlépés	*/

	double A_max, stepping;
	int i;

	A_max = -10000.0;
	
	for(i = 0; i < NGRID; i++) {
		if(Coeff_1(rvec[i]) > A_max) {
			A_max = Coeff_1(rvec[i]);
		}
	}

	stepping = DD * DD / (2.0 * A_max);

	return stepping;

}


/*	Runge-Kutta4 integrator	*/
// prad bemenet: AU-ban!
void int_step(double time, double prad, double *pressvec, double *dpressvec, double *sigmavec, double *sigmad, double *rdvec, double *rvec, double *ugvec, double step, double y, double *ynew, double *pradnew) {
	double dy1,dy2,dy3,dy4;
	double ytemp, ytemp2;
	double sigma, dpress, ugas; 
	double pdens, p;
	double pradtemp;
	int opt = 0;
	double sigmadd = 0.0;
	
/*	Mivel a kulongozo parametereket csak a megadott gridcella pontokban ismerjuk, de ez nem feltetlen egyezik meg a reszecskek pozicijaval, ezert minden fontos parametert interpolalunk a reszecskek tavolsagara	*/
	interpol(sigmavec,rvec,y,&sigma,DD,opt);
	interpol(dpressvec,rvec,y,&dpress,DD,opt);
	interpol(ugvec,rvec,y,&ugas,DD,opt);

	double dd = (RMAX - RMIN) / (PARTICLE_NUMBER-1);
	int dker = (int)(1./dd);//
	dker = dker * KEREK;
	double ddker = (double) dker;
	int temp;

	temp = (int)floor(y * ddker+0.5);
	ytemp2 = (double)temp / ddker;
	
	int i;
	for(i=0;i<PARTICLE_NUMBER;i++) {
		if(ytemp2 == rdvec[i]) {
			sigmadd = sigmad[i];
			break;
		}
	
	}
	

	if(optgr == 1.) {		// ha van reszecskenovekedes

		if(time != 0.) {	// ha nem t0 idopontban vagyunk
	
			pradtemp = prad;
			interpol(pressvec,rvec,y,&p,DD,opt);
			pdens = PDENSITY; 
			pradtemp = getSize(prad,pdens,sigma,sigmadd,y,p,dpress,step);	// itt szamolja a reszecskenovekedest
			prad = pradtemp;

		}

	}

	*pradnew = prad;

/*	Itt szamolja a reszecske poziciojat	*/
	eqrhs(prad, dpress, sigma, ugas, y, &dy1);

	ytemp = y + 0.5 * step * dy1;

	eqrhs(prad, dpress, sigma, ugas, ytemp, &dy2);
		
	ytemp = y + 0.5 * step * dy2;
	eqrhs(prad, dpress, sigma, ugas, ytemp, &dy3);
	
	ytemp = y + step * dy3;
	eqrhs(prad, dpress, sigma, ugas, ytemp, &dy4);

	*ynew = y + step * (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;

}


/*	Itt vegzi el az integralast, ha szukseg van ra	*/
void tIntegrate(char *nev, double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double *ugvec) {

   	int linesout;
   	PARTICLE_NUMBER = 0;
	linesout = 0;

/*	A reszecskek adatait tartalmazo file sorainak szama (azaz a reszecskek szama) elmentve egy integerbe, amennyiben a futas soran szamol a program driftet	*/
	if(optdr == 1.) {
		PARTICLE_NUMBER = reszecskek_szama(linesout,filenev1);  		
	} else {
		PARTICLE_NUMBER = 0;
	}

	char mv[2048], dens_name[1024],size_name[1024], porout[1024], poroutmicr[1024], massout[1024], dust_name[1024], dust_name2[1024];
   	double radius[PARTICLE_NUMBER][2],radiusmicr[PARTICLE_NUMBER][2],radiussec[4*PARTICLE_NUMBER][2],radius_rec[PARTICLE_NUMBER][2];
	double massvec[PARTICLE_NUMBER], massmicrvec[PARTICLE_NUMBER], masssecvec[4*PARTICLE_NUMBER];
	double max = 0.0, min = 0.0, max2 = 0.0, min2 = 0.0;
	int i;

/*	A reszecskek adatainak beolvasasa file-bol, ha az optdr erteke 1, egyebkent nem szamol a program driftet		*/
	if(optdr == 1.) {

		por_be(radius,radiusmicr,massvec,massmicrvec);			/*	porreszecskek adatainak beolvasasa	*/ 
		int dummy;
		sprintf(mv,"cp %s %s/",filenev1,nev);				/*	az adott helyre, ahova eppen menti a futast eredmenyeit a program, a porreszecskek adatait tartalmazo file atmasolasa, igy kesobb is visszanezheto, hogy mik voltak a kezdeti adatok	*/
		dummy = system(mv);						/*	itt masolja at a file-t			*/	
		(void)dummy;
/*	az aktualis mappaban a pormozgas.dat file letrehozasa: ebbe kerul be a porreszecske tavolsaga es indexe, valamint az adott idolepes	*/
		snprintf(porout,1024,"%s/pormozgas.dat",nev);			
/*	ha 2 populacios a futas, akkor a mikronos pornak is letrehoz egy pormozgas file-t, ebbe kerul be a tavolsag, index es az ido */
		if(opttwopop == 1.) snprintf(poroutmicr,1024,"%s/pormozgasmic.dat",nev);
/*	tomegnovekedesi file letrehozasa az aktualis mappaba - ez lehet, hogy egy kulon opcio lesz a kimeneti adatok meretenek csokkentesere	*/
		snprintf(massout,1024,"%s/mass.dat",nev);
   		fout = fopen(porout,"w");
		if(opttwopop == 1.) foutmicr = fopen(poroutmicr,"w");
   		massfil = fopen(massout,"w");

	}

  	double t = 0.0;
   	double t_integration = TMAX * 2.0 * M_PI; 						/*	numerikus integralas idotartama		*/
	double deltat = time_step(rvec)/5.;							/*	idolepes	*/

	if(DT <= deltat && DT != 0) deltat = DT;	

	double partmassind[PARTICLE_NUMBER][4], partmassmicrind[PARTICLE_NUMBER][4], partmasssecind[4*PARTICLE_NUMBER][4];
	double L = 0.;
	double sigmad[PARTICLE_NUMBER], sigmadm[PARTICLE_NUMBER], sigmads[4*PARTICLE_NUMBER], rdvec[PARTICLE_NUMBER], rmicvec[PARTICLE_NUMBER], rsvec[4*PARTICLE_NUMBER];
	int num = PARTICLE_NUMBER;
	
	double masstempiin = 0, massmtempiin = 0, masstempoin = 0, massmtempoin = 0; // a kulso es belso dze-n felgyulemlett por mennyisege -- bemeneti adat (Print_Mass fuggvenybe)
	double masstempiout = 0, masstempoout = 0, massmtempiout = 0, massmtempoout = 0; // a kulso es belso dze-n felgyulemlett por mennyisege -- kimeneti adat (Print_Mass fuggvenybol)
	double tavin,tavout;

	for(i = 0; i < 4*PARTICLE_NUMBER; i++) {
		radiussec[i][0] = 0;
		radiussec[i][1] = 0;
		partmasssecind[i][0] = 0;
		partmasssecind[i][1] = 0;
		masssecvec[i] = 0;
		sigmads[i] = 0;
		rsvec[i] = 0;
	}

	if(opttwopop == 0) {
		for(i = 0; i < PARTICLE_NUMBER; i++) {
			radiusmicr[i][0] = 0;
			radiusmicr[i][1] = 0;
			partmassmicrind[i][0] = 0;
			partmassmicrind[i][1] = 0;
			massmicrvec[i] = 0;
		}
	}

   	do {

/*	Ha van drift:	*/	

		if(optdr == 1.) {

			if(opttwopop == 1) secondaryGrowth(radius,radiusmicr,radiussec,partmassmicrind,partmasssecind,massmicrvec,masssecvec);
/*	A minimum kereseshez letrehozza a cm-es reszecskek tavolsaganak reciprokat	*/
			for (i=0; i < PARTICLE_NUMBER; i++) {
				if (radius[i][0] > 0. && radius[i][0] > RMIN) {
					radius_rec[i][0] = 1. / radius[i][0];
				} else {
					radius_rec[i][0] = 0.;
				}
			}
		
			max = find_max(radius,PARTICLE_NUMBER);					/*	Megkeresi, hogy melyik a legtavolabbi cm-es reszecske a kozponti csillagtol	*/
			min = find_max(radius_rec,PARTICLE_NUMBER);				/*	Megkeresi a tavolsag reciprokanak maximumat, azaz a legkisebb tavolsagra levo cm-es reszecsket	*/
			min = 1. / min;

			double mint, maxt;

/*	ha 2 populacios a szimulacio, a fentihez hasonloan megkeresi a legnagyobb es a legkisebb tavolsagra levo mikronos reszecsket	*/
			if(opttwopop == 1) {
				for (i=0; i < PARTICLE_NUMBER; i++) {
					if (radiusmicr[i][0] > 0. && radiusmicr[i][0] > RMIN) {
						radius_rec[i][0] = 1. / radiusmicr[i][0];
					} else {
						radius_rec[i][0] = 0.;
					}
				}

				max2 = find_max(radiusmicr,PARTICLE_NUMBER);					/*	Megkeresi, hogy melyik a legtavolabbi reszecske a kozponti csillagtol	*/
				min2 = find_max(radius_rec,PARTICLE_NUMBER);
				min2 = 1. / min2;

/*	megnezi, hogy mely reszecske van a legkozelebb, illetve legtavolabb (mikronos, vagy cm-es)	*/
				mint = find_min(min,min2,HUGE_VAL);
				maxt = find_min(1. / max, 1./max2,HUGE_VAL);
				maxt =  1./ maxt;
			} else {
				mint = min;
				maxt = max;
			}

/*	Ha a legtavolabbi reszecske tavolsaga nagyobb, mint RMIN (azaz meg a szimulacio tartomanyaban van), es a min es a max nem egyenlo (azaz nem gyult pl. ossze az osszes porreszecske 1 helyen), akkor a program tovabb szamol, egyebkent a futas leall	-- ez persze 1 dze eseten mukodik, meg kell oldani, hogy 2 dze eseten is lealljon akkor, ha az osszes por osszegyult a nyomasi maximumokban - meg kell persze csinalni, hogy ez is opcionalis legyen	*/
			if(maxt >= RMIN && mint != maxt) {

				double time = t / 2.0 / M_PI;

				if((fmod(time, (TMAX/WO)) < deltat || time == 0) && L-time < deltat){

/*	Az adatok kiirasahoz szukseges file-ok neveinek elmentese	*/
					if (t==0) {
						snprintf(dens_name,1024,"%s/surface.dat",nev);
					} else {
						if(optev == 1) {
							snprintf(dens_name,1024,"%s/dens.%d.dat",nev,(int)L);
						}
					}

					snprintf(dust_name,1024,"%s/dust.%i.dat",nev,(int)L);
					snprintf(dust_name2,1024,"%s/dustmic.%i.dat",nev,(int)L);
					snprintf(size_name,1024,"%s/size.%d.dat",nev,(int)L);
	
					if(t==0) {
/*	A reszecskek tomeget tartalmazo tomb inicializalasa		*/
						Count_Mass(radius,partmassind,massvec,t,PARTICLE_NUMBER);
						if(opttwopop==1) Count_Mass(radiusmicr,partmassmicrind,massmicrvec,t,PARTICLE_NUMBER);
						if(opttwopop==1) Count_Mass(radiussec,partmasssecind,masssecvec,t,4*PARTICLE_NUMBER);

/*	Ha van tomegnovekedes, akkor a por feluletisurusegenek kiszamolasa itt tortenik	*/
						if(optgr == 1.) {			
							Get_Sigmad(L,max,min,radius,radiusmicr,radiussec,sigmad,sigmadm,sigmads,massvec,massmicrvec,masssecvec,rdvec,rmicvec,rsvec);
						}
					}
		
/*	A sigma, p, dp kiirasa egy file-ba	*/
					if(optev == 1 || time == 0) {
						Print_Sigma(dens_name, rvec, sigmavec, pressvec, dpressvec);
					}

/*	Ha szamol a futas driftet, itt irja ki a reszecskek tavolsagat es meretet	*/
					if(optdr == 1) {
						Print_Pormozg_Size(size_name,(int)L,radius,radiusmicr);
					}

/*	A tomegnovekedesi file-ba az adatok kiirasa	*/
					masstempiout = 0; 
					massmtempiout = 0;
					masstempoout = 0;
					massmtempoout = 0;

					Print_Mass(L,rvec,partmassind,partmassmicrind,partmasssecind,dpressvec,masstempiin,masstempoin,massmtempiin,massmtempoin,&masstempiout,&masstempoout,&massmtempiout,&massmtempoout,&tavin,&tavout);

					if(r_dze_i != tavin) {
						masstempiin = masstempiout; 
						massmtempiin = massmtempiout;
					}
					if(r_dze_o != tavout) {
						masstempoin = masstempoout;
						massmtempoin = massmtempoout;
					}

/*	Ha van pornovekedes, kiirja a por felultisuruseget egy file-ba --> a pornovekedeshez szukseges egyaltalan ezt kiszamolni!	*/
					if(optgr == 1.) {
						Print_Sigmad(dust_name,dust_name2,rdvec,rmicvec,sigmad,sigmadm);
					}

					L = L+(double)(TMAX/WO);

				}

/*	Ha az optev erteke 1, akkor megoldja minden lepesben a sigmara vonatkozo diffuzios egyenletet	*/
				if(optev == 1.) {
					Get_Sigma_P_dP(rvec, sigmavec, pressvec, dpressvec, deltat);
				} 

/*	A reszecskek tomeget tartalmazo tomb adatainak frissitese	*/
				Count_Mass(radius,partmassind,massvec,t,PARTICLE_NUMBER);
				if(opttwopop==1) Count_Mass(radiusmicr,partmassmicrind,massmicrvec,t,PARTICLE_NUMBER);
				if(opttwopop==1) Count_Mass(radiussec,partmasssecind,masssecvec,t,4*PARTICLE_NUMBER);

/*	Ha van reszecskenovekedes, akkor kiszamolja a por feluletisuruseget	*/
				if(optgr == 1.) {
					Get_Sigmad(L,max,min,radius,radiusmicr,radiussec,sigmad,sigmadm,sigmads,massvec,massmicrvec,masssecvec,rdvec,rmicvec,rsvec);
				}

				int optsize = 0;		// ezt ki lehetne siman valtani opttwopop-pal!
/*	A cm-es reszecskek eseten az optsize erteke 0	*/
/*	Itt szamolja ki a program a cm-es reszecskek uj tavolsagat (es meretet, ha kell)	*/
				Get_Radius(nev,optsize,radius,pressvec,dpressvec,sigmavec,sigmad,rdvec,rvec,ugvec,deltat,t,PARTICLE_NUMBER);
 
/*	Ha a futas 2 populacios, akkor az optsize erteke 1	*/
/*	Itt szamolja ki a program a mikoronos reszecskek uj tavolsagat	*/
				if(opttwopop == 1.) {
					optsize = 1;
					Get_Radius(nev,optsize,radiusmicr,pressvec,dpressvec,sigmavec,sigmad,rdvec,rvec,ugvec,deltat,t,PARTICLE_NUMBER); 
					optsize = 2;
					Get_Radius(nev,optsize,radiussec,pressvec,dpressvec,sigmavec,sigmad,rdvec,rvec,ugvec,deltat,t,4*PARTICLE_NUMBER); 
				}

				t = t + deltat;						/*	Idoleptetes	*/

			} else {				/*	Ha a legmesszebbi reszecske tavolsaga mar nem nagyobb, vagy egyenlo, mint RMIN, vagy a legkisebb tavolsagra levo reszecske tavolsaga, akkor a program "figyelmezteto szoveg" mellett sikeresen kilep, nem fut "feleslegesen" tovabb.	*/
				printf("A program sikeresen lefutott az integralasi ido vege elott (t: %lg). \n\nNyomj ENTER-t a kilepeshez!\n",L);
				exit(EXIT_SUCCESS);
			}	
	
		} else {	/*	Ez az az eset, ha a program nem szamol driftet, azaz csak a gaz feluletisurusegenek fejlodesere vagyunk kivancsiak	*/

			double time = t / 2.0 / M_PI;

			if((fmod(time, (TMAX/WO)) < deltat || time == 0) && L-time < deltat){	

				if (t==0) {
					snprintf(dens_name,1024,"%s/surface.dat",nev);
				} else {
					snprintf(dens_name,1024,"%s/dens.%d.dat",nev,(int)L);
				}
	
				Print_Sigma(dens_name, rvec, sigmavec, pressvec, dpressvec);

				L = L+(double)(TMAX/WO);
			}

			Get_Sigma_P_dP(rvec, sigmavec, pressvec, dpressvec, deltat);
			t = t + deltat;						/*	Idoleptetes		*/
		}

   	} while (t <= t_integration);
  

/*	Az idoleptetes leteltevel a program sikeresen kilep	*/
	printf("\n\nA program sikeresen lefutott, azonban elkepzelheto, hogy az integralasi ido nem volt elegendo. A legtavolabbi reszecske %lg CsE tavolsagra van a kozponti csillagtol. \n\nNyomj ENTER-t a kilepeshez!\n",max); 


}


void secondaryGrowth(double rad[][2], double radmicr[][2], double radsec[][2], double partmicind[][4], double partsecind[][4], double *massmicvec, double *masssecvec) {

	int i, j; //, k=0;
	int step = 0;
	double temp[PARTICLE_NUMBER], temp2[4*PARTICLE_NUMBER], tempmic[PARTICLE_NUMBER], temp2ch1[PARTICLE_NUMBER], temp2ch2[PARTICLE_NUMBER], temp2ch3[PARTICLE_NUMBER], temp2ch4[PARTICLE_NUMBER];
	int hist[PARTICLE_NUMBER], histmic[PARTICLE_NUMBER];

/*	A porszemcse generalo program gridfelbontasa -- nem feltetlen egyezik meg a jelenlegi gridfelbontassal!	*/
	double dd = (RMAX - RMIN) / (PARTICLE_NUMBER-1), rtempvec[PARTICLE_NUMBER];
	double temptemp[4*PARTICLE_NUMBER], masstemp[4*PARTICLE_NUMBER][2], sizetemp[4*PARTICLE_NUMBER],mtp[4*PARTICLE_NUMBER], indtemp[4*PARTICLE_NUMBER][2];

	for(i=0; i < 4*PARTICLE_NUMBER; i++) {	
		temp2[i] = radsec[i][0];			// masodlagosan novesztett porszemcsek tavolsaga
		sizetemp[i] = radsec[i][1];
		temptemp[i] = 0;				// atmeneti vektor, amely elmeni, hogy hol noveszt uj reszecsket a program
		masstemp[i][0] = 0;				// atemeneti tomb, amely a novesztett reszecske tomeget tarolja el
		masstemp[i][1] = 0;				// atmeneti tomb, amely a novesztett reszecske tavolsagat tarolja el
		mtp[i] = masssecvec[i];
		if(temp2[i] != 0) step++;			// megszamolja, hogy hany novesztett reszecske van a fuggvenybe valo belepeskor
		indtemp[i][0] = partsecind[i][3];
	}

	for(i=0; i < PARTICLE_NUMBER; i++) {

		temp[i] = rad[i][0];				// cm-es porszemcsek tavolsaga
		tempmic[i] = radmicr[i][0];			// mikoronos porszemcsek tavolsaga
		rtempvec[i] = RMIN + i * dd;			// atmeneti r vektor
		rtempvec[i] = rtempvec[i] + dd / 2.0;		// ez adja a reszecskek eredeti helyet (a generalo programban ide helyeztuk el a reszecskeket	
		hist[i] = 0;					// cm-es porszemcsek hisztogramja
		histmic[i] = 0;					// mikoronos porszemcsek hisztogramja

		temp2ch1[i] = temp2[i];
		temp2ch2[i] = temp2[i+PARTICLE_NUMBER];
		temp2ch3[i] = temp2[i+2*PARTICLE_NUMBER];
		temp2ch4[i] = temp2[i+3*PARTICLE_NUMBER];
	}



/*	Megvizsgalja, hogy az adott tavolsagon hany reszecske talalhato	*/
/*	Valamiert csak ugy ad jo kepet, ha az if-es felteleket beirom.	*/
	for(i=0;i<PARTICLE_NUMBER;i++) {
		if(temp[i] != 0) histogram(temp[i],hist,dd);
		if(tempmic[i] != 0) histogram(tempmic[i],histmic,dd);
		if(temp2ch1[i] != 0) histogram(temp2ch1[i],hist,dd);
		if(temp2ch2[i] != 0) histogram(temp2ch2[i],hist,dd);
		if(temp2ch3[i] != 0) histogram(temp2ch3[i],hist,dd);
		if(temp2ch4[i] != 0) histogram(temp2ch4[i],hist,dd);
	}

	j=0;
	for(i=0;i<PARTICLE_NUMBER;i++) {
/*	Megvizsgalja, hogy kifogyott-e az adott tavolsagrol a cm-es por. Ha igen, viszont az adott tavolsagon van mikronos, akkor novel egy ujabb reszecsket, ha a mikronos por altal kepviselt tomeg az eredeti helyen levo tomeg szazadanal nagyobb	*/
		if(hist[i] == 0 && histmic[i] != 0) {
			if(rtempvec[i] != 0) {
				double sigcurr, siginit;
				sigcurr = partmicind[i][0]/(2.0 * M_PI * dd * (temp[i] - dd/2.));
				siginit = partmicind[i][2]/(2.0 * M_PI * dd * (temp[i] - dd/2.));
				if(sigcurr >= siginit / 100.) {
/*	Uj reszecske novelese --> az adott helyen levo mikronos reszecske tomegenek 75%-at kapja az uj reszecske, ezzel egyutt a mikronos reszecske altal kepviselt reszecske tomege 75%-kal csokken	*/
					massmicvec[i] = partmicind[i][0] * (1.-.75);	// a mikronos reszecske tomegenek csokkentese
					masstemp[j][1] = partmicind[i][0] * .75;	// novesztett reszecske altal kepviselt tomeg
					masstemp[j][0] = rtempvec[i];			// novesztett reszecske altal kepviselt tavolsag
					temptemp[j] = -1*rtempvec[i];			// atmeneti vektor, amely elmenti, hogy mely helyen noveszt uj reszecsket a program
					j++;
				}
			}
		}
	}

/*	A novesztett reszecske tavolsagat tartalmazo atmeneti vektor sorbarendezese	*/
	sort(temptemp,4*PARTICLE_NUMBER);	

	for(i=0; i < 4*PARTICLE_NUMBER-step; i++) {
		if((temptemp[i] != 0)) {
			radsec[i][0] = -1. * temptemp[i];
			temptemp[i] = 0;
		} else {
			radsec[i][0] = 0;
			radsec[i][1] = 0;
		}
	}

	for(i=0; i < 4*PARTICLE_NUMBER; i++) {
		temptemp[i] = radsec[i][0];
		radsec[i][0] = 0;
		radsec[i][1] = 0;
	}

	sort(temptemp,4*PARTICLE_NUMBER);
//	double index = find_max(indtemp,4*PARTICLE_NUMBER);

	for(i=0; i < 4*PARTICLE_NUMBER; i++) {
		if(temptemp[i] >= RMIN) {
			radsec[i][0] = temptemp[i];
			for(j = 0; j < 4*PARTICLE_NUMBER; j++) {
				if(masstemp[j][0] == radsec[i][0]) {
					masssecvec[i] = masstemp[j][1];
					partsecind[i][0] = masssecvec[i];
				}
				if(temp2[j] == radsec[i][0]) {
					masssecvec[i] = mtp[j];
					radsec[i][1] = sizetemp[j];
					partsecind[i][3] = indtemp[j][0];
					partsecind[i][0] = masssecvec[i];
				} else {
					radsec[i][1] = 1e-6/AU2CM;
				}
			}
		}
	}
}