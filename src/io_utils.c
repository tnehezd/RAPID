#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> // For mkdir
#include <unistd.h>   // For access (used in Mk_Dir)
#include <errno.h>    // For errno (used in Mk_Dir)

#include "io_utils.h"
#include "config.h"   // Most már tartalmazza a PARTICLE_NUMBER, AU2CM, filenev2 definíciókat
#include "dust_physics.h" // Add this if not already present



/*	A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut	*/
void find_r_annulus(double *rvec, double rin_val, int *ind_ii, int *ind_io, double rout_val, int *ind_oi, int *ind_oo) {

    int i;
    double rmid, rtemp;
    double roimH;
    double roipH;
    double roomH;
    double roopH;
    double riimH;
    double riipH;
    double riomH;
    double riopH;

    if(optdze == 0) { // optdze should be declared or extern somewhere (e.g., config.h)

        *ind_ii = 0;
        *ind_io = 0;

    }

    // Corrected: Use rin_val and rout_val instead of rin and rout
    riimH = (rin_val - scale_height(rin_val)) - DD / 2.0;    /*	A nyomasi maximum az rout pontban van, ettol rout - 1/2H - DD / 2 es rout + 1*2H -DD / 2 kozott van a korgyuru belso hatara (azert DD/2, hogy biztosan 1 cellat tudjunk kijelolni, ne pedig egy tartomanyt)	*/
    riipH = (rin_val - scale_height(rin_val)) + DD / 2.0;
    riomH = (rin_val + scale_height(rin_val)) - DD / 2.0;    /*	Az alabbi ketto pedig a kulso hatarat adja meg a korgyurunek	*/
    riopH = (rin_val + scale_height(rin_val)) + DD / 2.0;

    roimH = (rout_val - scale_height(rout_val)) - DD / 2.0;    /*	A nyomasi maximum az rout pontban van, ettol rout - 1/2H - DD / 2 es rout + 1*2H -DD / 2 kozott van a korgyuru belso hatara (azert DD/2, hogy biztosan 1 cellat tudjunk kijelolni, ne pedig egy tartomanyt)	*/
    roipH = (rout_val - scale_height(rout_val)) + DD / 2.0;
    roomH = (rout_val + scale_height(rout_val)) - DD / 2.0;    /*	Az alabbi ketto pedig a kulso hatarat adja meg a korgyurunek	*/
    roopH = (rout_val + scale_height(rout_val)) + DD / 2.0;

    for(i = 1; i <= NGRID; i++) { // NGRID and RMIN and DD should be declared or extern (e.g., config.h)

        if(optdze == 1) {
/*	Ha az r tavolsag a kijelolt hatarok kozott van, akkor az adott valtozo visszakapja r erteket	*/
            if(rvec[i] > riimH && rvec[i] < riipH) {
                 rmid = (rvec[i] - RMIN)/ DD;
                rtemp = (int) floor(rmid + 0.5);
                *ind_ii = rtemp;
            }

/*	Ha az r tavolsag a kijelolt hatarok kozott van, akkor az adott valtozo visszakapja r erteket	*/
            if(rvec[i] > riomH && rvec[i] < riopH) {
                 rmid = (rvec[i] - RMIN)/ DD;
                rtemp = (int) floor(rmid + 0.5);
                *ind_io = rtemp;
            }
        }


/*	Ha az r tavolsag a kijelolt hatarok kozott van, akkor az adott valtozo visszakapja r erteket	*/
        if(rvec[i] > roimH && rvec[i] < roipH) {
             rmid = (rvec[i] - RMIN)/ DD;
            rtemp = (int) floor(rmid + 0.5);
            *ind_oi = rtemp;
        }

/*	Ha az r tavolsag a kijelolt hatarok kozott van, akkor az adott valtozo visszakapja r erteket	*/
        if(rvec[i] > roomH && rvec[i] < roopH) {
             rmid = (rvec[i] - RMIN)/ DD;
            rtemp = (int) floor(rmid + 0.5);
            *ind_oo = rtemp;
        }

        if(rvec[i] > roopH) break;

    }

}


// --- FÜGGVÉNY DEFINÍCIÓK ---

// reszecskek_szama függvény implementáció
int reszecskek_szama(int numout, const char *filenev){

	char c;
	fin1 = fopen(filenev,"r+");		
	numout = 0;

/*	A porreszecskeket tartalmazo file megnyitasa es a sorok szamanak kiolvasasa while ciklussal				*/

		while((c = fgetc(fin1)) != EOF)				
/*	a file vegeig (EOF) keresse a c karakternek megadott '\n' sortorest: 							*/
			if(c == '\n')
				numout++;					
/*	amig talal sortorest, leptesse a lines integert, ezzel beolvastuk, hogy hany soros a file				*/
	fclose(fin1);	

	return numout;

}



/*	A porreszecskek adatainak beolvasasa	*/
void por_be(double radius[][2], double radiusmicr[][2], double *mass, double *massmicr) {

	int i, dummy;
	double distance, particle_radius, reprmass, reprmassmicr,radmicr;
	
   	fin1 = fopen(filenev1,"r");
 
/*	Beolvassa a file-ból a részecskék adatait: sorszámukat - ezt később nem használjuk; távolságukat; sugaruk méretét; a reprezentatív tömegüket - egyelőre ezt sem használjuk	*/  	
	for (i = 0; i < PARTICLE_NUMBER; i++) {			
            	if(fscanf(fin1,"%d %lg %lg %lg %lg %lg",&dummy,&distance,&reprmass,&reprmassmicr,&particle_radius,&radmicr) == 6) {	
/*	A beolvasás sikeres, ha az fscanf visszatérési értéke 6, mert 6 oszlopot szeretnénk beolvasni. Ekkor elmentjük a részecske távolságát (distance) és méretét (particle_radius) a megfelelő tömbbe	*/

/*	A cm-es porreszecskek adatainak beolvasasa!		*/
           		radius[i][0] = distance;			/*	a reszecske tavolsaga AU-ban		*/
	   		radius[i][1] = particle_radius / AU2CM;		/* 	a részecske mérete AU-ban		*/
			mass[i] = reprmass;				/*	a porreszecske altal kepviselt reprezentativ tomeg dimenziotlan egysegekben					*/


/*	A mikronos reszecskek adatainak beolvasasa!		*/
            		radiusmicr[i][0] = distance;			/*	a mikronos reszecske tavolsaga AU-ban --> kezdetben ugyanolyan messze van az osszes, mint a cm-es reszecskek!!	*/
	   		radiusmicr[i][1] = radmicr / AU2CM;		/* 	a mikronos részecske mérete AU-ban	*/
			massmicr[i] = reprmassmicr;			/*	a porreszecske altal kepviselt reprezentativ tomeg dimenziotlan egysegekben					*/

	    	} else {

/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/					
			printf("\n\n*******************     ERROR!     *********************\n\n  Nem sikerult a %i-ik sort beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n",dummy);
//			getchar();
			exit(EXIT_FAILURE);
   	        }
	}

	fclose(fin1);
	
	printf("\n\n *******   A file beolvasasa sikerult!   ******* \n ******* Uss egy ENTER-t a folytatashoz! ******* \n\n ");	

}


/*	A sigmat tartalmazo file parametereinek beolvasasa	*/
void sigIn(double *sigvec, double *rvec) {

	double sig,r;
	int i;

	FILE *densin;
	densin = fopen(inputsig,"r");
	
	for(i = 0; i < NGRID; i++) {
           	if(fscanf(densin,"%lg  %lg",&r,&sig) == 2) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 3, mert 3 oszlopot szeretnénk beolvasni.	*/
			rvec[i+1] = r;				/*	r vektor	*/
			sigvec[i+1] = sig; 			/*	adott lepeskozonkent irja ki a file-t, megvan, hogy hany evente, tudjuk, hogy meddig fusson a kod, igy egy egyszeru osztassal meg lehet adni, hogy az mindig hanyadik idolepes	*/
			
	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
			exit(EXIT_FAILURE);
   	        }
	}

	RMIN = rvec[1];
	RMAX = rvec[NGRID];

	fclose(densin);
	
}


/*	Fuggveny az adott futashoz mappa letrehozasara, igy egy adott futas mindig kulon mappaba kerul es nem kavarodnak ossze az adatok	*/
void Mk_Dir(char *nev) {

	int i, step, mappa;
	char comm[2048];

	step = 0;

/*	A kimeneti file-ok szamara egy tarolo mappa letrehozasa (system(mkdir ...)) paranccsal	*/
	mappa=system("mkdir output");

/*	A ciklus ellenorzi, hogy letezik-e mar ez adott mappa, ha nem, letrehozza output neven...	*/
		if (mappa == 0) {
			printf("... A kimeneti file-ok tarolasahoz az \"output\" mappa elkeszult ...\n\n\n");
					sprintf(nev,"output");				/*	A mappa nevenek eltarolasa azert, hogy a file-okat az eppen aktualis mappaba tudja kesobb kiirni a program	    */

/*	...ha letezik az output mappa, akkor egy do-while ciklussal szamozott mappat hoz letre (pl. output.0), es addig fut a ciklus, mig aztan nem talalja mar a soron kovetkezo szamozott mappat. Ekkor letre hozza azt	*/
		} else {
			printf("... A kimeneti file-ok tarolasahoz az \"output\" mappa letrehozasa meghiusult ...\n");
			printf("... A mappa valoszinuleg mar letezik, uj mappa letrehozasa ...\n\n\n");

			do{
				for(i=0;i<=step;i++){
					sprintf(nev,"output.%i",i);			/*	A mappa nevenek eltarolasa azert, hogy a file-okat az eppen aktualis mappaba tudja kesobb kiirni a program	    */
					sprintf(comm,"mkdir output.%i",i);	
				}

				mappa=system(comm);
				step++;					

			} while (mappa!=0);

			printf("... A kimeneti file-ok tarolasahoz a(z) %s mappa elkeszult ...\n\n\n",nev);
		}

}

/*	Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok	*/
void infoCurrent(char *nev) {

	char out[1024];

/*	A TCURR tartalmazza a a szimulacio inditasanak pillanataban az idot: honap, nap, ora, perc, mp formatumban	 */
	sprintf(out,"run_%i.dat",(int)TCURR);
	jelfut = fopen(out,"w");
	
/*	Abba a mappaba hoz letre egy file-t, ahonnan inditottuk a programot. A file-ban leolvashatjuk, hogy a szimulacio kimenete mely mappaban talalhato, illetve a szimulaciorol nehany informacio -- ezt ki fogom meg egeszitani	*/
	fprintf(jelfut,"A jelenlegi futás a %s mappaban taláható!\n",nev);
	fprintf(jelfut,"\n\nA korong paraméterei:\nRMIN: %lg, RMAX: %lg\nSIGMA0: %lg, SIGMA_EXP: %lg, flaring index: %lg\nALPHA_VISC: %lg, ALPHA_MOD: %lg\nR_DZE_I: %lg, R_DZE_O: %lg, DR_DZEI: %lg, DR_DZE_O: %lg  (*** R_DZE_I/O = 0, akkor azt a DZE-t nem szimulálja a futás! ***)\n\n\n",RMIN,RMAX,SIGMA0,SIGMAP_EXP,FLIND,alpha_visc,a_mod,r_dze_i,r_dze_o,Dr_dze_i,Dr_dze_o);
	fprintf(jelfut,"A központi csillag tömege: %lg M_Sun\n",STAR);
	fclose(jelfut);

}



/*	Fuggveny a tomegfile kiiratasara	*/
void Print_Mass(double step, double *rvec, double partmassind[][4], double partmassmicrind[][4], double partmasssecind[][4], double t, double *dpressvec, double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, double *tavin, double *tavout) {

	double ind_ii, ind_io, ind_oi, ind_oo, tav, tav2;	

	tav = r_dze_o;	
	tav2 = r_dze_i;	

	int dim = find_num_zero(dpressvec);			// megnezi, hogy hany, nyomasi maximumbol szarmazo, 0 pont van a derivaltban
	double r_count[dim];					// a nullpontnak megfelelo elemu tombot letrehozza (pl. ha van dze_i es dze_o is, akkor kulon igy lehet elmenteni azok helyet)
	double temp_new = 0.;
	double temp = 0.;
	double rout = r_dze_o;
	double rin = r_dze_i;
	double rin_new = 0.0;
	double rout_new = 0.0;

	int j, i;
	j=0;
	i=0;
	
	if(dim != 0) {						// ha van nullpont, akkor megkeresi, hogy hol
		for(i = 0; i < NGRID; i++) {
			temp_new = find_zero(i,rvec,dpressvec);	// ha van nullpont, akkor a temp_new valtozoba tarolja el --> ehhez vegig megy az egesz r-en 

			if(temp != temp_new && i > 3 && temp_new != 0.0) {
				r_count[j] = temp_new;		// a temp_new-t itt tarolja el, azaz a nullpontok szamanak megfeleloen itt tartolodnak el a nyomasi maximumok
				j++;
			}

/*	Ha csak kulso dze van, akkor a dze uj helyet itt tarolja el	*/
			if(optdze == 0) {
				if(temp_new > 0.) {			
					temp = temp_new;
					rout_new = temp;
				} 
			}
		}
	}


/*	Ha van belso dze is, akkor itt menti el a kuslo es belso dze helyet	*/		
	if(optdze == 1) {
		if(dim > 0) {
			if (dim == 1) {
				rin_new = r_count[0];
				rout_new = rout;
			} else {
				rin_new = r_count[0];
				rout_new = r_count[1];
			}
		} 
		if(dim == 0) {	// ha nincs nyomasi maximum meg, akkor a regi valtozok erteket (azaz a nyomasi dze_i es dze_o helyeket) menti el
			rin_new = rin;
			rout_new = rout;
		}
	}

	rin = rin_new;
	if(optdze == 0) rin = 0;	// ha nincs belso dze, akkor annak a helye 0 (ez vegulis ebben az esetben nem lenyeges)
	rout = rout_new;
	tav2 = rin;
	tav = rout;

/* 	MEG KELL OLDANI, HOGY AKKOR IS TUDJON TOMEGNOVEKEDEST SZAMOLNI, HA CSAK BELSO DZE VAN!	*/

	find_r_annulus(rvec,tav2,&ind_ii,&ind_io,tav,&ind_oi,&ind_oo);		/*	A belso es kuslo nyomasi maximum korul 2H tavolsagban keres korgyurut, a fuggveny visszaadja a cellak indexet	*/
	
	double masst0i = 0, massii = 0, massoi = 0;
	double masst0im = 0, massiim = 0,massoim = 0;
	double massis = 0, massos = 0;

	GetMass(PARTICLE_NUMBER,partmassind,(int)ind_ii,(int)ind_io,tav2,r_dze_i,&massii,(int)ind_oi,(int)ind_oo,tav,r_dze_o,&massoi);
	
	if(opttwopop == 1) {
		GetMass(4*PARTICLE_NUMBER,partmasssecind,(int)ind_ii,(int)ind_io,tav2,r_dze_i,&massis,(int)ind_oi,(int)ind_oo,tav,r_dze_o,&massos);
		GetMass(PARTICLE_NUMBER,partmassmicrind,(int)ind_ii,(int)ind_io,tav2,r_dze_i,&massiim,(int)ind_oi,(int)ind_oo,tav,r_dze_o,&massoim);
	} 

	double massi, massim, masso, massom;

	if(tav2 != r_dze_i) {
		massi = massii + massbtempii + massis;
		massim = massiim + massmtempii;
	} else {
		massi = massii + massis;
		massim = massiim;
	}
	if(tav != r_dze_o) {
		masso = massoi + massbtempoi + massos;
		massom = massoim + massmtempoi;
	} else {
		masso = massoi + massos;
		massom = massoim;
	}

	*massbtempio = massi;
	*massbtempoo = masso;
	*massmtempio = massim;
	*massmtempoo = massom;

	*tavin = tav2;
	*tavout = tav;

	fprintf(massfil,"%lg %lg %lg %lg %lg\n",step,tav2,massi+massim,tav,masso+massom);
	fflush(massfil);

}

/*	Fuggveny a sigma, p, dp kiiratasara	*/
void Print_Sigma(char *dens_name, double *rvec, double *sigmavec, double *pressvec, double *dpressvec) {

	int i;
	fmo = fopen(dens_name,"w");				

 	for(i = 1; i <= NGRID; i++) {
   		fprintf(fmo, "%lg   %lg   %lg   %lg\n", rvec[i],sigmavec[i],pressvec[i],dpressvec[i]);
	}

	fclose(fmo);

}


/*	Fuggveny a por feluletisurusegenek kiiratasara	*/
void Print_Sigmad(char *dust_name, char *dust_name2, double min, double *r, double *rm, double *sigmad, double *sigmadm) {

	int i;
	double rtempvec[PARTICLE_NUMBER];

	FILE *sid;	

	fil = fopen(dust_name,"w");
	if(opttwopop == 1) sid = fopen(dust_name2,"w");		/*	Ha 2pop --> megnyit egy kulon file-t a mikronos por feluletisurusegenek kiiratasara	*/

	for(i=0;i<PARTICLE_NUMBER;i++){

		if (r[i] >= RMIN) {			/*	a cm-es por feluletisurusege	*/
			fprintf(fil,"%.11lg  %lg \n",r[i],sigmad[i]);
		} 

		if(opttwopop == 1) {				/*	itt irja ki a mikronos por feluletisuruseget	*/

			if (rm[i] >= RMIN) {
				fprintf(sid,"%lg  %lg \n",rm[i],sigmadm[i]);
			}
		}
	}

	fclose(fil);
	if(opttwopop == 1) fclose(sid);

}


/*	Fuggveny a pormozgas kiiratasara	*/
void Print_Pormozg_Size(char *size_name, int step, double rad[][2], double radmicr[][2], double *rvec, double t){

	int i;
	if(optgr == 1) fout3 = fopen(size_name,"w");	/*	ha van pornovekedes	*/

	for(i=0;i<PARTICLE_NUMBER;i++){

/*	pormozgas.dat kiiratasa --> adott idokozonkent a cm-es reszecskek tavolsaga, indexe es az ido kiiaratasa egy file-ba	*/
		if(rad[i][0] >= RMIN) fprintf(fout,"%lg %d %lg\n",(double)step,i,rad[i][0]);
/*	ha a szimulacio 2 populacios, akkor a fentihez hasonlo file letrehozasa a mikronos reszecskekre	*/		
		if(opttwopop == 1.) if(radmicr[i][0] >= RMIN) fprintf(foutmicr,"%lg %d %lg\n",(double)step,i,radmicr[i][0]);
/*	ha van reszecske novekedes, akkor letrehoz egy file-t minden adott idolepesben, ami a centimeteres porreszecskek meretet tartalmazza	*/
		if(optgr == 1) {
			if (rad[i][0] >= RMIN) fprintf(fout3,"%lg  %lg  %lg \n",(double)step, rad[i][0], rad[i][1]*AU2CM);
		}
	}

	fflush(fout);
	if(opttwopop == 1.) fflush(foutmicr);
	if(optgr == 1) fclose(fout3);

}




/*	Az idot tartalmazo file parametereinek beolvasasa	*/
void timePar(double *tMax, double *step, double *current) {

	double tmax,stepping,curr;

	fin2 = fopen(filenev3,"r");

           	if(fscanf(fin2,"%lg  %lg %lg",&tmax,&stepping,&curr) == 3) { /*	A beolvasás sikeres, ha az fscanf visszatérési értéke 3, mert 3 oszlopot szeretnénk beolvasni.	*/
 			printf("\n\n *********** A korong parameterei sikeresen beolvasva!  *********** \n                  tmax: %lg, a program %lg evenkent irja ki a file-okat\n\n\n",tmax,stepping);  
			*tMax = tmax;				/*	meddig fusson a program	*/
			stepping = tmax/stepping;
			*step = stepping; 			/*	adott lepeskozonkent irja ki a file-t, megvan, hogy hany evente, tudjuk, hogy meddig fusson a kod, igy egy egyszeru osztassal meg lehet adni, hogy az mindig hanyadik idolepes	*/
			*current = curr;			/*	beolvassa, hogy mennyi volt az ido eppen a futas inditasanak pillanataban -- ez a kiiratasnal lenyeges	*/
			
	    	} else {					/*	Ha a beolvasás valamiért nem sikeres, akkor figyelmeztetés mellett a program kilép (megmondja, hogy melyik sort nem tudta kezelni)	*/
			printf("\n\n*******************     ERROR!     *********************\n\n  A file-t nem sikertult beolvasni, a program kilep!\n \t  A kilepeshez nyomj egy ENTER-t!\n");
			exit(EXIT_FAILURE);
   	        }

	fclose(fin2);
	
}
