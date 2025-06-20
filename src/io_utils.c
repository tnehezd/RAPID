#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> // For mkdir
#include <unistd.h>   // For access (used in Mk_Dir)
#include <errno.h>    // For errno (used in Mk_Dir)

#include "io_utils.h"
#include "config.h"   // Most már tartalmazza a PARTICLE_NUMBER, AU2CM, filenev2 definíciókat
#include "dust_physics.h" // Add this if not already present


// --- IDEIGLENES FORWARD DEKLARÁCIÓK ---
// Ezekre azért van szükség, mert az io_utils.c-ben lévő függvények hívnak olyan függvényeket,
// amelyek definíciói még más modulokban (vagy még az eredeti nagy fájlban) vannak.
// Amikor a többi modult is létrehozzuk és include-oljuk a header fájljaikat, ezeket törölni fogjuk.
int find_num_zero(double *rvec, double *dpressvec);
double find_zero(int i, double *rvec, double *dpressvec);

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



// por_be függvény implementáció
void por_be() {
    FILE *fin1; // Lokális fájlpointer
    int i; // Lokális ciklusváltozó

    fin1 = fopen(filenev1, "r"); // filenev1 globális a config.h-ból
    if (fin1 == NULL) {
        fprintf(stderr, "Error: Could not open file %s for por_be.\n", filenev1);
        exit(EXIT_FAILURE); // Fontos, hogy kilépjünk, ha kritikus fájlt nem tudunk megnyitni
    }

    // Itt kellene lennie az eredeti kódnak, ami beolvassa a por adatokat.
    // Feltételezzük, hogy a 'radius' globális tömb, vagy valahogy máshogy elérhető.
    // Ha nem, akkor paraméterként kellene kapnia por_be-nek.
    double particle_radius; // Lokális változó

    for (i = 0; i < PARTICLE_NUMBER; i++) { // PARTICLE_NUMBER globális a config.h-ból
        // Példa: ha az eredeti kód valahogy így olvasta be:
        // fscanf(fin1, "%lg %lg", &valami_radius_r, &particle_radius);
        // ... majd használná AU2CM-et:
        // valami_radius_au = particle_radius / AU2CM;
    }

    fclose(fin1);
}

// sigIn függvény implementáció (az 'inputsig' problémájának javításával)
// Most feltételezi, hogy a 'filenev2' (globális config.h-ból) a bemeneti szigma fájl neve.
void sigIn(double sigmavec[], double rvec[]) {
    FILE *densin; // Lokális fájlpointer
    int i; // Lokális ciklusváltozó

    densin = fopen(filenev2, "r"); // filenev2 globális a config.h-ból
    if (densin == NULL) {
        fprintf(stderr, "Error: Could not open file %s for sigIn.\n", filenev2);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < NGRID + 2; i++) { // NGRID globális a config.h-ból
        // Itt kellene lennie az eredeti kódnak, ami beolvassa a szigma adatokat.
        // fscanf(densin, "%lg %lg", &rvec[i], &sigmavec[i]);
    }
    fclose(densin);
}

// Mk_Dir függvény implementáció
void Mk_Dir(char *nev) {
    struct stat st = {0};
    if (stat(nev, &st) == -1) { // Ellenőrzi, hogy létezik-e a mappa
        if (mkdir(nev, 0777) == -1) { // Létrehozza, ha nem létezik
            perror("Error creating directory");
            exit(EXIT_FAILURE);
        }
    }
}

// infoCurrent függvény implementáció
void infoCurrent(char *nev) {
    char filename[1024];
    snprintf(filename, sizeof(filename), "%s/info.dat", nev);
    
    // jelfut globális fájlpointer a config.h-ból
    jelfut = fopen(filename, "w");
    if (jelfut == NULL) {
        fprintf(stderr, "Error: Could not open info file %s for writing.\n", filename);
        // Nem lépünk ki, csak kiírjuk a hibát és folytatjuk a program futását,
        // de jelfut NULL marad, így a további írási kísérletek hibát fognak okozni.
        return;
    }

    // Itt jön a rengeteg fprintf, ami a globális változókat használja a config.h-ból
    fprintf(jelfut, "RMIN: %lg\n", RMIN);
    fprintf(jelfut, "RMAX: %lg\n", RMAX);
    fprintf(jelfut, "NGRID: %d\n", NGRID);
    fprintf(jelfut, "DD: %lg\n", DD);
    fprintf(jelfut, "SIGMA0: %lg\n", SIGMA0);
    fprintf(jelfut, "SIGMAP_EXP: %lg\n", SIGMAP_EXP);
    fprintf(jelfut, "FLIND: %lg\n", FLIND);
    fprintf(jelfut, "alpha_visc: %lg\n", alpha_visc);
    fprintf(jelfut, "a_mod: %lg\n", a_mod);
    fprintf(jelfut, "STAR: %lg\n", STAR);
    fprintf(jelfut, "PDENSITY: %lg\n", PDENSITY);
    fprintf(jelfut, "PDENSITYDIMLESS: %lg\n", PDENSITYDIMLESS);
    fprintf(jelfut, "r_dze_i: %lg\n", r_dze_i);
    fprintf(jelfut, "r_dze_o: %lg\n", r_dze_o);
    fprintf(jelfut, "Dr_dze_i: %lg\n", Dr_dze_i);
    fprintf(jelfut, "Dr_dze_o: %lg\n", Dr_dze_o);
    fprintf(jelfut, "optdze: %d\n", optdze);
    fprintf(jelfut, "optev: %d\n", optev);
    fprintf(jelfut, "optdr: %d\n", optdr);
    fprintf(jelfut, "optgr: %d\n", optgr);
    fprintf(jelfut, "opttwopop: %d\n", opttwopop);
    fprintf(jelfut, "fFrag: %d\n", fFrag);
    fprintf(jelfut, "uFrag: %d\n", uFrag);
    fprintf(jelfut, "inputsig: %s\n", inputsig);
    fprintf(jelfut, "DT: %lg\n", DT);
    fprintf(jelfut, "TMAX: %lg\n", TMAX);
    fprintf(jelfut, "WO: %lg\n", WO);
    fprintf(jelfut, "TCURR: %lg\n", TCURR);
    fprintf(jelfut, "PARTICLE_NUMBER: %d\n", PARTICLE_NUMBER); // Globális konstans
    fprintf(jelfut, "AU2CM: %lg\n", AU2CM);                     // Globális konstans


    fclose(jelfut);
    jelfut = NULL; // Fontos, hogy NULL-ra állítsuk, miután bezártuk
}

/*	Fuggveny a tomegfile kiiratasara	*/
void Print_Mass(double step, double *rvec, double partmassind[][4], double partmassmicrind[][4], double partmasssecind[][4], double t, double *dpressvec, double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, double *tavin, double *tavout) {

	int ind_ii, ind_io, ind_oi, ind_oo, tav, tav2;	

	tav = r_dze_o;	
	tav2 = r_dze_i;	

	int dim = find_num_zero(rvec,dpressvec);		// megnezi, hogy hany, nyomasi maximumbol szarmazo, 0 pont van a derivaltban
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

// Print_Sigma függvény implementáció
void Print_Sigma(char *filename, double rvec[], double sigmavec[], double pressvec[], double dpressvec[]) {
    FILE *fout_local; // Lokális fájlpointer
    int i; // Lokális ciklusváltozó

    fout_local = fopen(filename, "w");
    if (fout_local == NULL) {
        fprintf(stderr, "Error: Could not open sigma file %s for writing.\n", filename);
        return;
    }

    for (i = 0; i < NGRID + 2; i++) { // NGRID globális a config.h-ból
        fprintf(fout_local, "%lg %lg %lg %lg\n", rvec[i], sigmavec[i], pressvec[i], dpressvec[i]);
    }
    fclose(fout_local);
}

/*	Fuggveny a por feluletisurusegenek kiiratasara	*/
void Print_Sigmad(char *dust_name, char *dust_name2, double min, double *r, double *rm, double *sigmad, double *sigmadm) {

	int i;
	double rtempvec[PARTICLE_NUMBER];
	double dd = (RMAX - RMIN) / (PARTICLE_NUMBER-1);	/*	Mivel a jelen futas gridfelbontasa nem feltetlen egyezik meg a porreszecskeket generalo program gridfelbontasaval - ez a feluletisuruseg miatt lenyeges! - ezert itt szamolja ki a program	*/

//	printf("dd: %lg  1/dd: %i\n",dd,(int)(1./dd)*100);

	FILE *sid = NULL;	

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


// Print_Pormozg_Size függvény implementáció
void Print_Pormozg_Size(char *size_name, int step, double rad[][2], double radmicr[][2], double *rvec, double t){
    // Itt is feltételezzük, hogy az 'fout3' globális fájlpointert használjuk a config.h-ból.
    // Ha nem, akkor deklarálni kell itt: FILE *fout3; és meg kell nyitni.
    int i; // Lokális ciklusváltozó

    // Ha az fout3 globális és nem itt nyitódik, akkor feltételezzük, hogy már nyitva van.
    // Ha itt kellene nyitni:
    // fout3 = fopen(size_name, "w");
    // if (fout3 == NULL) { /* handle error */ return; }

    for(i=0; i<PARTICLE_NUMBER; i++){ // PARTICLE_NUMBER globális a config.h-ból
        if (rad[i][0] >= RMIN) { // RMIN globális a config.h-ból
            // AU2CM globális a config.h-ból
            if (fout3 != NULL) { // Ellenőrizzük, hogy a fájlpointer érvényes
                fprintf(fout3,"%lg %lg %lg \n",(double)step, rad[i][0], rad[i][1]*AU2CM);
            } else {
                 fprintf(stderr, "Warning: fout3 is NULL in Print_Pormozg_Size. Cannot write.\n");
            }
        }
    }
    // Ha itt nyitottuk az fout3-at, akkor itt be is kell zárni: fclose(fout3); fout3 = NULL;
    // Ha globális és máshol záródik, akkor nincs itt dolgunk vele.
    // Az előző config.h-s definíciód szerint a fout3 globális, de NULL-ra van inicializálva,
    // tehát valahol a fő programban kellene megnyitni a használat előtt.
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
