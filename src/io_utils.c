#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> // For mkdir
#include <unistd.h>   // For access (used in Mk_Dir)
#include <errno.h>    // For errno (used in Mk_Dir)

#include "io_utils.h"
#include "config.h"   // Most már tartalmazza a PARTICLE_NUMBER, AU2CM, filenev2 definíciókat
#include "dust_physics.h" // Add this if not already present
#include "utils.h"




// --- FÜGGVÉNY DEFINÍCIÓK ---

/*	Visszaadja, hogy hany sora van a beolvasando file-nak, ez jelen esetben megadja a beolvasando reszecskek szamat!	*/
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


//*	A porreszecskek adatainak beolvasasa	*/
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
void sigIn(double *sigma_arr, double *r_arr, const disk_t *disk_params, const char *filename) {
    // ideiglenesen a globális filenev1-et használjuk, amíg azt is refaktoráljuk
    const char *input_filename = filenev1; // filenev1 a config.h-ból jön

    FILE *fp = fopen(input_filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR [sigIn]: Could not open input file '%s'.\n", input_filename);
        perror("Reason"); // Kiírja a rendszerhiba üzenetet
        exit(EXIT_FAILURE);
    }

    // Most jön a beolvasás logikája. Ügyelj a formátumra és az indexelésre!
    // Az init_tool kimenete: i, r, prad, reppmass
    int i_read;
    double r_read, prad_read, reppmass_read;

    // Feltételezve, hogy a fájl elején van egy fejlécsor, amit át kell ugrani, ha van
    // (ha nincs, akkor ezt a sort ki kell hagyni)
    char line[256];
    if (fgets(line, sizeof(line), fp) == NULL) {
        fprintf(stderr, "ERROR [sigIn]: Failed to read header line or empty file: %s\n", input_filename);
        fclose(fp);
        exit(EXIT_FAILURE);
    }
    // Ha van Press ENTER to continue! sor:
    if (fgets(line, sizeof(line), fp) == NULL) {
        fprintf(stderr, "ERROR [sigIn]: Failed to read second header line or empty file: %s\n", input_filename);
        fclose(fp);
        exit(EXIT_FAILURE);
    }
    // ... és így tovább, ahány fejléc sor van.

    // Valószínűleg 1-től NGRID-ig indexeljük, de a tömb NGRID+2 méretű
    // és a 0 és NGRID+1 indexek a peremfeltételekre vannak fenntartva.
    // Ezt ellenőrizd!
    for (int i = 1; i <= disk_params->NGRID; i++) { // Feltételezve, hogy a belső pontok NGRID számúak
        if (fscanf(fp, "%d %lf %lf %lf", &i_read, &r_read, &prad_read, &reppmass_read) != 4) {
            fprintf(stderr, "ERROR [sigIn]: Failed to read data for particle %d from file '%s'.\n", i, input_filename);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
        // Ellenőrizd, hogy az r_arr és sigma_arr indexelése helyes-e
        // Például ha az i_read 0-tól indul, de te 1-től tárolsz:
        r_arr[i] = r_read; // Vagy r_arr[i_read+1] ha i_read 0-tól van
        sigma_arr[i] = reppmass_read; // Vagy sigma_arr[i_read+1]
        // Lehet, hogy a sigma_arr-ba valamilyen gáz sűrűség érték kell, nem a reppmass_read?
        // Ezt tisztázni kell az init_tool és sigIn közötti konvenciókból.
    }

    fclose(fp);
    printf("DEBUG [sigIn]: Successfully loaded profile from %s.\n", input_filename);
}



// Change this line (the function signature)
void Mk_Dir(char *nev) { // 'nev' is now a modifiable buffer passed in

    int i, step, mappa;
    char comm[2048];

    step = 0;

    /*	A kimeneti file-ok szamara egy tarolo mappa letrehozasa (system(mkdir ...)) paranccsal	*/
    mappa=system("mkdir output");

    /*	A ciklus ellenorzi, hogy letezik-e mar ez adott mappa, ha nem, letrehozza output neven...	*/
    if (mappa == 0) {
        printf("... A kimeneti file-ok tarolasahoz az \"output\" mappa elkeszult ...\n\n\n");
        sprintf(nev,"output"); // Now 'nev' is a modifiable buffer, so this is safe

    /*	...ha letezik az output mappa, akkor egy do-while ciklussal szamozott mappat hoz letre (pl. output.0), es addig fut a ciklus, mig aztan nem talalja mar a soron kovetkezo szamozott mappat. Ekkor letre hozza azt	*/
    } else {
        printf("... A kimeneti file-ok tarolasahoz az \"output\" mappa letrehozasa meghiusult ...\n");
        printf("... A mappa valoszinuleg mar letezik, uj mappa letrehozasa ...\n\n\n");

        do{
            // The loop structure here is a bit unusual.
            // It seems like you intend to try "output.0", "output.1", etc.
            // However, the `for(i=0;i<=step;i++){ sprintf(nev,"output.%i",i); sprintf(comm,"mkdir output.%i",i); }`
            // will only process the LAST value of `i` (which is `step`) in each iteration of the `do-while`.
            // It should probably just be `sprintf(nev, "output.%i", step);` before the `mkdir` call.
            // Let's fix the logic while fixing the type.

            sprintf(nev, "output.%i", step); // Prepare the name for the current step
            sprintf(comm, "mkdir %s", nev);  // Prepare the mkdir command

            mappa=system(comm);
            step++;

        } while (mappa!=0);

        printf("... A kimeneti file-ok tarolasahoz a(z) %s mappa elkeszult ...\n\n\n",nev);
    }

}

/*	Elkeszit egy file-t, ami tartalmazza a jelenlegi futas parametereit, es hogy melyik mappaban talalhatoak a file-ok	*/
void infoCurrent(const char *nev) {

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
void Print_Mass(double step, double *rvec, double partmassind[][4], double partmassmicrind[][4], double partmasssecind[][4], double *dpressvec, double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, double *tavin, double *tavout) {

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
	
	double massii = 0, massoi = 0;
	double massiim = 0,massoim = 0;
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
void Print_Sigmad(char *dust_name, char *dust_name2, double *r, double *rm, double *sigmad, double *sigmadm) {

	int i;

	FILE *fil = NULL; // Érdemes ezt is inicializálni
	FILE *sid = NULL; // Ezt már inicializáltad, ez így jó

	fil = fopen(dust_name,"w");
    // MINDIG ellenőrizd az fopen sikerességét!
    if (fil == NULL) {
        perror("Hiba a dust_name fájl megnyitásakor írásra");
        return; // Kilépés a függvényből, ha nem tudja megnyitni az első fájlt
    }
	
	// Csak akkor próbáld megnyitni a sid fájlt, ha opttwopop == 1
	if(opttwopop == 1) {
		sid = fopen(dust_name2,"w");
        // Ha opttwopop == 1, de a sid fájlt mégsem tudja megnyitni, kezeld ezt!
        if (sid == NULL) {
            perror("Hiba a dust_name2 fájl megnyitásakor írásra (mikronos por)");
            fclose(fil); // Zárjuk be az első fájlt, ha a második nem megy
            return; // Kilépés a függvényből
        }
	}

	for(i=0;i<PARTICLE_NUMBER;i++){

		if (r[i] >= RMIN) {			/*	a cm-es por feluletisurusege	*/
			fprintf(fil,"%.11lg  %lg \n",r[i],sigmad[i]);
		}

		if(opttwopop == 1 && sid != NULL) { // Ide már bekerült a sid != NULL ellenőrzés is
                                            // Bár a fenti logikával ez már felesleges lehet,
                                            // de extra biztonságot ad.
			if (rm[i] >= RMIN) {
				fprintf(sid,"%lg  %lg \n",rm[i],sigmadm[i]);
			}
		}
	}

	fclose(fil);
	// Csak akkor zárd be a sid fájlt, ha sikeresen meg lett nyitva (azaz nem NULL)
	if(opttwopop == 1 && sid != NULL) { // Fontos ellenőrizni, hogy meg van-e nyitva!
        fclose(sid);
    }
}


/*	Fuggveny a pormozgas kiiratasara	*/
void Print_Pormozg_Size(char *size_name, int step, double rad[][2], double radmicr[][2]){

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
// A függvény definícióját megváltoztatjuk, hogy ne FILE* mutatót várjon, hanem közvetlenül az értékeket.
// Ezeket az értékeket majd a main.c fogja átadni a YAML-ból származó adatok alapján.
void timePar(double tMax_val, double stepping_val, double current_val) { // A paramétereket érték szerint kapja
    // Nincs szükség filenev3-ra és fin2-re, eltávolítjuk a fájlbeolvasást.
    // Nincs szükség fscanf-re sem.

    // A printf és a változók beállítása megmarad, de most a kapott értékeket használjuk.
    printf("\n\n *********** A korong parameterei sikeresen beolvasva!  *********** \n        tmax: %lg, a program %lg evenkent irja ki a file-okat\n\n\n", tMax_val, stepping_val);

    // Globális változók feltöltése (feltételezve, hogy TMAX, WO, TCURR globálisak és ezek a tényleges változónevek)
    // A C kódodból vettem át őket, kérlek ellenőrizd, hogy ezek a helyes globális változónevek
    // és hogy elérhetők az io_utils.c-ből (pl. "config.h" include-olásával).
    extern double TMAX, WO, TCURR; // Deklaráljuk, hogy globális változókat használunk

    TMAX = tMax_val;
    // Az eredeti kódodban a 'stepping' azzal volt egyenlő, hogy tmax/stepping.
    // Most már a Python adja át a kívánt 'lépésköz' értéket, tehát nem kell újraosztani.
    // Vagy ha az a cél, hogy WO az 'output frequency' legyen, akkor az osztás marad:
    WO = tMax_val / stepping_val; // 'stepping_val' mostantól az output frequency a YAML-ból

    TCURR = current_val; // Az aktuális idő, amit majd a Python ad át

    // Nincs exit(EXIT_FAILURE), mert a fájlbeolvasási hibát kivettük.
    // Nincs fclose(fin2), mert nincs fin2.
}