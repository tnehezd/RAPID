#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h> // For TCURR

#include "config.h"       // For global variables
#include "io_utils.h"    // Will be included when we move IO functions
#include "disk_model.h"  // Will be included later
#include "dust_physics.h" // Will be included later
#include "simulation_core.h" // Will be included later
#include "utils.h"         // Will be included later

// Forward declaration of the options_t struct and parsing functions
// This will eventually go into a separate 'options' module if it gets complex enough,
// or remain here if the parsing logic is tightly coupled to main.
/*	A futasok soran valtoztathato opciok letrehozasa egy strukturaval	*/
typedef struct options {

	// Option for dust drift
	double drift;
	// Option for dust growth
	double growth;
	// Option for solving the duffusion equation of the surface density
	double evol;
	// Option fot two population model
	double twopop;
	// Fragmentation barrier in cm/s
	double ufrag;
	// Fragmentation barrier constant
	double ffrag;
	// Number of grid cells
	int ngrid;
	// Input sigma file
	const char *input; 
	// Option for time stepping
	double tStep;

} options_t;

void create_default_options(options_t *def);
int parse_options(int argc, const char **argv, options_t *def);

// --- The rest of your main function will go here ---
// We will paste the original content of your main() function here later,
// and replace direct function calls with calls to our new modules.

int main(int argc, const char **argv) {
   
	options_t def;				/*	A letrehozott struktura elemeire .def-el lehet majd hivatkozni	*/
	create_default_options(&def);		/*	A struktura elemeinek feltoltese alapertelmezett parameterekkel	*/
/*	A terminalbol keresi a kapcsolok segitsegevel a struktura elemeire vonatkozo adatokat. A beolvasas kimenetet eltarolja egy integerbe	*/
	int retCode = parse_options(argc, argv, &def);
/*	Ha a parameterek beolvasasa sikertelen volt, akkor a program kilep	*/	
	if (0 != retCode) {
		exit(retCode);
	}

/*	A globalis valtozok feltoltese a parameterek ertekevel	*/
	optev = def.evol;
	optdr = def.drift;
	optgr = def.growth;
	opttwopop = def.twopop;
	fFrag = def.ffrag;
	uFrag = def.ufrag;
	inputsig = def.input;
	DT = def.tStep;

	if(inputsig == 0) {

		int lout, nout;
		lout = 0, nout = 0;

/*	A sigmat tartalmazo file sorainak szama elmentve egy integerbe	*/
		nout = reszecskek_szama(lout,inputsig); 
		NGRID = nout;

	} else {

		NGRID = def.ngrid;

	}

   	double sigmavec[NGRID+2], rvec[NGRID+2], pressvec[NGRID+2], dpressvec[NGRID+2], ugvec[NGRID+2];
	char dens_name[1024], nev[1024], mv[1024];

/*	A korong parametereinek beolvassa				*/
	disk_param_be(&SIGMA0, &SIGMAP_EXP, &RMIN, &RMAX, &r_dze_i, &r_dze_o, &Dr_dze_i, &Dr_dze_o, &a_mod, &PDENSITY, &PDENSITYDIMLESS, &alpha_visc,&STAR,&FLIND);
	DD = (RMAX - RMIN) / (NGRID - 1);						/*	rácsfelbontás	*/

/*	Ha _van_ bemeneti sigma file	*/
	if(inputsig == 0) {

		sigIn(sigmavec,rvec);
		Perem(rvec);
		Perem(sigmavec);
		
	} else {

/*	Kezdeti profilok betoltese a megadott vektorokba	*/
		load_R(rvec);
		Initial_Profile(sigmavec,rvec);

	}

	Initial_Press(pressvec,sigmavec,rvec);
	Initial_dPress(dpressvec,pressvec);
	Initial_Ugas(sigmavec,rvec,ugvec);

	timePar(&TMAX,&WO,&TCURR);

/*	az optdze globalis valtozo ertekenek megadasa: ha van belso nyomasi maximum is, akkor az erteke 1, egyebkent nulla. Ennek a tomegnovekedes szamolasanak szempontjabol van lenyege	*/ 
	optdze = 0;
	if(r_dze_i != 0) optdze = 1.;

/*	Mappa letrehozasa az adatok eltarolasahoz	*/
	Mk_Dir(nev);								

/*	Az aktualis mappaba a kezdeti adatokat tartalmazo file-ok atmasolasa - kesobb ezek az adatok visszanezhetok igy		*/
	int dummy;
	snprintf(mv,"cp %s %s/",filenev2,nev);
	dummy = system(mv);	
	snprintf(mv,"cp %s %s/",filenev3,nev);
	dummy = system(mv);

/*	Sigma file beolvasasa, ha szukseges	*/
	if(inputsig == 0) {
 	   	snprintf(mv, sizeof(mv), "cp %s %s/", filenev2, nev); // Itt a filenev2-t használjuk
		dummy = system(mv);
	}

/*	Abban a mappaban, ahol a futast inditottuk, egy file letrehozasa, amely az aktualis futasrol infokat ir ki (hol talalhatoak a kimeneti file-ok, es milyen parameterei vannak pl. a korongnak az adott futas eseten	*/
	infoCurrent(nev);

/*	Ha nincs sigma fejlodes es drift, akkor a kezdeti profilt kiirja egy file-ba es kilep "figyelmeztetes" mellett	*/
	if(optev == 0. && optdr == 0.) {
		printf("A megadott opciok szerint nem szamol sem sigmat, sem driftet, ezert a progam kilep!\n\nA kezdeti file-ok a %s mappaban talalhatoak!\n",nev);
/*	t=0-ban kiirja a sigma-t, a nyomast es a nyomasderivaltjat	*/
		snprintf(dens_name,1024,"%s/surface.dat",nev);
		Print_Sigma(dens_name,rvec,sigmavec,pressvec,dpressvec);
	} else {
		tIntegrate(nev,rvec,sigmavec,pressvec,dpressvec,ugvec);
	}

	return 0;

}



/*	A struktura elemeinek feltoltese alapertelmezett ertekekkel	*/
void create_default_options(options_t *opt) {

	opt->drift		 = 1.;		//Dust drift is included
	opt->growth		 = 1.;		//Particle growth is included
	opt->evol		 = 1.;		//The evolution of surface density is included
	opt->twopop		 = 1.;		//Two population simulation is included
	opt->ufrag		 = 1000.0;	//Fragmentation velocity in CGS -- Birnstiel et al 2012
	opt->ffrag		 = 0.37;	//?
	opt->ngrid		 = 2000;	//Number of the grid cells
	opt->input		 = "";
	opt->tStep		 = 0.;
	
}


/*	A valtoztathato parameterek beolvasa terminalbol	*/
int parse_options(int argc, const char **argv, options_t *opt){
	int i = 1;
	double temp = 1;

	while (i < argc) {
		char *p = malloc(sizeof(argv[i]));	/*	Ahhoz, hogy "szoveget" (char, mert a C nem tud stringet kezelni) tudjunk beolvasni, es ossze tudjuk vetni a lentebb megadott kapcsolokkal (pl. -drift), le kell foglalni a memoriateruletet p-nek. Ennek viszont elore nem tudjuk a meretet, ezert dinamikusan foglaljuk le azt	*/
   		strcpy(p, argv[i]);			/*	A p-be elmentjuk az adott argumentumot az strcpy paranccsal	*/	

		if(strcmp(p, "-drift") == 0) {		/*	Az strcmp parancs kepes osszevetni a karaktersorozatunkat egy masikkal: ha azok megegyeznek, akkor a kovetkezo beolvasott argumentumot elmenti a megadott struktura elembe	*/
			strcpy(p, argv[i]);
			i++;
			opt->drift = atof(argv[i]);
		}
		else if (strcmp(p, "-growth") == 0) {
			i++;
			opt->growth = atof(argv[i]);
		}
		else if (strcmp(p, "-evol") == 0) {
			i++;
			opt->evol = atof(argv[i]);
		}
		else if (strcmp(p, "-twopop") == 0) {
			i++;
			opt->twopop = atof(argv[i]);
		}
		else if (strcmp(p, "-tStep") == 0) {
			i++;
			opt->tStep = atof(argv[i]);
		}
		else if (strcmp(p, "-n") == 0) {
			i++;
			opt->ngrid = atoi(argv[i]);
		}
		else if (strcmp(p, "-i") == 0) {
			i++;
			opt->input = argv[i];
			temp = 0;
		}
		else {				/*	ha nem sikerul a parancssorbol a beolvasas, akkor a program kilep es feldobja, hogy mely kapcsolo nem volt jo, helyette mit probaljon meg	*/
			printf("\n\n**** Invalid switch on command-line: %s! ****\n\n\n",p);
			printf("**** Try following parameters: ****\n\n-drift\nif set to 0: dust drift is not included, if set to 1: dust drift is included (default)\n\n-growth\nif set to 0: partcle growth is not included, if set to 1: particle growth is included (default)\n\n-evol\nif set to 0: the evolution of the surface density is not included, if set to 1: evolution of the surface density is included (default)\n\n-n\nnumber of grid cells (2000 by default)\n");
			return 1;
		}
		i++;
		free(p);			/*	memoria foglalas miatt fel is kell szabaditanunk a memoria teruletet	*/
	}

	optinp = temp;
	
	return 0;
}

