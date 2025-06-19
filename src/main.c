#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h> // For TCURR

#include "config.h"       // For global variables
#include "io_utils.h"    // Will be included when we move IO functions
// #include "disk_model.h"  // Will be included later
// #include "dust_physics.h" // Will be included later
// #include "simulation_core.h" // Will be included later
#include "utils.h"         // Will be included later

// Forward declaration of the options_t struct and parsing functions
// This will eventually go into a separate 'options' module if it gets complex enough,
// or remain here if the parsing logic is tightly coupled to main.
typedef struct {
    int evol;
    int drift;
    int growth;
    int twopop;
    int ffrag;
    int ufrag;
    int input;
    double tStep;
    int ngrid; // Assuming ngrid is part of options for optinp == 1
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

	if(optinp == 0) {

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
	if(optinp == 0) {

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
	sprintf(mv,"cp %s %s/",filenev2,nev);
	dummy = system(mv);	
	sprintf(mv,"cp %s %s/",filenev3,nev);
	dummy = system(mv);

/*	Sigma file beolvasasa, ha szukseges	*/
	if(optinp == 0) {
		sprintf(mv,"cp %s %s/",inputsig,nev);
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


// --- Placeholder for create_default_options and parse_options ---
// These will be moved from your original large .c file here for now.
// If they become complex, they might get their own options.h/c module.

void create_default_options(options_t *def) {
    // Implement your default options here based on your original code
    def->evol = 0; // Example default
    def->drift = 0;
    def->growth = 0;
    def->twopop = 0;
    def->ffrag = 0;
    def->ufrag = 0;
    def->input = 1; // Assuming default is generating initial profile
    def->tStep = 1.0; // Example
    def->ngrid = 100; // Example
}

int parse_options(int argc, const char **argv, options_t *def) {
    // Implement your option parsing logic here based on your original code
    // This part involves parsing argc, argv and setting values in def
    // Return 0 on success, non-zero on error

    // Example of a simple loop, replace with your actual getopt/parsing logic
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-e") == 0 && i + 1 < argc) {
            def->evol = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-d") == 0 && i + 1 < argc) {
            def->drift = atoi(argv[++i]);
        }
        // ... add more parsing for other options ...
    }
    return 0; // Success
}
