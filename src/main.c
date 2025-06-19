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
    // Ezek az alapértelmezett értékek, ha nem adunk meg semmilyen kapcsolót.
    // Fontos, hogy ezeket az eredeti kódod alapértelmezett beállításaival egyeztesd!
    def->evol = 0;       // -e
    def->drift = 0;      // -d
    def->growth = 0;     // -g
    def->twopop = 0;     // -t
    def->ffrag = 0;      // -f
    def->ufrag = 0;      // -u
    def->input = 1;      // -i (0-ha fájlból, 1-ha generált)
    def->tStep = 0.0;    // -s (simulation time step)
    def->ngrid = 100;    // -n (NGRID, ha input=1)
    // Ha vannak további opcióid, add hozzá itt is az alapértelmezett értéküket.
}

int parse_options(int argc, const char **argv, options_t *def) {
    // Iterálunk a parancssori argumentumokon. Az első (argv[0]) a program neve,
    // ezért 1-től indulunk.
    for (int i = 1; i < argc; i++) {
        // -e: Gas evolution
        if (strcmp(argv[i], "-e") == 0 && i + 1 < argc) {
            def->evol = atoi(argv[++i]);
        }
        // -d: Particle drift
        else if (strcmp(argv[i], "-d") == 0 && i + 1 < argc) {
            def->drift = atoi(argv[++i]);
        }
        // -g: Particle growth
        else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
            def->growth = atoi(argv[++i]);
        }
        // -t: Two-population model
        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            def->twopop = atoi(argv[++i]);
        }
        // -f: Fragmentation flag
        else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            def->ffrag = atoi(argv[++i]);
        }
        // -u: Fragmentation update flag
        else if (strcmp(argv[i], "-u") == 0 && i + 1 < argc) {
            def->ufrag = atoi(argv[++i]);
        }
        // -i: Input sigma file option (0 for file, 1 for generated)
        // Ha ezt 0-ra állítjuk, akkor a következő argumentum a fájlnév lesz.
        // Ezt a fájlnevet külön kell majd kezelni, valószínűleg a globális filenev2-be kell másolni.
        else if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            def->input = atoi(argv[++i]);
            if (def->input == 0) {
                // Ha az input 0, akkor a következő argumentum a fájlnév.
                // Ezt a globális filenev2-be másoljuk át a config.h-ból.
                // A parse_options függvényből közvetlenül globális változókat írni elfogadható.
                if (i + 1 < argc) {
                    strncpy(filenev2, argv[++i], sizeof(filenev2) - 1);
                    filenev2[sizeof(filenev2) - 1] = '\0'; // Null-terminátor biztosítása
                } else {
                    fprintf(stderr, "Error: -i option requires a filename when set to 0.\n");
                    return 1; // Hiba kód
                }
            }
        }
        // -s: Time step (DT)
        else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            def->tStep = atof(argv[++i]);
        }
        // -n: NGRID (if inputsig is 1, i.e., generated profile)
        else if (strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            def->ngrid = atoi(argv[++i]);
        }
        // Kezeletlen vagy érvénytelen opció
        else {
            fprintf(stderr, "Warning: Unknown or incomplete option '%s'\n", argv[i]);
            // return 1; // Dönthetünk úgy, hogy hiba esetén kilépünk, vagy csak figyelmeztetünk.
        }
    }
    return 0; // Siker
}
