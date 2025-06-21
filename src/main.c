#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Header files
#include "config.h"           // Declarations of global variables and constants
#include "io_utils.h"         // Functions from io_utils.c
#include "disk_model.h"       // Functions from disk_model.c
#include "dust_physics.h"     // Functions from dust_physics.c
#include "simulation_core.h"  // Functions from simulation_core.c
#include "utils.h"            // Functions from utils.c
#include "init_tool_module.h" // run_init_tool and init_tool_options_t

// options_t structure: Gathers all command-line options.
// This structure also contains parameters used by init_tool
// for initial profile generation.
typedef struct options {
    double drift;
    double growth;
    double evol;
    double twopop;
    double ufrag;
    double ffrag;

    // Core simulation and init_tool parameters
    int ngrid_val;      // NGRID
    double rmin_val;    // RMIN
    double rmax_val;    // RMAX
    double sigma0_val;  // SIGMA0
    double sigmap_exp_val; // SIGMAP_EXP (Note: this is the positive exponent, init_tool will negate if needed)
    double alpha_visc_val; // alpha_visc
    double star_val;    // STAR
    double hasp_val;    // HASP
    double flind_val;    // FLIND
    double r_dze_i_val;  // r_dze_i
    double r_dze_o_val;  // r_dze_o
    double dr_dze_i_val; // Dr_dze_i (This is the multiplier for the transition width)
    double dr_dze_o_val; // Dr_dze_o (This is the multiplier for the transition width)
    double a_mod_val;    // a_mod

    const char *input_file; // Filename for initialization/loading. NULL if not specified.

    double tStep;
    double totalTime;
    double outputFrequency;
    double startTime;

    // Init tool specific options needed for profile generation
    // (these are copied into the init_tool_options_t structure)
    double md_val;
    long double eps_val;
    long double ratio_val;
    long double mic_val;
    long double onesize_val;

} options_t;

// Function declarations
void create_default_options(options_t *def);
// Corrected declaration for parse_options to match its usage
int parse_options(int argc, const char **argv, options_t *def);

// --- Global variables: These come from config.h, no need for further 'extern' declarations here. ---
// The compiler sees them by including 'config.h'.

int main(int argc, const char **argv) {
    // DEBUG: Program start
    printf("DEBUG [main]: Program started.\n");

    // Local structure to store command-line options
    options_t def;
    // Set default values
    create_default_options(&def);
    printf("DEBUG [main]: Default options created.\n");

    // Local structure to store init_tool parameters
    init_tool_options_t init_def;
    // Set default values for init_tool
    create_default_init_tool_options(&init_def);
    printf("DEBUG [main]: Default init_tool options created.\n");

    // Parancssori opciók értelmezése és a 'def' struktúra feltöltése
    int retCode = parse_options(argc, argv, &def);
    
    if (0 != retCode) {
        // Hiba esetén kilépünk
        printf("DEBUG [main]: Error parsing command-line options. Exiting with code %d.\n", retCode);
        return retCode;
    }
    printf("DEBUG [main]: Command-line options parsed successfully.\n");

    /* A globális változók feltöltése a 'def' struktúra értékeivel */
    // Ezek a változók a config.c-ben vannak definiálva, és a config.h-n keresztül érhetők el.
    optev = def.evol;
    optdr = def.drift;
    optgr = def.growth;
    opttwopop = def.twopop;
    fFrag = def.ffrag;
    uFrag = def.ufrag;
    printf("DEBUG [main]: Global simulation control variables set:\n");
    printf("DEBUG [main]:   evol=%.2f, drift=%.2f, growth=%.2f, twopop=%.2f, fFrag=%.2f, uFrag=%.2f\n",
           optev, optdr, optgr, opttwopop, fFrag, uFrag);

    // --- Input fájl kezelésének logikája ---
    if (def.input_file != NULL && strcmp(def.input_file, "") != 0) {
        // Ha van bemeneti fájl megadva (-i kapcsolóval)
        inputsig = def.input_file; // A globális inputsig pointer a fájlnévre mutat
        printf("DEBUG [main]: Input file specified: %s. Attempting to read initial profile.\n", inputsig);
        int lout = 0, nout = 0;
        nout = reszecskek_szama(lout, inputsig); // Meghatározzuk a részecskék számát a fájlból
        NGRID = nout; // A globális NGRID-et beállítjuk
        initialize_derived_config_variables(); // Inicializáljuk a származtatott konfigurációs változókat (pl. DD)
        printf("DEBUG [main]: NGRID set from input file: %d. DD calculated as %.4e.\n", NGRID, DD);
    } else {
        // Ha nincs bemeneti fájl megadva (alapértelmezett profil generálása)
        inputsig = NULL; // Jelzés, hogy nincs explicit input fájl
        printf("DEBUG [main]: No input file specified (-i flag not used). Generating default grid and profile.\n");
        NGRID = def.ngrid_val; // A globális NGRID-et az options_t-ből vesszük
        initialize_derived_config_variables(); // Inicializáljuk a származtatott konfigurációs változókat
        printf("DEBUG [main]: NGRID set to default/command-line value: %d. DD calculated as %.4e.\n", NGRID, DD);
    }
    // Időparaméterek beállítása
    timePar(def.totalTime, def.outputFrequency, def.startTime);
    printf("DEBUG [main]: Time parameters set: TMAX=%.2e, WO=%.2e, TCURR=%.2e\n", TMAX, WO, TCURR);


    // A többi globális konfigurációs változó beállítása az 'options_t'-ből
    RMIN = def.rmin_val;
    RMAX = def.rmax_val;
    SIGMA0 = def.sigma0_val;
    SIGMAP_EXP = def.sigmap_exp_val; // This is the positive exponent, as used by the main simulation
    alpha_visc = def.alpha_visc_val;
    STAR = def.star_val;
    HASP = def.hasp_val;
    FLIND = def.flind_val;
    r_dze_i = def.r_dze_i_val;
    r_dze_o = def.r_dze_o_val;
    Dr_dze_i = def.dr_dze_i_val;
    Dr_dze_o = def.dr_dze_o_val;
    a_mod = def.a_mod_val;
    printf("DEBUG [main]: Global disk parameters set:\n");
    printf("DEBUG [main]:   RMIN=%.2f, RMAX=%.2f, SIGMA0=%.2e, SIGMAP_EXP=%.2f\n", RMIN, RMAX, SIGMA0, SIGMAP_EXP);
    // Corrected printf format string: added HASP
    printf("DEBUG [main]:   alpha_visc=%.2e, STAR=%.2f, HASP=%.2f, FLIND=%.2f\n", alpha_visc, STAR, HASP, FLIND);
    printf("DEBUG [main]:   r_dze_i=%.2f, r_dze_o=%.2f, Dr_dze_i=%.2f, Dr_dze_o=%.2f, a_mod=%.2f\n",
           r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod);
    
    // Az 'options_t' paramétereinek átadása az 'init_tool_options_t' struktúrába
    // Ezt azért tesszük, hogy az init_tool minden szükséges paramétert megkapjon egyetlen struktúrában,
    // függetlenül attól, hogy a fő szimuláció mit használ.
    init_def.n = def.ngrid_val;
    init_def.ri = def.rmin_val;
    init_def.ro = def.rmax_val;
    init_def.sigma0 = def.sigma0_val;
    init_def.sigma0cgs = def.sigma0_val / SDCONV; // SDCONV a config.h-ból
    init_def.index = def.sigmap_exp_val; // Az index az init_tool-ban negatív kitevőként értelmezendő
    init_def.rdze_i = def.r_dze_i_val;
    init_def.rdze_o = def.r_dze_o_val;
    init_def.drdze_i = def.dr_dze_i_val;
    init_def.drdze_o = def.dr_dze_o_val;
    init_def.alphaParam = def.alpha_visc_val;
    init_def.amod = def.a_mod_val;
    init_def.h = def.hasp_val;
    init_def.flind = def.flind_val;
    init_def.m0 = def.star_val;
    init_def.md = def.md_val;
    init_def.eps = def.eps_val;
    init_def.ratio = def.ratio_val;
    init_def.mic = def.mic_val;
    init_def.onesize = def.onesize_val;
    printf("DEBUG [main]: init_tool_options_t (init_def) structure populated.\n");
    printf("DEBUG [main]:   init_def.n=%d, init_def.ri=%.2f, init_def.ro=%.2f, init_def.sigma0=%.2e\n",
           init_def.n, init_def.ri, init_def.ro, init_def.sigma0);
    printf("DEBUG [main]:   init_def.index=%.2f, init_def.alphaParam=%.2e, init_def.m0=%.2f\n",
           init_def.index, init_def.alphaParam, init_def.m0);


    // Tömbök deklarálása a diszk profilokhoz
    // A méret ellenőrzése fontos, NGRID+2 helyett NGRID (ha 0-tól NGRID-1-ig megy)
    // De ha NGRID pontot akarsz, és a 0. és NGRID+1. indexeket is használod, akkor rendben van.
    // feltételezve, hogy a NGRID már be van állítva a fentebbi if-else blokkban.
    double sigmavec[NGRID+2], rvec[NGRID+2], pressvec[NGRID+2], dpressvec[NGRID+2], ugvec[NGRID+2];
    char dens_name[4096], nev[4096], mv[4096];
    printf("DEBUG [main]: Disk profile arrays declared with size NGRID+2 = %d.\n", NGRID+2);


    // Diszk paraméterek beállítása (ez egy meglévő függvény, valószínűleg a disk_model.c-ből)
    printf("DEBUG [main]: Calling disk_param_be to set global disk parameters.\n");
    disk_param_be(&SIGMA0, &SIGMAP_EXP, &RMIN, &RMAX, &r_dze_i, &r_dze_o, &Dr_dze_i, &Dr_dze_o, &a_mod, &PDENSITY, &PDENSITYDIMLESS, &alpha_visc,&STAR,&FLIND);
    printf("DEBUG [main]: disk_param_be completed. PDENSITY=%.2e, PDENSITYDIMLESS=%.2e.\n", PDENSITY, PDENSITYDIMLESS);
    
    // DD globális változó beállítása, ha még nincs inicializálva a initialize_derived_config_variables() által
    // (Bár ez utóbbi már beállítja, így ez redundáns lehet, de nem árt.)
    DD = (RMAX - RMIN) / (NGRID - 1.0); // Ha NGRID pont van, akkor NGRID-1 intervallum
    printf("DEBUG [main]: DD re-calculated: %.4e (RMAX-RMIN)/(NGRID-1) = (%.2f-%.2f)/(%d-1)\n", DD, RMAX, RMIN, NGRID);


    // --- Kezdeti profil betöltése vagy generálása ---
    if(inputsig != NULL) {
        // Ha inputsig NEM NULL (azaz vagy a felhasználó adta meg, vagy az init_tool generálta és mi beállítottuk)
        printf("DEBUG [main]: inputsig is NOT NULL. Loading initial profile from %s.\n", inputsig);
        printf("DEBUG [main]: Calling sigIn(sigmavec,rvec)...\n");
        sigIn(sigmavec,rvec); // Fájlból olvassa be a sigma és r értékeket
        printf("DEBUG [main]: sigIn completed. Calling Perem for rvec and sigmavec...\n");
        Perem(rvec);           // Peremfeltételek beállítása az rvec-re
        Perem(sigmavec);       // Peremfeltételek beállítása a sigmavec-re
        printf("DEBUG [main]: Perem calls completed for initial profile.\n");
    } else {
        // Ha inputsig NULL (nincs bemeneti fájl, generálnunk kell)
        printf("DEBUG [main]: inputsig is NULL. Generating initial profile with init_tool...\n");
        printf("DEBUG [main]: Calling run_init_tool(&init_def)...\n");
        run_init_tool(&init_def); // Meghívjuk az init_tool-t a profil generálására
        printf("DEBUG [main]: run_init_tool completed.\n");

        // Az init_tool generált fájljának betöltése
        inputsig = filenev1; // A globális filenev1 (ami "init_data.dat") lesz az inputsig
        printf("DEBUG [main]: Generated profile set to be loaded from %s.\n", inputsig);

        // Mivel az NGRID esetleg változhatott az init_tool által generált fájl miatt,
        // újra meg kell határoznunk a részecskék számát és inicializálni a konfigurációt.
        printf("DEBUG [main]: Re-calculating NGRID from generated file %s.\n", inputsig);
        int lout = 0, nout = 0;
        nout = reszecskek_szama(lout, inputsig);
        NGRID = nout;
        initialize_derived_config_variables();
        printf("DEBUG [main]: NGRID re-set to %d. DD re-calculated as %.4e.\n", NGRID, DD);

        printf("DEBUG [main]: Calling sigIn(sigmavec,rvec) for generated profile...\n");
        sigIn(sigmavec, rvec); // Betöltjük a frissen generált fájlból
        printf("DEBUG [main]: sigIn completed for generated profile. Calling Perem...\n");
        Perem(rvec);
        Perem(sigmavec);
        printf("DEBUG [main]: Perem calls completed for generated profile.\n");
    }
    // --- Vége a kezdeti profil kezelésének ---

    // Nyomás- és gázsebesség profilok inicializálása
    printf("DEBUG [main]: Initializing pressure and gas velocity profiles...\n");
    Initial_Press(pressvec,sigmavec,rvec);
    Initial_dPress(dpressvec,pressvec);
    Initial_Ugas(sigmavec,rvec,ugvec);
    printf("DEBUG [main]: Pressure and gas velocity profiles initialized.\n");

    // Dead Zone opció beállítása
    optdze = 0; // Default to inactive
    if(r_dze_i != 0.0 || r_dze_o != 0.0) { // If either inner or outer dead zone radius is non-zero
        optdze = 1.;
        printf("DEBUG [main]: Dead Zone active (r_dze_i != 0.0 or r_dze_o != 0.0).\n");
    } else {
        printf("DEBUG [main]: Dead Zone inactive (both r_dze_i and r_dze_o are 0.0).\n");
    }

    // Kimeneti könyvtár létrehozása
    printf("DEBUG [main]: Creating output directory...\n");
    Mk_Dir(nev);
    printf("DEBUG [main]: Output directory created: %s\n", nev);

    int dummy; // Dummy változó a system() hívások visszatérési értékének figyelmen kívül hagyására

    // Konfigurációs és időzítő fájlok másolása a kimeneti könyvtárba
    printf("DEBUG [main]: Copying config and time files to output directory.\n");
    snprintf(mv, sizeof(mv), "cp %s %s/", filenev2, nev); // filenev2 ("disk_param.dat")
    dummy = system(mv); (void)dummy; // (void)dummy a warning elkerülésére
    snprintf(mv, sizeof(mv), "cp %s %s/", filenev3, nev); // filenev3 ("time.dat")
    dummy = system(mv); (void)dummy;
    printf("DEBUG [main]: disk_param.dat and time.dat copied.\n");

    // Az input profil fájl másolása a kimeneti könyvtárba (akár felhasználói, akár generált)
    if(inputsig != NULL) { // Ha volt bemeneti fájl (akár generált, akár külső)
        printf("DEBUG [main]: Copying initial profile file %s to output directory.\n", inputsig);
        snprintf(mv, sizeof(mv), "cp %s %s/", inputsig, nev);
        dummy = system(mv); (void)dummy;
        printf("DEBUG [main]: Initial profile file copied.\n");
    } else {
        printf("DEBUG [main]: No specific input file to copy (initial profile was generated, it's %s already).\n", filenev1);
    }

    // Aktuális információk kiírása
    printf("DEBUG [main]: Calling infoCurrent...\n");
    infoCurrent(nev);
    printf("DEBUG [main]: infoCurrent completed.\n");

    // Szimuláció futtatása vagy kilépés opciók alapján
    if(optev == 0. && optdr == 0.) {
        printf("DEBUG [main]: Evolution (optev=%.2f) and drift (optdr=%.2f) are OFF.\n", optev, optdr);
        printf("A megadott opciok szerint nem szamol sem sigmat, sem driftet, ezert a program kilep!\n\nA kezdeti file-ok a %s mappaban talalhatoak!\n",nev);
        snprintf(dens_name,sizeof(dens_name),"%s/surface.dat",nev);
        printf("DEBUG [main]: Printing initial surface density to %s.\n", dens_name);
        Print_Sigma(dens_name,rvec,sigmavec,pressvec,dpressvec);
        printf("DEBUG [main]: Print_Sigma completed. Program exiting.\n");
    } else {
        printf("DEBUG [main]: Evolution (optev=%.2f) or drift (optdr=%.2f) is ON. Starting main simulation loop.\n", optev, optdr);
        printf("DEBUG [main]: Calling tIntegrate...\n");
        tIntegrate(nev,rvec,sigmavec,pressvec,dpressvec,ugvec);
        printf("DEBUG [main]: tIntegrate completed. Program finished normally.\n");
    }

    printf("DEBUG [main]: Program exiting normally.\n");
    return 0;
}

/* --- create_default_options: A 'options_t' struktúra alapértelmezett értékeinek beállítása --- */
void create_default_options(options_t *opt) {
    printf("DEBUG [create_default_options]: Setting default values for options_t.\n");
    // Szimuláció vezérlő opciók
    opt->drift             = 1.;
    opt->growth            = 1.;
    opt->evol              = 1.;
    opt->twopop            = 1.;
    opt->ufrag             = 1000.0;
    opt->ffrag             = 0.37;
    
    // Alapvető diszk paraméterek (ezek az init_tool alapértelmezettjeivel is megegyeznek)
    opt->ngrid_val         = 2000;
    opt->rmin_val          = 1.0;
    opt->rmax_val          = 100.0;
    opt->sigma0_val        = 1.0;
    opt->sigmap_exp_val    = 1.0;    // e.g., 1.0 for r^-1 profile. This is the positive exponent.
    opt->alpha_visc_val    = 0.01;
    opt->star_val          = 1.0;
    opt->hasp_val          = 0.05;
    opt->flind_val         = 0.5;
    opt->r_dze_i_val       = 0.0;
    opt->r_dze_o_val       = 0.0;
    opt->dr_dze_i_val      = 0.0; // Multiplier for transition width (input parameter for init_tool)
    opt->dr_dze_o_val      = 0.0; // Multiplier for transition width (input parameter for init_tool)
    opt->a_mod_val         = 0.0;

    // Fájl input
    opt->input_file        = ""; // Default to empty string if no -i flag

    // Időparaméterek
    opt->tStep             = 0.; // This seems unused or implicitly handled by timePar
    opt->totalTime         = 1.0e6;
    opt->outputFrequency   = 1000.0;
    opt->startTime         = 0.0;

    // Init tool specific parameters' defaults (only if not set by the above parameter acquisition)
    // Note: The `init_tool_options_t` default values are more meaningful.
    // These `opt->md_val`, `opt->eps_val`, etc. here are just the initial state of the `options_t` structure.
    // They will be populated either by command-line arguments or by the `create_default_init_tool_options` call.
    opt->md_val = 0.01; // Matches init_tool default
    opt->eps_val = 0.01; // Matches init_tool default
    opt->ratio_val = 0.85; // Matches init_tool default
    opt->mic_val = 1e-4; // Matches init_tool default
    opt->onesize_val = 1.0; // Matches init_tool default
    printf("DEBUG [create_default_options]: Default options setting complete.\n");
}

/* --- parse_options: A parancssori argumentumok értelmezése --- */
// Corrected function definition to match the declaration and use the 'opt' parameter
int parse_options(int argc, const char **argv, options_t *opt){
    printf("DEBUG [parse_options]: Parsing command-line arguments (%d total).\n", argc);
    int i = 1; // Start from 1, as argv[0] is the program name

    while (i < argc) {
        printf("DEBUG [parse_options]: Processing argument %d: %s\n", i, argv[i]);
        if(strcmp(argv[i], "-drift") == 0) {
            i++;
            if (i < argc) opt->drift = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drift.\n"); return 1; }
            printf("DEBUG [parse_options]:   -drift set to %.2f\n", opt->drift);
        }
        else if (strcmp(argv[i], "-growth") == 0) {
            i++;
            if (i < argc) opt->growth = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -growth.\n"); return 1; }
            printf("DEBUG [parse_options]:   -growth set to %.2f\n", opt->growth);
        }
        else if (strcmp(argv[i], "-evol") == 0) {
            i++;
            if (i < argc) opt->evol = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -evol.\n"); return 1; }
            printf("DEBUG [parse_options]:   -evol set to %.2f\n", opt->evol);
        }
        else if (strcmp(argv[i], "-twopop") == 0) {
            i++;
            if (i < argc) opt->twopop = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -twopop.\n"); return 1; }
            printf("DEBUG [parse_options]:   -twopop set to %.2f\n", opt->twopop);
        }
        else if (strcmp(argv[i], "-ufrag") == 0) {
            i++;
            if (i < argc) opt->ufrag = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ufrag.\n"); return 1; }
            printf("DEBUG [parse_options]:   -ufrag set to %.2f\n", opt->ufrag);
        }
        else if (strcmp(argv[i], "-ffrag") == 0) {
            i++;
            if (i < argc) opt->ffrag = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ffrag.\n"); return 1; }
            printf("DEBUG [parse_options]:   -ffrag set to %.2f\n", opt->ffrag);
        }
        else if (strcmp(argv[i], "-tStep") == 0) {
            i++;
            if (i < argc) opt->tStep = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -tStep.\n"); return 1; }
            printf("DEBUG [parse_options]:   -tStep set to %.2e\n", opt->tStep);
        }
        else if (strcmp(argv[i], "-n") == 0) { // Main simulation NGRID
            i++;
            if (i < argc) opt->ngrid_val = atoi(argv[i]); else { fprintf(stderr, "Error: Missing value for -n.\n"); return 1; }
            printf("DEBUG [parse_options]:   -n (NGRID) set to %d\n", opt->ngrid_val);
        }
        else if (strcmp(argv[i], "-i") == 0) {
            i++;
            if (i < argc) opt->input_file = argv[i]; else { fprintf(stderr, "Error: Missing value for -i.\n"); return 1; }
            printf("DEBUG [parse_options]:   -i (input_file) set to '%s'\n", opt->input_file);
        }
        else if (strcmp(argv[i], "-tmax") == 0) {
            i++;
            if (i < argc) opt->totalTime = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -tmax.\n"); return 1; }
            printf("DEBUG [parse_options]:   -tmax set to %.2e\n", opt->totalTime);
        }
        else if (strcmp(argv[i], "-outfreq") == 0) {
            i++;
            if (i < argc) opt->outputFrequency = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -outfreq.\n"); return 1; }
            printf("DEBUG [parse_options]:   -outfreq set to %.2e\n", opt->outputFrequency);
        }
        else if (strcmp(argv[i], "-curr") == 0) {
            i++;
            if (i < argc) opt->startTime = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -curr.\n"); return 1; }
            printf("DEBUG [parse_options]:   -curr (startTime) set to %.2e\n", opt->startTime);
        }
        // --- Init_tool specific opciók feldolgozása a fő parserben ---
        else if (strcmp(argv[i], "-n_init") == 0) { i++; if (i < argc) opt->ngrid_val = atoi(argv[i]); else { fprintf(stderr, "Error: Missing value for -n_init.\n"); return 1; } printf("DEBUG [parse_options]:   -n_init (NGRID for init) set to %d\n", opt->ngrid_val); }
        else if (strcmp(argv[i], "-ri") == 0) { i++; if (i < argc) opt->rmin_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ri.\n"); return 1; } printf("DEBUG [parse_options]:   -ri (RMIN for init) set to %.2f\n", opt->rmin_val); }
        else if (strcmp(argv[i], "-ro") == 0) { i++; if (i < argc) opt->rmax_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ro.\n"); return 1; } printf("DEBUG [parse_options]:   -ro (RMAX for init) set to %.2f\n", opt->rmax_val); }
        else if (strcmp(argv[i], "-sigma0_init") == 0) { i++; if (i < argc) opt->sigma0_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -sigma0_init.\n"); return 1; } printf("DEBUG [parse_options]:   -sigma0_init set to %.2e\n", opt->sigma0_val); }
        else if (strcmp(argv[i], "-index_init") == 0) { i++; if (i < argc) opt->sigmap_exp_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -index_init.\n"); return 1; } printf("DEBUG [parse_options]:   -index_init set to %.2f\n", opt->sigmap_exp_val); }
        else if (strcmp(argv[i], "-rdzei") == 0) { i++; if (i < argc) opt->r_dze_i_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -rdzei.\n"); return 1; } printf("DEBUG [parse_options]:   -rdzei set to %.2f\n", opt->r_dze_i_val); }
        else if (strcmp(argv[i], "-rdzeo") == 0) { i++; if (i < argc) opt->r_dze_o_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -rdzeo.\n"); return 1; } printf("DEBUG [parse_options]:   -rdzeo set to %.2f\n", opt->r_dze_o_val); }
        else if (strcmp(argv[i], "-drdzei") == 0) { i++; if (i < argc) opt->dr_dze_i_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drdzei.\n"); return 1; } printf("DEBUG [parse_options]:   -drdzei set to %.2f\n", opt->dr_dze_i_val); }
        else if (strcmp(argv[i], "-drdzeo") == 0) { i++; if (i < argc) opt->dr_dze_o_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -drdzeo.\n"); return 1; } printf("DEBUG [parse_options]:   -drdzeo set to %.2f\n", opt->dr_dze_o_val); }
        else if (strcmp(argv[i], "-alpha_init") == 0) { i++; if (i < argc) opt->alpha_visc_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -alpha_init.\n"); return 1; } printf("DEBUG [parse_options]:   -alpha_init set to %.2e\n", opt->alpha_visc_val); }
        else if (strcmp(argv[i], "-amod") == 0) { i++; if (i < argc) opt->a_mod_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -amod.\n"); return 1; } printf("DEBUG [parse_options]:   -amod set to %.2f\n", opt->a_mod_val); }
        else if (strcmp(argv[i], "-h_init") == 0) { i++; if (i < argc) opt->hasp_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -h_init.\n"); return 1; } printf("DEBUG [parse_options]:   -h_init set to %.2f\n", opt->hasp_val); }
        else if (strcmp(argv[i], "-flind_init") == 0) { i++; if (i < argc) opt->flind_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -flind_init.\n"); return 1; } printf("DEBUG [parse_options]:   -flind_init set to %.2f\n", opt->flind_val); }
        else if (strcmp(argv[i], "-m0_init") == 0) { i++; if (i < argc) opt->star_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -m0_init.\n"); return 1; } printf("DEBUG [parse_options]:   -m0_init set to %.2f\n", opt->star_val); }
        else if (strcmp(argv[i], "-md_init") == 0) { i++; if (i < argc) opt->md_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -md_init.\n"); return 1; } printf("DEBUG [parse_options]:   -md_init set to %.2e\n", opt->md_val); }
        else if (strcmp(argv[i], "-eps_init") == 0) { i++; if (i < argc) opt->eps_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -eps_init.\n"); return 1; } printf("DEBUG [parse_options]:   -eps_init set to %Le\n", opt->eps_val); }
        else if (strcmp(argv[i], "-ratio_init") == 0) { i++; if (i < argc) opt->ratio_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -ratio_init.\n"); return 1; } printf("DEBUG [parse_options]:   -ratio_init set to %Le\n", opt->ratio_val); }
        else if (strcmp(argv[i], "-onesize_init") == 0) { i++; if (i < argc) opt->onesize_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -onesize_init.\n"); return 1; } printf("DEBUG [parse_options]:   -onesize_init set to %Le\n", opt->onesize_val); }
        else if (strcmp(argv[i], "-mic_init") == 0) { i++; if (i < argc) opt->mic_val = atof(argv[i]); else { fprintf(stderr, "Error: Missing value for -mic_init.\n"); return 1; } printf("DEBUG [parse_options]:   -mic_init set to %Le\n", opt->mic_val); }
        // --- Vége az init_tool opciók feldolgozásának ---
        else {
            // Ismeretlen kapcsoló esetén hibaüzenet és kilépés
            fprintf(stderr, "ERROR [parse_options]: Invalid switch on command-line: %s!\n", argv[i]);
            printf("\n\n**** Invalid switch on command-line: %s! ****\n\n", argv[i]);
            printf("**** Try following parameters: ****\n\n-drift <val>\n-growth <val>\n-evol <val>\n-twopop <val>\n-ufrag <val>\n-ffrag <val>\n-tStep <val>\n-n <val>\n-i <filename>\n-tmax <val>\n-outfreq <val>\n-curr <val>\n");
            printf("\n**** Initial profile parameters (used if -i is not provided): ****\n-n_init <val>\n-ri <val>\n-ro <val>\n-sigma0_init <val>\n-index_init <val>\n-rdzei <val>\n-rdzeo <val>\n-drdzei <val>\n-drdzeo <val>\n-alpha_init <val>\n-amod <val>\n-h_init <val>\n-flind_init <val>\n-m0_init <val>\n-md_init <val>\n-eps_init <val>\n-ratio_init <val>\n-onesize_init <val>\n-mic_init <val>\n");
            return 1; // Hiba kód visszaadása
        }
        i++;
    }

    printf("DEBUG [parse_options]: Command-line parsing complete.\n");
    return 0; // Sikeres értelmezés
}