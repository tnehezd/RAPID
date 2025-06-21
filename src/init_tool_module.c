#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "init_tool_module.h"
#include "config.h"

// Alapértelmezett opciók beállítása
void create_default_init_tool_options(init_tool_options_t *opt) {
    opt->n               = 1000;
    opt->ri              = 0.1;
    opt->ro              = 5.0;
    opt->sigma0          = 0.0001;
    opt->sigma0cgs       = 0.0001 / SDCONV; // Feltételezve, hogy SDCONV elérhető itt
    opt->index           = 0.5; // Ez az alapértelmezett pozitív érték
    opt->rdze_i          = 0.0;
    opt->rdze_o          = 0.0;
    opt->drdze_i         = 0.0; // Ez a szorzó, nem a kiszámolt vastagság
    opt->drdze_o         = 0.0; // Ez a szorzó, nem a kiszámolt vastagság
    opt->alphaParam      = 1.0e-2;
    opt->amod            = 0.01;
    opt->h               = 5.0e-2;
    opt->flind           = 0.0;
    opt->m0              = 1.0;
    opt->md              = 0.01; // Az eredeti kód alapértelmezett MD értéke
    opt->eps             = 0.01;
    opt->ratio           = 0.85;
    opt->mic             = 1e-4;
    opt->onesize         = 1.0;
}

// alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel
// Ezek a függvények mostantól közvetlenül az `init_opts` paramétereit használják.
static double alpha_turb_it(double r, init_tool_options_t *init_opts) {
    // drdze_i és drdze_o számítása a hívó függvényben történik, és átadódik az init_opts-ban
    // Az eredeti kód logikája alapján itt újra kiszámoljuk, ha a fő függvény nem tette meg
    // Fontos, hogy a drdze_i/o ne legyen 0 a nevezőben, ezért a 1e-6 biztonsági érték.
    double drdze_i_val = pow(init_opts->rdze_i, 1.0 + init_opts->flind) * init_opts->h * ((init_opts->drdze_i == 0.0) ? 1e-6 : init_opts->drdze_i);
    double drdze_o_val = pow(init_opts->rdze_o, 1.0 + init_opts->flind) * init_opts->h * ((init_opts->drdze_o == 0.0) ? 1e-6 : init_opts->drdze_o);

    double alpha_r = 1.0 - 0.5 * (1.0 - init_opts->amod) *
                     (tanh((r - init_opts->rdze_i) / drdze_i_val) +
                      tanh((init_opts->rdze_o - r) / drdze_o_val));
    return alpha_r * init_opts->alphaParam;
}

// sigma0 kiszámolása (illetve megadása) M_NAP / AU / AU-ban
static double sigma_null_it(init_tool_options_t *init_opts) {
    // Az eredeti `SIGMAP_EXP` az `index` paraméter negáltja volt.
    // Itt `(-init_opts->index)`-et használunk, mert az `init_opts->index` pozitív, ahogy bejön.
    double alpha2 = (-init_opts->index) + 2;
    double denominator = pow(init_opts->ro, alpha2) - pow(init_opts->ri, alpha2);
    if (fabs(denominator) < 1e-12) {
        fprintf(stderr, "Error: Denominator is zero or too small in Sigma0 calculation! Check RMAX, RMIN, and SIGMAP_EXP values.\n");
        return 0.0;
    }
    return (alpha2 / (2.0 * M_PI)) * init_opts->md / denominator;
}

// sigma kiszámolása r helyen M_NAP / AU / AU-ban
static double sigma_gas_it(double r, init_tool_options_t *init_opts, long double current_sigma0) {
    // Az eredeti kód `pow(r, SIGMAP_EXP)`-et használt, ahol `SIGMAP_EXP` már negálva volt.
    // Itt az `init_opts->index` pozitív, ezért `-init_opts->index`-et használunk.
    return current_sigma0 * pow(r, -init_opts->index);
}

// a por feluletisurusegenek kiszamitasa, a hohatar figyelembevetelevel M_SUN / AU / AU-ban
static long double sigma_dust_it(double r, init_tool_options_t *init_opts, long double current_sigma0) {
    long double sig_dust;
    // Az eredeti kód `pow(r, SIGMAP_EXP)`-et használt, ahol `SIGMAP_EXP` már negálva volt.
    // Itt az `init_opts->index` pozitív, ezért `-init_opts->index`-et használunk.
    sig_dust = current_sigma0 * pow(r, -init_opts->index) * init_opts->eps;

    // IDE JÖHETNE AZ ICEFACTOR ÉS SNOWLINE KEZELÉS, HA SZÜKSÉGES:
    // if (r >= SNOWLINE) {
    //     sig_dust *= ICEFACTOR;
    // }
    return sig_dust;
}

// minimum megkeresese harom elem kozul
static double find_min_it(double s1, double s2, double s3) {
    double min_val = s1;
    if (s2 < min_val) {
        min_val = s2;
    }
    if (s3 < min_val) {
        min_val = s3;
    }
    return min_val;
}

// Az init_tool "main" függvénye, átnevezve és paraméterezve
int run_init_tool(init_tool_options_t *init_opts) {
    FILE *fout = NULL;
    FILE *fout2 = NULL;

    double r, DD;
    long double current_sigma0; // Lokális változó a kiszámolt/átvett sigma0 számára

    // A drdze_i és drdze_o kiszámítása, ahogy az eredeti `main` függvényben volt.
    // Ezeket az értékeket fogjuk kiírni a disk_param.dat fájlba.
    double drdze_i_calculated = pow(init_opts->rdze_i, 1.0 + init_opts->flind) * init_opts->h *
                               ((init_opts->drdze_i == 0.0) ? 1e-6 : init_opts->drdze_i);
    double drdze_o_calculated = pow(init_opts->rdze_o, 1.0 + init_opts->flind) * init_opts->h *
                               ((init_opts->drdze_o == 0.0) ? 1e-6 : init_opts->drdze_o);


    // Döntés sigma0 számításáról: az eredeti kód szerint, ha az `md` paramétert megadták (az alapértelmezett 0.01-től eltér),
    // akkor abból számolunk, egyébként a direkt `sigma0` értéket használjuk.
    // Az `fabs(init_opts->md - 0.01) > 1e-9` ellenőrzés arra szolgál, hogy megnézzük, az `md` eltér-e az alapértelmezettől.
    if (fabs(init_opts->md - 0.01) > 1e-9) {
        current_sigma0 = sigma_null_it(init_opts);
        printf("Sigma0 calculated from Disk Mass (Md): %Lg M_Sun/AU^2\n", current_sigma0);
    } else {
        current_sigma0 = init_opts->sigma0;
        printf("Using explicit Sigma0: %Lg M_Sun/AU^2\n", current_sigma0);
    }

    // Az eredeti kód az `onesize` paraméter megléte esetén felülírta a `ratio`-t 1.0-ra.
    // `init_opts->onesize != 1.0` azt jelenti, hogy a felhasználó valamilyen értékkel beállította, ami nem az alapértelmezett 1.0.
    if (init_opts->onesize != 1.0) {
        init_opts->ratio = 1.0;
    }

    printf("Surface density exponent (index): %lg\n", -init_opts->index); // Kiíratáskor negált formában, ahogy a képletekben használjuk

    fout = fopen("init_data.dat", "w");
    if (fout == NULL) {
        perror("Error opening init_data.dat in init_tool_module");
        return 1;
    }

    DD = (init_opts->ro - init_opts->ri) / (init_opts->n - 1.0);

    printf("\n--- Simulation Parameters ---\n");
    printf("Disk mass (Solar Mass): %lg\n", init_opts->md);
    printf("Inner disk edge (AU): %lg\n", init_opts->ri);
    printf("Outer disk edge (AU): %lg\n", init_opts->ro);
    printf("Surface density profile exponent: %lg\n", -init_opts->index);
    printf("Snowline position (AU): %lg\n", SNOWLINE);
    printf("Ice factor beyond snowline: %lg\n", ICEFACTOR);
    printf("Gas surface density at 1 AU (Solar Mass/AU^2): %Lg\n", current_sigma0);
    printf("Dust to gas ratio: %lg\n", init_opts->eps);
    printf("Number of representative particles: %d\n", (int)init_opts->n);
    printf("------------------------------\n\n");

    double f_drift = 0.55;
    double f_frag = 0.37;

    double rho_p = 1.6; // Por sűrűsége (g/cm^3), ez fix
    double u_frag = 1000.0; // Ütközési sebesség (cm/s)
    u_frag = u_frag * CMPSECTOAUPYRP2PI; // Konverzió AU / (yr/2pi) egységbe
    double u_frag2 = u_frag * u_frag;

    for (int i = 0; i < init_opts->n; i++) {
        r = init_opts->ri + i * DD;
        double reval = r + DD / 2.0; // A cella középpontja
        double reval2 = reval * reval;

        long double reppmass = 2.0 * M_PI * r * DD * sigma_dust_it(reval, init_opts, current_sigma0);

        double s_max;
        if (init_opts->onesize != 1.0) {
            s_max = init_opts->onesize; // Ha explicit méret van megadva, használja azt
        } else {
            double H = pow(reval, 1.0 + init_opts->flind) * init_opts->h;
            double v_kep = sqrt(G_GRAV_CONST2 * init_opts->m0 / reval);
            double omega = v_kep / reval;

            double c_s = omega * H;
            double v_kep2 = v_kep * v_kep;
            double c_s2 = c_s * c_s;

            double Sigma = sigma_gas_it(reval, init_opts, current_sigma0);
            double Sigma_cgs = Sigma / SDCONV;
            long double Sigmad_cgs = sigma_dust_it(reval, init_opts, current_sigma0) / SDCONV;

            double rho_mp = 1.0 / sqrt(2.0 * M_PI) * Sigma / H;
            double P = rho_mp * c_s * c_s;

            // A dPdr számításához `current_sigma0`-t és `init_opts`-ot használunk
            // Az eredeti képletben `SIGMAP_EXP` volt, ami már negatív volt.
            // Itt `-init_opts->index`-et használunk, mert `init_opts->index` pozitív.
            double dPdr = (init_opts->flind + (-init_opts->index) - 2.0) * pow(reval, (init_opts->flind + (-init_opts->index) - 3.0)) * init_opts->h * G_GRAV_CONST2 * init_opts->m0 * current_sigma0 / sqrt(2.0 * M_PI);

            double dlnPdlnr;
            if (P == 0.0) {
                fprintf(stderr, "Error: P is zero in dlnPdlnr calculation at r = %lg. Check input parameters.\n", reval);
                dlnPdlnr = 0.0;
            } else {
                dlnPdlnr = reval / P * dPdr;
            }

            double s_drift = f_drift * 2.0 / M_PI * Sigmad_cgs / rho_p * v_kep2 / c_s2 * fabs(1.0 / dlnPdlnr);
            double s_frag = f_frag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * alpha_turb_it(reval, init_opts)) * u_frag2 / c_s2;

            double dlnPdlnr_abs_cs2_half = fabs(dlnPdlnr * c_s2 * 0.5);
            double s_df;
            if (dlnPdlnr_abs_cs2_half == 0.0) {
                fprintf(stderr, "Error: Denominator is zero in s_df calculation at r = %lg. Check dlnPdlnr value.\n", reval);
                s_df = 1e99; // Nagyon nagy érték, hogy ne ez legyen a minimum
            } else {
                s_df = u_frag * v_kep / dlnPdlnr_abs_cs2_half * 2.0 * Sigma_cgs / (M_PI * rho_p);
            }

            s_max = find_min_it(s_drift, s_frag, s_df);
        }

        if (s_max <= 0) {
            fprintf(stderr, "Warning: s_max <= 0 at r = %lg. This might indicate problematic physical parameters. Setting to a small positive value.\n", reval);
            s_max = 1e-10; // Biztosítsuk, hogy pozitív legyen
        }

        // Kiírás az init_data.dat fájlba
        // reppmass * init_opts->ratio az első populáció tömege
        // reppmass * (1.0 - init_opts->ratio) a második populáció tömege
        // s_max a maximális részecskeméret
        // init_opts->mic a mikroméretű részecskék mérete (ha a TWOPOP modell van használva)
        fprintf(fout, "%d %lg  %Lg %Lg %lg %Lg\n", i, reval,
                reppmass * init_opts->ratio, reppmass * (1.0 - init_opts->ratio),
                s_max, init_opts->mic);
    }

    fclose(fout);

    printf("Particle data file created (init_data.dat). Writing disk parameters file!\n\nFile content:\n i (particle index), r (distance from star), prad (particle size in cm), reppmass (representative particle mass)\n");
    printf("Press ENTER to continue!\n");
    // getchar(); // Kikommentelve, hogy a program ne álljon meg a felhasználói bevitelre várva

    fout2 = fopen("disk_param.dat", "w");
    if (fout2 == NULL) {
        perror("Error opening disk_param.dat in init_tool_module");
        return 1;
    }

    // Paraméterek kiírása a disk_param.dat fájlba
    fprintf(fout2, "%lg %lg %d %lg %Lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
            init_opts->ri, init_opts->ro, (int)init_opts->n, init_opts->index, current_sigma0,
            G_GRAV_CONST, init_opts->rdze_i, init_opts->rdze_o,
            drdze_i_calculated, drdze_o_calculated, // Kiszámolt Dr_dze értékek
            init_opts->amod, rho_p, init_opts->alphaParam, init_opts->m0, init_opts->flind);

    fclose(fout2);

    printf("Disk parameters file created (disk_param.dat).\n\nFile content:\n RMIN, RMAX, NGRID (number of particles), profile exponent (SIGMAP_EXP), and sigma0 (M_Sun / AU / AU) - sigma at r=1AU, G (gravitational constant), inner deadzone edge (AU), transition width (AU), outer deadzone edge (AU), transition width (AU), viscosity reduction factor, average particle density, alpha, and central star mass!\n\n");

    return 0;
}