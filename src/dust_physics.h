// src/dust_physics.h

#ifndef DUST_PHYSICS_H
#define DUST_PHYSICS_H

#include "config.h" // Szükséges a globális konstansokhoz (G2, STAR, HASP, NGRID, M_PI, DD, optdze stb.)
#include <math.h>   // Szükséges matematikai függvényekhez (pow, sqrt)

// -----------------------------------------------------------------------------
// NYILVÁNOS FÜGGVÉNYEK (Más modulok is hívhatják őket)
// -----------------------------------------------------------------------------

// pressure scale-height
double scale_height(double r);


// Korongnyomás számítása
double press(double sigma, double r);

// Nyomás gradiensének számítása
void dpress(double *dp, double *p);

// Gáz energia fluxusának számítása (feltételezem, hogy ez az "u_gas" jelentése)
void u_gas(double *sigmavec, double *rvec, double *ug);

void Get_Sigmad(double L, double max, double min, double rad[][2], double radmicr[][2], double radsec[][2], double *sigma_d, double *sigma_dm, double *sigma_ds, double *massvec, double *massmicrvec, double *masssecvec, double *rd, double *rmic, double *rs);

// Részecsketömegek összesítése gyűrűkben
void GetMass(int n, double partmassind[][4],int indii, int indio, double tavi, double dzei, double *massiout,int indoi, int indoo, double tavo, double dzeo, double *massoout);

// -----------------------------------------------------------------------------
// BELSŐ SEGÉDFÜGGVÉNYEK (Csak a dust_physics.c-n belülről hívják őket)
//
// Ezeket a függvényeket technikailag lehetne "static" deklarációval is ellátni a .c fájlban
// (pl. `static double scale_height(double r);`), ha biztosan tudjuk, hogy soha nem lesznek
// más .c fájlokból hívva. Ez jó gyakorlat a modularitásra.
// DE most, hogy a fordítási hibákat orvosoljuk, a legegyszerűbb, ha deklaráljuk őket itt.
// Később optimalizálhatjuk a hatókört.
// -----------------------------------------------------------------------------

// Stokes number
double Stokes_Number(double pradius, double sigma);

// viscosity
double visc(double r);

// Lokális skála magasság számítása
double scale_height(double r);

// Lokális kepleri sebesség számítása
double v_kep(double r);

// Lokális kepleri körfrekvencia számítása
double kep_freq(double r);

// Lokális hangsebesség számítása
double c_sound(double r);

// alpha turb param
double alpha_turb(double r);

// Sűrűség a középsíkon
double rho_mp(double sigma, double r);

// Együttható az u_gas számításához
double Coeff_3(double sigma, double r);




//reprezentativ reszecske kezdeti meretenek meghatarozasa
// 1. radialis drift altal meghatarozott maximalis meret            --> kimenet cm-ben!
double a_drift(double sigmad, double r, double p, double dp, double rho_p);

// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret --> kimenet cm-ben!
double a_turb(double sigma, double r, double rho_p);

// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret        --> kimenet cm-ben!
double a_df(double sigma, double r, double p, double dp, double rho_p);

//  a reszecskek novekedesenek idoskalaja
double tauGr(double r, double eps);

//    kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKKK
double getSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress, double dt);


void Get_Sigma_P_dP(double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double deltat);

void Get_Radius(char *nev, int opt, double radius[][2], double *pressvec, double *dpressvec, double *sigmavec, double *sigmad, double *rdvec, double *rvec, double *ugvec, double deltat, double t, int n);
#endif // DUST_PHYSICS_H
