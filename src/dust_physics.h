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

// Részecsketömegek összesítése gyűrűkben
void GetMass(int n, double partmassind[][4], int indii, int indio, double tavi, double dzei,
             double *massiout, int indoi, int indoo, double tavo, double dzeo, double *massoout);

// Nulla pontok számolása a nyomás gradiensében (általában nyomás maximumok)
int find_num_zero(double *rvec, double *dp);

// Nulla pont meghatározása lineáris interpolációval
double find_zero(int i, double *rvec, double *dp);


// -----------------------------------------------------------------------------
// BELSŐ SEGÉDFÜGGVÉNYEK (Csak a dust_physics.c-n belülről hívják őket)
//
// Ezeket a függvényeket technikailag lehetne "static" deklarációval is ellátni a .c fájlban
// (pl. `static double scale_height(double r);`), ha biztosan tudjuk, hogy soha nem lesznek
// más .c fájlokból hívva. Ez jó gyakorlat a modularitásra.
// DE most, hogy a fordítási hibákat orvosoljuk, a legegyszerűbb, ha deklaráljuk őket itt.
// Később optimalizálhatjuk a hatókört.
// -----------------------------------------------------------------------------

// Lokális skála magasság számítása
double scale_height(double r);

// Lokális kepleri sebesség számítása
double v_kep(double r);

// Lokális kepleri körfrekvencia számítása
double kep_freq(double r);

// Lokális hangsebesség számítása
double c_sound(double r);

// Sűrűség a középsíkon
double rho_mp(double sigma, double r);

// Együttható az u_gas számításához
double Coeff_3(double sigma, double r);

// Nulla sugarának meghatározása lineáris illesztéssel (find_zero hívja)
double find_r_zero(double r1, double r2, double dp1, double dp2);


#endif // DUST_PHYSICS_H
