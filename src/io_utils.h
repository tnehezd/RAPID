#ifndef IO_UTILS_H
#define IO_UTILS_H

// Include necessary headers if functions here depend on types/definitions
// from other modules (e.g., if Print_Mass uses a Particle struct).
// For now, minimal includes:
#include <stdio.h> // For FILE*

// --- Function Declarations ---
int reszecskek_szama(int lout, int inputsig);
void por_be(); // Add parameters if it has any
void sigIn(double sigmavec[], double rvec[]); // Add parameters
void Mk_Dir(char *nev);
void infoCurrent(char *nev);
void Print_Mass(double L, double rvec[], double partmassind[][2], double partmassmicrind[][2], double partmasssecind[][2], double t, double dpressvec[], double masstempiin, double masstempoin, double massmtempiin, double massmtempoin, double *masstempiout, double *masstempoout, double *massmtempiout, double *massmtempoout, double *tavin, double *tavout);
void Print_Sigma(char *filename, double rvec[], double sigmavec[], double pressvec[], double dpressvec[]);
void Print_Sigmad(char *dust_name, char *dust_name2, double mint, double rdvec[], double rmicvec[], double sigmad[], double sigmadm[]);
void Print_Pormozg_Size(char *size_name, int L, double radius[][2], double radiusmicr[][2], double rvec[], double t);


#endif // IO_UTILS_H
