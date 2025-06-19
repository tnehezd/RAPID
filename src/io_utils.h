#ifndef IO_UTILS_H
#define IO_UTILS_H

// Include necessary headers if functions here depend on types/definitions
// from other modules (e.g., if Print_Mass uses a Particle struct).
// For now, minimal includes:
#include <stdio.h> // For FILE*

// --- Function Declarations ---
// --- Function Declarations ---
int reszecskek_szama(int lout, int inputsig); // Will re-evaluate parameters after content is moved
void por_be();
void sigIn(double sigmavec[], double rvec[]); // You had `inputsig` here as an int. If it's a filename, it should be const char*
                                              // Based on the error, sigIn(const char* filename, double sigmavec[], double rvec[]) is more likely.
                                              // Please correct this signature based on how it's used in your original code.
                                              // For now, I'll assume the original 'inputsig' was actually meant to be a filename string.
                                              // Let's assume it should be:
// void sigIn(const char *filename, double sigmavec[], double rvec[]); // Corrected
// If 'inputsig' in main is an int controlling *whether* to read a file, and the filename itself is a global or passed separately:
// void sigIn(double sigmavec[], double rvec[]); // As it was, but the error suggests its first param is used as a filename.
// Let's go with the error's implication for now:

void Mk_Dir(char *nev);
void infoCurrent(char *nev);

// Correct Print_Mass signature: Ensure the array dimensions match ([][4])
void Print_Mass(double step, double *rvec, double partmassind[][4], double partmassmicrind[][4], double partmasssecind[][4], double t, double *dpressvec, double massbtempii, double massbtempoi, double massmtempii, double massmtempoi, double *massbtempio, double *massbtempoo, double *massmtempio, double *massmtempoo, double *tavin, double *tavout);

void Print_Sigma(char *filename, double rvec[], double sigmavec[], double pressvec[], double dpressvec[]);

// Correct Print_Sigmad signature: The error message indicated 'min' was an unused parameter. Let's adjust or keep based on actual use.
// Assuming 'min' was a placeholder and 'r' and 'rm' are double arrays:
void Print_Sigmad(char *dust_name, char *dust_name2, double *r, double *rm, double *sigmad, double *sigmadm);

void Print_Pormozg_Size(char *size_name, int step, double rad[][2], double radmicr[][2], double *rvec, double t);



#endif // IO_UTILS_H
