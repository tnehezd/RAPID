#ifndef DISK_MODEL_H
#define DISK_MODEL_H

// Necessary includes for types used in function prototypes
#include "config.h" // For global variables like RMIN, DD, NGRID, filenev2, etc.
#include <stdio.h>  // For FILE* (if needed, though not strictly in prototypes)

/*
 * @brief Reads disk parameters from a file and updates global variables.
 * @param sigma0 Pointer to SIGMA0 (surface density at 1AU)
 * @param sdexp Pointer to SIGMAP_EXP (surface density profile exponent)
 * @param Rmin Pointer to RMIN (inner disk boundary)
 * @param Rmax Pointer to RMAX (outer disk boundary)
 * @param r_dzei Pointer to r_dze_i (inner dead zone radius)
 * @param r_dzeo Pointer to r_dze_o (outer dead zone radius)
 * @param dr_dzei Pointer to Dr_dze_i (inner dead zone transition width)
 * @param dr_dzeo Pointer to Dr_dze_o (outer dead zone transition width)
 * @param alph_mod Pointer to a_mod (viscosity reduction factor)
 * @param rho_p Pointer to PDENSITY (average particle density in cgs)
 * @param rho_p_dimless Pointer to PDENSITYDIMLESS (dimensionless particle density)
 * @param alphav Pointer to alpha_visc (alpha viscosity parameter)
 * @param mStar Pointer to STAR (central star mass)
 * @param gamma Pointer to FLIND (flaring index)
 */
void disk_param_be(double *sigma0, double *sdexp, double *Rmin, double *Rmax,
                   double *r_dzei, double *r_dzeo, double *dr_dzei, double *dr_dzeo,
                   double *alph_mod, double *rho_p, double *rho_p_dimless,
                   double *alphav, double *mStar, double *gamma);

/*
 * @brief Fits a parabola to boundary conditions for a given vector.
 * @param vec The array to apply boundary conditions to.
 * @param i1 Index for the first point.
 * @param i2 Index for the second point.
 * @param i3 Index for the third point.
 * @param a Pointer to the 'a' coefficient of the parabola.
 * @param b Pointer to the 'b' coefficient of the parabola.
 * @param c Pointer to the 'c' coefficient of the parabola.
 * @param dd Grid spacing.
 */
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd);

/*
 * @brief Initializes the radial grid (rvec).
 * @param rvec Array to store radial grid points.
 */
void load_R(double *rvec);

/*
 * @brief Initializes the gas surface density profile.
 * @param sigmavec Array to store gas surface density values.
 * @param r Array of radial grid points.
 */
void Initial_Profile(double *sigmavec, double *r);

/*
 * @brief Initializes the gas pressure profile.
 * @param pressvec Array to store gas pressure values.
 * @param sigmavec Array of gas surface density values.
 * @param rvec Array of radial grid points.
 */
void Initial_Press(double *pressvec, double *sigmavec, double *rvec);

/*
 * @brief Initializes the radial pressure gradient profile.
 * @param dpressvec Array to store radial pressure gradient values.
 * @param pressvec Array of gas pressure values.
 */
void Initial_dPress(double *dpressvec, double *pressvec);

/*
 * @brief Initializes the gas velocity perturbation profile.
 * @param sigmavec Array of gas surface density values.
 * @param rvec Array of radial grid points.
 * @param ug Array to store gas velocity perturbation values.
 */
void Initial_Ugas(double *sigmavec, double *rvec, double *ug);

#endif // DISK_MODEL_H
