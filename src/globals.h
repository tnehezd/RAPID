#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdio.h> 
#include <math.h> 

/**
 * @file globals.h
 * @brief This file contains global constants (macros) used throughout the simulation.
 *
 * These constants include physical conversion factors, astronomical values,
 * and other simulation-specific parameters defined as preprocessor macros.
 * All constants are defined using #define directives.
 */

// --- Physical Constants (Macros) ---
/** @def GAS_SD_CONV_RATE
 * @brief Surface density conversion factor from g/cm^2 to M_Sun/AU^2.
 */
#define GAS_SD_CONV_RATE             1.12521e-7                   

/** @def G_GRAV_CONST
 * @brief Gravitational constant (G=1) in the chosen simulation units.
 */
#define G_GRAV_CONST 1.0                             

/** @def CM_PER_SEC_TO_AU_PER_YEAR_OVER_2PI
 * @brief Conversion factor from cm/s to AU/(year / 2*pi).
 */
#define CM_PER_SEC_TO_AU_PER_YEAR_OVER_2PI   3.35725e-07                                        

/** @def SUN2GR
 * @brief Conversion rate between Solar Mass and grams.
 */
#define SUN2GR 1.989e33                              

/** @def AU2CM
 * @brief Conversion rate between Astronomical Unit and centimeters (AU -> cm).
 */
#define AU2CM 1.496e13                               

/** @def TWOPI
 * @brief Represents 2 * PI (two times Pi).
 */
#define TWOPI (2.0 * M_PI)                           

/** @def ROUND
 * @brief Precision factor used for rounding floating-point numbers.
 * Set to 1.0: rounding performed by adding 0.5 before the floor operation (for reaching the nearest interger),
 * without additional scaling.
 */
#define ROUND 1.0                                    

#endif // GLOBALS_H