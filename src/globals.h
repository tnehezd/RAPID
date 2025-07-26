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
/** @def SDCONV
 * @brief Surface density conversion factor.
 */
#define SDCONV              1.12521e-7                   

/** @def ICEFACTOR
 * @brief Factor for dust density beyond the snowline.
 * This typically accounts for the condensation of ice beyond a certain radial distance.
 */
#define ICEFACTOR           3.0                          

/** @def SNOWLINE
 * @brief Snowline radius in Astronomical Units (AU).
 * This defines the radial distance where ice condensation becomes significant.
 */
#define SNOWLINE            2.7                          

/** @def G_GRAV_CONST
 * @brief Gravitational constant (G=1) in the chosen simulation units.
 */
#define G_GRAV_CONST 1.0                             

/** @def G_GRAV_CONST2
 * @brief Square of the gravitational constant (G^2).
 */
#define G_GRAV_CONST2 (G_GRAV_CONST * G_GRAV_CONST)   

/** @def AUPDAY2CMPSEC
 * @brief Conversion factor from Astronomical Units per Day to centimeters per second.
 */
#define AUPDAY2CMPSEC       1.7314568e8                  

/** @def CMPSECTOAUPYRP2PI
 * @brief Conversion factor from centimeters per second to AU per (year / 2*pi).
 * This often relates to angular velocity or orbital periods.
 */
#define CMPSECTOAUPYRP2PI   3.35725e-07                    

/** @def GRPCM32MSUNAU3
 * @brief Conversion factor from grams per cubic centimeter to solar masses per cubic AU.
 */
#define GRPCM32MSUNAU3      1.68329e6                    

/** @def SUN2GR
 * @brief Solar Mass in grams (M_solar -> g).
 * @warning PLEASE VERIFY THIS VALUE! It's critical for accurate mass conversions.
 */
#define SUN2GR 1.989e33                              

/** @def AU2CM
 * @brief Astronomical Unit in centimeters (AU -> cm).
 * @warning PLEASE VERIFY THIS VALUE! It's critical for accurate distance conversions.
 */
#define AU2CM 1.496e13                               

/** @def TWOPI
 * @brief Represents 2 * PI (two times Pi).
 * Parentheses are used for safety in expressions.
 */
#define TWOPI (2.0 * M_PI)                           

/** @def ROUND
 * @brief Precision factor used for rounding floating-point numbers.
 * This is typically used to align continuous values to a discrete grid,
 * such as mapping 'y' values to specific 'rdvec' elements.
 * A value of 1.0 might indicate no explicit scaling for rounding,
 * or that rounding logic is handled by adding 0.5 before floor.
 */
#define ROUND 1.0                                    

#endif // GLOBALS_H