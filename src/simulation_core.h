// src/simulation_core.h

#ifndef SIMULATION_CORE_H
#define SIMULATION_CORE_H

// Include any headers needed *by this header itself* for types used in declarations.
// For tIntegrate, it uses 'char' and 'double', which are built-in, so no specific
// includes are strictly necessary *just for the prototype*. However, if this header
// were to define structs or other types that are part of the public interface,
// you'd include their definitions here.

// Function Declarations
void tIntegrate(char *nev, double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double *ugvec);

// If you add other functions to simulation_core.c that are called elsewhere,
// declare them here as well.

#endif // SIMULATION_CORE_H
