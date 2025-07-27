#ifndef SIMULATION_TYPES_H
#define SIMULATION_TYPES_H

#include <stdio.h>

/**
 * @file simulation_types.h
 * @brief Contains the definitions for core data structures used in the simulation.
 *
 * This header file defines structures for particles, dust populations,
 * disk parameters, simulation options, and output file pointers.
 * These structures encapsulate related data, promoting code organization
 * and reusability throughout the RAPID simulation project.
 */

#define MAX_PATH_LEN 16384 // Define a maximum length for file paths

/**
 * @brief Represents a single dust particle in the simulation.
 */
typedef struct {
    double r;        ///< Radial distance from star (AU)
    double prad;     ///< Particle radius (cm)
    double reppmass; ///< Representative particle mass (mass of the material this particle represents)
    // Add any other per-particle properties here, e.g.:
    // double velocity_r; ///< Radial velocity if you track it explicitly
    // int    index_in_grid; ///< Or calculated on the fly
} particle_t;

/**
 * @brief Represents a collection of particles of a certain type (e.g., cm-sized, micron-sized).
 */
typedef struct {
    particle_t *particles;   ///< Pointer to an array of particle_t structures
    int num_particles;       ///< Number of particles in this specific population

    // If needed, specific properties for the population itself:
    // double total_mass; ///< Total mass represented by this population
    // double initial_mass; ///< Initial total mass for comparison
} dust_population_t;

/**
 * @brief Encapsulates all global disk parameters and dynamically allocated arrays.
 *
 * This structure holds all physical constants, geometrical parameters,
 * grid properties, and pointers to arrays representing the disk's state.
 */
typedef struct disk_t {
    // Geometrical and grid parameters
    double RMIN;       ///< Inner boundary of the disk simulation domain (AU)
    double RMAX;       ///< Outer boundary of the disk simulation domain (AU)
    int NGRID;         ///< Number of grid cells
    double DD;         ///< Grid cell width (calculated based on RMIN, RMAX, NGRID)

    // Gas disk parameters
    double SIGMA0;     ///< Initial reference gas surface density (e.g., at 1 AU)
    double SIGMAP_EXP; ///< Exponent of the initial gas surface density power-law profile
    double alpha_visc; ///< Alpha viscosity parameter (Shakura-Sunyaev alpha)
    double STAR_MASS;  ///< Mass of the central star
    double DISK_MASS;  ///< Total mass of the disk
    double HASP;       ///< Aspect ratio (H/R) of the disk (scale height parameter)
    double FLIND;      ///< Flaring index of the disk (exponent of H ~ r^FLIND)

    // Dead zone (DZE) parameters
    double r_dze_i;    ///< Inner radius of the dead zone transition region (AU)
    double r_dze_o;    ///< Outer radius of the dead zone transition region (AU)
    double Dr_dze_i;   ///< Width of the inner dead zone transition region (AU)
    double Dr_dze_o;   ///< Width of the outer dead zone transition region (AU)
    double a_mod;      ///< Alpha modification factor within the dead zone

    // Dust parameters
    double PDENSITY;       ///< Physical density of dust particle material (e.g., g/cm^3)
    double PDENSITYDIMLESS;///< Dimensionless density of dust particle material (M_Sun/AU^3)



    // Dynamically allocated arrays (these are expected to be allocated in main or dedicated init function)
    double *rvec;           ///< Radial grid vector (centers of grid cells)
    double *sigmavec;       ///< Gas surface density at grid points
    double *pressvec;       ///< Gas pressure at grid points
    double *dpressvec;      ///< Gas pressure gradient (radial derivative) at grid points
    double *ugvec;          ///< Gas radial velocity at grid points
    double *sigmadustvec;   ///< Dust surface density at grid points
    // Note: The `rhopvec` from your original comments might refer to local midplane dust density,
    // which would be derived from `sigmadustvec` and disk scale height.
    // If it's a separate array of values, ensure its purpose is distinct.

    // Dust coagulation/fragmentation parameters
    double fFrag;           ///< Fragmentation efficiency factor (f_frag in Birnstiel 2012)
    double uFrag;           ///< Fragmentation velocity threshold (u_frag in Birnstiel 2012)
    double fDrift;          ///< Drift efficiency factor (f_D in Birnstiel 2012, usually ~0.55)

} disk_t;

/**
 * @brief Groups all simulation control options and flags.
 */
typedef struct {
    double evol;           ///< Gas evolution flag (1.0 for enabled, 0.0 for disabled)
    double drift;          ///< Particle drift flag (1.0 for enabled, 0.0 for disabled)
    double growth;         ///< Particle growth flag (1.0 for enabled, 0.0 for disabled)
    double twopop;         ///< Two-population simulation flag (1.0 for enabled, 0.0 for disabled)
    double dzone;          ///< Dead zone dynamics flag (1.0 for dynamic DZE, 0.0 for fixed DZE)
    double DT;             ///< User-defined fixed time step
    double TMAX;           ///< Maximum simulation time
    double WO;             ///< Write-out interval (output frequency, e.g., TMAX/WO)
    double TCURR;          ///< Current simulation time

    int num_dust_particles; ///< Total number of dust particles simulated

    char input_filename[MAX_PATH_LEN];    ///< Path to the general input parameter file
    char output_dir_name[MAX_PATH_LEN];   ///< Path to the main output directory
    char dust_input_filename[MAX_PATH_LEN]; ///< Path to the input file for dust particle initial properties.
} simulation_options_t;

/**
 * @brief Stores pointers to all open output files.
 *
 * This structure helps manage file I/O by centralizing all FILE pointers,
 * ensuring they can be accessed and properly closed across different modules.
 */
typedef struct {
    FILE *por_motion_file;    ///< File pointer for particle motion data (e.g., pormozgas.dat)
    FILE *micron_motion_file; ///< File pointer for micron-sized particle motion data (e.g., pormozgasmic.dat)
    FILE *mass_file;          ///< File pointer for mass output (e.g., mass.dat)
    FILE *surface_file;       ///< File pointer for gas surface density output (e.g., surface.dat)
    FILE *dust_file;          ///< File pointer for dust surface density output (e.g., dust.dat)
    FILE *micron_dust_file;   ///< File pointer for micron dust surface density output (e.g., dustmic.dat)
    FILE *time_scale_file;    ///< File pointer for various timescale outputs.
    // Add other file pointers here if you need more output files
} output_files_t;

#endif // SIMULATION_TYPES_H    