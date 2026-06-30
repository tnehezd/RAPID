/**
 * @file simulation_types.h
 * @brief Core data structures describing disk physics and simulation settings.
 *
 * Defines the fundamental types used across the simulation, including
 * disk parameters, user‑configurable options, output file handles,
 * snapshot modes, and the PressureTrap structure used to identify
 * dust‑accumulation regions in the disk.
 */

#ifndef SIMULATION_TYPES_H
#define SIMULATION_TYPES_H

#include <stdio.h>
#include <stdbool.h>
#define MAX_PATH_LEN 16384 

/**
 * @brief Physical and numerical parameters describing the protoplanetary disk.
 *
 * This structure stores all disk‑related quantities used throughout the
 * simulation, including radial grid information, gas surface density,
 * pressure gradients, and viscosity parameters. Many of these values
 * are initialized from user input or configuration files.
 */
typedef struct {
    double r_min;                         /**< Inner radius of the simulation domain [AU]. */
    double r_max;                         /**< Outer radius of the simulation domain [AU]. */
    int    grid_number;                   /**< Number of radial grid points. */
    double delta_r;                       /**< Radial grid spacing [AU]. */
    double sigma_0;                       /**< Gas surface density normalization at 1 AU. */
    double sigma_power_law_index;         /**< Exponent of the gas surface density power law. */
    double alpha_parameter;               /**< Shakura–Sunyaev viscosity parameter α. */
    double stellar_mass;                  /**< Stellar mass [M_sun]. */
    double disk_mass;                     /**< Total disk mass [M_sun]. */
    double h_aspect_ratio;                /**< Disk aspect ratio H/R at reference radius. */
    double flaring_index;                 /**< Disk flaring index (temperature gradient). */
    double r_dze_i;                       /**< Inner dead‑zone edge radius [AU]. */
    double r_dze_o;                       /**< Outer dead‑zone edge radius [AU]. */
    double dr_dze_i;                      /**< Width of inner dead‑zone transition [AU]. */
    double dr_dze_o;                      /**< Width of outer dead‑zone transition [AU]. */
    double alpha_parameter_modification;  /**< α‑viscosity reduction factor inside dead zones. */
    double particle_density;              /**< Dust particle material density [g/cm³]. */
    double particle_density_dimensionless;/**< Dimensionless dust density (code units). */
    double *radial_grid;                  /**< Radial grid array [AU]. */
    double *gas_surface_density_vector;   /**< Gas surface density at each grid point. */
    double *gas_pressure_vector;          /**< Gas pressure at each grid point. */
    double *gas_pressure_gradient_vector; /**< Radial pressure gradient dP/dr. */
    double *gas_velocity_vector;          /**< Radial gas velocity profile. */
    double fragmentation_factor;          /**< Mass reduction factor after fragmentation. */
    double fragmentation_velocity;        /**< Fragmentation threshold velocity [m/s]. */
    double drift_factor;                  /**< Scaling factor for dust drift speed. */        
    bool cutoff;             // Replaces or adds the cutoff flag
    double r_cutoff;         // The characteristic tapering radius (Rc)
    double n_cutoff;         // The exponent factor (n) for sharpness
    double total_disk_mass;  // The clean name for your total disk mass variable    
    double *sigma_dot_photoevap;           /**< Sink term for photoevaporation */
    double *photoevap_sink;
    double *dr_array;
    int    hole_flag;
    double density_floor;   /**< Minimum allowed gas surface density */
    double r_hole;
    double *sink_photoevap;
    bool   enable_photoevaporation;       /**< Global switch to turn photoevaporation ON/OFF. */
    char   photoevaporation_mode_string[32]; /**< Model selection: "Owen" or "Picogna" */
    double xray_luminosity;               /**< Stellar X-ray luminosity [erg/s]. */
    double gaussian_sigma;
    double gaussian_cutoff;

} DiskParameters;


/**
 * @enum OutputFormat
 * @brief Supported output formats for simulation snapshots.
 *
 * This enumeration specifies the available output backends used when
 * writing simulation data to disk. The selected format determines how
 * gas and dust fields, particle properties, and metadata are stored.
 *
 * The naming follows the project's coding conventions
 * (see Coding Standard).
 */
typedef enum {

    /**
     * @brief Plain-text ASCII output.
     *
     * Human-readable column-based text files.  
     * Useful for quick inspection, debugging, and lightweight post-processing.
     */
    OUTPUT_ASCII = 0,

    /**
     * @brief HDF5 output format.
     *
     * Binary, hierarchical, self-describing file format suitable for
     * large datasets and high-performance workflows.  
     * Enables structured storage of gas fields, dust populations,
     * metadata, and multi-snapshot datasets.
     */
    OUTPUT_HDF5 = 1

} OutputFormat;


/**
 * @brief User‑defined simulation options and runtime configuration.
 *
 * This structure contains all high‑level simulation settings, including
 * which physical modules are enabled, time‑stepping parameters, and
 * input/output file paths.
 */
typedef struct {
    double option_for_evolution;              /**< Enable/disable gas evolution. */
    double option_for_dust_drift;             /**< Enable/disable dust drift module. */
    double option_for_dust_growth;            /**< Enable/disable dust growth module. */
    double option_for_dust_secondary_population; /**< Enable micron‑dust population. */
    double flag_for_deadzone;                 /**< Enable/disable dead‑zone physics. */
    double user_defined_time_step;            /**< Fixed time step (if > 0). */
    double maximum_simulation_time;           /**< Maximum simulation time [years]. */
    double output_frequency;                  /**< Snapshot output frequency. */
    int number_of_dust_particles;             /**< Number of dust particles to simulate. */
    char input_filename[MAX_PATH_LEN];        /**< Path to gas input file. */
    char output_dir_name[MAX_PATH_LEN];       /**< Output directory for results. */
    char dust_input_filename[MAX_PATH_LEN];   /**< Path to dust input file. */
    OutputFormat output_format;               /**< Selects ASCII or HDF5 output backend. */
    int dust_smoothing_mode;                  /**< Dust smoothing method for mapping Lagrangian to Euler grid */
    double gaussian_sigma;   /**< Gaussian kernel width in grid-cell units (Δr) */
    double gaussian_cutoff;       /**< Gaussian kernel cutoff in multiples of sigma */
} SimulationOptions;


/**
 * @brief File handles and metadata for simulation output.
 *
 * This structure centralizes all output-related resources used during
 * the simulation. It stores file pointers for ASCII output streams
 * as well as placeholders for HDF5 output handling.
 *
 * The structure allows unified initialization, writing, and cleanup
 * of output resources independent of the selected output format.
 *
 * @note
 * The codebase currently supports ASCII output and is being extended
 * with optional HDF5 support. To avoid forcing a hard dependency on
 * the HDF5 library in public headers, the HDF5 file handle is stored
 * as a generic pointer (`void *`). This will later be mapped to the
 * native `hid_t` type internally when HDF5 support is enabled.
 *
 * The `snapshot_index` member is format-independent and provides a
 * consistent mechanism for indexing time-dependent outputs. While
 * ASCII output typically appends data sequentially to files, HDF5
 * output will use this index to write into structured datasets.
 */
typedef struct {

    /* =========================
     * ASCII output format
     * ========================= */

    FILE *dust_motion_file;   /**< Drift velocities of dust particles. */
    FILE *micron_motion_file; /**< Drift velocities of micron-sized dust. */
    FILE *mass_file;          /**< Dust mass evolution output. */
    FILE *surface_file;       /**< Gas surface density and pressure output. */
    FILE *dust_file;          /**< Dust surface density output. */
    FILE *micron_dust_file;   /**< Micron dust surface density output. */
    FILE *size_file;          /**< Dust particle size distribution output. */

    /* =========================
     * HDF5 output handling
     * ========================= */

    void *hdf5_file;          /**< Opaque handle to the HDF5 file object.
                                   Currently stored as void* to avoid
                                   direct dependency on the HDF5 library
                                   in headers. Will later map to hid_t. */

    int snapshot_index;       /**< Sequential snapshot counter used for
                                   indexing time-dependent datasets. */

} OutputFiles;


/**
 * @brief Enumeration of available snapshot output modes.
 *
 * These modes determine which physical modules produce output during
 * the simulation (gas, drift, growth, two‑population models, etc.).
 */
typedef enum {
    SnapshotNonevolving     = 0,    /**< No snapshots written. */
    SnapshotGas             = 1,    /**< Gas‑only snapshots. */
    SnapshotDrift           = 2,    /**< Dust drift snapshots. */
    SnapshotGrowth          = 3,    /**< Dust growth snapshots. */
    SnapshotDriftTwoPop     = 4,    /**< Drift snapshots for two dust populations. */
    SnapshotGrowthTwoPop    = 5     /**< Growth snapshots for two dust populations. */

} SnapshotMode;

/**
 * @brief Determines the snapshot mode based on simulation options.
 *
 * @param sim_opts Pointer to SimulationOptions.
 * @return The selected SnapshotMode.
 */
SnapshotMode determineSnapshotMode(const SimulationOptions *sim_opts);

/**
 * @brief Converts a SnapshotMode value to a human‑readable string.
 *
 * @param mode Snapshot mode enumeration value.
 * @return String representation of the mode.
 */
const char* snapshotModeToString(SnapshotMode mode);

/**
 * @struct PressureTrap
 * @brief Data structure representing a gas pressure maximum and its dust accumulation.
 *
 * Pressure maxima in protoplanetary disks act as "dust traps" by halting
 * inward drift of solids. This structure stores the trap location, its
 * integration boundaries, and the mass of dust accumulated within it.
 */
typedef struct {
    double radial_position;      /**< Radial location of the pressure maximum [AU]. */
    double inner_boundary;       /**< Inner integration boundary [AU]. */
    double outer_boundary;       /**< Outer integration boundary [AU]. */
    double primary_dust_mass;    /**< Mass of cm‑sized dust in the trap [M_sun]. */
    double secondary_dust_mass;  /**< Mass of micron‑sized dust in the trap [M_sun]. */
    double total_dust_mass;      /**< Total dust mass (primary + secondary) [M_sun]. */
    int trap_id;                 /**< Unique identifier for the trap. */                     
} PressureTrap;


typedef enum {
    SMOOTHING_CIC = 0,
    SMOOTHING_NGP = 1,
    SMOOTHING_TOPHAT = 2,
    SMOOTHING_GAUSSIAN = 3
} DustSmoothingMode;




#endif // SIMULATION_TYPES_H