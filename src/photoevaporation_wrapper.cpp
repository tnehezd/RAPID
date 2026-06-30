/**
 * @file photoevaporation_wrapper.cpp
 * @brief Multi-mode wrapper with on-the-fly logarithmic transformation to support
 * log-space photoevaporation models on a linear simulation grid.
 * @date 2026-06-26
 */

#include "simulation_types.h"
#include "photoevaporation_wrapper.h"
#include "../extern/include/photoevaporation.hpp"
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>

#ifdef __cplusplus
extern "C" {
#endif

#undef ACTUAL_SIGMA_CRIT
#define ACTUAL_SIGMA_CRIT 1.1e-15

void computePhotoevaporationSink(void *disk_opaque)
{
    DiskParameters *disk = static_cast<DiskParameters*>(disk_opaque);
    int N = disk->grid_number; // Linear grid size

    // 0. Early exit if photoevaporation is globally disabled
    if (!disk->enable_photoevaporation) {
        if (disk->sigma_dot_photoevap != nullptr) {
            std::fill(disk->sigma_dot_photoevap, disk->sigma_dot_photoevap + N + 2, 0.0);
        }
        return;
    }

    // =========================================================================
    // 1. GENERATE ANNA'S EQUIVALENT LOGARITHMIC GRID
    // =========================================================================
    double r_min = disk->r_min;
    double r_max = disk->r_max;
    
    std::vector<double> log_r(N + 1, 0.0);
    std::vector<double> log_dr(N, 0.0);
    
    double log_ratio = std::log(r_max / r_min);
    
    // Generate radius coordinates exactly like Anna's discinit/main log-grid
    for (int i = 0; i <= N; i++) {
        log_r[i] = r_min * std::exp((double)i / (double)(N + 1) * log_ratio);
    }
    // Generate the spatial step sizes (dr) for the log grid
    for (int i = 0; i < N; i++) {
        log_dr[i] = log_r[i+1] - log_r[i];
    }

    // =========================================================================
    // 2. INTERPOLATE LINEAR SIGMA TO THE LOGARITHMIC GRID
    // =========================================================================
    std::vector<double> log_sigma(N, 0.0);
    double *linear_r = disk->radial_grid + 1;
    double *linear_sigma = disk->gas_surface_density_vector + 1;

    for (int i = 0; i < N; i++) {
        double target_r = log_r[i];
        
        // Simple linear interpolation from your linear grid to the log grid point
        if (target_r <= linear_r[0]) {
            log_sigma[i] = linear_sigma[0];
        } else if (target_r >= linear_r[N-1]) {
            log_sigma[i] = linear_sigma[N-1];
        } else {
            // Find cell in linear grid
            int idx = (int)((target_r - disk->r_min) / disk->delta_r);
            if (idx < 0) idx = 0;
            if (idx >= N - 1) idx = N - 2;
            
            double r_low = linear_r[idx];
            double r_high = linear_r[idx+1];
            double weight = (target_r - r_low) / (r_high - r_low);
            log_sigma[i] = (1.0 - weight) * linear_sigma[idx] + weight * linear_sigma[idx+1];
        }
    }

    // =========================================================================
    // 3. EXECUTE PHOTOEVAPORATION IN LOGARITHMIC SPACE
    // =========================================================================
    std::vector<double> log_evap(N, 0.0);
    static int gap = 0;
    static int hole = 0;
    static double r_hole = 0.0;

    // Search hole using the log-spaced coordinates and interpolated density
    // 1. Monitor gap/hole formation dynamically
    Search_hole(log_r.data(), log_sigma.data(), gap, hole, r_hole, disk);

    // DEBUG LOG TO CATCH THE 20,000 YR CRASH
/*    if (gap == 1 || hole == 1) {
        std::cerr << "[WRAPPER DETECTED HOLE]: gap=" << gap 
                  << ", hole=" << hole 
                  << ", r_hole=" << r_hole << " AU" << std::endl;
        if (hole == 1 && r_hole == 0.0) {
            std::cerr << "--> CRITICAL ERROR: Hole detected but r_hole is 0.0! Outer radius interpolation will explode." << std::endl;
        }
    }
*/
    std::string model_mode(disk->photoevaporation_mode_string);
    std::transform(model_mode.begin(), model_mode.end(), model_mode.begin(), [](unsigned char c){ return std::tolower(c); });

    if (model_mode == "owen") 
    {
        double norm = Norm(log_r.data(), hole, r_hole, log_dr.data(), disk);
        Photoevaporation_2012(log_evap.data(), log_r.data(), norm, hole, r_hole, log_dr.data(), disk);
    } 
    else if (model_mode == "picogna") 
    {
        double lx_cgs = disk->xray_luminosity; 
        New_Photoevaporation(log_evap.data(), log_r.data(), lx_cgs, hole, r_hole, log_dr.data(), disk);
    }
    else 
    {
        for (int i = 0; i < N; i++) log_evap[i] = 0.0;
    }

    // =========================================================================
    // 4. INTERPOLATE CALCULATED RATES BACK TO YOUR LINEAR RECIPIENT GRID
    // =========================================================================
    std::vector<double> local_evap_linear(N, 0.0);
    const double days_per_year = 365.242199;

    for (int i = 0; i < N; i++) {
        double target_linear_r = linear_r[i];
        
        // Find where your linear R falls inside Anna's log r array
        if (target_linear_r <= log_r[0]) {
            local_evap_linear[i] = log_evap[0];
        } else if (target_linear_r >= log_r[N-1]) {
            local_evap_linear[i] = log_evap[N-1];
        } else {
            // Find bounding indices in the log grid via inverse function
            double fractional_index = std::log(target_linear_r / r_min) / log_ratio * (N + 1);
            int idx = (int)std::floor(fractional_index);
            if (idx < 0) idx = 0;
            if (idx >= N - 1) idx = N - 2;
            
            double r_low = log_r[idx];
            double r_high = log_r[idx+1];
            double weight = (target_linear_r - r_low) / (r_high - r_low);
            local_evap_linear[i] = (1.0 - weight) * log_evap[idx] + weight * log_evap[idx+1];
        }
    }

    // 5. Map the final linear-spaced profile back to your tracking array
    if (disk->sigma_dot_photoevap != nullptr) {
        for (int i = 0; i < N; i++) {
            disk->sigma_dot_photoevap[i + 1] = local_evap_linear[i];
        }
        
        disk->sigma_dot_photoevap[0] = disk->sigma_dot_photoevap[1];
        disk->sigma_dot_photoevap[N + 1] = disk->sigma_dot_photoevap[N];
    }
}

#ifdef __cplusplus
}
#endif