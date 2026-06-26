#include "simulation_types.h"
#include "photoevaporation_wrapper.h"
#include "../extern/include/photoevaporation.hpp"
#include <vector>
#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

void computePhotoevaporationSink(void *disk_opaque)
{
    DiskParameters *disk = static_cast<DiskParameters*>(disk_opaque);

    int N = disk->grid_number; // Physical grid size (excluding ghost cells)
    
    // Shift pointers by +1 to bypass the inner ghost cell and point to the actual physical disk
    double *physical_r = disk->radial_grid + 1;
    double *physical_sigma = disk->gas_surface_density_vector + 1;
    
    // Generate dr array using the simulation's spatial step size
    std::vector<double> dr_array(N, disk->delta_r);
    std::vector<double> local_evap(N, 0.0);

    // Dynamic hole tracking variables inside the wrapper (or you can add them to DiskParameters later)
    static int gap = 0;
    static int hole = 0;
    static double r_hole = 0.0;

    // 1. Search for a gap or an inner hole based on the current surface density profile
    Search_hole(physical_r, physical_sigma, gap, hole, r_hole, disk);

    // 2. Calculate the normalization factor using the active physical grid
    double norm = Norm(physical_r, hole, r_hole, dr_array.data(), disk);
//    std::cout << "WRAPPER CHECK: hole=" << hole << ", r_hole=" << r_hole << ", m_star=" << disk->stellar_mass << ", norm=" << norm << std::endl;

    // 3. Compute the photoevaporation profile into the local array
    Photoevaporation_2012(local_evap.data(), physical_r, norm, hole, r_hole, dr_array.data(), disk);

    // 4. Map the calculated sink profile back to your simulation's tracking array
    if (disk->sigma_dot_photoevap != nullptr) {
        for (int i = 0; i < N; i++) {
            // Synchronize index with your inner grid system (i+1 to account for the inner ghost cell)
            disk->sigma_dot_photoevap[i + 1] = local_evap[i];
        }
        
        // Extrapolate to boundary ghost cells to maintain smooth boundaries
        disk->sigma_dot_photoevap[0] = disk->sigma_dot_photoevap[1];
        disk->sigma_dot_photoevap[N + 1] = disk->sigma_dot_photoevap[N];
    }
}

#ifdef __cplusplus
}
#endif