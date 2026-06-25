#include "simulation_types.h"
#include "photoevaporation_wrapper.h"
#include "../extern/include/photoevaporation.hpp"
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

void computePhotoevaporationSink(void *disk_opaque)
{
    DiskParameters *disk = static_cast<DiskParameters*>(disk_opaque);

    int N = disk->grid_number; // A TE rácspontszámod
    double *r = disk->radial_grid;
    
    // Generálunk egy valódi dr tömböt a Delta_R-ből, mivel a szimulációd fix lépésközzel megy
    std::vector<double> dr_array(N, disk->delta_r);
    std::vector<double> local_evap(N, 0.0);

    // KIKÉNYSZERÍTJÜK, hogy az extern kód a te csillagtömegedet használja, 
    // ha a globális makró felülírható, vagy átadjuk a függvénynek:
    // (Ha az M_STAR makró, akkor a Func_C-ben nem tudod átírni, de a hole_flag-et igen!)

    // A TE szimulációdból származó lyuk adatok (ha bevezeted őket a simulation_types.h-ba):
    int hole_flag = 0;   // disk->hole_flag; 
    double r_hole = 0.0; // disk->r_hole;

    // 1. Kiszámoljuk a normát a TE sugárrácsoddal és a valódi dr_array-jel
    double norm = Norm(r, hole_flag, r_hole, dr_array.data());

    // 2. Futtatjuk a lefolyást
    Photoevaporation_2012(local_evap.data(), r, norm, hole_flag, r_hole, dr_array.data());

    // 3. Beírjuk a kiszámolt profilt a TE struktúrádba
    if (disk->sigma_dot_photoevap != nullptr) {
        for (int i = 0; i < N; i++) {
            // Beindexeljük a belső rácsra (i+1 a ghost cellák miatt, ha az init_tool-ból indulunk ki)
            disk->sigma_dot_photoevap[i + 1] = local_evap[i];
        }
    }
}

#ifdef __cplusplus
}
#endif