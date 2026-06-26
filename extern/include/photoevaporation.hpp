#pragma once // if everything goes well, this header will be included once with this line

// photoevaporation.hpp

#include "config.h"

#define EVAPORATION_ON 0      // 0: no photoevaporation, 1: Owen et al. 2012 X-ray photoevaporation, 2: Picogna et al. 2019 X-ray photoevaporation
#define WIND_CONV 1.0     // Conversion factor for wind profile output [from M_sol/AU^2/day to desired unit]
#define WIND_CRIT 1e-20  // Critical value for wind profile (below this value, wind is considered zero) [M_sol/AU^2/day]
//-------------------Global parameters----------------//
constexpr double L_X = 1e30;     // X-ray luminosity of star [erg/sec]
constexpr double AMAX = -10000.0; // Maximum value for Timestep_calculation

//-------------------X-ray photoevaporation (Owen et al 2012)----------------//
namespace owen2012 {
    constexpr double a1 = 0.15138;
    constexpr double a2 = -0.438226;
    constexpr double b1 = -1.2182;
    constexpr double b2 = -0.10658387;
    constexpr double c1 = 3.4046;
    constexpr double c2 = 0.5699464;
    constexpr double d1 = -3.5717;
    constexpr double d2 = 0.010732277;
    constexpr double e1 = -0.32762;
    constexpr double e2 = -0.131809597;
    constexpr double f1 = 3.6064;
    constexpr double f2 = -1.32285709;
    constexpr double g1 = -2.4918;
}

//-------------------New photoevaporation (Picogna et al 2019)----------------//
namespace picogna2019 {
    constexpr double a = -0.5885;
    constexpr double b = 4.313;
    constexpr double c = -12.1214;
    constexpr double d = 16.3587;
    constexpr double e = -11.4721;
    constexpr double f = 5.7248;
    constexpr double g = -2.8562;
    constexpr double aa = 0.11843;
    constexpr double bb = 0.99695;
    constexpr double cc = 0.48835;
    constexpr double A_L = -2.7326;
    constexpr double B_L = 3.3307;
    constexpr double C_L = -0.0029868;
    constexpr double D_L = -7.258;
}


//----------------------------------  Photoevaporation 2012 Owen et al.  -----------------------------------------//
/// @brief For photoevaporation you have to norm Func_C ---> int(2*pi*r*sigma_dot*dr) = M_dot must be fulfilled !!
/// @param radius_array Radius array of the disk [AU]
/// @param hole      0 -> no hole, 1 -> hole exists
/// @param r_hole    Radius of the hole [AU]
/// @return norm     Unitless M_dot/sum, Sigma_dot -> Sigma_dot*norm
/// @def L_X     X-ray luminosity of star [1e30 erg/sec] ---- from include/photoevaporation.hpp
/// @def M_STAR  Mass of the star [M_sol] ---- from include/defines.hpp
/// @def NGRID    Number of grid points ---- from include/defines.hpp
/// @def DR       Radial step size [AU] ---- from include/defines.hpp
/// @def a1,b1,c1,d1,e1,f1,g1  Coefficients for no hole case ---- from include/photoevaporation.hpp
/// @def a2,b2,c2,d2,e2,f2  Coefficients for hole case ---- from include/photoevaporation.hpp
double Norm(double *radius_array, int hole, double r_hole, double *dr_array, const DiskParameters *disk_params);

/// @brief Counting sigma_dot for X-ray photoevaporation
/// @param r         Radius where sigma_dot is calculated [AU]
/// @param norm      Normalization factor from Norm function
/// @param hole      0 -> no hole, 1 -> hole exists
/// @param r_hole    The hole extends up until this radius [AU]
/// @return func_c   Sigma_dot at radius r [M_sol/AU^2/day]
double Func_C2012(double r, double norm, int hole, double r_hole, const DiskParameters *disk_params);

/// @brief Creating array for disk wind values for easier calculasion
/// @param evap_array    disk wind values --- Astro unit: [M_sol/AU^2/yr]
/// @param radius_array  Spatial grid --- Astro unit: [AU]
/// @param norm          Normalizing factor, used by "Func_C2012" --- dimensionless
/// @param b_hole        Bool for hole, used by "Func_C2012"
/// @param r_hole        Radius of hole, used by "Func_C2012" --- The hole extends up until this radius in [AU]
void Photoevaporation_2012(double *evap_array, double *radius_array, double lx, int b_hole, double r_hole, double *dr_array, const DiskParameters *disk_params);

void Debug_Photoevaporation_2012(double *evap_array, double *radius_array, double norm, int b_hole, double r_hole, double *dr_array, const DiskParameters *disk_params);
//----------------------------------  Photoevaporation 2019 Picogna et al.  -------------------------------------//
/// @brief Counting sigma_dot for new photoevaporation profile from Picogna et al. 2019
/// @param r         Radius where sigma_dot is calculated [AU]
/// @param L_x       X-ray luminosity of the star [erg/sec]--- can be constant from defines.hpp or time dependent
/// @param b_hole    bool of hole existence
/// @param r_hole    The hole extends up until this radius [AU]
/// @return func_c   Sigma_dot at radius r [M_sol/AU^2/day
double Func_C2019(double r, double L_x, bool b_hole, double r_hole);

/// @brief Calculates new photoevaporation profile into an array based on Picogna et al. 2019, usint Func_C2019
/// @param evap_array Photoevaporation array to be filled
/// @param radius_array Radius array of the disk [AU] 
/// @param lx x-ray luminosity of the star [erg/sec] --- can be constant from defines.hpp or time dependent
/// @param b_hole bool of hole existence
/// @param r_hole  The hole extends up until this radius [AU]
/// @param dr Radial step size [AU] --- can be constant from defines.hpp or time dependent
void New_Photoevaporation(double *a_evap, double *a_radius, double lx, bool b_hole, double r_hole, double *dr_array, const DiskParameters *disk_params);


void Search_hole(double *radius_array, double *sigma_array, int &gap, int &hole, double &r_hole, const DiskParameters *disk_params);