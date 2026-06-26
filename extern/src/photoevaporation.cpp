#ifndef __APPLE__
#define _GNU_SOURCE
#endif

struct WindParams{
  double norm;
  int hole;
  double Rhole;
};

#include "../include/photoevaporation.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

// CRITICAL FIX: Changed from 1e-10 to a realistic gas surface density threshold 
// matched to your code's M_sun/AU^2 unit system (e.g., 1e-8 or your WIND_CRIT).
#define ACTUAL_SIGMA_CRIT 1e-8  
#define DEBUG 0

using namespace owen2012;   
using namespace picogna2019; 


//----------------------------------  Photoevaporation 2012 Owen et al.  -----------------------------------------//
/// @brief Calculates the normalization factor for photoevaporation mass-loss rate  ---> int(2*pi*r*sigma_dot*dr) = M_dot must be fulfilled !!
/// @param radius_array Pointer to the array of disk radii [AU].
/// @param hole      Integer flag indicating if a hole exists: 0 = no hole, 1 = hole exists
/// @param r_hole    Radius of the hole [AU]
/// @return norm     The normalization factor (unitless): M_dot/sum, Sigma_dot -> Sigma_dot*norm
double Norm(double *radius_array, int hole, double r_hole, double *dr_array, const DiskParameters *disk_params){  
  double func_c = 0.0;
  double sum = 0.0;
  double norm = 0.0;  
  double M_dot = 0.0;

  if(hole == 0){  //Primordial disk - no hole
    
    M_dot = 6.25e-9*pow(disk_params->stellar_mass,-0.068)*pow(L_X/disk_params->xray_luminosity, 1.14);  // Total Mass-loss rate for primordial disk [M_sol/yr]

    for(int i = 0; i < disk_params->grid_number; i++ ){
      double x = 0.85*radius_array[i]/disk_params->stellar_mass;     // Radius must be in AU, M_star in M_sol. 
      
      if(x > 0.7){
        func_c = pow(10, a1 * pow(log10(x), 6) + b1 * pow(log10(x), 5) + c1 * pow(log10(x), 4))
            * pow(10 , d1 * pow(log10(x), 3) + e1 * pow(log10(x), 2) + f1 * log10(x)) * pow (10, g1)
            * (6*a1*pow(log(x),5)/(x*x*pow(log(10),7)) + 5*b1*pow(log(x),4)/(x*x*pow(log(10),6)) + 4*c1*pow(log(x),3)/(x*x*pow(log(10),5))
            + 3*d1*pow(log(x),2)/(x*x*pow(log(10),4)) + 2*e1*log(x)/(x*x*pow(log(10),3)) + f1/(x*x*pow(log(10),2)))
            * exp((-1)*pow(x/100,10)) ;
      }
      else func_c = 0.0;

      sum += 2*M_PI*radius_array[i]*func_c*dr_array[i]; // Counting integral [M_sol/yr]
    }

    norm = M_dot/sum; // Unitless
  } 

  else if(hole == 1){  //Transition disk - Hole exists
    M_dot = 4.8e-9*pow(disk_params->stellar_mass,-0.148)*pow(L_X/disk_params->xray_luminosity, 1.14);  // Total Mass-loss rate for transition disk [M_sol/yr]

    for(int i = 0; i < disk_params->grid_number; i++ ){
      double y = 0.95*(radius_array[i]-r_hole)/disk_params->stellar_mass;

      if(y >= 0) func_c = (a2*b2*exp(b2*y) + c2*d2*exp(d2*y) + e2*f2*exp(f2*y))/radius_array[i]* exp((-1)*pow(y/57,10));
      else       func_c = 0.0;

      sum += 2*M_PI*radius_array[i]*func_c*dr_array[i]; // Counting integral [M_sol/yr]
    }
    norm = M_dot/sum; 
  }  
 
  return norm; 
}

/// @brief Calculates sigma_dot for X-ray photoevaporation based on Owen et al. 2012
/// @param r         Radius where sigma_dot is calculated [AU]
/// @param norm      Normalization factor from Norm function
/// @param hole      0 -> no hole, 1 -> hole exists
/// @param r_hole    The hole extends up until this radius [AU]
/// @return func_c   Sigma_dot at radius r [M_sol/AU^2/day]
double Func_C2012(double r, double norm, int hole, double r_hole, const DiskParameters *disk_params) 
{
  double func_c = 0.0;  // Unit: [M_sol/AU^2/yr]

  if(hole == 0) 
  {
    double x = 0.85*r/disk_params->stellar_mass;

    if(x > 0.7)
    {
      func_c = pow(10, a1 * pow(log10(x), 6) + b1 * pow(log10(x), 5) + c1 * pow(log10(x), 4))
      * pow(10 , d1 * pow(log10(x), 3) + e1 * pow(log10(x), 2) + f1 * log10(x)) * pow (10, g1)
      *(6*a1*pow(log(x),5)/(x*x*pow(log(10),7)) + 5*b1*pow(log(x),4)/(x*x*pow(log(10),6)) + 4*c1*pow(log(x),3)/(x*x*pow(log(10),5))
      + 3*d1*pow(log(x),2)/(x*x*pow(log(10),4)) + 2*e1*log(x)/(x*x*pow(log(10),3)) +  f1/(x*x*pow(log(10),2)))
      * exp((-1)*pow(x/100,10)) ;
    }
    else func_c = 0.0;
  }
  
  else if(hole == 1) 
  {
    double y = 0.95*(r-r_hole)/disk_params->stellar_mass;

    if(y >= 0) func_c = (a2*b2*exp(b2*y) + c2*d2*exp(d2*y) + e2*f2*exp(f2*y))/r * exp((-1)*pow(y/57,10));
    else func_c = 0.0;
  }
 return(func_c*norm);//*DAYS_PER_YEAR_CONVERSION_FACTOR);   // Convert to astro unit [M_sol/AU^2/day]
}

/// @brief Creating array for disk wind values for easier calculation
void Photoevaporation_2012(double *evap_array, double *radius_array, double norm, int b_hole, double r_hole, double *dr_array, const DiskParameters *disk_params) {
  for(int i = 0; i < disk_params->grid_number; ++i) {
    evap_array[i] = Func_C2012(radius_array[i], norm, b_hole, r_hole, disk_params);
  }
  if(DEBUG == 1)Debug_Photoevaporation_2012(evap_array, radius_array, norm, b_hole, r_hole, dr_array, disk_params);
}

void Debug_Photoevaporation_2012(double *evap_array, double *radius_array, double norm, int b_hole, double r_hole, double *dr_array, const DiskParameters *disk_params) {
  double sum = 0.0;
  double M_dot = 6.25e-9*pow(disk_params->stellar_mass,-0.068)*pow(L_X/disk_params->xray_luminosity, 1.14);
  for(int i = 0; i < disk_params->grid_number; ++i) sum += 2*M_PI*radius_array[i]* evap_array[i]*dr_array[i];   
  double norm_check = 0.0;
  for(int i = 0; i < disk_params->grid_number; ++i) norm_check += 2*M_PI*radius_array[i]* Func_C2012(radius_array[i], norm, b_hole, r_hole,disk_params)*dr_array[i];
  std::cout << "Sum [M_sol/yr]: " << sum*YEARS_PER_DAY_CONVERSION_FACTOR <<  std::endl;
  std::cout << "M_dot [M_sol/yr]: " << M_dot << std::endl;
  std::cout << "Norm check [M_sol/day]: " << norm_check << std::endl;
}

//----------------------------------  Photoevaporation 2019 Picogna et al.  -------------------------------------//
/// @brief Calculates sigma_dot for the new photoevaporation profile from Picogna et al. 2019
/// @return func_c   Sigma_dot at radius r [M_sol/AU^2/day]
double Func_C2019(double r, double L_x, bool b_hole, double r_hole)
{
  double f_c = 0.0;  // Unit: [M_sol/AU^2/yr]
  double M_W_dot = 0.0;
  double m_lx_dot = pow( 10.0, A_L * exp( pow(log(log10(L_x))-B_L, 2) / C_L ) + D_L );
  double log_formula = ((6.0*a*pow(log(r), 5))/(pow(log(10), 5))) + ((5.0*b*pow(log(r), 4))/(pow(log(10), 4))) +
                       ((4.0*c*pow(log(r), 3))/(pow(log(10), 3))) + ((3.0*d*pow(log(r), 2))/(pow(log(10), 2))) +
                       ((2.0*e* log(r))    /     log(10))      +       f;

  // Primordial disk case
  if (!b_hole) {
    M_W_dot = m_lx_dot * pow(10.0, a*pow(log10(r),6) + b*pow(log10(r),5) + c*pow(log10(r),4) + d*pow(log10(r),3) + e*pow(log10(r),2) + f*log10(r) + g);
    f_c = log_formula * M_W_dot/(2.0*M_PI*r*r);
  }

  // Transition disk case (hole exists)
  else if(b_hole){
    double x = r-r_hole;
    f_c = aa*pow(bb,x)*pow(x,cc-1.0)*(x*log(bb) + cc)*1.12*m_lx_dot/(2.0*M_PI*r);
  }

  if(f_c < 1e-20) f_c = 0.0; 
  return(f_c*DAYS_PER_YEAR_CONVERSION_FACTOR); 
}

/// @brief Calculates new photoevaporation profile into an array based on Picogna et al. 2019
void New_Photoevaporation(double *evap_array, double *radius_array, double lx, bool b_hole, double r_hole, double *dr_array, const DiskParameters *disk_params)
{
  double sum = 0.0;
  for(int i = 0; i < disk_params->grid_number; ++i) evap_array[i] = Func_C2019(radius_array[i], lx, b_hole, r_hole);
  for(int i = 0; i < disk_params->grid_number; ++i) sum += 2*M_PI*radius_array[i]* evap_array[i]*dr_array[i];
  if(DEBUG == 1) std::cout << "Sum: " << sum << std::endl;
}

/// @brief Dynamically monitors and updates gap and inner cavity (hole) radii
void Search_hole(double *radius_array,
                 double *sigma_array,
                 int &gap,
                 int &hole,
                 double &r_hole,
                 const DiskParameters *disk_params)
{
  const int N = disk_params->grid_number;

  // ============================================================
  // 1. INITIAL GAP DETECTION
  // ============================================================
  if (gap == 0) {
    for (int i = 1; i < N - 1; i++) {

      if (sigma_array[i] <= ACTUAL_SIGMA_CRIT) {

        printf("\n====================================\n");
        printf(" GAP CREATED\n");
        printf(" r = %.3f AU\n", radius_array[i]);
        printf(" sigma = %.3e\n", sigma_array[i]);
        printf(" crit = %.3e\n", ACTUAL_SIGMA_CRIT);
        printf(" i = %d\n", i);
        printf("====================================\n");

        gap = 1;
        break;
      }
    }
  }

  // ============================================================
  // 2. ACTIVE CAVITY TRACKING
  // ============================================================
  else if (gap == 1) {

    // ------------------------------------------------------------
    // 2a. HOLE ACTIVATION
    // ------------------------------------------------------------
    if (hole == 0) {

      if (sigma_array[1] < ACTUAL_SIGMA_CRIT ||
          sigma_array[2] < ACTUAL_SIGMA_CRIT) {
        hole = 1;
      }

      int alive_cells = 0;
      for (int i = 1; i < N - 1; i++) {
        if (sigma_array[i] > ACTUAL_SIGMA_CRIT) {
          alive_cells++;
        }
      }

      if (alive_cells == 0) {
        hole = 1;
      }
    }

    // ------------------------------------------------------------
    // 2b. HOLE RADIUS TRACKING
    // ------------------------------------------------------------
    if (hole == 1) {

      double current_r_hole = 0.0;
      int found = 0;

      // --------------------------------------------------------
      // Algorithm A: sign-change boundary detection
      // --------------------------------------------------------
      for (int i = 1; i < N - 1; i++) {

        double sigma_left  = sigma_array[i - 1] - ACTUAL_SIGMA_CRIT;
        double sigma_right = sigma_array[i]     - ACTUAL_SIGMA_CRIT;

        if (sigma_left * sigma_right <= 0.0) {
          current_r_hole = 0.5 * (radius_array[i - 1] + radius_array[i]);
          found = 1;
          break;
        }
      }

      // --------------------------------------------------------
      // Algorithm B: fallback (first surviving gas cell)
      // --------------------------------------------------------
      if (!found) {
        for (int i = 0; i < N; i++) {
          if (sigma_array[i] > ACTUAL_SIGMA_CRIT) {
            current_r_hole = radius_array[i];
            found = 1;
            break;
          }
        }
      }

      // --------------------------------------------------------
      // VALIDATION (CRITICAL — AFTER COMPUTATION)
      // --------------------------------------------------------
      if (found) {
        if (current_r_hole <= radius_array[0] ||
            current_r_hole >= radius_array[N - 1]) {
          found = 0; // reject unphysical result
        }
      }

      // --------------------------------------------------------
      // COMMIT ONLY VALID VALUES
      // --------------------------------------------------------
      if (found) {
        r_hole = current_r_hole;
      } else {
        // safety: do NOT allow silent 0 AU hole in active state
        if (r_hole <= 0.0) {
          hole = 0; // rollback to avoid corrupt Picogna/Owen profile
        }
      }
    }
  }


}
