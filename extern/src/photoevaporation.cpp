#ifndef __APPLE__
#define _GNU_SOURCE
#endif

struct WindParams{
  double norm;
  int hole;
  double Rhole;
};

#include "../include/photoevaporation.hpp"
#include "../include/defines.hpp"
#include <cmath>
#include <iostream>
#include <fstream>


using namespace owen2012;   
using namespace picogna2019; 


//----------------------------------  Photoevaporation 2012 Owen et al.  -----------------------------------------//
/// @brief Calculates the normalization factor for photoevaporation mass-loss rate  ---> int(2*pi*r*sigma_dot*dr) = M_dot must be fulfilled !!
/// @param radius_array Pointer to the array of disk radii [AU].
/// @param hole      Integer flag indicating if a hole exists: 0 = no hole, 1 = hole exists
/// @param r_hole    Radius of the hole [AU]
/// @return norm     The normalization factor (unitless): M_dot/sum, Sigma_dot -> Sigma_dot*norm
/// @def L_X     X-ray luminosity of star [1e30 erg/sec] ---- from include/photoevaporation.hpp
/// @def M_STAR  Mass of the star [M_sol] ---- from include/defines.hpp
/// @def NGRID    Number of grid points ---- from include/defines.hpp
/// @def DR       Radial step size [AU] ---- from include/defines.hpp
/// @def a1,b1,c1,d1,e1,f1,g1  Coefficients for no hole case ---- from include/photoevaporation.hpp
/// @def a2,b2,c2,d2,e2,f2  Coefficients for hole case ---- from include/photoevaporation.hpp
double Norm(double *radius_array, int hole, double r_hole, double *dr_array){  //Unitless M_dot/sum = norm, Sigma_dot -> Sigma_dot*norm
  double func_c = 0.0;
  double sum = 0.0;
  double norm = 0.0;  
  double M_dot = 0.0;

  if(hole == 0){  //Primordial disk - no hole
    
    M_dot = 6.25e-9*pow(M_STAR,-0.068)*pow(L_X/1e30, 1.14);  //total Mass-loss rate for primordial disk [M_sol/yr]

    for(int i = 0; i < NGRID; i++ ){
      double x = 0.85*radius_array[i]/M_STAR;     //Radius must be in AU, M_star in M_sol. 
      
      if(x > 0.7){
        func_c = pow(10, a1 * pow(log10(x), 6) + b1 * pow(log10(x), 5) + c1 * pow(log10(x), 4))

            * pow(10 , d1 * pow(log10(x), 3) + e1 * pow(log10(x), 2) + f1 * log10(x)) * pow (10, g1)

            * (6*a1*pow(log(x),5)/(x*x*pow(log(10),7)) + 5*b1*pow(log(x),4)/(x*x*pow(log(10),6)) + 4*c1*pow(log(x),3)/(x*x*pow(log(10),5))

            + 3*d1*pow(log(x),2)/(x*x*pow(log(10),4)) + 2*e1*log(x)/(x*x*pow(log(10),3)) + f1/(x*x*pow(log(10),2)))

            * exp((-1)*pow(x/100,10)) ;
      }
      else func_c = 0.0;

      sum += 2*M_PI*radius_array[i]*func_c*dr_array[i]; //counting integral [M_sol/yr]
    }

    norm = M_dot/sum; //Unitless
  } //end if(hole == 0)

  else if(hole == 1){  //Transition disk - Hole exists
    M_dot = 4.8e-9*pow(M_STAR,-0.148)*pow(L_X/1e30, 1.14);  //total Mass-loss rate for transition disk [M_sol/yr]

    for(int i = 0; i < NGRID; i++ ){
      double y = 0.95*(radius_array[i]-r_hole)/M_STAR;

      if(y >= 0) func_c = (a2*b2*exp(b2*y) + c2*d2*exp(d2*y) + e2*f2*exp(f2*y))/radius_array[i]* exp((-1)*pow(y/57,10));
      else       func_c = 0.0;

      sum += 2*M_PI*radius_array[i]*func_c*dr_array[i]; //counting integral [M_sol/yr]
    }
    norm = M_dot/sum; //counting norm from integral
  }  //end if(hole == 1)
 
  //std::cout << "M_dot: " << M_dot << std::endl;
  //std::cout << "Sum * norm: " << sum*norm << std::endl;
  return norm; //unitless
}

/// @brief Calculates sigma_dot for X-ray photoevaporation based on Owen et al. 2012
/// @param r         Radius where sigma_dot is calculated [AU]
/// @param norm      Normalization factor from Norm function
/// @param hole      0 -> no hole, 1 -> hole exists
/// @param r_hole    The hole extends up until this radius [AU]
/// @return func_c   Sigma_dot at radius r [M_sol/AU^2/day]
double Func_C2012(double r, double norm, int hole, double r_hole) //Sigma_dot is in [M_sol/AU^2/yr] in the paper -> Astro unit: [M_sol/AU^2/day]
{
  double func_c = 0.0;  //unit : [M_sol/AU^2/yr]

  if(hole == 0) //no hole
  {
    double x = 0.85*r/M_STAR;

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
  
  else if(hole == 1) //hole exists
  {
    double y = 0.95*(r-r_hole)/M_STAR;

    if(y >= 0) func_c = (a2*b2*exp(b2*y) + c2*d2*exp(d2*y) + e2*f2*exp(f2*y))/r * exp((-1)*pow(y/57,10));
    else func_c = 0.0;
  }
  //if(func_c < 1e-20) func_c = 0.0; //cutoff for very small values
 return(func_c*norm*DAY2YEAR);   //Astro unit [M_sol/AU^2/day]
}

/// @brief Creating array for disk wind values for easier calculasion
/// @param evap_array    disk wind values --- Astro unit: M_sol/AU^2/yr
/// @param radius_array  Spatial grid --- Astro unit: AU
/// @param norm          Normalizing factor, used by "Func_C2012"
/// @param b_hole        Bool for hole, used by "Func_C2012"
/// @param r_hole        Radius of hole, used by "Func_C2012"
void Photoevaporation_2012(double *evap_array, double *radius_array, double norm, int b_hole, double r_hole, double *dr_array) {
  for(int i = 0; i < NGRID; ++i) {

    evap_array[i] = Func_C2012(radius_array[i], norm, b_hole, r_hole);

   // if(evap_array[i]<1e-20) evap_array[i] = 0.0;
  }
  if(DEBUG == 1)Debug_Photoevaporation_2012(evap_array, radius_array, norm, b_hole, r_hole, dr_array);
}

void Debug_Photoevaporation_2012(double *evap_array, double *radius_array, double norm, int b_hole, double r_hole, double *dr_array) {
  double sum = 0.0;
  double M_dot = 6.25e-9*pow(M_STAR,-0.068)*pow(L_X/1e30, 1.14);
  for(int i = 0; i < NGRID; ++i) sum += 2*M_PI*radius_array[i]* evap_array[i]*dr_array[i];   //Astro unit [M_sol/day]
  double norm_check = 0.0;
  for(int i = 0; i < NGRID; ++i) norm_check += 2*M_PI*radius_array[i]* Func_C2012(radius_array[i], norm, b_hole, r_hole)*dr_array[i];
  std::cout << "Sum [M_sol/yr]: " << sum*YEAR2DAY <<  std::endl;
  std::cout << "M_dot [M_sol/yr]: " << M_dot << std::endl;
  std::cout << "Norm check [M_sol/day]: " << norm_check << std::endl;
}

//----------------------------------  Photoevaporation 2019 Picogna et al.  -------------------------------------//
/// @brief Calculates sigma_dot for the new photoevaporation profile from Picogna et al. 2019
/// @param r         Radius where sigma_dot is calculated [AU]
/// @param L_x       X-ray luminosity of the star [erg/sec]--- can be constant from defines.hpp or time dependent
/// @param b_hole    bool of hole existence
/// @param r_hole    The hole extends up until this radius [AU]
/// @return func_c   Sigma_dot at radius r [M_sol/AU^2/day
double Func_C2019(double r, double L_x, bool b_hole, double r_hole)
{
  double f_c = 0.0;  //in [M_sol/AU^2/yr] in the paper -> Astro unit: [M_sol/AU^2/day]
  double M_W_dot = 0.0;
  double m_lx_dot = pow( 10.0, A_L * exp( pow(log(log10(L_x))-B_L, 2) / C_L ) + D_L );
  double log_formula = ((6.0*a*pow(log(r), 5))/(pow(log(10), 5))) + ((5.0*b*pow(log(r), 4))/(pow(log(10), 4))) +
                       ((4.0*c*pow(log(r), 3))/(pow(log(10), 3))) + ((3.0*d*pow(log(r), 2))/(pow(log(10), 2))) +
                       ((2.0*e*    log(r))    /     log(10))      +       f;

  //double norm = 2 / log_formula;
  //double norm = 1.0; 

  //no hole/////
  if (!b_hole) {
    M_W_dot = m_lx_dot * pow(10.0, a*pow(log10(r),6) + b*pow(log10(r),5) + c*pow(log10(r),4) + d*pow(log10(r),3) + e*pow(log10(r),2) + f*log10(r) + g);
    f_c = log_formula * M_W_dot/(2.0*M_PI*r*r);
  }

  //hole exists////
  else if(b_hole){
    double x = r-r_hole;
    f_c = aa*pow(bb,x)*pow(x,cc-1.0)*(x*log(bb) + cc)*1.12*m_lx_dot/(2.0*M_PI*r);
  }

    //return(M_W_dot); //for plotting article figure
    if(f_c < 1e-20) f_c = 0.0; 
    return(f_c*DAY2YEAR); //for plotting article figure
    //return(f_c*DAY2YEAR*norm*r*r*M_PI); //for ???? time evolution
}

/// @brief Calculates new photoevaporation profile into an array based on Picogna et al. 2019, usint Func_C2019
/// @param evap_array Photoevaporation array to be filled
/// @param radius_array Radius array of the disk [AU] 
/// @param lx x-ray luminosity of the star [erg/sec] --- can be constant from defines.hpp or time dependent
/// @param b_hole bool of hole existence
/// @param r_hole  The hole extends up until this radius [AU]
/// @param dr Radial step size [AU] --- can be constant from defines.hpp or time dependent
void New_Photoevaporation(double *evap_array, double *radius_array, double lx, bool b_hole, double r_hole, double *dr_array)
{
  double sum = 0.0;
  for(int i = 0; i < NGRID; ++i) evap_array[i] = Func_C2019(radius_array[i], lx, b_hole, r_hole);
  for(int i = 0; i < NGRID; ++i) sum += 2*M_PI*radius_array[i]* evap_array[i]*dr_array[i];
  if(DEBUG == 1) std::cout << "Sum: " << sum << std::endl;
}

void Search_hole(double *radius_array, double *sigma_array, int gap, int hole, double r_hole){
  int gap_index = 0;
  if(gap == 0){
    for(int i = 1; i < NGRID-1; ++i){
      if(sigma_array[i] <= SIGMA_CRIT){
        gap = 1;
        gap_index = i;
      } 
    }
  }
  else if(gap == 1){
    if(hole == 0){
      if(gap == 1 && sigma_array[2] < SIGMA_CRIT ) hole = 1;
      //-------------------------------------------------------
      int disk = 0;
      for(int i = 1; i < gap_index; ++i) {
          if(sigma_array[i] > SIGMA_CRIT) disk ++;
      }
      if(disk == 0) hole = 1;
    }
    else if (hole == 1){
      r_hole = 0.0;
      double sigma_temp = 0.0;
      int temp_for_hole = 0;
      for(int i = 2; i < NGRID-1; i++){
        sigma_temp = (sigma_array[i-1]-SIGMA_CRIT)*(sigma_array[i]-SIGMA_CRIT);
        if( sigma_temp < 0.0 && temp_for_hole == 0) {
          r_hole = (radius_array[i]+radius_array[i+1])/2.0;
          temp_for_hole = 1;
          break;
        }
      } 
    }
  }

}

