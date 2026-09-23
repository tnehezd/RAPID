#include <math.h>
#include <string.h>
#include "photoevaporation.h"
#include "photoevaporation_coeffs.h"
#include "simulation_types.h"
#include "logger.h"

static const double DAYS_PER_YEAR = 365.25;
static const double ACTUAL_SIGMA_CRIT = 1e-8;

/* ============================================================================
   HOLE / GAP DETECTION
   ============================================================================ */
void PhotoevapSearchHole(const double *radius_array,
                          const double *sigma_array,
                          DiskParameters *disk)
{
    int N = disk->grid_number;

    int *gap  = &disk->gap_flag;
    int *hole = &disk->hole_flag;
    double *r_hole = &disk->r_hole;

    if (*gap == 0) {
        for (int i = 1; i < N - 1; i++) {
            if (sigma_array[i] <= ACTUAL_SIGMA_CRIT) {
                *gap = 1;
                break;
            }
        }
    }

    if (*gap == 1) {

        if (*hole == 0) {

            if (sigma_array[1] < ACTUAL_SIGMA_CRIT ||
                sigma_array[2] < ACTUAL_SIGMA_CRIT) {
                *hole = 1;
            }

            int alive_cells = 0;
            for (int i = 1; i < N - 1; i++) {
                if (sigma_array[i] > ACTUAL_SIGMA_CRIT) alive_cells++;
            }
            if (alive_cells == 0) *hole = 1;
        }

        if (*hole == 1) {

            double current_r_hole = 0.0;
            int found = 0;

            for (int i = 1; i < N - 1; i++) {
                double L = sigma_array[i - 1] - ACTUAL_SIGMA_CRIT;
                double R = sigma_array[i]     - ACTUAL_SIGMA_CRIT;

                if (L * R <= 0.0) {
                    current_r_hole = 0.5 * (radius_array[i - 1] + radius_array[i]);
                    found = 1;
                    break;
                }
            }

            if (!found) {
                for (int i = 0; i < N; i++) {
                    if (sigma_array[i] > ACTUAL_SIGMA_CRIT) {
                        current_r_hole = radius_array[i];
                        found = 1;
                        break;
                    }
                }
            }

            if (found &&
                current_r_hole > radius_array[0] &&
                current_r_hole < radius_array[N - 1]) {
                *r_hole = current_r_hole;
            } else {
                if (*r_hole <= 0.0) *hole = 0;
            }
        }
    }
}

/* ============================================================================
   OWEN 2012 – NORMALIZATION
   ============================================================================ */
double PhotoevapNormOwen2012(const double *radius_array,
                               int hole,
                               double r_hole,
                               const double *dr_array,
                               const DiskParameters *disk)
{
    int N = disk->grid_number;
    double M_star = disk->stellar_mass;
    double Lx     = disk->xray_luminosity;

    double sum = 0.0;
    double M_dot = 0.0;

    if (hole == 0) {

        M_dot = 6.25e-9 * pow(M_star, -0.068) * pow(Lx / 1e30, 1.14);

        for (int i = 1; i <= N; i++) {
            double r = radius_array[i];
            double x = 0.85 * r / M_star;

            double func_c = 0.0;
            if (x > 0.7) {
                double logx10 = log10(x);

                double P1 = OWEN_A1 * pow(logx10, 6) +
                            OWEN_B1 * pow(logx10, 5) +
                            OWEN_C1 * pow(logx10, 4);

                double P2 = OWEN_D1 * pow(logx10, 3) +
                            OWEN_E1 * pow(logx10, 2) +
                            OWEN_F1 * logx10 +
                            OWEN_G1;

                double F = pow(10.0, P1) * pow(10.0, P2);

                double lx   = log(x);
                double ln10 = log(10.0);

                double dFdx =
                    (6*OWEN_A1*pow(lx,5)/(x*x*pow(ln10,7))) +
                    (5*OWEN_B1*pow(lx,4)/(x*x*pow(ln10,6))) +
                    (4*OWEN_C1*pow(lx,3)/(x*x*pow(ln10,5))) +
                    (3*OWEN_D1*pow(lx,2)/(x*x*pow(ln10,4))) +
                    (2*OWEN_E1*lx       /(x*x*pow(ln10,3))) +
                    (OWEN_F1            /(x*x*pow(ln10,2)));

                func_c = F * dFdx * exp(-pow(x / 100.0, 10.0));
            }

            sum += 2.0 * M_PI * r * func_c * dr_array[i];
        }
    } else {

        M_dot = 4.8e-9 * pow(M_star, -0.148) * pow(Lx / 1e30, 1.14);

        for (int i = 1; i <= N; i++) {
            double r = radius_array[i];
            double y = 0.95 * (r - r_hole) / M_star;

            double func_c = 0.0;
            if (y >= 0.0) {
                func_c = (OWEN_A2*OWEN_B2*exp(OWEN_B2*y) +
                          OWEN_C2*OWEN_D2*exp(OWEN_D2*y) +
                          OWEN_E2*OWEN_F2*exp(OWEN_F2*y)) / r
                         * exp(-pow(y / 57.0, 10.0));
            }

            sum += 2.0 * M_PI * r * func_c * dr_array[i];
        }
    }

    if (sum <= 0.0) return 0.0;
    return M_dot / sum;
}

/* ============================================================================
   OWEN 2012 – LOCAL σ̇(r)
   ============================================================================ */
double PhotoevapFuncOwen2012(double r,
                            double norm,
                            int hole,
                            double r_hole,
                            const DiskParameters *disk)
{
    double M_star = disk->stellar_mass;
    double func_c = 0.0;

    double a1 = OWEN_A1;
    double b1 = OWEN_B1;
    double c1 = OWEN_C1;
    double d1 = OWEN_D1;
    double e1 = OWEN_E1;
    double f1 = OWEN_F1;
    double g1 = OWEN_G1;

    double a2 = OWEN_A2;
    double b2 = OWEN_B2;
    double c2 = OWEN_C2;
    double d2 = OWEN_D2;
    double e2 = OWEN_E2;
    double f2 = OWEN_F2;

    if(hole == 0) {
        
        double x = 0.85*r/M_star;

        if(x > 0.7) {
                func_c = pow(10, a1 * pow(log10(x), 6) + b1 * pow(log10(x), 5) + c1 * pow(log10(x), 4))

                * pow(10 , d1 * pow(log10(x), 3) + e1 * pow(log10(x), 2) + f1 * log10(x)) * pow (10, g1)

                *(6*a1*pow(log(x),5)/(x*x*pow(log(10),7)) + 5*b1*pow(log(x),4)/(x*x*pow(log(10),6)) + 4*c1*pow(log(x),3)/(x*x*pow(log(10),5))

                + 3*d1*pow(log(x),2)/(x*x*pow(log(10),4)) + 2*e1*log(x)/(x*x*pow(log(10),3)) +  f1/(x*x*pow(log(10),2)))

                * exp((-1)*pow(x/100,10)) ;
        }
        else func_c = 0.0;

    } else if(hole == 1) {
        double y = 0.95*(r-r_hole)/M_star;

        if(y >= 0) func_c = (a2*b2*exp(b2*y) + c2*d2*exp(d2*y) + e2*f2*exp(f2*y))/r * exp((-1)*pow(y/57,10));
        else func_c = 0.0;
    }
  
    return(func_c*norm*DAYS_PER_YEAR);   //Astro unit [M_sol/AU^2/day]
}

/* ============================================================================
   OWEN 2012 – FULL PROFILE
   ============================================================================ */
void PhotoevaporationOwen2012(double *evap_array,
                               const double *radius_array,
                               double norm,
                               int hole,
                               double r_hole,
                               const double *dr_array,
                               const DiskParameters *disk)
{
    int N = disk->grid_number;
    (void)dr_array;

    evap_array[0] = 0.0;
    evap_array[N + 1] = 0.0;

    for (int i = 1; i <= N; i++) {
        double sigma_dot_yr = PhotoevapFuncOwen2012(radius_array[i],
                                                      norm,
                                                      hole,
                                                      r_hole,
                                                      disk);
        evap_array[i] = sigma_dot_yr / DAYS_PER_YEAR;
    }
}

/* ============================================================================
   PICOGNA 2019 – LOCAL σ̇(r)
   ============================================================================ */
static double calculatePicognaMassLossRate(double L_x)
{
    double log_lx = log10(L_x);
    if (log_lx <= 0.0) return 0.0;

    return pow(10.0,
               PICO_A_L * exp(pow(log(log_lx) - PICO_B_L, 2.0) /
                              PICO_C_L) +
               PICO_D_L);
}

double PhotoevapFuncPicogna2019(double r, double L_x, int hole, double r_hole)
{
    if (r <= 0.0) return 0.0;

    double m_lx_dot = calculatePicognaMassLossRate(L_x);
    double flux = 0.0;

    if (hole == 0) {
        double x = log10(r);
        double P = PICO_A * pow(x, 6) + PICO_B * pow(x, 5) +
                   PICO_C * pow(x, 4) + PICO_D * pow(x, 3) +
                   PICO_E * pow(x, 2) + PICO_F * x + PICO_G;
        double log_formula =
            6.0 * PICO_A * pow(log(r), 5) / pow(log(10.0), 5) +
            5.0 * PICO_B * pow(log(r), 4) / pow(log(10.0), 4) +
            4.0 * PICO_C * pow(log(r), 3) / pow(log(10.0), 3) +
            3.0 * PICO_D * pow(log(r), 2) / pow(log(10.0), 2) +
            2.0 * PICO_E * log(r) / log(10.0) + PICO_F;
        double wind_mass_profile = m_lx_dot * pow(10.0, P);
        flux = log_formula * wind_mass_profile /
               (2.0 * M_PI * r * r);
    } else {
        double x = r - r_hole;
        if (x > 0.0) {
            flux = PICO_AA * pow(PICO_BB, x) * pow(x, PICO_CC - 1.0) *
                   (x * log(PICO_BB) + PICO_CC) * 1.12 * m_lx_dot /
                   (2.0 * M_PI * r);
        }
    }

    return (flux < 1e-20) ? 0.0 : flux;
}

/* ============================================================================
   PICOGNA 2019 – FULL PROFILE
   ============================================================================ */
void PhotoevaporationPicogna2019(double *evap_array,
                                  const double *radius_array,
                                  double L_x,
                                  const double *dr_array,
                                  const DiskParameters *disk)
{
    int N = disk->grid_number;

    evap_array[0] = 0.0;
    evap_array[N + 1] = 0.0;

    for (int i = 1; i <= N; i++) {
        evap_array[i] = PhotoevapFuncPicogna2019(radius_array[i], L_x,
                                                 disk->hole_flag,
                                                 disk->r_hole);
    }

    (void)dr_array;
}

/* ============================================================================
   MAIN DISPATCHER – called from dust_physics.c
   ============================================================================ */
void computePhotoevaporationSink(DiskParameters *disk)
{
    int N = disk->grid_number;
 
    if (!disk->enable_photoevaporation || disk->sigma_dot_photoevap == NULL) {
        if (disk->sigma_dot_photoevap != NULL) {
            for (int i = 0; i < N + 2; i++) {
                disk->sigma_dot_photoevap[i] = 0.0;
            }
        }
        return;
    } 

    PhotoevapSearchHole(disk->radial_grid,
                         disk->gas_surface_density_vector,
                         disk);

    const char *mode = disk->photoevaporation_mode_string;
    

    if (mode == NULL || mode[0] == '\0' ||
        strcasecmp(mode, "none") == 0) {

        for (int i = 0; i < N + 2; i++) {
            disk->sigma_dot_photoevap[i] = 0.0;
        }
        return;
    }

    if (strcasecmp(mode, "owen") == 0) {

        double norm = PhotoevapNormOwen2012(disk->radial_grid,
                                              disk->hole_flag,
                                              disk->r_hole,
                                              disk->delta_r_array,
                                              disk);

        PhotoevaporationOwen2012(disk->sigma_dot_photoevap,
                                  disk->radial_grid,
                                  norm,
                                  disk->hole_flag,
                                  disk->r_hole,
                                  disk->delta_r_array,
                                  disk);

    } else if (strcasecmp(mode, "picogna") == 0) {

        double Lx_cgs = disk->xray_luminosity;

        PhotoevaporationPicogna2019(disk->sigma_dot_photoevap,
                                     disk->radial_grid,
                                     Lx_cgs,
                                     disk->delta_r_array,
                                     disk);

    } else {

        for (int i = 0; i < N + 2; i++) {
            disk->sigma_dot_photoevap[i] = 0.0;
        }
    }

    disk->sigma_dot_photoevap[0]     = disk->sigma_dot_photoevap[1];
    disk->sigma_dot_photoevap[N + 1] = disk->sigma_dot_photoevap[N];
}