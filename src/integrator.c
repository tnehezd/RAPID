
// src/integrator.h

// Standard C Library Includes
#include <stdio.h>    // For printf, fopen, fclose, fscanf, snprintf, sprintf
#include <stdlib.h>   // For exit, EXIT_FAILURE, EXIT_SUCCESS, system
#include <math.h>     // For M_PI, fmod, HUGE_VAL (and pow if used by other functions)
#include <string.h>   // For snprintf, sprintf

#include <omp.h>

// Your Project Header Includes
#include "config.h"       // For particle_number, TMAX, WO, r_min, DT, optdr, sim_opts->twopop, sim_opts->growth, optev, r_dze_i, r_dze_o
#include "io_utils.h"     
#include "disk_model.h"   
#include "dust_physics.h" 
#include "utils.h"        
#include "simulation_core.h"
#include "particle_data.h" // Új include
#include "gas_physics.h"
#include "boundary_conditions.h"
#include "integrator.h"



/*	Runge-Kutta4 integrator	*/
// prad bemenet: AU-ban!
void integrateParticleRungeKutta4(double time, double prad, const double *sigmad, const double *rdvec, double step, double y, double *ynew, double *pradnew, const DiskParameters *disk_params, const simulation_options_t *sim_opts){
    double dy1,dy2,dy3,dy4;
    double ytemp, ytemp2;
    double sigma, dpress, ugas; 
    double pdens, p;
    double pradtemp;
    int opt = 0;
    double sigmadd = 0.0;
    
/*	Mivel a kulongozo parametereket csak a megadott gridcella pontokban ismerjuk, de ez nem feltetlen egyezik meg a reszecskek pozicijaval, ezert minden fontos parametert linearInterpolationalunk a reszecskek tavolsagara	*/
    linearInterpolation(disk_params->sigmavec,disk_params->rvec,y,&sigma,disk_params->delta_r,opt,disk_params);
    linearInterpolation(disk_params->dpressvec,disk_params->rvec,y,&dpress,disk_params->delta_r,opt,disk_params);
    linearInterpolation(disk_params->ugvec,disk_params->rvec,y,&ugas,disk_params->delta_r,opt,disk_params);

    double dd = (disk_params->r_max - disk_params->r_min) / (particle_number-1);
    int dker = (int)(1./dd);//
    dker = dker * ROUNDING_FACTOR;
    double ddker = (double) dker;
    int temp;

    temp = (int)floor(y * ddker+0.5);
    ytemp2 = (double)temp / ddker;
    
    int i;
    for(i=0;i<particle_number;i++) {
        if(ytemp2 == rdvec[i]) {
            sigmadd = sigmad[i];
            break;
        }
    }

    if(sim_opts->growth == 1.) {		// ha van reszecskenovekedes
        if(time != 0.) {	// ha nem t0 idopontban vagyunk
            pradtemp = prad;
            linearInterpolation(disk_params->pressvec,disk_params->rvec,y,&p,disk_params->delta_r,opt,disk_params);
            pdens = disk_params->particle_density; 
            pradtemp = calculateDustParticleSize(prad,pdens,sigma,sigmadd,y,p,dpress,step,disk_params);	// itt szamolja a reszecskenovekedest
            prad = pradtemp;
        }
    }

    *pradnew = prad;

/*	Itt szamolja a reszecske poziciojat	*/
    calculate1DDustDrift(prad, dpress, sigma, ugas, y, &dy1,disk_params);

    ytemp = y + 0.5 * step * dy1;
    calculate1DDustDrift(prad, dpress, sigma, ugas, ytemp, &dy2,disk_params);
        
    ytemp = y + 0.5 * step * dy2;
    calculate1DDustDrift(prad, dpress, sigma, ugas, ytemp, &dy3,disk_params);
    
    ytemp = y + step * dy3;
    calculate1DDustDrift(prad, dpress, sigma, ugas, ytemp, &dy4,disk_params);

    *ynew = y + step * (dy1 + 2.0 * dy2 + 2.0 * dy3 + dy4) / 6.0;

}



