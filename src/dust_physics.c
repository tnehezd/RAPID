// src/dust_physics.c
#include "dust_physics.h" // A saját headerjét mindig includolni kell
#include "config.h"       // Szükséges lehet a globális konstansokhoz
#include "simulation_types.h"

#include "simulation_core.h"
#include "utils.h"
#include <stdio.h>
#include <math.h>

/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
double alpha_turb(double r) {
	double alpha_r;
	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));
	return alpha_r * alpha_visc;
}


/*	kiszamolja az adott reszecskehez tartozo Stokes szamot	*/
/*	St = rho_particle * radius_particle * PI / (2 * sigma)	*/
double Stokes_Number(double pradius, double sigma) {		/*	in the Epstein drag regime	*/

	return PDENSITYDIMLESS * pradius * M_PI / (2.0 * sigma);

} 

/*	Lokalis viszkozitas erteke	*/
double visc(double r){
	double nu;
	double cs, H;
	
	H = scale_height(r);
	cs = c_sound(r);

	nu = alpha_turb(r) * cs * H;
	return nu;
}

/*	local scale height	*/
double scale_height(double r) {

	return pow(r,1.+FLIND) * HASP;

}

/*	lokális kepleri sebesség	*/
double v_kep(double r) {
	return sqrt(G_GRAV_CONST * STAR / r); // Uses G_GRAV_CONST now
}

/*	lokalis kepleri korfrekvencia	*/
double kep_freq(double r) {
	return sqrt(G_GRAV_CONST * STAR / r / r / r); // Uses G_GRAV_CONST now
}

/*	local sound speed		*/
double c_sound(double r) {
	return kep_freq(r) * scale_height(r);
}


/*	Suruseg a midplane-ben	*/
double rho_mp(double sigma, double r) {
	return 1. / sqrt(2.0 * M_PI) * sigma / scale_height(r);
}

/* local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
double press(double sigma, double r) {
	return rho_mp(sigma,r) * c_sound(r) * c_sound(r);
}

/*	a nyomas derivaltja	*/
void dpress(double *dp, double *p) {
	int i;
	double ptemp, pvec[NGRID+1];

	for(i = 1; i <= NGRID; i++) {
		ptemp = (p[i+1] - p[i-1]) / (2.0 * DD);
		pvec[i] = ptemp;
	}
	
	for(i = 1; i <= NGRID; i++) {
		dp[i] = pvec[i];
	}
}	

/*	u_gas kiszamolasahoz eltarolt koefficiens	*/
double Coeff_3(double sigma, double r){
	return -1.0 * (3.0 / (sigma * sqrt(r)));
}

/*	u_gas = -3/(Sigma*R^0.5)*(d/dR)(nu*Sigma*R^0.5) kiszamolasa	*/
void u_gas(double *sigmavec, double *rvec, double *ug) {
	double tempug, ugvec[NGRID+2],ugvectemp[NGRID+1];
	int i;

	for(i = 0; i <= NGRID+1; i++) {
		ugvec[i] = sigmavec[i] * visc(rvec[i]) * sqrt(rvec[i]);
	}

	for(i = 1; i <= NGRID; i++) {
		tempug = (ugvec[i+1] - ugvec[i-1]) / (2.0 * DD);
		ugvectemp[i] = Coeff_3(sigmavec[i],rvec[i]) * tempug;
	}
	
	for(i = 1; i <= NGRID; i++) {
		ug[i] = ugvectemp[i];
	}
}

void GetMass(int n, double partmassind[][4],int indii, int indio, double tavi, double dzei, double *massiout,int indoi, int indoo, double tavo, double dzeo, double *massoout) {

	int i;
	double index_i;
	double massitemp = 0., massotemp = 0;

	if(optdze == 1) {
		for (i = 0; i < n; i++) {

/*	A partmassind tombben van eltarolva a reszecske indexe (ez az integralas soran allando, ezzel lehet abrazolni szinesben), az altala kepviselt reprezentativ tomeg (ez az integralas soran allando) es a grid cella indexe, amelyben eppen tartozkodik (ezt minden lepesben kiszamolja a program	*/
			index_i = partmassind[i][1];			/*	a gridcella indexe, amelyben a reszecske tartozkodik eppen	*/

/*	Ha a reszecse a kiszamolt korgyuruben tartozkodik, akkor a mass valtozohoz hozzaadjuk a tomeget --> igy tudjuk a tomegnovekedest kiszamolni	*/

			if(tavi != dzei) {

				if ((index_i >= (int)indii) && (index_i <= (int)indio) && (partmassind[i][3] == 0)) {

					partmassind[i][3] = 1;
					massitemp = massitemp + partmassind[i][0];
				}
			
			} else {
	
				if ((index_i >= (int)indii) && (index_i <= (int)indio)) {
					massitemp = partmassind[i][0] + massitemp;
				}
			}

			if(tavo != dzeo) {

				if ((index_i >= (int)indoi) && (index_i <= (int)indoo) && (partmassind[i][3] == 0))  {
					partmassind[i][3] = 1;
					massotemp = massotemp + partmassind[i][0];
				}

			} else {
	
				if ((index_i >= (int)indoi) && (index_i <= (int)indoo)) {
					massotemp = partmassind[i][0] + massotemp;
				}
			}

		}
	} else {

		for (i = 0; i < n; i++) {

/*	A partmassind tombben van eltarolva a reszecske indexe (ez az integralas soran allando, ezzel lehet abrazolni szinesben), az altala kepviselt reprezentativ tomeg (ez az integralas soran allando) es a grid cella indexe, amelyben eppen tartozkodik (ezt minden lepesben kiszamolja a program	*/
			index_i = partmassind[i][1];			/*	a gridcella indexe, amelyben a reszecske tartozkodik eppen	*/

/*	Ha a reszecse a kiszamolt korgyuruben tartozkodik, akkor a mass valtozohoz hozzaadjuk a tomeget --> igy tudjuk a tomegnovekedest kiszamolni	*/
		
			if(tavo != dzeo) {

				if ((index_i >= (int)indoi) && (index_i <= (int)indoo) && (partmassind[i][3] == 0))  {
					partmassind[i][3] = 1;
					massotemp = massotemp + partmassind[i][0];
				}

			} else {
	
				if ((index_i >= (int)indoi) && (index_i <= (int)indoo)) {
					massotemp = partmassind[i][0] + massotemp;
				}
			}
		}

	}

	*massiout = massitemp;
	*massoout = massotemp;
}



/*			BIRNSTIEL EL AL 2012			*/

//reprezentativ reszecske kezdeti meretenek meghatarozasa
// 1. radialis drift altal meghatarozott maximalis meret			--> kimenet cm-ben!
double a_drift(double sigmad, double r, double p, double dp, double rho_p) {

	double f_drift = 0.55;
	double Sigmad_cgs = sigmad / SDCONV;

	double vkep = v_kep(r);
	double vkep2 = vkep * vkep;
	double c_s = c_sound(r);
	double c_s2 = c_s * c_s;
	double dlnPdlnr = r / p * dp;
	double s_drift = f_drift * 2.0 / M_PI * Sigmad_cgs/rho_p * vkep2 / c_s2 * fabs(1.0 / dlnPdlnr);
	return s_drift;

}



// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret	--> kimenet cm-ben!
double a_turb(double sigma, double r, double rho_p) {

	double s_frag, u_frag, u_frag2, Sigma_cgs, c_s, c_s2;

	u_frag = uFrag * CMPSECTOAUPYRP2PI; 	/*	cm/sec --> AU / (yr/2pi)	*/
  	u_frag2 = u_frag * u_frag;
	Sigma_cgs = sigma / SDCONV;
	c_s = c_sound(r); // / CMPSECTOAUPYRP2PI;
	c_s2 = c_s * c_s;

	s_frag = fFrag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * alpha_turb(r)) * u_frag2 / c_s2;
	
	return s_frag;

}


// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret		--> kimenet cm-ben!
double a_df(double sigma, double r, double p, double dp, double rho_p) {

	double u_frag, vkep, dlnPdlnr, c_s, c_s2, s_df, Sigma_cgs;

	u_frag = uFrag * CMPSECTOAUPYRP2PI; 	/*	cm/sec --> AU / (yr/2pi)	*/
	Sigma_cgs = sigma / SDCONV;
	c_s = c_sound(r);
	c_s2 = c_s * c_s;
	dlnPdlnr = r / p * dp;
	vkep = v_kep(r);

	s_df = u_frag * vkep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);

	return s_df;

}


/*	a reszecskek novekedesenek idoskalaja	*/
double tauGr(double r, double eps) {

	double omega = kep_freq(r);
	double taugr = eps / omega;
	
	return taugr;

}

/*	kiszamolja az adott helyen a reszecske meretet --> BIRNSTIEL CIKK	*/
double getSize(double prad, double pdens, double sigma, double sigmad, double y, double p, double dpress, double dt) {

	double sturb = a_turb(sigma,y,pdens);			// cm-ben
	double sdf = a_df(sigma,y,p,dpress,pdens);		// cm-ben
	double srdf = a_drift(sigmad,y,p,dpress,pdens);		// cm-ben
	double smin = find_min(sturb,sdf,srdf);		// cm-ben -- megadja, hogy a fenti ket reszecske korlatbol melyik ad kisebb meretet (az a reszecskenovekedes felso korlatja	
//	double eps = sigma / 100.;
	double eps = sigma / sigmad;
	double tau_gr = tauGr(y,eps);
	double rt = 0.0;

	smin = smin / AU2CM;					// AU-ban

/*	kiszamolja, hogy a fenti smin, vagy a novekedesi idoskalabol szarmazo meret korlatozza a reszecske meretet	*/
	if (prad < smin) rt = find_min(prad*exp(dt/tau_gr),smin,HUGE_VAL);
	if (prad >= smin) rt = smin;

	return rt;

}

void Get_Sigmad(double L, double max, double min, double rad[][2], double radmicr[][2], double radsec[][2], double *sigma_d, double *sigma_dm, double *sigma_ds, double *massvec, double *massmicrvec, double *masssecvec, double *rd, double *rmic, double *rs) {

	double dd = (RMAX - RMIN) / (PARTICLE_NUMBER-1);	/*	Mivel a jelen futas gridfelbontasa nem feltetlen egyezik meg a porreszecskeket generalo program gridfelbontasaval - ez a feluletisuruseg miatt lenyeges! - ezert itt szamolja ki a program	*/
	int i,j,k;
	double sigdtemp[PARTICLE_NUMBER][3], sigdmicrtemp[PARTICLE_NUMBER][3], sigdsectemp[4*PARTICLE_NUMBER][3], rtempvec[PARTICLE_NUMBER][2], sigtempvec[PARTICLE_NUMBER];
	double radtemp[PARTICLE_NUMBER], sdtemp[PARTICLE_NUMBER], radmicrtemp[PARTICLE_NUMBER], sdmicrtemp[PARTICLE_NUMBER], radsectemp[4*PARTICLE_NUMBER], sdsectemp[4*PARTICLE_NUMBER];
	double intsig[PARTICLE_NUMBER], intsigmicr[PARTICLE_NUMBER];
	double rtemp[PARTICLE_NUMBER];
	int rtempi[PARTICLE_NUMBER];

	for(i=0;i<PARTICLE_NUMBER;i++){

		rtempvec[i][0] = RMIN + i * dd;
		rtempvec[i][0] = rtempvec[i][0] + dd / 2.0;
		rtempvec[i][1] = i;
		radtemp[i] = 0;
		radmicrtemp[i] = 0;
		sdtemp[i] = 0;
		sdmicrtemp[i] = 0;
		intsig[i] = 0;
		intsigmicr[i] = 0;
		rd[i] = 0;
		rmic[i] = 0;
		rtemp[i] = 0;
		rtempi[i] = 0;
		sigtempvec[i] = 0;
	}


	loadSigDust(rad,massvec,sigdtemp,dd,PARTICLE_NUMBER);
	if(opttwopop == 1) {
		loadSigDust(radmicr,massmicrvec,sigdmicrtemp,dd,PARTICLE_NUMBER);
	}


//	for(i=0; i<PARTICLE_NUMBER;i++) printf("r: %lg  sd: %lg\n",rad[i][0],sigdtemp[i][0]);


	contract(L,sigdtemp,dd,PARTICLE_NUMBER);
	if(opttwopop == 1) {
		contract(L,sigdmicrtemp,dd,PARTICLE_NUMBER);
	}


	for(i=0; i < PARTICLE_NUMBER; i++) {
		rd[i] = sigdtemp[i][1];
		rmic[i] = 0;
		sigma_d[i] = sigdtemp[i][0];
		sigma_dm[i] = 0;
	}


}


/*	Fuggveny a sigma, p, dp kiszamolasara	*/
void Get_Sigma_P_dP(double *rvec, double *sigmavec, double *pressvec, double *dpressvec, double deltat) {

	double u, u_bi, u_fi, sigma_temp[NGRID+2], uvec[NGRID+2];
	int i;

	sigma_temp[0] = sigmavec[0];
	sigma_temp[NGRID+1] = sigmavec[NGRID+1];
	uvec[0] = sigmavec[0] * visc(rvec[0]);
	uvec[NGRID+1] = sigmavec[NGRID+1] * visc(rvec[NGRID+1]);

	for(i = 1; i <= NGRID; i++) {					/* 	solving d(sigma*nu)/dt = 3*nu*d2(sigma*nu)/dr2 + 9*nu/(2*r)*dsigma/dr	*/
		u = sigmavec[i] * visc(rvec[i]);			/*	az elozo lepesbol szarmazo sigma*nu elmentese	*/
		u_bi = sigmavec[i-1] * visc(rvec[i-1]);
		u_fi = sigmavec[i+1] * visc(rvec[i+1]);
		uvec[i] = u;
/*	az adott idolepesben az uj sigma*nu kiszamolasa a diffuzios egyenletbol	*/
		double temp = Coeff_1(rvec[i]) * (u_fi - 2.0 * u + u_bi) / (DD * DD) + Coeff_2(rvec[i]) * (u_fi - u_bi) / (2.0 * DD);
		sigma_temp[i] = uvec[i] + deltat * temp;		/*	regi + dt*uj	*/

    
	}

	for(i = 1; i <= NGRID; i++) {
		sigmavec[i] = sigma_temp[i]/visc(rvec[i]);		/*	nekunk sigma kell, nem sigma*nu		*/
		pressvec[i] = press(sigmavec[i],rvec[i]);		/*	nyomas kiszamolasa az uj sigmaval	*/
	}

	Perem(sigmavec);						/*	loading boundary condition in each timestep		*/
	dpress(dpressvec,pressvec);					/*	loading boundary condition in each timestep		*/
	Perem(pressvec);						/*	loading boundary condition in each timestep		*/
	Perem(dpressvec);						/*	loading boundary condition in each timestep		*/
	
}

/*	Fuggveny a porszemcsek uj tavolsaganak elraktarozasara		*/
void Get_Radius(char *nev, int opt, double radius[][2], double *pressvec, double *dpressvec, double *sigmavec, double *sigmad, double *rdvec, double *rvec, double *ugvec, double deltat, double t, int n) {		/*	a bemenet dimenziótlan egységekben van!!!	*/

	int i;
	double y, y_out, prad_new, particle_radius;
	char scout[1024];

/*	A futas kezdeten kiirja egy file-ba a drift idoskalajat!	*/
	if(t==0) {
		sprintf(scout,"%s/timescale.dat",nev);
		fout2 = fopen(scout,"w");
	}

	for (i = 0; i < n; i++) {	

		double drdt[n];

		if (radius[i][0] > RMIN && radius[i][0] < RMAX) { 

			y = radius[i][0];					// a reszecske tavolsaga
     			particle_radius = radius[i][1]; 			// a reszecske merete	

/*	a reszecske tavolsaganak kiszamolasa	*/	
     			int_step(t,particle_radius,pressvec,dpressvec,sigmavec,sigmad,rdvec,rvec,ugvec,deltat,y,&y_out,&prad_new);

/*	itt az opt a 2 populacios modellre utal, ha az 1, akkor a mikronos reszecskekre szamol - ehelyett egyebkent valoszinu lehetne az opttwopop-ot is hasznalni	*/
			if(t == 0) {
				if(opt == 0) {
					drdt[i] = (fabs(y_out - y)/(deltat));	/*	csak a cm-es esetben irja ki a drift idoskalajat	*/
						
/*	timescale.dat, ez tartalmazza a reszecskek bearamlasanak idoskalajat (evben)				*/
/*	A timescale.dat-hoz a t_drift = r / v(r,peak) kepletet a kovetkezo modon oldja meg a program:		*/
/*	Kiszamolja, hogy a kovetkezo idopillanatban hol lenne minden reszecske, majd az uj koordinatabol (x(1))	*/
/*	a regit (x(0)) levonva, majd az egeszet elosztva dt-vel megkapjuk a bearamlas sebesseget (drdt)		*/
/*	ez utan kiszamolhato az eredeti t_drift keplet, ezt iratom ki a file-ba evben (ezert osztok le 2PI-vel	*/

					fprintf(fout2,"%lg %lg\n",radius[i][0], (radius[i][0] / drdt[i])/2.0/M_PI);
				}

			}
	
/*	cm-es porokra az uj meret es pozicio elmentese				*/
			if (opt != 1){
				radius[i][1] = prad_new;
		     		radius[i][0] = y_out;
			} 

/*	mikronos esetben a az uj pozicio elmentese, azonban nincs novekedes	*/	
			if (opt == 1){
		     		radius[i][0] = y_out;
			}

/*	ha a reszecske RMIN-en belulre kerul, akkor a tavolsaga 0 lesz	*/		
		} else {
			radius[i][0] = 0.0;
		}
	}

	if(t==0) fclose(fout2);
}

