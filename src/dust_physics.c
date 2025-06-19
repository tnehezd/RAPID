// src/dust_physics.c
#include "dust_physics.h" // A saját headerjét mindig includolni kell
#include "config.h" // Szükséges lehet a globális konstansokhoz
#include <stdio.h>
#include <math.h>

/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
double alpha_turb(double r) {

	double alpha_r;
	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));
   	return alpha_r * alpha_visc;

}


/*	Lokalis viszkozitas erteke	*/
double visc(double r){
 
  	double nu, alpha_r;
	double cs, H;
	
	H = scale_height(r);
	cs = c_sound(r);

/*	alpha_r: a redukcio merteke	*/
 	alpha_r = alpha_turb(r);

    	nu = alpha_r * cs * H;
  
  	return nu;
  
}

/*	local scale height	*/
double scale_height(double r) {

	return pow(r,1.+FLIND) * HASP;

}


/*	lokális kepleri sebesség	*/
double v_kep(double r) {

	return sqrt(G2 * STAR / r);
	
}


/*	lokalis kepleri korfrekvencia	*/
double kep_freq(double r) {

	return sqrt(G2 * STAR / r / r / r);			/*	v_kepler in AU / (yr/2pi)	*/

}


/*	local sound speed		*/
double c_sound(double r) {

	return kep_freq(r) * scale_height(r);

}

/*	Suruseg a midplane-ben	*/
double rho_mp(double sigma, double r) {

	return 1. / sqrt(2.0 * M_PI) * sigma / scale_height(r);

}

/* 	local pressure of the gas p = rho_gas * cs * cs kepletbol!!	*/
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

				if ((index_i >= (int)indii) && (index_i <= (int)indio) & partmassind[i][3] == 0) {
					partmassind[i][3] = 1;
					massitemp = massitemp + partmassind[i][0];
				}
			
			} else {
	
				if ((index_i >= (int)indii) && (index_i <= (int)indio)) {
					massitemp = partmassind[i][0] + massitemp;
				}
			}

			if(tavo != dzeo) {

				if ((index_i >= (int)indoi) && (index_i <= (int)indoo) & partmassind[i][3] == 0) {
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

				if ((index_i >= (int)indoi) && (index_i <= (int)indoo) & partmassind[i][3] == 0) {
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

/*	counting the number of zero points of the pressure gradient function	*/
int find_num_zero(double *rvec, double *dp) {

	int i,count;
	count = 0;

	for(i = 0; i < NGRID-1; i++) {
		if(((dp[i] * dp[i+1]) <= 0.)  && (dp[i] > dp[i+1])) {	/*	Osszeszorozza a ket ponton a nyomas derivaltjanak erteket, ahol a szorzat negativ, ott elojelvaltas tortenik --> negativbol pozitivba, vagy pozitivbol negativba valt --> nyomasi maximum. Maximum pedig ott talalhato, ahol a fuggveny pozitivbol negativba valt at (ezt keresi a masodik feltetel).	*/
			count++;
		} 

	}

	return count;
}

/*	mivel a dp csak diszkret pontokban ismert, ezert 2 pontra illeszt egy egyenest, es megnezi, hogy annak hol lenne 0 az erteke	*/
/*	solving a*x + b = y (here a = r, y = dp)	*/
double find_r_zero(double r1, double r2, double dp1, double dp2) {

	double a, b, r_zero;
	a = (dp2 - dp1) / (r2 - r1);
	b = dp1 - a * r1;
	r_zero = - b / a;

	return r_zero;

}

/*	this function counts where (which r) the pressure maximum is	*/
double find_zero(int i, double *rvec, double *dp) {

	double r;
	
	if(((dp[i] * dp[i+1]) <= 0.) && (dp[i] > dp[i+1])) {		/*	Ha a ket pont szorzata negativ --> elojel valtas a dp-ben, nyomasi min/max. Maximum hely ott van, ahol pozitivbol negativba valt az ertek	*/
		r = find_r_zero(rvec[i],rvec[i+1],dp[i],dp[i+1]);	/*	Ha elojel valtas tortenik es nyomasi maximum van, akkor kiszamolja a ket pont kozott, hogy hol lenne a zerus hely pontosan	*/
	} else {
		r = 0.0;
	}

	return r;

}
