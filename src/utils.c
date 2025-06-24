#include "utils.h"    // Ezt kell includolni, mert ebben lesz a Parabola deklarációja
#include "config.h"   // Szükséges a RMIN és DD makrók miatt, amiket a Parabola használ
#include <math.h>     // Bár a Parabola most nem használ math.h függvényt,
                      // más utility függvényeknek szüksége lehet rá.
                      // Jó gyakorlat ide tenni.

#include "simulation_types.h" 

/*	Parabola illesztés a peremen	*/
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd, const disk_t *disk_params) {

	double x1, x2, x3;	/*	meghatározott x pontok, ahol illesztünk					*/
	double y1, y2, y3;	/*	amit illesztünk a meghatározott pontokban				*/
	double av, bv, cv;	/*	illesztéshez szükséges együtthatók --> ezt adja vissza a függvény	*/

	x1 = disk_params->RMIN + (i1-1) * dd;
	x2 = disk_params->RMIN + (i2-1) * dd;
	x3 = disk_params->RMIN + (i3-1) * dd;
 
	y1 = vec[i1];
	y2 = vec[i2];
	y3 = vec[i3];

	av = ((y1 - y3) / (x1 - x3) - (y1 - y2) / (x1 - x2)) / (x3 - x2);
	bv = (y1 - y2) / (x1 - x2) - av * (x1 + x2);
	cv = y1 - av * x1 * x1 - bv * x1;

	*a = av;
	*b = bv;
	*c = cv;

}


/*	A peremen parabolat illeszt	*/
void Perem(double *vec, const disk_t *disk_params) {					/*	boundary condition for sigma, p, dp...	*/

	double a, b, c; 

//	Parabola(vec, 1, 2, 3, &a, &b, &c, disk_params->DD,disk_params);
//	vec[0] =  a * (disk_params->RMIN - disk_params->DD) * (disk_params->RMIN - disk_params->DD) + b * (disk_params->RMIN - disk_params->DD) + c;
	vec[0] = vec[1];
//	Parabola(vec, disk_params->NGRID - 2, disk_params->NGRID - 1, disk_params->NGRID, &a, &b, &c, disk_params->DD,disk_params);
//	vec[disk_params->NGRID+1] = a * (disk_params->RMAX + disk_params->DD) * (disk_params->RMAX + disk_params->DD) + b * (disk_params->RMAX + disk_params->DD) + c;
	vec[disk_params->NGRID+1] = vec[disk_params->NGRID];
}


/*	egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere	*/
void interpol(double *invec, double *rvec, double pos, double *out, double rd, int opt, const disk_t *disk_params) {

	double rmid, rindex, coef1, temp;
	int index; 

    rmid = pos - disk_params->RMIN;
	rmid = rmid / rd;     						/* 	the integer part of this gives at which index is the body	*/
	index = (int) floor(rmid);					/* 	ez az rmid egesz resze	(kerekites 0.5-tol)			*/
	rindex = rvec[index];       					/* 	the corresponding r, e.g rd[ind] < r < rd[ind+1]		*/

 	coef1 = (invec[index + 1] - invec[index]) / rd; 		/*	ez az alabbi ket sor a linearis interpolacio - remelem, jo!!!	*/
	temp = invec[index] + coef1 * (pos - rindex);          		/*	a beerkezo dimenzionak megfelelo mertekegysegben		*/

	if(opt == 1) if(temp < 0) temp = -1.*temp;

	*out = temp;

}


/*	megkeresi egy tomb maximumat	*/
double find_max(double r[][2], int n) {

	int i;
	double maxim = -1.;

	for(i = 0; i < n; i++) {

		if (r[i][0] > maxim) {
			maxim = r[i][0];
		}
	
	}

	return maxim;
}

/*	minimum megkeresese harom elem kozul	*/
double find_min(double s1, double s2, double s3) {

	double min;
	
	if (s1 < s2) {
		if (s1 < s3) {
			min = s1;
		} else {
			min = s3;
		}
	} else {
		if (s2 < s3) {
			min = s2;	
		} else {
			min = s3;
		}
		
	}  

	return min;
}



/*	counting the number of zero points of the pressure gradient function	*/
int find_num_zero(const disk_t *disk_params) {

	int i,count;
	count = 0;

	for(i = 0; i < disk_params->NGRID-1; i++) {
		if(((disk_params->dpressvec[i] * disk_params->dpressvec[i+1]) <= 0.)  && (disk_params->dpressvec[i] > disk_params->dpressvec[i+1])) {	/*	Osszeszorozza a ket ponton a nyomas derivaltjanak erteket, ahol a szorzat negativ, ott elojelvaltas tortenik --> negativbol pozitivba, vagy pozitivbol negativba valt --> nyomasi maximum. Maximum pedig ott talalhato, ahol a fuggveny pozitivbol negativba valt at (ezt keresi a masodik feltetel).	*/
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
double find_zero(int i, const double *rvec, const double *dp) {

	double r;
	
	if(((dp[i] * dp[i+1]) <= 0.) && (dp[i] > dp[i+1])) {		/*	Ha a ket pont szorzata negativ --> elojel valtas a dp-ben, nyomasi min/max. Maximum hely ott van, ahol pozitivbol negativba valt az ertek	*/
		r = find_r_zero(rvec[i],rvec[i+1],dp[i],dp[i+1]);	/*	Ha elojel valtas tortenik es nyomasi maximum van, akkor kiszamolja a ket pont kozott, hogy hol lenne a zerus hely pontosan	*/
	} else {
		r = 0.0;
	}

	return r;

}


/*	A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut	*/
void find_r_annulus(const double *rvec, double rin, double *ind_ii, double *ind_io, double rout, double *ind_oi, double *ind_oo, const simulation_options_t *sim_opts, const disk_t *disk_params) {

    // Debug kiírás a függvény belépéskor, hogy ellenőrizzük a disk_params érvényességét
    fprintf(stderr, "DEBUG [find_r_annulus]: Belépés. disk_params címe=%p, FLIND=%.2f, HASP=%.2f\n",
            (void*)disk_params, disk_params->FLIND, disk_params->HASP);

    int i;
    double rmid, rtemp;
    double roimH;
    double roipH;
    double roomH;
    double roopH;
    double riimH;
    double riipH;
    double riomH;
    double riopH;

    if(sim_opts->dzone == 0) {
        *ind_ii = 0;
        *ind_io = 0;
    }

    // --- KRITIKUS JAVÍTÁS: Add hozzá a disk_params-t a scale_height hívásokhoz ---
    riimH = (rin - scale_height(rin, disk_params)) - disk_params->DD / 2.0;
    riipH = (rin - scale_height(rin, disk_params)) + disk_params->DD / 2.0;
    riomH = (rin + scale_height(rin, disk_params)) - disk_params->DD / 2.0;
    riopH = (rin + scale_height(rin, disk_params)) + disk_params->DD / 2.0;

    roimH = (rout - scale_height(rout, disk_params)) - disk_params->DD / 2.0;
    roipH = (rout - scale_height(rout, disk_params)) + disk_params->DD / 2.0;
    roomH = (rout + scale_height(rout, disk_params)) - disk_params->DD / 2.0;
    roopH = (rout + scale_height(rout, disk_params)) + disk_params->DD / 2.0;
    // --- KRITIKUS JAVÍTÁS VÉGE ---

    // Debug kiírások a számított értékekről (opcionális, de jó az ellenőrzéshez)
    fprintf(stderr, "DEBUG [find_r_annulus]: rin=%.2e, rout=%.2e\n", rin, rout);
    fprintf(stderr, "DEBUG [find_r_annulus]: riimH=%.2e, riipH=%.2e, riomH=%.2e, riopH=%.2e\n",
            riimH, riipH, riomH, riopH);
    fprintf(stderr, "DEBUG [find_r_annulus]: roimH=%.2e, roipH=%.2e, roomH=%.2e, roopH=%.2e\n",
            roimH, roipH, roomH, roopH);


    for(i = 1; i <= disk_params->NGRID; i++) {

        if(sim_opts->dzone == 1) {
            /* Ha az r távolság a kijelölt határok között van, akkor az adott változó visszakapja r értékét */
            if(rvec[i] > riimH && rvec[i] < riipH) {
                rmid = (disk_params->rvec[i] - disk_params->RMIN)/ disk_params->DD;
                rtemp = (int) floor(rmid + 0.5);
                *ind_ii = rtemp;
            }
                
            /* Ha az r távolság a kijelölt határok között van, akkor az adott változó visszakapja r értékét */
            if(rvec[i] > riomH && rvec[i] < riopH) {
                rmid = (rvec[i] - disk_params->RMIN)/ disk_params->DD;
                rtemp = (int) floor(rmid + 0.5);
                *ind_io = rtemp;
            }
        }


        /* Ha az r távolság a kijelölt határok között van, akkor az adott változó visszakapja r értékét */
        if(rvec[i] > roimH && rvec[i] < roipH) {
            rmid = (rvec[i] - disk_params->RMIN)/ disk_params->DD;
            rtemp = (int) floor(rmid + 0.5);
            *ind_oi = rtemp;
        }
                
        /* Ha az r távolság a kijelölt határok között van, akkor az adott változó visszakapja r értékét */
        if(rvec[i] > roomH && rvec[i] < roopH) {
            rmid = (rvec[i] - disk_params->RMIN)/ disk_params->DD;
            rtemp = (int) floor(rmid + 0.5);
            *ind_oo = rtemp;
        }

        if(rvec[i] > roopH) break;

    }
    fprintf(stderr, "DEBUG [find_r_annulus]: Kilépés. disk_params címe=%p, FLIND=%.2f, HASP=%.2f\n",
            (void*)disk_params, disk_params->FLIND, disk_params->HASP);
}



/*	fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	*/
void sort(double *rv,int n) {

	double temp, temp2;
	int i, step;

	for(step = 1; step <= (n-1); step++) {

		for(i = 0; i <= (n-2); i++) {

			if(rv[i] > rv[i + 1]) {

				temp = rv[i];
				rv[i] = rv[i + 1];
				rv[i + 1] = temp;

			}
		}
	}
}
// If NGRID is not directly available, you might need to pass the array size
// void histogram(double r, int *hist, double dd, int hist_size) {
void histogram(double r, int *hist, double dd, disk_t *disk_params) {
    int index;
    double rmid; // hist_i is no longer needed as a separate variable

    // 1. Clamp 'r' to ensure it's within the valid range [RMIN, RMAX]
    // This prevents negative indices or indices that are too large.
    if (r < disk_params->RMIN) {
        r = disk_params->RMIN;
    } else if (r > disk_params->RMAX) {
        r = disk_params->RMAX;
    }

    // Calculate the potential index
    rmid = (r - disk_params->RMIN) / dd;
    index = (int) floor(rmid);

    // 2. Explicitly check and clamp the index to the array bounds
    // Assuming 'hist' is an array of size NGRID, valid indices are 0 to NGRID - 1.
    // Replace NGRID with PARTICLE_NUMBER if that's the actual array size used for hist.
    if (index < 0) {
        index = 0; // Ensure index is not negative
        // Optionally, you could print a debug message if this happens unexpectedly:
        // fprintf(stderr, "DEBUG WARNING: histogram index became negative. Clamped to 0. r=%.10f, RMIN=%.10f, dd=%.10e, rmid=%.10f\n", r, RMIN, dd, rmid);
    }
    // Make sure NGRID is the correct size of the array 'hist'
    // If hist is int hist[PARTICLE_NUMBER], then upper bound is PARTICLE_NUMBER-1
    if (index >= disk_params->NGRID) { // NGRID should be defined and accessible here
        index = disk_params->NGRID - 1; // Ensure index does not exceed the array's upper bound
        // Optionally print a debug message:
        // fprintf(stderr, "DEBUG WARNING: histogram index exceeded NGRID. Clamped to NGRID-1. r=%.10f, RMAX=%.10f, dd=%.10e, rmid=%.10f\n", r, RMAX, dd, rmid);
    }

    // 3. Increment the counter directly (since hist is an int array)
    hist[index]++;
}

/*	fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	*/
void sortarray(double rv[][3],int n) {

	double temp, temp2, temp3;
	int i, step;

	for(step = 1; step <= (n-1); step++) {

		for(i = 0; i <= (n-2); i++) {

			if(rv[i][1] > rv[i + 1][1]) {

				temp = rv[i][1];
				rv[i][1] = rv[i + 1][1];
				rv[i + 1][1] = temp;

				temp2 = rv[i][0];
				rv[i][0] = rv[i + 1][0];
				rv[i + 1][0] = temp2;

				temp3 = rv[i][2];
				rv[i][2] = rv[i + 1][2];
				rv[i + 1][2] = temp3;
			}
		}
	}
}


void kerekit(double in[][3], int n, const disk_t *disk_params) {

	double dd = (disk_params->RMAX - disk_params->RMIN) / (PARTICLE_NUMBER-1);
	int dker = (int)(1./dd);//
	dker = dker * KEREK;
	double ddker = (double) dker;
	int i;
	int temp;

	for(i = 0; i<n; i++) {

		temp = (int)floor(in[i][1] * ddker+0.5);
		in[i][1] = (double)temp / ddker;
	
	}

}




void contract(double in[][3], double dd, int n, const disk_t *disk_params) {

	int i;
	int j;
	int k;
	double sig = 0, radout[n], sigout[n];

	for(i = 0; i < n; i++) {
		radout[i] = 0;
		sigout[i] = 0;
		in[i][2] = 0;
	}

	i = 0;
	j = 0;
	k = 0;

	kerekit(in,n,disk_params);
	sortarray(in,n);


	do {
		if(in[i][1] != in[i+1][1]) {

			radout[j] = in[i][1];
			sigout[j] = in[i][0];
			sig = 0;
			k = 0;
			j++;
			i++;
		} else {

			do {

				radout[j] = in[i][1];
				sig = sig + in[i+k][0];
				sigout[j] = sig;
				k++;

			} while (in[i][1] == in[i+k][1]);
			i = i+k;
			k = 0;
			j++;
		}

	} while (i < n);

	for(i = 0; i < n; i++) {
	
		in[i][0] = sigout[i];
		in[i][1] = radout[i];
  		double rmid = (in[i][1] - disk_params->RMIN) / dd;     						/* 	The integer part of this gives at which index is the body			*/
		int rindex = (int) floor(rmid);							/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
		in[i][2] = (double)rindex;
	}

}

/*	fuggveny a reszecskek tomegenek, indexenek es r_indexenek (a reszecske tavolsagaban levo gridcella indexe) frissitesere	*/
void Count_Mass(double radin[][2], double partmassindin[][4], double *massvecin, double t, int n, const disk_t *disk_params) {

	int i, rindex;
	double rmid;	

	for (i = 0; i < n; i++) {	

  		rmid = (radin[i][0] - disk_params->RMIN) / disk_params->DD;     						/* 	The integer part of this gives at which index is the body			*/
		rindex = (int) floor(rmid+0.5);							/* 	Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/
		if(rmid < 0) rindex = 0;
		if(isnan(rmid)) rindex = 0;	/*	neha az indexre - ha az mar RMIN-en belul van, nan-t ad, ezert abban az esetben is 0 lesz a rindex	 */

		if(n == PARTICLE_NUMBER) partmassindin[i][0] = massvecin[i];							/*	mass of the particles				*/
		partmassindin[i][1] = rindex;							/*	initial distance of the particles				*/
		if(t == 0) {
			partmassindin[i][2] = partmassindin[i][0];							/*	initial distance of the particles				*/
			partmassindin[i][3] = 0;
		}
 	}
}