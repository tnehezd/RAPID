#include "utils.h"    // Ezt kell includolni, mert ebben lesz a Parabola deklarációja
#include "config.h"   // Szükséges a RMIN és DD makrók miatt, amiket a Parabola használ
#include <math.h>     // Bár a Parabola most nem használ math.h függvényt,
#include <stdlib.h>                      // más utility függvényeknek szüksége lehet rá.
                      // Jó gyakorlat ide tenni.

#include "simulation_types.h" 
#include "dust_physics.h"

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
void calculate_boundary(double *vec, const disk_t *disk_params) {					/*	boundary condition for sigma, p, dp...	*/

	double a, b, c; 

//	Parabola(vec, 1, 2, 3, &a, &b, &c, disk_params->DD,disk_params);
//	vec[0] =  a * (disk_params->RMIN - disk_params->DD) * (disk_params->RMIN - disk_params->DD) + b * (disk_params->RMIN - disk_params->DD) + c;
	vec[0] = vec[1];
//	Parabola(vec, disk_params->NGRID - 2, disk_params->NGRID - 1, disk_params->NGRID, &a, &b, &c, disk_params->DD,disk_params);
//	vec[disk_params->NGRID+1] = a * (disk_params->RMAX + disk_params->DD) * (disk_params->RMAX + disk_params->DD) + b * (disk_params->RMAX + disk_params->DD) + c;
	vec[disk_params->NGRID+1] = vec[disk_params->NGRID];
}


/* egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere */
void interpol(double *invec, double *rvec, double pos, double *out, double rd, int opt, const disk_t *disk_params) {
    // Védőfeltétel a pos értékére
    if (isnan(pos) || isinf(pos) || pos < disk_params->RMIN || pos > disk_params->RMAX) {
        fprintf(stderr, "DEBUG [interpol]: Invalid position %.4e AU. Returning 0.\n", pos);
        *out = 0.0;
        return;
    }

    // Ellenőrizd a bemeneti tömbök érvényességét
    if (invec == NULL || rvec == NULL) {
        fprintf(stderr, "ERROR [interpol]: Input vectors are NULL. Returning 0.\n");
        *out = 0.0;
        return;
    }
    
    int index;
    double rmid, rindex, coef1, temp;

    // A rács indexének kiszámítása
    rmid = pos - disk_params->RMIN;
    rmid = rmid / rd;
    index = (int) floor(rmid);

    // Biztonsági ellenőrzés a tömb határainál
    if (index < 0) {
        index = 0;
    }
    if (index >= disk_params->NGRID - 1) {
        index = disk_params->NGRID - 2;
    }

    rindex = rvec[index];
    
    // Védőfeltételek a nullával való osztás ellen
    if (fabs(rvec[index + 1] - rvec[index]) < 1e-12) {
        fprintf(stderr, "ERROR [interpol]: Grid spacing is too small or zero at index %d. Returning value from current point.\n", index);
        *out = invec[index];
        return;
    }
    
    coef1 = (invec[index + 1] - invec[index]) / (rvec[index + 1] - rvec[index]);
    temp = invec[index] + coef1 * (pos - rindex);

    // Fizikailag értelmetlen értékek ellenőrzése
    // Ezeket a feltételeket csak fizikai mennyiségekre kell alkalmazni, amelyeknek pozitívnak kell lenniük.
    // Tegyük fel, hogy opt=1 a gáznyomásra (gas_pressure) vonatkozik.
    if (opt == 1 && temp < 0.0) {
        // Ha az érték negatív, de fizikailag pozitívnak kellene lennie, akkor beállítjuk egy kis pozitív számra.
        fprintf(stderr, "WARNING [interpol]: Interpolated value is negative (%.4e) for a positive quantity. Clamping to 0.0.\n", temp);
        temp = 0.0; // Vagy egy nagyon kicsi pozitív érték, pl. 1e-30
    }

    // Végleges ellenőrzés NaN/inf értékekre
    if (isnan(temp) || isinf(temp)) {
        fprintf(stderr, "CRITICAL ERROR [interpol]: Interpolated value is NaN/INF (%.4e) at pos=%.4e. Clamping to 0.0.\n", temp, pos);
        *out = 0.0;
    } else {
        *out = temp;
    }
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

// calculate_index_from_radius függvény (melyet korábban megbeszéltünk, valahol globálisan)
double calculate_index_from_radius(double r_coord, disk_t *disk_params) {
    if (r_coord < disk_params->RMIN) return 0.0;
    return fmax(0.0, fmin((double)(disk_params->NGRID - 1), floor((r_coord - disk_params->RMIN) / disk_params->DD + 0.5)));
}

// Az átalakított find_r_annulus függvény
// Paraméterek módosítva a struktúrák átadására és a konstansok használatára
/*	A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut	*/
void find_r_annulus(double rin, double *ind_ii, double *ind_io,
                            double rout, double *ind_oi, double *ind_oo,
                            const simulation_options_t *sim_opts, disk_t *disk_params) {

	    volatile int debug_marker = 0; // Adj hozzá ezt a sort


    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [find_r_annulus]: disk_params is NULL!\n");
        exit(1); // Program leállítása
    }


    // Lokális változók deklarálása
    int i;
    double rmid, rtemp;
    double roimH, roipH, roomH, roopH;
    double riimH, riipH, riomH, riopH;

    // --- FONTOS: Inicializáljuk az összes kimeneti indexet a ciklus ELŐTT! ---
    *ind_ii = 0.0;
    *ind_io = 0.0;
    *ind_oi = 0.0;
    *ind_oo = 0.0;

    // --- ITT HÍVJUK MEG A calculate_scale_height-et EGYSZER, ÉS MENTSÜK EL AZ EREDMÉNYT ---
    double h_rin = calculate_scale_height(rin, disk_params); // Első hívás, eredmény mentése

    double h_rout = calculate_scale_height(rout, disk_params); // Rout-ra is számoljuk ki egyszer

    // Számítsuk ki a határokhoz szükséges "rin +/- h_rin" és "rout +/- h_rout" értékeket
    // Ezeket a változókat használjuk majd a riimH, roimH stb. számításoknál
    double rin_minus_h_rin = rin - h_rin;
    double rin_plus_h_rin = rin + h_rin;
    double rout_minus_h_rout = rout - h_rout;
    double rout_plus_h_rout = rout + h_rout;


    // Határok kiszámítása: HASZNÁLJUK A MENTETT h_rin ÉS h_rout VÁLTOZÓKAT!
    // Ez kritikus, az eredeti elírásokat javítja.
    riimH = rin_minus_h_rin - disk_params->DD / 2.0;
    riipH = rin_minus_h_rin + disk_params->DD / 2.0;
    riomH = rin_plus_h_rin - disk_params->DD / 2.0;
    riopH = rin_plus_h_rin + disk_params->DD / 2.0;

    roimH = rout_minus_h_rout - disk_params->DD / 2.0;
    roipH = rout_minus_h_rout + disk_params->DD / 2.0;
    roomH = rout_plus_h_rout - disk_params->DD / 2.0;
    roopH = rout_plus_h_rout + disk_params->DD / 2.0;


    // Iteráció az rvec tömbön
    for (i = 0; i < disk_params->NGRID; i++) {
        // Ezen a ponton érdemes ellenőrizni disk_params->rvec[i] értékét
        // fprintf(stderr, "DEBUG_FIRA_LOOP: i=%d, rvec[i]=%.10lg\n", i, disk_params->rvec[i]);

        // Ez az if blokk csak akkor aktív, ha sim_opts->dzone == 1
        if (sim_opts->dzone == 1) {
            // INNER (RIN) határok
            if (disk_params->rvec[i] > riimH && disk_params->rvec[i] < riipH) {
                rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
                rtemp = floor(rmid + 0.5);
                *ind_ii = rtemp;
            }

            if (disk_params->rvec[i] > riomH && disk_params->rvec[i] < riopH) {
                rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
                rtemp = floor(rmid + 0.5);
                *ind_io = rtemp;
            }
        }

        // OUTER (ROUT) határok
        if (disk_params->rvec[i] > roimH && disk_params->rvec[i] < roipH) {
            rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
            rtemp = floor(rmid + 0.5);
            *ind_oi = rtemp;
        }

        if (disk_params->rvec[i] > roomH && disk_params->rvec[i] < roopH) {
            rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
            rtemp = floor(rmid + 0.5);
            *ind_oo = rtemp;
        }

        // KILÉPÉS feltétele
        if (disk_params->rvec[i] > roopH) break;
    }


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
	dker = dker * ROUND_PRECISION_FACTOR;
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

void Count_Mass(double radin[][2], double partmassindin[][5], double *massvecin, double t, int n, const disk_t *disk_params) {

    int i, rindex;
    double rmid;  


    for (i = 0; i < n; i++) {   
        // A részecske aktuális sugara radin[i][0]-ban van
        rmid = (radin[i][0] - disk_params->RMIN) / disk_params->DD; 
        rindex = (int) floor(rmid+0.5);
        if(rmid < 0) rindex = 0;
        if(isnan(rmid)) rindex = 0;


		if(n == PARTICLE_NUMBER) partmassindin[i][0] = massvecin[i];							/*	mass of the particles				*/
		partmassindin[i][1] = rindex;							/*	initial distance of the particles				*/
		if(t == 0) {
			partmassindin[i][2] = partmassindin[i][0];							/*	initial distance of the particles				*/
			partmassindin[i][3] = 0;
		}
 	
	}

}
