#include "utils.h"      // Ezt kell includolni, mert ebben lesz a Parabola deklarációja
#include "config.h"     // Szükséges a RMIN és DD makrók miatt, amiket a Parabola használ
#include <math.h>       // Bár a Parabola most nem használ math.h függvényt,
#include <stdlib.h>     // más utility függvényeknek szüksége lehet rá.
                        // Jó gyakorlat ide tenni.
#include <stdio.h>      // For fprintf and stderr

#include "simulation_types.h"
#include "dust_physics.h"

/*	Parabola illesztés a peremen	*/
void Parabola(long double *vec, int i1, int i2, int i3, long double *a, long double *b, long double *c, long double dd, const disk_t *disk_params) { // Changed to long double

    long double x1, x2, x3;	/*	meghatározott x pontok, ahol illesztünk					*/ // Changed to long double
    long double y1, y2, y3;	/*	amit illesztünk a meghatározott pontokban				*/ // Changed to long double
    long double av, bv, cv;	/*	illesztéshez szükséges együtthatók --> ezt adja vissza a függvény	*/ // Changed to long double

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
void Perem(long double *vec, const disk_t *disk_params) {					/*	boundary condition for sigma, p, dp...	*/ // Changed to long double

    long double a, b, c; // Changed to long double

//	Parabola(vec, 1, 2, 3, &a, &b, &c, disk_params->DD,disk_params);
//	vec[0] =  a * (disk_params->RMIN - disk_params->DD) * (disk_params->RMIN - disk_params->DD) + b * (disk_params->RMIN - disk_params->DD) + c;
    vec[0] = vec[1];
//	Parabola(vec, disk_params->NGRID - 2, disk_params->NGRID - 1, disk_params->NGRID, &a, &b, &c, disk_params->DD,disk_params);
//	vec[disk_params->NGRID+1] = a * (disk_params->RMAX + disk_params->DD) * (disk_params->RMAX + disk_params->DD) + b * (disk_params->RMAX + disk_params->DD) + c;
    vec[disk_params->NGRID+1] = vec[disk_params->NGRID];
}


/*	egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere	*/
void interpol(long double *invec, long double *rvec, long double pos, long double *out, long double rd, int opt, const disk_t *disk_params) { // Changed to long double

    long double rmid, rindex, coef1, temp; // Changed to long double
    int index;

    rmid = pos - disk_params->RMIN;
    rmid = rmid / rd;                           /* the integer part of this gives at which index is the body	*/
    index = (int) floorl(rmid);                 /* ez az rmid egesz resze	(kerekites 0.5-tol)			*/ // Changed to floorl
    rindex = rvec[index];                       /* the corresponding r, e.g rd[ind] < r < rd[ind+1]		*/

    coef1 = (invec[index + 1] - invec[index]) / rd;         /*	ez az alabbi ket sor a linearis interpolacio - remelem, jo!!!	*/
    temp = invec[index] + coef1 * (pos - rindex);                  /*	a beerkezo dimenzionak megfelelo mertekegysegben		*/

    if(opt == 1) if(temp < 0) temp = -1.L*temp; // Changed to long double literal

    *out = temp;

}


/*	megkeresi egy tomb maximumat	*/
long double find_max(long double r[][2], int n) { // Changed to long double

    int i;
    long double maxim = -1.L; // Changed to long double literal

    for(i = 0; i < n; i++) {

        if (r[i][0] > maxim) {
            maxim = r[i][0];
        }

    }

    return maxim;
}

/*	minimum megkeresese harom elem kozul	*/
long double find_min(long double s1, long double s2, long double s3) { // Changed to long double

    long double min; // Changed to long double

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
        if(((disk_params->dpressvec[i] * disk_params->dpressvec[i+1]) <= 0.L)  && (disk_params->dpressvec[i] > disk_params->dpressvec[i+1])) {	/*	Osszeszorozza a ket ponton a nyomas derivaltjanak erteket, ahol a szorzat negativ, ott elojelvaltas tortenik --> negativbol pozitivba, vagy pozitivbol negativba valt --> nyomasi maximum. Maximum pedig ott talalhato, ahol a fuggveny pozitivbol negativba valt at (ezt keresi a masodik feltetel).	*/ // Changed to long double literal
            count++;
        }

    }

    return count;
}


/*	mivel a dp csak diszkret pontokban ismert, ezert 2 pontra illeszt egy egyenest, es megnezi, hogy annak hol lenne 0 az erteke	*/
/*	solving a*x + b = y (here a = r, y = dp)	*/
long double find_r_zero(long double r1, long double r2, long double dp1, long double dp2) { // Changed to long double

    long double a, b, r_zero; // Changed to long double
    a = (dp2 - dp1) / (r2 - r1);
    b = dp1 - a * r1;
    r_zero = - b / a;

    return r_zero;

}



/*	this function counts where (which r) the pressure maximum is	*/
long double find_zero(int i, const long double *rvec, const long double *dp) { // Changed to long double

    long double r; // Changed to long double

    if(((dp[i] * dp[i+1]) <= 0.L) && (dp[i] > dp[i+1])) {		/*	Ha a ket pont szorzata negativ --> elojel valtas a dp-ben, nyomasi min/max. Maximum hely ott van, ahol pozitivbol negativba valt az ertek	*/ // Changed to long double literal
        r = find_r_zero(rvec[i],rvec[i+1],dp[i],dp[i+1]);	/*	Ha elojel valtas tortenik es nyomasi maximum van, akkor kiszamolja a ket pont kozott, hogy hol lenne a zerus hely pontosan	*/

    } else {
        r = 0.0L; // Changed to long double literal
    }

    return r;

}

// calculate_index_from_radius függvény (melyet korábban megbeszéltünk, valahol globálisan)
long double calculate_index_from_radius(long double r_coord, disk_t *disk_params) { // Changed to long double
    if (r_coord < disk_params->RMIN) return 0.0L; // Changed to long double literal
    return fmaxl(0.0L, fminl((long double)(disk_params->NGRID - 1), floorl((r_coord - disk_params->RMIN) / disk_params->DD + 0.5L))); // Changed to long double literals and functions
}

// Az átalakított find_r_annulus függvény
// Paraméterek módosítva a struktúrák átadására és a konstansok használatára
/*	A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut	*/
void find_r_annulus(long double rin, long double *ind_ii, long double *ind_io, // Changed to long double
                            long double rout, long double *ind_oi, long double *ind_oo, long double r_buff, long double *ind_oi_buff, long double *ind_oo_buff, // Changed to long double
                            const simulation_options_t *sim_opts, disk_t *disk_params) {

    volatile int debug_marker = 0; // Adj hozzá ezt a sort


    if (disk_params == NULL) {
        fprintf(stderr, "ERROR [find_r_annulus]: disk_params is NULL!\n");
        exit(1); // Program leállítása
    }


    // Lokális változók deklarálása
    int i;
    long double rmid, rtemp; // Changed to long double
    long double roimH, roipH, roomH, roopH; // Changed to long double
    long double riimH, riipH, riomH, riopH; // Changed to long double
    long double roimH_buff, roipH_buff, roomH_buff, roopH_buff; // Changed to long double

    // --- FONTOS: Inicializáljuk az összes kimeneti indexet a ciklus ELŐTT! ---
    *ind_ii = 0.0L; // Changed to long double literal
    *ind_io = 0.0L; // Changed to long double literal
    *ind_oi = 0.0L; // Changed to long double literal
    *ind_oo = 0.0L; // Changed to long double literal
    *ind_oi_buff = 0.0L; // Changed to long double literal
    *ind_oo_buff = 0.0L; // Changed to long double literal

    // --- ITT HÍVJUK MEG A scale_height-et EGYSZER, ÉS MENTSÜK EL AZ EREDMÉNYT ---
    long double h_rin = scale_height(rin, disk_params); // Első hívás, eredmény mentése // Changed to long double

    long double h_rout = scale_height(rout, disk_params); // Rout-ra is számoljuk ki egyszer // Changed to long double


    // Számítsuk ki a határokhoz szükséges "rin +/- h_rin" és "rout +/- h_rout" értékeket
    // Ezeket a változókat használjuk majd a riimH, roimH stb. számításoknál
    long double rin_minus_h_rin = rin - h_rin; // Changed to long double
    long double rin_plus_h_rin = rin + h_rin; // Changed to long double
    long double rout_minus_h_rout = rout - h_rout; // Changed to long double
    long double rout_plus_h_rout = rout + h_rout; // Changed to long double



    // Határok kiszámítása: HASZNÁLJUK A MENTETT h_rin ÉS h_rout VÁLTOZÓKAT!
    // Ez kritikus, az eredeti elírásokat javítja.
    riimH = rin_minus_h_rin - disk_params->DD / 2.0L; // Changed to long double literal
    riipH = rin_minus_h_rin + disk_params->DD / 2.0L; // Changed to long double literal
    riomH = rin_plus_h_rin - disk_params->DD / 2.0L; // Changed to long double literal
    riopH = rin_plus_h_rin + disk_params->DD / 2.0L; // Changed to long double literal

    roimH = rout_minus_h_rout - disk_params->DD / 2.0L; // Changed to long double literal
    roipH = rout_minus_h_rout + disk_params->DD / 2.0L; // Changed to long double literal
    roomH = rout_plus_h_rout - disk_params->DD / 2.0L; // Changed to long double literal
    roopH = rout_plus_h_rout + disk_params->DD / 2.0L; // Changed to long double literal


    long double h_rout_buff = 0.0L; // Changed to long double literal
    long double rout_buff_minus_h_rout = 0.0L; // Changed to long double literal
    long double rout_buff_plus_h_rout = 0.0L; // Changed to long double literal

    if (r_buff != 0.0L) { // Changed to long double literal
        h_rout_buff = scale_height(r_buff, disk_params);
        rout_buff_minus_h_rout = r_buff - h_rout_buff;
        rout_buff_plus_h_rout = r_buff + h_rout_buff;
        roimH_buff = rout_buff_minus_h_rout - disk_params->DD / 2.0L; // Changed to long double literal
        roipH_buff = rout_buff_minus_h_rout + disk_params->DD / 2.0L; // Changed to long double literal
        roomH_buff = rout_buff_plus_h_rout - disk_params->DD / 2.0L; // Changed to long double literal
        roopH_buff = rout_buff_plus_h_rout + disk_params->DD / 2.0L; // Changed to long double literal
    } else {
        roimH_buff = 0.0L; // Changed to long double literal
        roipH_buff = 0.0L; // Changed to long double literal
        roomH_buff = 0.0L; // Changed to long double literal
        roopH_buff = 0.0L; // Changed to long double literal
    }



    // Iteráció az rvec tömbön
    for (i = 0; i < disk_params->NGRID; i++) {
        // Ezen a ponton érdemes ellenőrizni disk_params->rvec[i] értékét
        // fprintf(stderr, "DEBUG_FIRA_LOOP: i=%d, rvec[i]=%.10Lg\n", i, disk_params->rvec[i]); // Changed to %Lg

        // Ez az if blokk csak akkor aktív, ha sim_opts->dzone == 1
        if (sim_opts->dzone == 1.0L) { // Changed to long double literal
            // INNER (RIN) határok
            if (disk_params->rvec[i] > riimH && disk_params->rvec[i] < riipH) {
                rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
                rtemp = floorl(rmid + 0.5L); // Changed to floorl and long double literal
                *ind_ii = rtemp;
            }

            if (disk_params->rvec[i] > riomH && disk_params->rvec[i] < riopH) {
                rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
                rtemp = floorl(rmid + 0.5L); // Changed to floorl and long double literal
                *ind_io = rtemp;
            }
        }

        // OUTER (ROUT) határok
        if (disk_params->rvec[i] > roimH && disk_params->rvec[i] < roipH) {
            rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
            rtemp = floorl(rmid + 0.5L); // Changed to floorl and long double literal
            *ind_oi = rtemp;
        }

        if (disk_params->rvec[i] > roomH && disk_params->rvec[i] < roopH) {
            rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
            rtemp = floorl(rmid + 0.5L); // Changed to floorl and long double literal
            *ind_oo = rtemp;
        }


        // OUTER BUFFER (ROUT) határok
        if (disk_params->rvec[i] > roimH_buff && disk_params->rvec[i] < roipH_buff) {
            rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
            rtemp = floorl(rmid + 0.5L); // Changed to floorl and long double literal
            *ind_oi_buff = rtemp;
        }

        if (disk_params->rvec[i] > roomH_buff && disk_params->rvec[i] < roopH_buff) {
            rmid = (disk_params->rvec[i] - disk_params->RMIN) / disk_params->DD;
            rtemp = floorl(rmid + 0.5L); // Changed to floorl and long double literal
            *ind_oo_buff = rtemp;
        }

        // KILÉPÉS feltétele
        if (disk_params->rvec[i] > roopH_buff) break;
    }


}

/*	fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	*/
void sort(long double *rv,int n) { // Changed to long double

    long double temp; // Changed to long double
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
// void histogram(long double r, int *hist, long double dd, int hist_size) { // Changed to long double
void histogram(long double r, int *hist, long double dd, disk_t *disk_params) { // Changed to long double
    int index;
    long double rmid; // hist_i is no longer needed as a separate variable // Changed to long double

    // 1. Clamp 'r' to ensure it's within the valid range [RMIN, RMAX]
    // This prevents negative indices or indices that are too large.
    if (r < disk_params->RMIN) {
        r = disk_params->RMIN;
    } else if (r > disk_params->RMAX) {
        r = disk_params->RMAX;
    }

    // Calculate the potential index
    rmid = (r - disk_params->RMIN) / dd;
    index = (int) floorl(rmid); // Changed to floorl

    // 2. Explicitly check and clamp the index to the array bounds
    // Assuming 'hist' is an array of size NGRID, valid indices are 0 to NGRID - 1.
    // Replace NGRID with PARTICLE_NUMBER if that's the actual array size used for hist.
    if (index < 0) {
        index = 0; // Ensure index is not negative
        // Optionally, you could print a debug message if this happens unexpectedly:
        // fprintf(stderr, "DEBUG WARNING: histogram index became negative. Clamped to 0. r=%.10Lf, RMIN=%.10Lf, dd=%.10Le, rmid=%.10Lf\n", r, RMIN, dd, rmid); // Changed to %Lf, %Le
    }
    // Make sure NGRID is the correct size of the array 'hist'
    // If hist is int hist[PARTICLE_NUMBER], then upper bound is PARTICLE_NUMBER-1
    if (index >= disk_params->NGRID) { // NGRID should be defined and accessible here
        index = disk_params->NGRID - 1; // Ensure index does not exceed the array's upper bound
        // Optionally print a debug message:
        // fprintf(stderr, "DEBUG WARNING: histogram index exceeded NGRID. Clamped to NGRID-1. r=%.10Lf, RMAX=%.10Lf, dd=%.10Le, rmid=%.10Lf\n", r, RMAX, dd, rmid); // Changed to %Lf, %Le
    }

    // 3. Increment the counter directly (since hist is an int array)
    hist[index]++;
}

/*	fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	*/
void sortarray(long double rv[][3],int n) { // Changed to long double

    long double temp, temp2, temp3; // Changed to long double
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


void kerekit(long double in[][3], int n, const disk_t *disk_params) { // Changed to long double

    long double dd = (disk_params->RMAX - disk_params->RMIN) / (PARTICLE_NUMBER-1); // Changed to long double
    int dker = (int)(1.L/dd); // Changed to long double literal
    dker = dker * KEREK;
    long double ddker = (long double) dker; // Changed to long double
    int i;
    int temp;

    for(i = 0; i<n; i++) {

        temp = (int)floorl(in[i][1] * ddker+0.5L); // Changed to floorl and long double literal
        in[i][1] = (long double)temp / ddker; // Changed to long double
    }

}


void contract(long double in[][3], long double dd, int n, const disk_t *disk_params) { // Changed to long double

    int i;
    int j;
    int k;
    long double sig = 0, radout[n], sigout[n]; // Changed to long double

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
        long double rmid = (in[i][1] - disk_params->RMIN) / dd;                            /* The integer part of this gives at which index is the body			*/ // Changed to long double
        int rindex = (int) floorl(rmid);                           /* Ez az rmid egesz resze --> floor egeszreszre kerekit lefele, a +0.5-el elerheto, hogy .5 felett felfele, .5 alatt lefele kerekitsen						*/ // Changed to floorl
        in[i][2] = (long double)rindex; // Changed to long double
    }

}

void Count_Mass(long double radin[][2], long double partmassindin[][5], long double *massvecin, long double t, int n, const disk_t *disk_params) { // Changed to long double

    int i, rindex;
    long double rmid;   // Changed to long double


    for (i = 0; i < n; i++) {
        // A részecske aktuális sugara radin[i][0]-ban van
        rmid = (radin[i][0] - disk_params->RMIN) / disk_params->DD;
        rindex = (int) floorl(rmid+0.5L); // Changed to floorl and long double literal
        if(rmid < 0.L) rindex = 0; // Changed to long double literal
        if(isnan(rmid)) rindex = 0;


        if(n == PARTICLE_NUMBER) partmassindin[i][0] = massvecin[i];                           /*	mass of the particles				*/
        partmassindin[i][1] = rindex;                           /*	initial distance of the particles				*/
        if(t == 0.L) { // Changed to long double literal
            partmassindin[i][2] = partmassindin[i][0];                           /*	initial distance of the particles				*/
            partmassindin[i][3] = 0;
        }

    }

}