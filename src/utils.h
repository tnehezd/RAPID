#ifndef UTILS_H
#define UTILS_H

// Szükséges include-ok, ha a deklarációk más headerekben lévő típusokra vagy makrókra hivatkoznak.
// A Parabola függvénynek nincsenek komplex függőségei,
// de ha más utility függvények lennének itt, azoknak lehetnének.
// Jelenleg nem szükséges ide tenni semmit a Parabola miatt.

// Függvény deklarációk (prototípusok)
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd);
void Perem(double *vec);
/*	egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere	*/
void interpol(double *invec, double *rvec, double pos, double *out, double rd, int opt);

// megkeresi egy tomb maximumat
double find_max(double r[][2], int n);
double find_min(double s1, double s2, double s3);

//counting the number of zero points of the pressure gradient function	
int find_num_zero(double *dp);

//mivel a dp csak diszkret pontokban ismert, ezert 2 pontra illeszt egy egyenest, es megnezi, hogy annak hol lenne 0 az erteke
//	solving a*x + b = y (here a = r, y = dp)	
double find_r_zero(double r1, double r2, double dp1, double dp2);


// Nulla pont meghatározása lineáris interpolációval
double find_zero(int i, double *rvec, double *dp);

/* A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut */
// Feltehetően a disk_model.h tartalmazza a scale_height deklarációját.
void find_r_annulus(double *rvec, double rin_val, double *ind_ii, double *ind_io, double rout_val, double *ind_oi, double *ind_oo);


//fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	
void sort(double *rv,int n);

void histogram(double r, int *hist, double dd);



/*	fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	*/
void sortarray(double rv[][3],int n);

void kerekit(double in[][3], int n);

void contract(double L, double in[][3], double dd, int n); 

void Count_Mass(double radin[][2], double partmassindin[][4], double *massvecin, double t, int n);

// Ide kerülhetnek majd a jövőbeni egyéb általános segédfüggvények deklarációi is.
// Például:
// double calculate_average(double *data, int count);
// void print_timestamp(FILE *f);

#endif // UTILS_H
