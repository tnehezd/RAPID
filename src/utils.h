#ifndef UTILS_H
#define UTILS_H

#include "simulation_types.h" 

// Szükséges include-ok, ha a deklarációk más headerekben lévő típusokra vagy makrókra hivatkoznak.
// A Parabola függvénynek nincsenek komplex függőségei,
// de ha más utility függvények lennének itt, azoknak lehetnének.
// Jelenleg nem szükséges ide tenni semmit a Parabola miatt.

// Függvény deklarációk (prototípusok)
void Parabola(double *vec, int i1, int i2, int i3, double *a, double *b, double *c, double dd, const disk_t *disk_params);
void Perem(double *vec,const disk_t *disk_params);
/*	egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere	*/
void interpol(double *invec, double *rvec, double pos, double *out, double rd, int opt, const disk_t *disk_params);

// megkeresi egy tomb maximumat
double find_max(double r[][2], int n);
double find_min(double s1, double s2, double s3);

//counting the number of zero points of the pressure gradient function	
int find_num_zero(const disk_t *disk_params);

//mivel a dp csak diszkret pontokban ismert, ezert 2 pontra illeszt egy egyenest, es megnezi, hogy annak hol lenne 0 az erteke
//	solving a*x + b = y (here a = r, y = dp)	
double find_r_zero(double r1, double r2, double dp1, double dp2);


// Nulla pont meghatározása lineáris interpolációval
double find_zero(int i, const double *rvec, const double *dp);

/* A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut */
// Feltehetően a disk_model.h tartalmazza a scale_height deklarációját.
void find_r_annulus(double rin_val, double *ind_ii, double *ind_io, double rout_val, double *ind_oi, double *ind_oo, const simulation_options_t *sim_opts, const disk_t *disk_params);

//fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	
void sort(double *rv,int n);

void histogram(double r, int *hist, double dd, disk_t *disk_params);



/*	fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	*/
void sortarray(double rv[][3],int n);

void kerekit(double in[][3], int n, const disk_t *disk_params);

void contract(double in[][3], double dd, int n, const disk_t *disk_params); 

void Count_Mass(double radin[][2], double partmassindin[][4], double *massvecin, double t, int n, const disk_t *disk_params);

// Ide kerülhetnek majd a jövőbeni egyéb általános segédfüggvények deklarációi is.
// Például:
// double calculate_average(double *data, int count);
// void print_timestamp(FILE *f);

#endif // UTILS_H
