#ifndef UTILS_H
#define UTILS_H

#include "simulation_types.h"

// Szükséges include-ok, ha a deklarációk más headerekben lévő típusokra vagy makrókra hivatkoznak.
// A Parabola függvénynek nincsenek komplex függőségei,
// de ha más utility függvények lennének itt, azoknak lehetnének.
// Jelenleg nem szükséges ide tenni semmit a Parabola miatt.

// Függvény deklarációk (prototípusok)
void Parabola(long double *vec, int i1, int i2, int i3, long double *a, long double *b, long double *c, long double dd, const disk_t *disk_params);
void Perem(long double *vec,const disk_t *disk_params);
/*	egy megadott, diszkret pontokban ismert fuggvenyt interpolal a reszecske aktualis helyere	*/
void interpol(long double *invec, long double *rvec, long double pos, long double *out, long double rd, int opt, const disk_t *disk_params);

// megkeresi egy tomb maximumat
long double find_max(long double r[][2], int n);
long double find_min(long double s1, long double s2, long double s3);

//counting the number of zero points of the pressure gradient function
int find_num_zero(const disk_t *disk_params);

//mivel a dp csak diszkret pontokban ismert, ezert 2 pontra illeszt egy egyenest, es megnezi, hogy annak hol lenne 0 az erteke
//	solving a*x + b = y (here a = r, y = dp)
long double find_r_zero(long double r1, long double r2, long double dp1, long double dp2);


// Nulla pont meghatározása lineáris interpolációval
long double find_zero(int i, const long double *rvec, const long double *dp);

long double calculate_index_from_radius(long double r_coord, disk_t *disk_params);

/* A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut */
// Feltehetően a disk_model.h tartalmazza a scale_height deklarációját.
void find_r_annulus(long double rin_val, long double *ind_ii, long double *ind_io, long double rout_val, long double *ind_oi, long double *ind_oo, long double r_buff, long double *ind_oi_buff, long double *ind_oo_buff, const simulation_options_t *sim_opts, disk_t *disk_params);

//fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program
void sort(long double *rv,int n);

void histogram(long double r, int *hist, long double dd, disk_t *disk_params);


/*	fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	*/
void sortarray(long double rv[][3],int n);

void kerekit(long double in[][3], int n, const disk_t *disk_params);

void contract(long double in[][3], long double dd, int n, const disk_t *disk_params);

void Count_Mass(long double radin[][2], long double partmassindin[][5], long double *massvecin, long double t, int n, const disk_t *disk_params);


// Ide kerülhetnek majd a jövőbeni egyéb általános segédfüggvények deklarációi is.
// Például:
// long double calculate_average(long double *data, int count);
// void print_timestamp(FILE *f);

#endif // UTILS_H