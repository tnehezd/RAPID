#ifndef UTILS_H
#define UTILS_H

#include "simulation_types.h" 

// Szükséges include-ok, ha a deklarációk más headerekben lévő típusokra vagy makrókra hivatkoznak.
// A Parabola függvénynek nincsenek komplex függőségei,
// de ha más utility függvények lennének itt, azoknak lehetnének.
// Jelenleg nem szükséges ide tenni semmit a Parabola miatt.


/*	egy megadott, diszkret pontokban ismert fuggvenyt linearInterpolational a reszecske aktualis helyere	*/
void linearInterpolation(double *invec, double *rvec, double pos, double *out, double rd, int opt, const DiskParameters *disk_params);

// megkeresi egy tomb maximumat
double findMaximumOfAnArray(double r[][2], int n);
double findMinimumOfAnArray(double s1, double s2, double s3);

//counting the number of zero points of the pressure gradient function	
int countZeroPoints(const DiskParameters *disk_params);

//mivel a dp csak diszkret pontokban ismert, ezert 2 pontra illeszt egy egyenest, es megnezi, hogy annak hol lenne 0 az erteke
//	solving a*x + b = y (here a = r, y = dp)	
double findZeroPointRadius(double r1, double r2, double dp1, double dp2);


// Nulla pont meghatározása lineáris linearInterpolationációval
double findZeroPoint(int i, const double *rvec, const double *dp);

//double calculateIndexFromRadius(double r_coord, DiskParameters *disk_params);

/* A nyomasi maximum korul 1H tavolsagban jeloli ki a korgyurut */
// Feltehetően a disk_model.h tartalmazza a scale_height deklarációját.
void findRAnnulusAroundDZE(double rin_val, double *ind_ii, double *ind_io, double rout_val, double *ind_oi, double *ind_oo, const SimulationOptions *sim_opts, DiskParameters *disk_params);

//fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	
void sortAnArray(double *rv,int n);

void histogram(double r, int *hist, double dd, DiskParameters *disk_params);



/*	fuggveny egy tomb elemeinek sorbarendezesere --> ezt jelenleg nem hasznalja sehol a program	*/
void sortAnArrayarray(double rv[][3],int n);

void roundParticleRadii(double in[][3], int n, const DiskParameters *disk_params);

void mergeParticlesByRadius(double in[][3], double dd, int n, const DiskParameters *disk_params); 

void updateParticleGridIndices(double radin[][2], double partmassindin[][5], double *massvecin, double t, int n, const DiskParameters *disk_params);


// Ide kerülhetnek majd a jövőbeni egyéb általános segédfüggvények deklarációi is.
// Például:
// double calculate_average(double *data, int count);
// void print_timestamp(FILE *f);

#endif // UTILS_H
