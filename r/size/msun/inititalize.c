//peldaprogram arra, hogyan kellene kiszamolni a reprezentativ reszecskek tomeget
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define SDCONV			1.12521e-7
#define ICEFACTOR       	3.0
#define SNOWLINE        	2.7
#define G			1.0
#define G2			G*G
#define AUPDAY2CMPSEC   	1.7314568e8    	//conversion between the velocities 
#define CMPSECTOAUPYRP2PI	3.35725e-07	//conversion between the velocities
#define GRPCM32MSUNAU3		1.68329e6

double ALPHA, STAR, DISKMASS, HASP,  FLIND, DUST2GASR;
double RMIN, RMAX;
double NGRID;
double r_dze_i;
double Dr_dze_i;
double r_dze_o;
double Dr_dze_o;
double a_mod;
double SIGMAP_EXP;
long double SIGMA0;
long double RATIO;
long double TWOPOP;
long double UNIFSIZE;

FILE *fin, *fout, *fout2;

typedef struct options {
	//! Number of particles
	double	n;
	//! Inner edge of the gaseous disk
	double	ri;
	//! Outer edge of the gaseous disk
	double	ro;
	//! Surface density of the gaseous disk
	long double	sigma0;
	//! Surface density of the gaseous disk in gr/cm/cm
	long double	sigma0cgs;
	//! The index of the profile of the initial surface density of the gaseous disk
	double	index;
	//! The inner deadzone edge
	double 	rdze_i;
	//! The outer deadzone edge
	double 	rdze_o;
	//! The width of the viscosity jump at the inner deadzone edge
	double 	drdze_i;
	//! The width of the viscosity jump at the outer deadzone edge
	double 	drdze_o;
	//! The alpha parameter (Shakura & Sunyaev (1973))
	double	alphaParam;
	//! Viscosity reduction parameter
	double 	amod;
	//! The aspect ratio of the disk (h = H/r)
	double	h;
	//! Flaring index
	double	flind;
	//! The mass of the central star (solar mass)
	double	m0;
	//! The mass of the disk (solar mass)
	double	md;
	//! Dust to gas ratio
	long double eps;
	//! Mass ratio between the population of dust
	long double ratio;
	//! Size of micronsized particle if twopopulation model is used
	long double mic;
	//! Size of the particles is uniform size is used
	long double onesize;

} options_t;


void create_default_options(options_t *opt) {



	opt->n				 = 1000;	
	opt->ri				 = 0.1;		
	opt->ro				 = 5.0;
	opt->sigma0			 = 0.0001;
	opt->sigma0cgs			 = 1.0e-4/SDCONV;
	opt->index			 = 0.5;
	opt->rdze_i			 = 0.0;		// rdze_i: nincs belso dze
	opt->rdze_o			 = 0.0;		// rdze_o: nincs kulso dze
	opt->drdze_i			 = 0.0;		
	opt->drdze_o			 = 0.0;
	opt->alphaParam		 	 = 1.0e-2;
	opt->amod			 = 0.01;
	opt->h				 = 5.0e-2;
	opt->flind			 = 0.0;		// flat korong!
	opt->m0				 = 1.0;		
	opt->md				 = 0.01;
	opt->eps			 = 0.01;
	opt->ratio			 = 0.85;
	opt->mic			 = 1e-4;
	opt->onesize			 = 1.0;
	
}


int parse_options(int argc, const char **argv, options_t *opt, double *par, double *dpar, double *exppar, double *sizepar){
	int i = 1;
	double temp1 = 0., temp2 = 0., temp3 = 0., temp4 = 0., temp5 = 0.;

	while (i < argc) {
		char *p = malloc(sizeof(argv[i]));
   		strcpy(p, argv[i]);

		// Number of intervals
		if(strcmp(p, "-n") == 0) {
			i++;
			opt->n = atof(argv[i]);
		}
		else if (strcmp(p, "-ri") == 0) {
			i++;
			opt->ri = atof(argv[i]);
		}
		else if (strcmp(p, "-ro") == 0) {
			i++;
			opt->ro = atof(argv[i]);
		}
		else if (strcmp(p, "-sigma0") == 0) {
			temp1 = 1.;
			i++;
			opt->sigma0 = atof(argv[i]);
		}
		else if (strcmp(p, "-sigma0cgs") == 0) {
			temp1 = 1.;
			i++;
			opt->sigma0 = atof(argv[i])/SDCONV;
			opt->sigma0cgs = atof(argv[i]);
		}
		else if (strcmp(p, "-index") == 0) {
			temp3 = 1.0;
			i++;
			opt->index = atof(argv[i]);
		}
		else if (strcmp(p, "-rdzei") == 0) {
			i++;
			opt->rdze_i = atof(argv[i]);
		}
		else if (strcmp(p, "-rdzeo") == 0) {
			i++;
			opt->rdze_o = atof(argv[i]);
		}
		else if (strcmp(p, "-drdzei") == 0) {
			i++;
			opt->drdze_i = atof(argv[i]);
		}
		else if (strcmp(p, "-drdzeo") == 0) {
			i++;
			opt->drdze_o = atof(argv[i]);
		}
		else if (strcmp(p, "-alpha") == 0) {
			i++;
			opt->alphaParam = atof(argv[i]);
		}	
		else if (strcmp(p, "-amod") == 0) {
			i++;
			opt->amod = atof(argv[i]);
		}
		else if (strcmp(p, "-h") == 0) {
			i++;
			opt->h = atof(argv[i]);
		}

		else if (strcmp(p, "-flind") == 0) {
			i++;
			opt->flind = atof(argv[i]);
		}
		else if (strcmp(p, "-m0") == 0) {
			i++;
			opt->m0 = atof(argv[i]);
		}
		else if (strcmp(p, "-md") == 0) {
			temp2= 1.;
			i++;
			opt->md = atof(argv[i]);
		}
		else if (strcmp(p, "-eps") == 0) {
			i++;
			opt->eps = atof(argv[i]);
		}
		else if (strcmp(p, "-ratio") == 0) {
			i++;
			opt->ratio = atof(argv[i]);
		}
		else if (strcmp(p, "-onesize") == 0) {
			temp4 = 1.;
			i++;
			opt->onesize = atof(argv[i]);
		}
		else if (strcmp(p, "-mic") == 0) {
			i++;
			opt->mic = atof(argv[i]);
		}
		else {
			printf("\n\n**** Invalid switch on command-line: %s! ****\n\n\n",p);
			printf("**** Try following parameters: ****\n\n-n: number of grids (particles)\n-ri\ninner edge of disk\n-ro\nother edge of disk\n-sigma0\nsigma at r=1AU\n-sigma0cgs\nsigma at r=1AU in cgs\n-index\nsigma exponent\n-rdzei\ninner DZE\n-rdzeo\nother DZE\n-drdzei\nthe width (in H, eg. 2H, 1.5H, etc) of the viscosity jump at the inner deadzone edge (if not given, but rdzei is given, set to 2*H)\n-drdzeo\nthe width (in H, eg. 2H, 1.5H, etc) of the viscosity jump at the outer deadzone edge (if not given, but rdzeo is given, set to 2*H)\n-alpha\nalpha viscosity\n-amod: viscosity reduction\n-h\naspect ratio\n-flind\nflaring index\n-m0\nmass of the central star (in solar masses)\n-md\nmass of the disk\n-eps\ndust to gas ratio\n\n\n");
			return 1;
		}
		i++;
		free(p);
	}
	
	*par = temp1;
	*dpar = temp2;
	*exppar = temp3;
	*sizepar = temp4;
	
	return 0;
}



/*	alpha turbulens paraméter kiszámolása --> alfa csökkentése alpha_r-rel	*/
double alpha_turb(double r) {
   
   	double alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));
   	return alpha_r * ALPHA;

}


/*	sigma0 kiszámolása (illetve megadása) M_NAP / AU / AU-ban	*/
double sigma_null() {
	
	double alpha2 = SIGMAP_EXP + 2;
	return (alpha2 / (2.0 * M_PI)) * DISKMASS / (pow(RMAX,alpha2) - pow(RMIN,alpha2));
	
}


/*	sigma kiszámolása r helyen M_NAP / AU / AU-ban	*/
double sigma_gas(double r) {

	return SIGMA0 * pow(r,SIGMAP_EXP);
	
}


/*	a por feluletisurusegenek kiszamitasa, a hohatar figyelembevetelevel M_SUN / AU / AU-ban	*/
long double sigma_dust(double r) {

	double sig_dust;
	
/*    	if (r <= SNOWLINE) {
	       	sig_dust = DUST2GASR * sigma0 * pow(r,SIGMAP_EXP);
    	} else {
        	sig_dust = DUST2GASR * sigma0 * ICEFACTOR * pow(r,SIGMAP_EXP);
	}
*/
	sig_dust = SIGMA0 * pow(r,SIGMAP_EXP) * DUST2GASR;

//	printf("gaz: %lg, por: %lg\n", sigma0*pow(r,SIGMAP_EXP), sig_dust);
		
        return sig_dust;

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


int main(int argc, const char **argv) {  

	options_t def;
	create_default_options(&def);
	double par, par2, par3, par4, par5;
	int retCode = parse_options(argc, argv, &def, &par, &par2, &par3, &par4);

	if (0 != retCode) {
		exit(retCode);
	}

	NGRID = def.n;

	SIGMAP_EXP = def.index;
	SIGMAP_EXP = SIGMAP_EXP*(-1.0);

	ALPHA = def.alphaParam;
	STAR = def.m0;
	DISKMASS = def.md; 
	HASP = def.h;
	FLIND = def.flind;
	DUST2GASR = def.eps;
	RMIN = def.ri;
	RMAX = def.ro;
	r_dze_i = def.rdze_i;
	Dr_dze_i = pow(r_dze_i,1.+FLIND)*HASP*def.drdze_i;
	r_dze_o = def.rdze_o;
	Dr_dze_o = pow(r_dze_o,1.+FLIND)*HASP*def.drdze_o;
	a_mod = def.amod;
	RATIO = def.ratio;
	UNIFSIZE = def.onesize;
	TWOPOP = def.mic;
	SIGMA0 = def.sigma0;

	if(par4 == 0) {
		RATIO = def.ratio;
	} else {
		RATIO = 1.;
	}

	if(par2 == 1) SIGMA0 = sigma_null();

  	int i; 		//ez egyben a reprezentativ reszecskek szamat is jelenti
  	double r, DD, reval, reval2;
	long double reppmass;
  	double f_drift, f_frag, u_frag, u_frag2, H;
	long double v_kep, c_s, v_kep2, c_s2, Sigma, Sigma_cgs, Sigmad_cgs, P, dPdr, dlnPdlnr;
	double rho_p, rho_mp, s_drift, s_frag, s_df;
  	double s_max;

	printf("sdexp: %lg\n",SIGMAP_EXP);
  
  	fout = fopen("init_data.dat","w");
  	fin = fopen("denspress.dat","r");

	DD = (RMAX-RMIN)/(NGRID - 1);
   
  	printf("a parameterek osszegzese\n");
  	printf("a korong tomege Naptomegben %lg\n",DISKMASS);
  	printf("korong belso hatara %lg CsE kulso hatara %lg CsE\n",RMIN,RMAX);
  	printf("a feluleti suruseg hatvanyfuggveny alakjanak kitevoje %lg\n",SIGMAP_EXP);
  	printf("a hohatar helyzete %lg CsE \n", SNOWLINE);
  	printf("a por feluleti suruseg ugrasanak tenyezoje a hohataron tul %lg\n",ICEFACTOR);
  	printf("a gaz feluleti surusege a csillagtol 1 CsE tavolsagban %Lg\n",SIGMA0);
  	printf("a por gaz arany %lg\n", DUST2GASR);
  	printf("a reprezentativ reszecskek szama %lg\n",NGRID);
  
  	f_drift = 0.55;
  	f_frag  = 0.37;
  
  	rho_p = 1.6;    	// g/cm^3
//	rho_p = rho_p * GRPCM32MSUNAU3; 
  	u_frag = 1000.0; 	// cm/s 
	//persze ez a reszecske osszeteteletol is fugg, valahogyan ezt is figyelembe lehetne/kellene venni
  
//  	u_frag = u_frag / AUPDAY2CMPSEC; //u_frag CsE/nap mertekegysegben
	u_frag = u_frag * CMPSECTOAUPYRP2PI; 	/*	cm/sec --> AU / (yr/2pi)	*/
  	u_frag2 = u_frag * u_frag;
  
  	for (i=0; i<NGRID; i++) {
       
       		r = RMIN + i * DD;
       		reval = r + DD / 2.0;
	   	reval2 = reval * reval;
      
       		reppmass = 2.0 * M_PI * r * DD * sigma_dust(reval); //ez itt a reprezentativ reszecske tomege naptomegben
       
		if(par4 == 0) {
 	    		H = pow(reval,FLIND+1.) * HASP;
       			v_kep = sqrt(G2 * STAR/reval); 
			double omega = v_kep / reval; 
                                       
		   	c_s = omega * H;
        	
       			v_kep2 = v_kep * v_kep;
       			c_s2 = c_s * c_s;

		   	Sigma = sigma_gas(reval);
			Sigma_cgs = sigma_gas(reval) / SDCONV;
			Sigmad_cgs = sigma_dust(reval) / SDCONV;
//		   	P = HASP * G2 * STAR * Sigma / (sqrt(2.0 * M_PI) * reval2);

			rho_mp = 1. / sqrt(2. * M_PI) * Sigma / H;
			P = rho_mp * c_s * c_s;
			dPdr = (FLIND + SIGMAP_EXP - 2) * pow(reval,(FLIND + SIGMAP_EXP - 3.0)) * HASP * G2 * STAR * SIGMA0 / sqrt(2.0 * M_PI);
	
		   	dlnPdlnr = reval/P * dPdr;
		  
	
	   		//reprezentativ reszecske kezdeti meretenek meghatarozasa
	       		// 1. radialis drift altal meghatarozott maximalis meret
	       		s_drift = f_drift * 2.0 / M_PI * Sigmad_cgs/rho_p * v_kep2 / c_s2 * fabs(1.0 / dlnPdlnr);
       
	       		// 2. a kis skalaju turbulencia altal okozott fragmentacio szerinti maximalis meret
	       		s_frag = f_frag * 2.0 / (3.0 * M_PI) * Sigma_cgs / (rho_p * alpha_turb(reval)) * u_frag2 / c_s2;
	       		//s_frag = f_frag * 2.0/(3.0*M_PI) * Sigma_cgs/(rho_p*ALPHA) * u_frag2/c_s2;
	
	       		// 3. radialis drift altal okozott fragmentacio szerinti maximalis meret
	       		s_df = u_frag * v_kep / fabs(dlnPdlnr * c_s2 * 0.5) * 2.0 * Sigma_cgs / (M_PI * rho_p);
	
	       		s_max = find_min(s_drift,s_frag,s_df);
	
// 	    		printf("s_drift: %lg, s_frag: %lg, s_df: %lg, smin: %lg\n", s_drift,s_frag,s_df, s_max);
       
 //     	 	printf("a reprezentativ reszecske sorszama, csillagtol valo tavolsaga es tomege %d %lg %lg %lg %lg %lg\n", i, reval, reppmass,s_drift,s_frag,s_df);

		} else {
			s_max = UNIFSIZE;
		}

       		fprintf(fout,"%d %lg  %Lg %Lg %lg %Lg\n",i,reval,reppmass*RATIO,reppmass*(1.-RATIO),s_max,TWOPOP);//,s_drift,s_frag,s_df);
//       		fprintf(fout,"%d %lg  %Lg %Lg %lg %Lg\n",i,reval,reppmass*RATIO,reppmass*(1.-RATIO),50.,TWOPOP);//,s_drift,s_frag,s_df);
  	}

  	printf("A reszecskeket tartalmazo file elkeszult (init_data.dat), a korong parametereit tartalmazo file kiirasa! \n\n A file tartalma: \n i (a reszecske sorszama), r (a reszecske tavolsaga), prad (a reszecske merete cm-ben), reppmass (a reszecske reprezentativ tomege)\n");
  	printf("A folytatashoz nyomj ENTER-t!\n");

  	fout2 = fopen("disk_param.dat","w"); 
       	fprintf(fout2,"%lg %lg %lg %lg %Lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",RMIN, RMAX, NGRID, SIGMAP_EXP, SIGMA0, G, r_dze_i, r_dze_o, Dr_dze_i, Dr_dze_o, a_mod, rho_p, ALPHA, STAR, FLIND);

  	fclose(fout2);   

  	printf("A korong parametereit tartalmazo file elkeszult (disk_param.dat) \n\n A file tartalma: \n RMIN, RMAX, NGRID (ez a reszecskek darabszama is), a profil kitevoje (SIGMAP_EXP), és sigma0 (M_Sun / AU / AU) - sigma az r=1AU-nál, G (gravit. konst), a belso deadzone hatara (AU), az atmenet vastagsaga (AU), a deadzone hatara (AU), az atmenet vastagsaga (AU), a viszkozitas csokkentes merteke, a reszecskek atlagsurusege, alfa es a kozponti csillag tomge! \n\n");    


  	return 0;

}
