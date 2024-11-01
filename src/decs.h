
/***********************************************************************************
    Copyright 2013 Joshua C. Dolence, Charles F. Gammie, Monika Mo\'scibrodzka,
                   and Po Kin Leung

                        GRMONTY  version 1.0   (released February 1, 2013)

    This file is part of GRMONTY.  GRMONTY v1.0 is a program that calculates the
    emergent spectrum from a model using a Monte Carlo technique.

    This version of GRMONTY is configured to use input files from the HARM code
    available on the same site.   It assumes that the source is a plasma near a
    black hole described by Kerr-Schild coordinates that radiates via thermal 
    synchrotron and inverse compton scattering.
    
    You are morally obligated to cite the following paper in any
    scientific literature that results from use of any part of GRMONTY:

    Dolence, J.C., Gammie, C.F., Mo\'scibrodzka, M., \& Leung, P.-K. 2009,
        Astrophysical Journal Supplement, 184, 387

    Further, we strongly encourage you to obtain the latest version of 
    GRMONTY directly from our distribution website:
    http://rainman.astro.illinois.edu/codelib/

    GRMONTY is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    GRMONTY is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GRMONTY; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include <omp.h>
#include "constants.h"

#define NDIM	4
#define NPRIM	8

/* Range of initial superphoton frequencies */
#define NUMIN 1.e9
#define NUMAX 1.e16

#define THETAE_MAX	1000.	/* Only used for harm3d models */
#define THETAE_MIN	0.3
#define TP_OVER_TE	(3.)

#define WEIGHT_MIN	(1e31)

#define HAMR (1)
#define HAMR3D (0) /* Leave it equal 0 to do HAMR2D */

/* mnemonics for primitive vars; conserved vars */
#define KRHO     0
#define UU      1
#define U1      2
#define U2      3
#define U3      4
#define B1      5
#define B2      6
#define B3      7

/* numerical convenience */
#define SMALL	1.e-40

/* physical parameters */
#define MMW	0.5		/* mean molecular weight, in units of mp */

/** data structures **/
struct of_photon {
	double X[NDIM];
	double K[NDIM];
	double dKdlam[NDIM];
	double w;
	double E;
	double L;
	double X1i;
	double X2i;
	double tau_abs;
	double tau_scatt;
	double ne0;
	double thetae0;
	double b0;
	double E0;
	double E0s;
	int nscatt;
};

struct of_geom {
	double gcon[NDIM][NDIM];
	double gcov[NDIM][NDIM];
	double g;
};

struct of_spectrum {
	double dNdlE;
	double dEdlE;
	double nph;
	double nscatt;
	double X1iav;
	double X2isq;
	double X3fsq;
	double tau_abs;
	double tau_scatt;
	double ne0;
	double thetae0;
	double b0;
	double E0;
};

#define N_ESAMP		200
#define N_EBINS		200
#define N_THBINS	6

struct of_grid {
	struct of_spectrum spec[N_EBINS];
	double th, phi;
	int nlist;
	int *in;
};

/** global variables **/
/** model independent */
extern gsl_rng *r;
extern int photon_scattered;
extern double F[N_ESAMP + 1], wgt[N_ESAMP + 1];

extern int Ns;
extern int N_superph_recorded, N_scatt;

/* HARM model globals */
extern struct of_geom *geom;
extern int N1, N2, N3;
extern int n_within_horizon;


/* some coordinate parameters */
extern double a;
extern double R0, Rin, Rh, Rout, Rms;
extern double hslope;
extern double startx[NDIM], stopx[NDIM], dx[NDIM];
extern double dlE, lE0;
extern double gam;
extern double dMsim;

extern double M_unit;
extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;
extern double Ne_unit;
extern double Thetae_unit;

extern double max_tau_scatt, Ladv, dMact, bias_norm;

/* some useful macros */
#define DLOOP  for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++)
#define SPATIAL_INDEX2D(i,j) ((j + N2 * i))
#define SPATIAL_INDEX3D(i,j,k) ((k) + N3 * (j + N2 * i))
#define INDEX(i,j,k)	(NPRIM*( (k) + N3*((j) + N2*(i))))
#define MYSIN(x,sx) 	{							\
			double _xp = (x)-M_PI; 					\
			double _yp = _xp*(FOUR_PI - FOUR_PISQ*fabs(_xp)); 	\
			sx = -_yp*(0.225*fabs(_yp)+0.775);			\
			}
#define MYCOS(x,cx) 	{							\
			double _xp = (x)-THREEPI_TWO; 					\
			_xp += (_xp<-M_PI)*2.*M_PI; 				\
			double _yp = _xp*(FOUR_PI - FOUR_PISQ*fabs(_xp));		\
			cx = _yp*(0.225*fabs(_yp)+0.775);			\
			}

/** model-independent subroutines **/
/* core monte carlo/radiative transport routines */
void track_super_photon(struct of_photon *ph, int round, int photon_index);
void record_super_photon(struct of_photon *ph);
void report_spectrum(int N_superph_made);
void scatter_super_photon(struct of_photon *ph, struct of_photon *php,
			  double Ne, double Thetae, double B,
			  double Ucon[NDIM], double Bcon[NDIM],
			  double Gcov[NDIM][NDIM]);

/* OpenMP specific functions */
void omp_reduce_spect(void);

/* MC/RT utilities */
void init_monty_rand(int seed);
double monty_rand(void);

/* geodesic integration */
void init_dKdlam(double X[], double Kcon[], double dK[]);
void push_photon_ham(double X[NDIM], double Kcon[][NDIM], double dl[]);
void push_photon(double X[NDIM], double Kcon[NDIM], double dKcon[NDIM],
		 double dl, double *E0, int n);
void push_photon4(double X[NDIM], double Kcon[NDIM], double dKcon[NDIM],
		  double dl);
void push_photon_cart(double X[NDIM], double Kcon[NDIM],
		      double dKcon[NDIM], double dl);
double stepsize(double X[NDIM], double K[NDIM]);
void push_photon_gsl(double X[NDIM], double Kcon[NDIM], double dl);
int geodesic_deriv(double t, const double y[], double dy[], void *params);
void interpolate_geodesic(double Xi[], double X[], double Ki[], double K[],
			  double frac, double del_l);

/* basic coordinate functions supplied by grmonty */
void boost(double k[NDIM], double p[NDIM], double ke[NDIM]);
void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov);
double gdet_func(double gcov[][NDIM]);  /* calculated numerically */
void coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM],
			  double K_tetrad[NDIM]);
void tetrad_to_coordinate(double Ecov[NDIM][NDIM], double K_tetrad[NDIM],
			  double K[NDIM]);
double delta(int i, int j);
void normalize(double Ucon[NDIM], double Gcov[NDIM][NDIM]);
void normalize_null(double Gcov[NDIM][NDIM], double K[NDIM]);
void make_tetrad(double Ucon[NDIM], double Bhatcon[NDIM],
		 double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
		 double Ecov[NDIM][NDIM]);

/* functions related to basic radiation functions & physics */
	/* physics-independent */
double get_fluid_nu(double X[4], double K[4], double Ucov[NDIM]);
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
		    double Bcov[NDIM], double B);
double alpha_inv_scatt(double nu, double thetae, double Ne);
double alpha_inv_abs(double nu, double thetae, double Ne, double B,
		     double theta);
double Bnu_inv(double nu, double thetae);
double jnu_inv(double nu, double thetae, double ne, double B,
	       double theta);

	/* thermal synchrotron */
double jnu_synch(double nu, double Ne, double Thetae, double B,
		 double theta);
double int_jnu(double Ne, double Thetae, double Bmag, double nu);
void init_emiss_tables(void);
double F_eval(double Thetae, double Bmag, double nu);
double K2_eval(double Thetae);

	/* compton scattering */
void init_hotcross(void);
double total_compton_cross_lkup(double nu, double theta);
double klein_nishina(double a, double ap);
double kappa_es(double nu, double theta);
void sample_electron_distr_p(double k[NDIM], double p[NDIM], double theta);
void sample_beta_distr(double theta, double *gamma_e, double *beta_e);
double sample_klein_nishina(double k0);
double sample_thomson(void);
double sample_mu_distr(double beta_e);
double sample_y_distr(double theta);
void sample_scattered_photon(double k[NDIM], double p[NDIM],
			     double kp[NDIM]);

/** model dependent functions required by code: these 
   basic interfaces define the model **/

/* physics related */
void init_model(char *args[]);
void make_super_photon(struct of_photon *ph, int *quit_flag);
double bias_func(double Te, double w);
void get_fluid_params(double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
		      double *Thetae, double *B, double Ucon[NDIM],
		      double Ucov[NDIM], double Bcon[NDIM],
		      double Bcov[NDIM]);
int stop_criterion(struct of_photon *ph);
int record_criterion(struct of_photon *ph);

/* coordinate related */
void get_connection(double *X, double lconn[][NDIM][NDIM]);
void gcov_func(double *X, double gcov[][NDIM]);
void gcon_func(double *X, double gcon[][NDIM]);


