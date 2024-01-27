
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


/*

	HARM model specification routines 

*/

#include "decs.h"
#define global
#include "harm_model.h"
#undef global
struct of_spectrum spect[N_THBINS][N_EBINS] = { };

#pragma omp threadprivate(spect)

/*

	encapsulates all initialization routines 
	
*/

void init_model(char *args[])
{
	/* find dimensional quantities from black hole
	   mass and its accretion rate */
	set_units(args[3]);

	fprintf(stderr, "getting simulation data...\n");
	#if(HAMR)
	init_hamr_data(args[2]); /*PEDRO EDIT -> file to read H-AMR data*/
	#else
	init_harm_data(args[2]);	/* read in HARM simulation data */
	#endif
	/* initialize the metric */
	fprintf(stderr, "initializing geometry...\n");
	fflush(stderr);
	init_geometry();
	fprintf(stderr, "done.\n\n");
	fflush(stderr);
	a = 0.9375;
	Rh = 1 + sqrt(1. - a * a);
	fprintf(stderr, "Rh = %le\n", Rh);

	/* make look-up table for hot cross sections */
	init_hotcross();

	/* make table for solid angle integrated emissivity and K2 */
	init_emiss_tables();

	/* make table for superphoton weights */
	init_weight_table();

	/* make table for quick evaluation of ns_zone */
	init_nint_table();

}

/*
	make super photon 
*/

int n2gen = -1;
double dnmax;
int zone_i, zone_j;

void make_super_photon(struct of_photon *ph, int *quit_flag)
{

	while (n2gen <= 0) {
		n2gen = get_zone(&zone_i, &zone_j, &dnmax);
	}

	n2gen--;
	//fprintf(stderr, "zone_i = %d, zone_j = %d\n", zone_i, zone_j);
	if (zone_i == N1)
		*quit_flag = 1;
	else
		*quit_flag = 0;

	if (*quit_flag != 1) {
		/* Initialize the superphoton energy, direction, weight, etc. */
		sample_zone_photon(zone_i, zone_j, dnmax, ph);
	}

	return;
}


double bias_func(double Te, double w)
{
	double bias, max, avg_num_scatt;

	max = 0.5 * w / WEIGHT_MIN;

	avg_num_scatt = N_scatt / (1. * N_superph_recorded + 1.);
	bias =
	    100. * Te * Te / (bias_norm * max_tau_scatt *
			      (avg_num_scatt + 2));

	if (bias < TP_OVER_TE)
		bias = TP_OVER_TE;
	if (bias > max)
		bias = max;

	return bias / TP_OVER_TE;
}

/* 

	these supply basic model data to grmonty

*/

void get_fluid_zone(int i, int j, double *Ne, double *Thetae, double *B,
		    double Ucon[NDIM], double Bcon[NDIM])
{

	int l, m;
	double Ucov[NDIM], Bcov[NDIM];
	double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;

	*Ne = p[KRHO][i][j] * Ne_unit;
	*Thetae = p[UU][i][j] / (*Ne) * Ne_unit * Thetae_unit;

	Bp[1] = p[B1][i][j];
	Bp[2] = p[B2][i][j];
	Bp[3] = p[B3][i][j];

	Vcon[1] = p[U1][i][j];
	Vcon[2] = p[U2][i][j];
	Vcon[3] = p[U3][i][j];

	/* Get Ucov */
	VdotV = 0.;
	for (l = 1; l < NDIM; l++)
		for (m = 1; m < NDIM; m++)
			VdotV += geom[i][j].gcov[l][m] * Vcon[l] * Vcon[m];
	Vfac = sqrt(-1. / geom[i][j].gcon[0][0] * (1. + fabs(VdotV)));
	Ucon[0] = -Vfac * geom[i][j].gcon[0][0];
	for (l = 1; l < NDIM; l++){
		Ucon[l] = Vcon[l] - Vfac * geom[i][j].gcon[0][l];
	}
	lower(Ucon, geom[i][j].gcov, Ucov);
	/* Get B and Bcov */
	UdotBp = 0.;
	for (l = 1; l < NDIM; l++)
		UdotBp += Ucov[l] * Bp[l];
	Bcon[0] = UdotBp;
	for (l = 1; l < NDIM; l++){
		Bcon[l] = (Bp[l] + Ucon[l] * UdotBp) / Ucon[0];
	}
	lower(Bcon, geom[i][j].gcov, Bcov);
	*B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
		  Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]) * B_unit;


}

void get_fluid_params(double X[NDIM], double gcov[NDIM][NDIM], double *Ne,
		      double *Thetae, double *B, double Ucon[NDIM],
		      double Ucov[NDIM], double Bcon[NDIM],
		      double Bcov[NDIM])
{
	int i, j;
	double del[NDIM];
	double rho, uu;
	double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
	double gcon[NDIM][NDIM], coeff[4];
	double interp_scalar(double **var, int i, int j, double del[4]);
	//fprintf(stderr, "startx[1] = %le, stopx[1] = %le, startx[2] = %le, stopx[2] = %le\n", startx[1], stopx[1], startx[2], stopx[2]);

	//checks if it's within the grid
	if (X[1] < startx[1] ||
	    X[1] > stopx[1] || X[2] < startx[2] || X[2] > stopx[2]) {

		*Ne = 0.;

		return;
	}
	// Finds out i and j index as well as fraction displacement del from the coordinates X[1], X[2], X[3]
	Xtoij(X, &i, &j, del);

	//Calculate the coeficient of displacement
	coeff[0] = (1. - del[1]) * (1. - del[2]);
	coeff[1] = (1. - del[1]) * del[2];
	coeff[2] = del[1] * (1. - del[2]);
	coeff[3] = del[1] * del[2];

	//interpolate based on the displacement
	rho = interp_scalar(p[KRHO], i, j, coeff);
	uu = interp_scalar(p[UU], i, j, coeff);
	// if (fabs(rho - p[KRHO][i][j]) > 1e-2){
	// 	fprintf(stderr, "X[1] = %le, X[2] = %le, i = %d, j = %d\n", X[1], X[2], i, j);
	// 	fprintf(stderr, "rho = %le, interp_rho = %le\n", p[KRHO][i][j], rho);
	// }
	*Ne = rho * Ne_unit;
	*Thetae = uu / rho * Thetae_unit;

	Bp[1] = interp_scalar(p[B1], i, j, coeff);
	Bp[2] = interp_scalar(p[B2], i, j, coeff);
	Bp[3] = interp_scalar(p[B3], i, j, coeff);

	Vcon[1] = interp_scalar(p[U1], i, j, coeff);
	Vcon[2] = interp_scalar(p[U2], i, j, coeff);
	Vcon[3] = interp_scalar(p[U3], i, j, coeff);

	#if(HAMR)
	gcon_func_hamr(gcov, gcon);
	//gcon_func(X, gcon);
	#else
	gcon_func(X, gcon);
	#endif
	
	/* Get Ucov */
	VdotV = 0.;
	for (i = 1; i < NDIM; i++)
		for (j = 1; j < NDIM; j++)
			VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
	Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
	Ucon[0] = -Vfac * gcon[0][0];
	for (i = 1; i < NDIM; i++){
		Ucon[i] = Vcon[i] - Vfac * gcon[0][i];
	}
	lower(Ucon, gcov, Ucov);

	/* Get B and Bcov */
	UdotBp = 0.;
	for (i = 1; i < NDIM; i++)
		UdotBp += Ucov[i] * Bp[i];
	Bcon[0] = UdotBp;
	for (i = 1; i < NDIM; i++)
		Bcon[i] = (Bp[i] + Ucon[i] * UdotBp) / Ucon[0];
	lower(Bcon, gcov, Bcov);

	*B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
		  Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]) * B_unit;

}


/* 
   Current metric: modified Kerr-Schild, squashed in theta
   to give higher resolution at the equator 
*/

/* mnemonics for dimensional indices */
#define TT      0
#define RR      1
#define TH      2
#define PH      3

void gcon_func(double *X, double gcon[][NDIM])
{

	int k, l;
	double sth, cth, irho2;
	double r, th, phi;
	double hfac;
	/* required by broken math.h */
	void sincos(double in, double *sth, double *cth);

	DLOOP gcon[k][l] = 0.;
	#if(HAMR)
	bl_coord_hamr(X, &r, &th, &phi);
	#else
	bl_coord(X, &r, &th);
	#endif


	sincos(th, &sth, &cth);
	sth = fabs(sth) + SMALL;

	irho2 = 1. / (r * r + a * a * cth * cth);

	// transformation for Kerr-Schild -> modified Kerr-Schild 
	hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);

	#if(HAMR)
	hfac = 1.;
	#endif

	gcon[TT][TT] = -1. - 2. * r * irho2;
	gcon[TT][1] = 2. * irho2;

	gcon[1][TT] = gcon[TT][1];
	gcon[1][1] = irho2 * (r * (r - 2.) + a * a) / (r * r);
	gcon[1][3] = a * irho2 / r;

	gcon[2][2] = irho2 / (hfac * hfac);

	gcon[3][1] = gcon[1][3];
	gcon[3][3] = irho2 / (sth * sth);
}

void gcov_func(double *X, double gcov[][NDIM])
{
	int k, l;
	double sth, cth, s2, rho2;
	double r, th, phi;
	double tfac, rfac, hfac, pfac;
	/* required by broken math.h */
	void sincos(double th, double *sth, double *cth);

	DLOOP gcov[k][l] = 0.;
	#if(HAMR)
	bl_coord_hamr(X, &r, &th, &phi);
	#else
	bl_coord(X, &r, &th);
	#endif
	sincos(th, &sth, &cth);
	sth = fabs(sth) + SMALL;
	s2 = sth * sth;
	rho2 = r * r + a * a * cth * cth;

	/* transformation for Kerr-Schild -> modified Kerr-Schild */
	tfac = 1.;
	rfac = r - R0;
	hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
	pfac = 1.;

	#if(HAMR)
	tfac = 1.;
	rfac = 1.;
	hfac = 1.;
	pfac = 1.;
	#endif

	gcov[TT][TT] = (-1. + 2. * r / rho2) * tfac * tfac;
	gcov[TT][1] = (2. * r / rho2) * tfac * rfac;
	gcov[TT][3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

	gcov[1][TT] = gcov[TT][1];
	gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
	gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

	gcov[2][2] = rho2 * hfac * hfac;

	gcov[3][TT] = gcov[TT][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][3] =
	    s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;

	//fprintf(stderr, "gcov[3][3]_harm = %lf\n", gcov[3][3]);

}




#undef TT
#undef RR
#undef TH
#undef PH

/* 

   connection calculated analytically for modified Kerr-Schild
   	coordinates 


   this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   where i = {1,2,3,4} corresponds to, e.g., {t,ln(r),theta,phi}
*/

void get_connection(double X[4], double lconn[4][4][4])
{
	double r1, r2, r3, r4, sx, cx;
	double th, dthdx2, dthdx22, d2thdx22, sth, cth, sth2, cth2, sth4,
	    cth4, s2th, c2th;
	double a2, a3, a4, rho2, irho2, rho22, irho22, rho23, irho23,
	    irho23_dthdx2;
	double fac1, fac1_rho23, fac2, fac3, a2cth2, a2sth2, r1sth2,
	    a4cth4;
	/* required by broken math.h */
	void sincos(double th, double *sth, double *cth);

	r1 = exp(X[1]);
	r2 = r1 * r1;
	r3 = r2 * r1;
	r4 = r3 * r1;

	sincos(2. * M_PI * X[2], &sx, &cx);

	/* HARM-2D MKS */
	#if(HAMR)
	double x2_mod;
	x2_mod = (X[2] + 1.)/2.;
	th = M_PI * x2_mod;
    dthdx2 = M_PI * (1./2.);
    d2thdx22 = 0;
	#else
	th = M_PI * X[2] + 0.5 * (1 - hslope) * sx;
	dthdx2 = M_PI * (1. + (1 - hslope) * cx);
	d2thdx22 = -2. * M_PI * M_PI * (1 - hslope) * sx;
	#endif
	dthdx22 = dthdx2 * dthdx2;

	sincos(th, &sth, &cth);
	sth2 = sth * sth;
	r1sth2 = r1 * sth2;
	sth4 = sth2 * sth2;
	cth2 = cth * cth;
	cth4 = cth2 * cth2;
	s2th = 2. * sth * cth;
	c2th = 2 * cth2 - 1.;

	a2 = a * a;
	a2sth2 = a2 * sth2;
	a2cth2 = a2 * cth2;
	a3 = a2 * a;
	a4 = a3 * a;
	a4cth4 = a4 * cth4;

	rho2 = r2 + a2cth2;                
	rho22 = rho2 * rho2;
	rho23 = rho22 * rho2;
	irho2 = 1. / rho2;
	irho22 = irho2 * irho2;
	irho23 = irho22 * irho2;
	irho23_dthdx2 = irho23 / dthdx2;

	fac1 = r2 - a2cth2;
	fac1_rho23 = fac1 * irho23;
	fac2 = a2 + 2 * r2 + a2 * c2th;
	fac3 = a2 + r1 * (-2. + r1);

	lconn[0][0][0] = 2. * r1 * fac1_rho23;
	lconn[0][0][1] = r1 * (2. * r1 + rho2) * fac1_rho23;
	lconn[0][0][2] = -a2 * r1 * s2th * dthdx2 * irho22;
	// fprintf(stderr, "a2 = %le\n", a2);
	// fprintf(stderr, "r1 = %le\n", r1);
	// fprintf(stderr, "irho22 = %le\n", irho22);
	// fprintf(stderr, "s2th = %le\n", s2th);
	// fprintf(stderr, "dthdx2 = %le\n", dthdx2);

	lconn[0][0][3] = -2. * a * r1sth2 * fac1_rho23;

	//lconn[0][1][0] = lconn[0][0][1];
	lconn[0][1][1] = 2. * r2 * (r4 + r1 * fac1 - a4cth4) * irho23;
	lconn[0][1][2] = -a2 * r2 * s2th * dthdx2 * irho22;
	lconn[0][1][3] =
	    a * r1 * (-r1 * (r3 + 2 * fac1) + a4cth4) * sth2 * irho23;

	//lconn[0][2][0] = lconn[0][0][2];
	//lconn[0][2][1] = lconn[0][1][2];
	lconn[0][2][2] = -2. * r2 * dthdx22 * irho2;
	lconn[0][2][3] = a3 * r1sth2 * s2th * dthdx2 * irho22;
	//lconn[0][3][0] = lconn[0][0][3];
	//lconn[0][3][1] = lconn[0][1][3];
	//lconn[0][3][2] = lconn[0][2][3];
	lconn[0][3][3] =
	    2. * r1sth2 * (-r1 * rho22 + a2sth2 * fac1) * irho23;

	lconn[1][0][0] = fac3 * fac1 / (r1 * rho23);
	lconn[1][0][1] = fac1 * (-2. * r1 + a2sth2) * irho23;
	lconn[1][0][2] = 0.;
	lconn[1][0][3] = -a * sth2 * fac3 * fac1 / (r1 * rho23);

	//lconn[1][1][0] = lconn[1][0][1];
	lconn[1][1][1] =
	    (r4 * (-2. + r1) * (1. + r1) +
	     a2 * (a2 * r1 * (1. + 3. * r1) * cth4 + a4cth4 * cth2 +
		   r3 * sth2 + r1 * cth2 * (2. * r1 + 3. * r3 -
					    a2sth2))) * irho23;
	lconn[1][1][2] = -a2 * dthdx2 * s2th / fac2;
	lconn[1][1][3] =
	    a * sth2 * (a4 * r1 * cth4 + r2 * (2 * r1 + r3 - a2sth2) +
			a2cth2 * (2. * r1 * (-1. + r2) + a2sth2)) * irho23;

	//lconn[1][2][0] = lconn[1][0][2];
	//lconn[1][2][1] = lconn[1][1][2];
	lconn[1][2][2] = -fac3 * dthdx22 * irho2;
	lconn[1][2][3] = 0.;

	//lconn[1][3][0] = lconn[1][0][3];
	//lconn[1][3][1] = lconn[1][1][3];
	//lconn[1][3][2] = lconn[1][2][3];
	lconn[1][3][3] =
	    -fac3 * sth2 * (r1 * rho22 - a2 * fac1 * sth2) / (r1 * rho23);

	lconn[2][0][0] = -a2 * r1 * s2th * irho23_dthdx2;
	lconn[2][0][1] = r1 * lconn[2][0][0];
	lconn[2][0][2] = 0.;
	lconn[2][0][3] = a * r1 * (a2 + r2) * s2th * irho23_dthdx2;

	//lconn[2][1][0] = lconn[2][0][1];
	lconn[2][1][1] = r2 * lconn[2][0][0];
	lconn[2][1][2] = r2 * irho2;
	lconn[2][1][3] =
	    (a * r1 * cth * sth *
	     (r3 * (2. + r1) +
	      a2 * (2. * r1 * (1. + r1) * cth2 + a2 * cth4 +
		    2 * r1sth2))) * irho23_dthdx2;

	//lconn[2][2][0] = lconn[2][0][2];
	//lconn[2][2][1] = lconn[2][1][2];
	lconn[2][2][2] =
	    -a2 * cth * sth * dthdx2 * irho2 + d2thdx22 / dthdx2;

	lconn[2][2][3] = 0.;

	//lconn[2][3][0] = lconn[2][0][3];
	//lconn[2][3][1] = lconn[2][1][3];
	//lconn[2][3][2] = lconn[2][2][3];
	lconn[2][3][3] =
	    -cth * sth * (rho23 +
			  a2sth2 * rho2 * (r1 * (4. + r1) + a2cth2) +
			  2. * r1 * a4 * sth4) * irho23_dthdx2;

	lconn[3][0][0] = a * fac1_rho23;
	lconn[3][0][1] = r1 * lconn[3][0][0];
	lconn[3][0][2] = -2. * a * r1 * cth * dthdx2 / (sth * rho22);
	lconn[3][0][3] = -a2sth2 * fac1_rho23;

	//lconn[3][1][0] = lconn[3][0][1];
	lconn[3][1][1] = a * r2 * fac1_rho23;
	lconn[3][1][2] =
	    -2 * a * r1 * (a2 + 2 * r1 * (2. + r1) +
			   a2 * c2th) * cth * dthdx2 / (sth * fac2 * fac2);
	lconn[3][1][3] = r1 * (r1 * rho22 - a2sth2 * fac1) * irho23;

	//lconn[3][2][0] = lconn[3][0][2];
	//lconn[3][2][1] = lconn[3][1][2];
	lconn[3][2][2] = -a * r1 * dthdx22 * irho2;
	lconn[3][2][3] =
	    dthdx2 * (0.25 * fac2 * fac2 * cth / sth +
		      a2 * r1 * s2th) * irho22;

	//lconn[3][3][0] = lconn[3][0][3];
	//lconn[3][3][1] = lconn[3][1][3];
	//lconn[3][3][2] = lconn[3][2][3];
	lconn[3][3][3] = (-a * r1sth2 * rho22 + a3 * sth4 * fac1) * irho23;

}

/* stopping criterion for geodesic integrator */
/* K not referenced intentionally */

#define RMAX	10.
#define ROULETTE	1.e4
int stop_criterion(struct of_photon *ph)
{
	double wmin, X1min, X1max;

	wmin = WEIGHT_MIN;	/* stop if weight is below minimum weight */
	
	X1min = log(Rh);	/* this is coordinate-specific; stop
				   at event horizon */
	X1max = log(RMAX);	/* this is coordinate and simulation
				   specific: stop at large distance */

	//fprintf(stderr, "w is: %le, wmin is: %le\n", ph->w, wmin);
	//fprintf(stderr, "X[1] = %le, X1min = %le\n", ph ->X[1], X1min);
	if (ph->X[1] < X1min)
		//fprintf(stderr, "it's getting here 1\n");
		return 1;

	if (ph->X[1] > X1max) {
		if (ph->w < wmin) {
			if (monty_rand() <= 1. / ROULETTE) {
				//fprintf(stderr, "it's getting here 2\n");
				ph->w *= ROULETTE;
			} else
				//fprintf(stderr, "it's getting here 3\n");
				ph->w = 0.;
		}
		return 1;
	}

	if (ph->w < wmin) {
		if (monty_rand() <= 1. / ROULETTE) {
			ph->w *= ROULETTE;
			//fprintf(stderr, "it's getting here 4\n");
		} else {
			ph->w = 0.;
			//fprintf(stderr, "it's getting here 5\n");
			return 1;
		}
	}

	return (0);
}

/* criterion for recording photon */

int record_criterion(struct of_photon *ph)
{
	const double X1max = log(RMAX);
	/* this is coordinate and simulation
	   specific: stop at large distance */
	//fprintf(stderr, "X[1] coord = %le, X1max = %le\n", ph->X[1], X1max);
	if (ph->X[1] > X1max)
		return (1);

	else
		return (0);

}


/* EPS really ought to be related to the number of
   zones in the simulation. */
//#define EPS	1e-5
#define EPS   0.04


double stepsize(double X[NDIM], double K[NDIM])
{
	double dl, dlx1, dlx2, dlx3;
	double idlx1, idlx2, idlx3;
	#if(HAMR)
		double x2_normal, stopx2_normal;
		x2_normal = (1 + X[2])/2;
		stopx2_normal = 1.; 
		dlx2 = EPS * GSL_MIN(x2_normal, stopx2_normal - x2_normal) / (fabs(K[2]) + SMALL);
	#else
		dlx2 = EPS * GSL_MIN(X[2], stopx[2] - X[2]) / (fabs(K[2]) + SMALL);
	#endif
	dlx1 = EPS * X[1] / (fabs(K[1]) + SMALL);
	dlx3 = EPS / (fabs(K[3]) + SMALL);

	idlx1 = 1. / (fabs(dlx1) + SMALL);
	idlx2 = 1. / (fabs(dlx2) + SMALL);
	idlx3 = 1. / (fabs(dlx3) + SMALL);

	dl = 1. / (idlx1 + idlx2 + idlx3);
	//fprintf(stderr, "dlx1 = %le, dlx2 = %le, dlx3 = %le, idlx1 = %le, idlx2 = %le, idlx3 = %le, dl = %le\n", dlx1, dlx2, dlx3, idlx1, idlx2, idlx3, dl);
	// for (int i = 0; i < NDIM; i++){
	// 	fprintf(stderr, "X[%d] = %le, k[%d] = %le\n", i, X[i], i, K[i]);
	// }
	//exit(1);

	return (dl);
}

/* 
	record contribution of super photon to spectrum.

	This routine should make minimal assumptions about the
	coordinate system.

*/
void record_super_photon(struct of_photon *ph)
{
	double lE, dx2;
	int iE, ix2;

	if (isnan(ph->w) || isnan(ph->E)) {
		fprintf(stderr, "record isnan: %g %g\n", ph->w, ph->E);
		return;
	}
#pragma omp critical (MAXTAU)
	{
		if (ph->tau_scatt > max_tau_scatt)
			max_tau_scatt = ph->tau_scatt;
	}
	/* currently, bin in x2 coordinate */

	/* get theta bin, while folding around equator */
	dx2 = (stopx[2] - startx[2]) / (2. * N_THBINS);
	if (ph->X[2] < 0.5 * (startx[2] + stopx[2]))
		ix2 = (int) (ph->X[2] / dx2);
	else
		ix2 = (int) ((stopx[2] - ph->X[2]) / dx2);

	/* check limits */
	if (ix2 < 0 || ix2 >= N_THBINS)
		return;

	/* get energy bin */
	lE = log(ph->E);
	iE = (int) ((lE - lE0) / dlE + 2.5) - 2;	/* bin is centered on iE*dlE + lE0 */

	/* check limits */
	if (iE < 0 || iE >= N_EBINS)
		return;

#pragma omp atomic
	N_superph_recorded++;
#pragma omp atomic
	N_scatt += ph->nscatt;

	/* sum in photon */
	spect[ix2][iE].dNdlE += ph->w;
	spect[ix2][iE].dEdlE += ph->w * ph->E;
	spect[ix2][iE].tau_abs += ph->w * ph->tau_abs;
	spect[ix2][iE].tau_scatt += ph->w * ph->tau_scatt;
	spect[ix2][iE].X1iav += ph->w * ph->X1i;
	spect[ix2][iE].X2isq += ph->w * (ph->X2i * ph->X2i);
	spect[ix2][iE].X3fsq += ph->w * (ph->X[3] * ph->X[3]);
	spect[ix2][iE].ne0 += ph->w * (ph->ne0);
	spect[ix2][iE].b0 += ph->w * (ph->b0);
	spect[ix2][iE].thetae0 += ph->w * (ph->thetae0);
	spect[ix2][iE].nscatt += ph->nscatt;
	spect[ix2][iE].nph += 1.;

}

struct of_spectrum shared_spect[N_THBINS][N_EBINS] = { };

void omp_reduce_spect()
{
/* Combine partial spectra from each OpenMP process		*
 * Inefficient, but only called once so doesn't matter	*/

	int i, j;

#pragma omp critical (UPDATE_SPECT)
	{
		for (i = 0; i < N_THBINS; i++) {
			for (j = 0; j < N_EBINS; j++) {
				shared_spect[i][j].dNdlE +=
				    spect[i][j].dNdlE;
				shared_spect[i][j].dEdlE +=
				    spect[i][j].dEdlE;
				shared_spect[i][j].tau_abs +=
				    spect[i][j].tau_abs;
				shared_spect[i][j].tau_scatt +=
				    spect[i][j].tau_scatt;
				shared_spect[i][j].X1iav +=
				    spect[i][j].X1iav;
				shared_spect[i][j].X2isq +=
				    spect[i][j].X2isq;
				shared_spect[i][j].X3fsq +=
				    spect[i][j].X3fsq;
				shared_spect[i][j].ne0 += spect[i][j].ne0;
				shared_spect[i][j].b0 += spect[i][j].b0;
				shared_spect[i][j].thetae0 +=
				    spect[i][j].thetae0;
				shared_spect[i][j].nscatt +=
				    spect[i][j].nscatt;
				shared_spect[i][j].nph += spect[i][j].nph;
			}
		}
	}
#pragma omp barrier
#pragma omp master
	{
		for (i = 0; i < N_THBINS; i++) {
			for (j = 0; j < N_EBINS; j++) {
				spect[i][j].dNdlE =
				    shared_spect[i][j].dNdlE;
				spect[i][j].dEdlE =
				    shared_spect[i][j].dEdlE;
				spect[i][j].tau_abs =
				    shared_spect[i][j].tau_abs;
				spect[i][j].tau_scatt =
				    shared_spect[i][j].tau_scatt;
				spect[i][j].X1iav =
				    shared_spect[i][j].X1iav;
				spect[i][j].X2isq =
				    shared_spect[i][j].X2isq;
				spect[i][j].X3fsq =
				    shared_spect[i][j].X3fsq;
				spect[i][j].ne0 = shared_spect[i][j].ne0;
				spect[i][j].b0 = shared_spect[i][j].b0;
				spect[i][j].thetae0 =
				    shared_spect[i][j].thetae0;
				spect[i][j].nscatt =
				    shared_spect[i][j].nscatt;
				spect[i][j].nph = shared_spect[i][j].nph;
			}
		}
	}
}

/* 

	output spectrum to file 

*/

#define SPECTRUM_FILE_NAME	"grmonty_hamr.spec"

void report_spectrum(int N_superph_made)
{
	int i, j;
	double dx2, dOmega, nuLnu, tau_scatt, L;
	FILE *fp;

	fp = fopen(SPECTRUM_FILE_NAME, "w");
	if (fp == NULL) {
		fprintf(stderr, "trouble opening spectrum file\n");
		exit(0);
	}

	/* output */
	max_tau_scatt = 0.;
	L = 0.;
	for (i = 0; i < N_EBINS; i++) {

		/* output log_10(photon energy/(me c^2)) */
		fprintf(fp, "%10.5g ", (i * dlE + lE0) / M_LN10);

		for (j = 0; j < N_THBINS; j++) {

			/* convert accumulated photon number in each bin 
			   to \nu L_\nu, in units of Lsun */
			dx2 = (stopx[2] - startx[2]) / (2. * N_THBINS);

			/* factor of 2 accounts for folding around equator */
			dOmega = 2. * dOmega_func(j * dx2, (j + 1) * dx2);

			nuLnu =
			    (ME * CL * CL) * (4. * M_PI / dOmega) * (1. /
								     dlE);

			nuLnu *= spect[j][i].dEdlE;
			nuLnu /= LSUN;

			tau_scatt =
			    spect[j][i].tau_scatt / (spect[j][i].dNdlE +
						     SMALL);
			fprintf(fp,
				"%10.5g %10.5g %10.5g %10.5g %10.5g %10.5g ",
				nuLnu,
				spect[j][i].tau_abs / (spect[j][i].dNdlE +
						       SMALL), tau_scatt,
				spect[j][i].X1iav / (spect[j][i].dNdlE +
						     SMALL),
				sqrt(fabs
				     (spect[j][i].X2isq /
				      (spect[j][i].dNdlE + SMALL))),
				sqrt(fabs
				     (spect[j][i].X3fsq /
				      (spect[j][i].dNdlE + SMALL)))
			    );

			if (tau_scatt > max_tau_scatt)
				max_tau_scatt = tau_scatt;

			L += nuLnu * dOmega * dlE;
		}
		fprintf(fp, "\n");
	}
	fprintf(stderr,
		"luminosity %g, dMact %g, efficiency %g, L/Ladv %g, max_tau_scatt %g\n",
		L, dMact * M_unit / T_unit / (MSUN / YEAR),
		L * LSUN / (dMact * M_unit * CL * CL / T_unit),
		L * LSUN / (Ladv * M_unit * CL * CL / T_unit),
		max_tau_scatt);
	fprintf(stderr, "\n");
	fprintf(stderr, "N_superph_made: %d\n", N_superph_made);
	fprintf(stderr, "N_superph_recorded: %d\n", N_superph_recorded);

	fclose(fp);

}
