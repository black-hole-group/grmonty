
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


#include "decs.h"
#include "harm_model.h"

/* Harm globals */
extern double *****econ;
extern double *****ecov;
extern double ****bcon;
extern double ****bcov;
extern double ****ucon;
extern double ****ucov;
extern double **p;
extern struct of_geom *geom;
extern double ***ne;
extern double ***thetae;
extern double ***b;

void Xtoij(double X[NDIM], int *i, int *j, double del[NDIM]);
void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);
void coord(int i, int j, double *X);
void get_fluid_zone(int i, int j, int k, double *Ne, double *Thetae, double *B,
		    double Ucon[NDIM], double Bcon[NDIM]);

/** HARM utilities **/

/********************************************************************

				Interpolation routines

 ********************************************************************/

// double interp_scalar(double ***var, int i, int j, double coeff[4])
// {

// 	double interp;

// 	interp =
// 	    var[i][j] * coeff[0] +
// 	    var[i][j + 1] * coeff[1] +
// 	    var[i + 1][j] * coeff[2] + var[i + 1][j + 1] * coeff[3];

// 	return interp;
// }

double interp_scalar(double *var, int i, int j, int k, double coeff[8]){
	double interp;

	interp = coeff[0] * var[SPATIAL_INDEX3D(i,j,k)] + coeff[5] * var[SPATIAL_INDEX3D(i+1,j,k)] +
	coeff[4] * var[SPATIAL_INDEX3D(i,j+1,k)] + coeff[7]  * var[SPATIAL_INDEX3D(i+1,j+1,k)] +
	coeff[1] * var[SPATIAL_INDEX3D(i,j,k+1)] + coeff[6] * var[SPATIAL_INDEX3D(i+1,j,k+1)] +
	coeff[2] * var[SPATIAL_INDEX3D(i,j+1,k+1)] + coeff[3] * var[SPATIAL_INDEX3D(i+1,j+1,k+1)];

	return interp;
}

double lnu_min, lnu_max, dlnu;

static void init_linear_interp_weight()
{

	lnu_min = log(NUMIN);
	lnu_max = log(NUMAX);
	dlnu = (lnu_max - lnu_min) / (N_ESAMP);
}

static double linear_interp_weight(double nu)
{

	int i;
	double di, lnu;

	lnu = log(nu);

	di = (lnu - lnu_min) / dlnu;
	i = (int) di;
	di = di - i;

	return exp((1. - di) * wgt[i] + di * wgt[i + 1]);

}


/***********************************************************************************

					End interpolation routines

 ***********************************************************************************/


#define JCST	(M_SQRT2*EE*EE*EE/(27*ME*CL*CL))
void init_weight_table(void)
{
	int k;
	int i, j, l, lstart, lend, myid, nthreads;
	double Ne, Thetae, B, K2;
	double sum[N_ESAMP + 1], nu[N_ESAMP + 1];
	double fac, sfac;
	double Ucon[NDIM], Bcon[NDIM];

	fprintf(stderr, "Building table for superphoton weights\n");
	fflush(stderr);

	/*      Set up interpolation */
	init_linear_interp_weight();

#pragma omp parallel for schedule(static) private(i)
	for (i = 0; i <= N_ESAMP; i++) {
		sum[i] = 0.;
		nu[i] = exp(i * dlnu + lnu_min);
	}

	sfac = dx[1] * dx[2] * dx[3] * L_unit * L_unit * L_unit;

#pragma omp parallel private(i,j,Thetae, K2, Ne, B, fac, l, lstart, lend,myid,nthreads,Ucon,Bcon)
	{
		nthreads = omp_get_num_threads();
		myid = omp_get_thread_num();
		lstart = myid * (N_ESAMP / nthreads);
		lend = (myid + 1) * (N_ESAMP / nthreads);
		if (myid == nthreads - 1)
			lend = N_ESAMP + 1;

		for (i = 0; i < N1; i++)
			for (j = 0; j < N2; j++)
				for (k = 0; k < N3; k++){
					get_fluid_zone(i, j, k, &Ne, &Thetae, &B,
							Ucon, Bcon);
					if (Ne == 0. || Thetae < THETAE_MIN)
						continue;
					K2 = K2_eval(Thetae);
					fac =
						(JCST * Ne * B * Thetae * Thetae /
						K2) * sfac * geom[SPATIAL_INDEX2D(i,j)].g;
					for (l = lstart; l < lend; l++){
						//fprintf(stderr, "external B = %le\n", B);
						sum[l] +=
							fac * F_eval(Thetae, B, nu[l]);	
				}
			}
#pragma omp barrier
	}
#pragma omp parallel for schedule(static) private(i)
	for (i = 0; i <= N_ESAMP; i++)
		wgt[i] = log(sum[i] / (HPL * Ns));

	fprintf(stderr, "done.\n\n");
	fflush(stderr);

	return;
}

#undef JCST

#define BTHSQMIN	(1.e-8)//Pedro edit from 1e-4 //1e-8 works for riaf
#define BTHSQMAX	(1.e14) //1.e8 Pedro edit  // 1e12 works for riaf
#define	NINT		(20000) //20000
double lb_min, dlb;
double nint[NINT + 1];
double dndlnu_max[NINT + 1];
void init_nint_table(void)
{

	int i, j;
	double Bmag, dn;
	static int firstc = 1;

	if (firstc) {
		lb_min = log(BTHSQMIN);
		dlb = log(BTHSQMAX / BTHSQMIN) / NINT;
		firstc = 0;
	}

	for (i = 0; i <= NINT; i++) {
		nint[i] = 0.;
		Bmag = exp(i * dlb + lb_min);
		dndlnu_max[i] = 0.;
		for (j = 0; j < N_ESAMP; j++) {
			dn = F_eval(1., Bmag,
				    exp(j * dlnu +
					lnu_min)) / (exp(wgt[j]) +
						     1.e-100);
			if (dn > dndlnu_max[i])
				dndlnu_max[i] = dn;
			nint[i] += dlnu * dn;
		}
		nint[i] *= dx[1] * dx[2] * dx[3] * L_unit * L_unit * L_unit
		    * M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL)
		    * 1. / HPL;
		nint[i] = log(nint[i]);
		dndlnu_max[i] = log(dndlnu_max[i]);
	}

	return;
}

static void init_zone(int i, int j, int k, double *nz, double *dnmax)
{
	int l;
	double Ne, Thetae, Bmag, lbth;
	double dl, dn, ninterp, K2;
	double Ucon[NDIM], Bcon[NDIM];

	get_fluid_zone(i, j, k, &Ne, &Thetae, &Bmag, Ucon, Bcon);

	if (Ne == 0. || Thetae < THETAE_MIN) {
		*nz = 0.;
		*dnmax = 0.;
		return;
	}

	lbth = log(Bmag * Thetae * Thetae);

	dl = (lbth - lb_min) / dlb;
	l = (int) dl;
	dl = dl - l;
	if (l < 0) {
		*dnmax = 0.;
		*nz = 0.;
		return;
	} else if (l >= NINT) {
		fprintf(stderr,
			"warning: outside of nint table range %g...change in harm_utils.c\n",
			Bmag * Thetae * Thetae);
		fprintf(stderr, "lbth = %le, lb_min = %le, dlb = %le l = %d, Bmag = %le, Thetae = %le, (%d, %d)\n", lbth, lb_min, dlb, l, Bmag, Thetae, i,j);
		ninterp = 0.;
		*dnmax = 0.;
		for (l = 0; l <= N_ESAMP; l++) {
			dn = F_eval(Thetae, Bmag,
				    exp(j * dlnu +
					lnu_min)) / (exp(wgt[l]) +
						     1.e-100);
			if (dn > *dnmax)
				*dnmax = dn;
			ninterp += dlnu * dn;
		}
		ninterp *= dx[1] * dx[2] * dx[3] * L_unit * L_unit * L_unit
		    * M_SQRT2 * EE * EE * EE / (27. * ME * CL * CL)
		    * 1. / HPL;
	} else {
		if (isinf(nint[l]) || isinf(nint[l + 1])) {
			ninterp = 0.;
			*dnmax = 0.;
		} else {
			ninterp =
			    exp((1. - dl) * nint[l] + dl * nint[l + 1]);
			*dnmax =
			    exp((1. - dl) * dndlnu_max[l] +
				dl * dndlnu_max[l + 1]);
		}
	}

	K2 = K2_eval(Thetae);
	if (K2 == 0.) {
		*nz = 0.;
		*dnmax = 0.;
		return;
	}

	*nz = geom[SPATIAL_INDEX2D(i,j)].g * Ne * Bmag * Thetae * Thetae * ninterp / K2;
	if (*nz > Ns * log(NUMAX / NUMIN)) {
		fprintf(stderr,
			"Something very wrong in zone %d %d: \nB=%g  Thetae=%g  K2=%g  ninterp=%g\n\n",
			i, j, Bmag, Thetae, K2, ninterp);
		*nz = 0.;
		*dnmax = 0.;
	}

	return;
}

int zone_flag;
int get_zone(int *i, int *j, int *k, double *dnmax)
{
/* Return the next zone and the number of superphotons that need to be		*
 * generated in it.								*/

	int in2gen;
	double n2gen;
	static int zi = 0;
	static int zj = 0;
	static int zk = -1;

	zone_flag = 1;
	zk++;
	if(zk >= N3){
		zk = 0;
		zj++;
		if (zj >= N2) {
			zj = 0;
			fprintf(stderr, "zi = %d\n", zi);
			zi++;
			if (zi >= N1) {
				in2gen = 1;
				*i = N1;
				return 1;
			}
		}
	}
	init_zone(zi, zj, zk, &n2gen, dnmax);
	if (fmod(n2gen, 1.) > monty_rand()) {
		in2gen = (int) n2gen + 1;
	} else {
		in2gen = (int) n2gen;
	}

	*i = zi;
	*j = zj;
	*k = zk;

	return in2gen;
}

void sample_zone_photon(int i, int j, int k, double dnmax, struct of_photon *ph)
{
/* Set all initial superphoton attributes */
	int l;
	int z = 0;
	double K_tetrad[NDIM], tmpK[NDIM], E, Nln;
	double nu, th, cth, sth, phi, sphi, cphi, jmax, weight;
	double Ne, Thetae, Bmag, Ucon[NDIM], Bcon[NDIM], bhat[NDIM];
	static double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
	#if(HAMR)
	coord_hamr(i, j, z, CENT, ph -> X);
	#else
	coord(i, j, ph->X);
	#endif
	Nln = lnu_max - lnu_min;

	get_fluid_zone(i, j, z, &Ne, &Thetae, &Bmag, Ucon, Bcon);

	/* Sample from superphoton distribution in current simulation zone */
	do {
		nu = exp(monty_rand() * Nln + lnu_min);
		weight = linear_interp_weight(nu);
	} while (monty_rand() >
		 (F_eval(Thetae, Bmag, nu) / (weight + 1.e-100)) / dnmax);

	ph->w = weight;
	jmax = jnu_synch(nu, Ne, Thetae, Bmag, M_PI / 2.);
	do {
		cth = 2. * monty_rand() - 1.;
		th = acos(cth);

	} while (monty_rand() >
		 jnu_synch(nu, Ne, Thetae, Bmag, th) / jmax);

	sth = sqrt(1. - cth * cth);
	phi = 2. * M_PI * monty_rand();
	cphi = cos(phi);
	sphi = sin(phi);

	E = nu * HPL / (ME * CL * CL);
	K_tetrad[0] = E;
	K_tetrad[1] = E * cth;
	K_tetrad[2] = E * cphi * sth;
	K_tetrad[3] = E * sphi * sth;

	if (zone_flag) {	/* first photon created in this zone, so make the tetrad */
		if (Bmag > 0.) {
			for (l = 0; l < NDIM; l++)
				bhat[l] = Bcon[l] * B_unit / Bmag;
		} else {
			for (l = 1; l < NDIM; l++)
				bhat[l] = 0.;
			bhat[1] = 1.;
		}
		make_tetrad(Ucon, bhat, geom[SPATIAL_INDEX2D(i,j)].gcov, Econ, Ecov);
		zone_flag = 0;
	}

	tetrad_to_coordinate(Econ, K_tetrad, ph->K);

	K_tetrad[0] *= -1.;
	tetrad_to_coordinate(Ecov, K_tetrad, tmpK);

	ph->E = ph->E0 = ph->E0s = -tmpK[0];
	ph->L = tmpK[3];
	ph->tau_scatt = 0.;
	ph->tau_abs = 0.;
	ph->X1i = ph->X[1];
	ph->X2i = ph->X[2];
	ph->nscatt = 0;
	ph->ne0 = Ne;
	ph->b0 = Bmag;
	ph->thetae0 = Thetae;

	return;
}

void Xtoij(double X[NDIM], int *i, int *j, double del[NDIM])
{

	*i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
	*j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;

	if (*i < 0) {
		*i = 0;
		del[1] = 0.;
	} else if (*i > N1 - 2) {
		*i = N1 - 2;
		del[1] = 1.;
	} else {
		del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
	}

	if (*j < 0) {
		*j = 0;
		del[2] = 0.;
	} else if (*j > N2 - 2) {
		*j = N2 - 2;
		del[2] = 1.;
	} else {
		del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2]; //fractional displacement of the center of the grid cel
	}

	return;
}


void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM])
{

	*i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
	*j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
	if (*i < 0) {
		*i = 0;
		del[1] = 0.;
	} else if (*i > N1 - 2) {
		*i = N1 - 2;
		del[1] = 1.;
	} else {
		del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
	}

	if (*j < 0) {
		*j = 0;
		del[2] = 0.;
	} else if (*j > N2 - 2) {
		*j = N2 - 2;
		del[2] = 1.;
	} else {
		del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2]; //fractional displacement of the center of the grid cel
	}
	*k = 0;
	del[3] = 0;
	#if(HAMR3D)
	*k= (int) ((X[3] - startx[3]) / dx[3] - 0.5 + 1000) - 1000;
	if (*k < 0) {
		*k = 0;
		del[3] = 0.;
	} else if (*k > N3 - 2) {
		*k = N3 - 2;
		del[3] = 1.;
	} else {
		del[3] = (X[3] - ((*k + 0.5) * dx[3] + startx[3])) / dx[3]; //fractional displacement of the center of the grid cel
	}
	#endif
	return;
}

/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th)
{
	//fprintf(stderr,"X[1] = %le, X[2] = %le, X[3] = %le \n", X[1], X[2], X[3]);
	*r = exp(X[1]) + R0;
	*th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);

	return;
}



void coord(int i, int j, double *X)
{
	/* returns zone-centered values for coordinates */
	X[0] = startx[0];
	X[1] = startx[1] + (i + 0.5) * dx[1];
	X[2] = startx[2] + (j + 0.5) * dx[2];
	X[3] = startx[3];

	return;
}

void set_units(char *munitstr)
{
	double MBH;

	/* set black hole mass */
	/** could be read in from file here, 
	    along with M_unit and other parameters **/
	//MBH = 4.e6;
	MBH= 10;
	sscanf(munitstr, "%lf", &M_unit);

	/** input parameters appropriate to Sgr A* **/
	MBH *= MSUN;

	/** from this, calculate units of length, time, mass,
	    and derivative units **/
	L_unit = GNEWT * MBH / (CL * CL);
	T_unit = L_unit / CL;

	fprintf(stderr, "\nUNITS\n");
	fprintf(stderr, "L,T,M: %g %g %g\n", L_unit, T_unit, M_unit);

	RHO_unit = M_unit / pow(L_unit, 3);
	U_unit = RHO_unit * CL * CL;
	B_unit = CL * sqrt(4. * M_PI * RHO_unit);

	fprintf(stderr, "rho,u,B: %g %g %g\n", RHO_unit, U_unit, B_unit);

	Ne_unit = RHO_unit / (MP + ME);

	max_tau_scatt = (6. * L_unit) * RHO_unit * 0.4;

	fprintf(stderr, "max_tau_scatt: %g\n", max_tau_scatt);

}

/* set up all grid functions */
void init_geometry()
{
	int i, j, k;
	double X[NDIM];

	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			for (k = 0; k < N3; k++) {

			/* zone-centered */
			#if(HAMR)
			coord_hamr(i, j, k, CENT, X);
			//fprintf(stderr, "i = %d, j = %d, k = %d, X= %lf %lf %lf %lf\n", i, j, k, X[0], X[1], X[2], X[3]);
			gcov_func_hamr(X, geom[SPATIAL_INDEX2D(i,j)].gcov);
			#else
			coord(i, j, X);
			gcov_func(X, geom[SPATIAL_INDEX2D(i,j)].gcov);
			#endif


			geom[SPATIAL_INDEX2D(i,j)].g = gdet_func(geom[SPATIAL_INDEX2D(i,j)].gcov);

			#if(HAMR)
			gcon_func_hamr(geom[SPATIAL_INDEX2D(i,j)].gcov, geom[SPATIAL_INDEX2D(i,j)].gcon);
			//gcon_func(X, geom[SPATIAL_INDEX2D(i,j)].gcon);
			//fprintf(stderr, "gcon = %lf\n", geom[SPATIAL_INDEX2D(i,j)].gcon[3][3]);
			#else
			gcon_func(X, geom[SPATIAL_INDEX2D(i,j)].gcon);
			#endif
			}
		}
	}

	/* done! */
}

/*

	return solid angle between points x2i, x2f
	and over all x3.

*/

double dOmega_func(double x2i, double x2f)
{
	double dO;

	dO = 2. * M_PI *
	    (-cos(M_PI * x2f + 0.5 * (1. - hslope) * sin(2 * M_PI * x2f))
	     + cos(M_PI * x2i + 0.5 * (1. - hslope) * sin(2 * M_PI * x2i))
	    );

	return (dO);
}

void init_storage(void)
{

    p = (double **) malloc(NPRIM * sizeof(double *));
    for (int i = 0; i < NPRIM; i++)
    {
		p[i] = (double *) malloc(N1 * N2 * N3 * sizeof(double));
    }

    geom = (struct of_geom *) malloc(N1* N2* sizeof(struct of_geom));
    return;
}

