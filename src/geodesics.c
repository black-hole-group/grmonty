
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

/*

this is the main photon orbit integrator 

*/

#define FAST_CPY(in,out) {out[0] = in[0]; out[1] = in[1]; out[2] = in[2]; out[3] = in[3];}

#define ETOL 1.e-3
#define MAX_ITER 2

void push_photon(double X[NDIM], double Kcon[NDIM], double dKcon[NDIM],
		 double dl, double *E0, int n)
{
	double lconn[NDIM][NDIM][NDIM];
	double Kcont[NDIM], K[NDIM], dK;
	double Xcpy[NDIM], Kcpy[NDIM], dKcpy[NDIM];
	double Gcov[NDIM][NDIM], E1;
	double dl_2, err, errE;
	int i, k, iter;

	if (X[1] < startx[1])
		return;

	FAST_CPY(X, Xcpy);
	FAST_CPY(Kcon, Kcpy);
	FAST_CPY(dKcon, dKcpy);

	dl_2 = 0.5 * dl;
	/* Step the position and estimate new wave vector */
	// fprintf(stderr, "iteration = %d\n", n);
	// fprintf(stderr, "Inside Push Photon before X[0] = %lf, X[1] = %lf, X[2] = %lf, X[3] = %lf\n", X[0], X[1], X[2], X[3]);
	// fprintf(stderr, "Inside Push Photon before k[0] = %lf, k[1] = %lf, k[2] = %lf, k[3] = %lf\n", K[0], K[1], K[2], K[3]);
	// fprintf(stderr, "dl = %lf, dKcon[0] = %lf, dKcon[1] = %lf, dKcon[2] = %lf, dKcon[3] = %lf\n", dl, dKcon[0], dKcon[1], dKcon[2], dKcon[3]);
	for (i = 0; i < NDIM; i++) {
		dK = dKcon[i] * dl_2;
		Kcon[i] += dK;
		K[i] = Kcon[i] + dK;
		X[i] += Kcon[i] * dl;
	}
	// fprintf(stderr, "Inside Push Photon after X[0] = %lf, X[1] = %lf, X[2] = %lf, X[3] = %lf\n", X[0], X[1], X[2], X[3]);
	// fprintf(stderr, "Inside Push Photon after k[0] = %lf, k[1] = %lf, k[2] = %lf, k[3] = %lf\n", K[0], K[1], K[2], K[3]);

	// if(omp_get_thread_num() == 0){
	// 	fprintf(stderr, "X1 after = %le, K1 after = %le\n", X[1], K[1]);
	// }
	get_connection(X, lconn);

	/* We're in a coordinate basis so take advantage of symmetry in the connection */
	iter = 0;
	do {
		iter++;
		FAST_CPY(K, Kcont);

		err = 0.;
		for (k = 0; k < 4; k++) {
			dKcon[k] =
			    -2. * (Kcont[0] *
				   (lconn[k][0][1] * Kcont[1] +
				    lconn[k][0][2] * Kcont[2] +
				    lconn[k][0][3] * Kcont[3])
				   +
				   Kcont[1] * (lconn[k][1][2] * Kcont[2] +
					       lconn[k][1][3] * Kcont[3])
				   + lconn[k][2][3] * Kcont[2] * Kcont[3]
			    );

			dKcon[k] -=
			    (lconn[k][0][0] * Kcont[0] * Kcont[0] +
			     lconn[k][1][1] * Kcont[1] * Kcont[1] +
			     lconn[k][2][2] * Kcont[2] * Kcont[2] +
			     lconn[k][3][3] * Kcont[3] * Kcont[3]
			    );

			K[k] = Kcon[k] + dl_2 * dKcon[k];
			err += fabs((Kcont[k] - K[k]) / (K[k] + SMALL));
		}
	} while (err > ETOL && iter < MAX_ITER);

	FAST_CPY(K, Kcon);

	#if(HAMR)
	gcov_func_hamr(X, Gcov);
	//gcov_func(X, Gcov);
	#else
	gcov_func(X, Gcov);
	#endif
	E1 = -(Kcon[0] * Gcov[0][0] + Kcon[1] * Gcov[0][1] +
	       Kcon[2] * Gcov[0][2] + Kcon[3] * Gcov[0][3]);
	errE = fabs((E1 - (*E0)) / (*E0));

	if (n < 7
	    && (errE > 1.e-4 || err > ETOL || isnan(err) || isinf(err))) {
		FAST_CPY(Xcpy, X);
		FAST_CPY(Kcpy, Kcon);
		FAST_CPY(dKcpy, dKcon);
		push_photon(X, Kcon, dKcon, 0.5 * dl, E0, n + 1);
		push_photon(X, Kcon, dKcon, 0.5 * dl, E0, n + 1);
		E1 = *E0;
	}

	*E0 = E1;

	/* done! */
}

/* spare photon integrator: 4th order Runge-Kutta */
void push_photon4(double X[], double K[], double dK[], double dl)
{

	int k;
	double lconn[NDIM][NDIM][NDIM];
	double Kt[NDIM], Xt[NDIM];
	double f1x[NDIM], f2x[NDIM], f3x[NDIM], f4x[NDIM];
	double f1k[NDIM], f2k[NDIM], f3k[NDIM], f4k[NDIM];
	double dl_2 = 0.5 * dl;

	for (k = 0; k < NDIM; k++)
		f1x[k] = K[k];

	get_connection(X, lconn);


	for (k = 0; k < NDIM; k++) {
		f1k[k] =
		    -2. * (K[0] *
			   (lconn[k][0][1] * K[1] + lconn[k][0][2] * K[2] +
			    lconn[k][0][3] * K[3]) +
			   K[1] * (lconn[k][1][2] * K[2] +
				   lconn[k][1][3] * K[3]) +
			   lconn[k][2][3] * K[2] * K[3]
		    );

		f1k[k] -=
		    (lconn[k][0][0] * K[0] * K[0] +
		     lconn[k][1][1] * K[1] * K[1] +
		     lconn[k][2][2] * K[2] * K[2] +
		     lconn[k][3][3] * K[3] * K[3]
		    );
	}

	for (k = 0; k < NDIM; k++) {
		Kt[k] = K[k] + dl_2 * f1k[k];
		f2x[k] = Kt[k];
		Xt[k] = X[k] + dl_2 * f1x[k];
	}

	get_connection(Xt, lconn);

	for (k = 0; k < NDIM; k++) {
		f2k[k] =
		    -2. * (Kt[0] *
			   (lconn[k][0][1] * Kt[1] +
			    lconn[k][0][2] * Kt[2] +
			    lconn[k][0][3] * Kt[3]) +
			   Kt[1] * (lconn[k][1][2] * Kt[2] +
				    lconn[k][1][3] * Kt[3]) +
			   lconn[k][2][3] * Kt[2] * Kt[3]
		    );

		f2k[k] -=
		    (lconn[k][0][0] * Kt[0] * Kt[0] +
		     lconn[k][1][1] * Kt[1] * Kt[1] +
		     lconn[k][2][2] * Kt[2] * Kt[2] +
		     lconn[k][3][3] * Kt[3] * Kt[3]
		    );
	}

	for (k = 0; k < NDIM; k++) {
		Kt[k] = K[k] + dl_2 * f2k[k];
		f3x[k] = Kt[k];
		Xt[k] = X[k] + dl_2 * f2x[k];
	}

	get_connection(Xt, lconn);

	for (k = 0; k < NDIM; k++) {
		f3k[k] =
		    -2. * (Kt[0] *
			   (lconn[k][0][1] * Kt[1] +
			    lconn[k][0][2] * Kt[2] +
			    lconn[k][0][3] * Kt[3]) +
			   Kt[1] * (lconn[k][1][2] * Kt[2] +
				    lconn[k][1][3] * Kt[3]) +
			   lconn[k][2][3] * Kt[2] * Kt[3]
		    );

		f3k[k] -=
		    (lconn[k][0][0] * Kt[0] * Kt[0] +
		     lconn[k][1][1] * Kt[1] * Kt[1] +
		     lconn[k][2][2] * Kt[2] * Kt[2] +
		     lconn[k][3][3] * Kt[3] * Kt[3]
		    );
	}

	for (k = 0; k < NDIM; k++) {
		Kt[k] = K[k] + dl * f3k[k];
		f4x[k] = Kt[k];
		Xt[k] = X[k] + dl * f3x[k];
	}

	get_connection(Xt, lconn);

	for (k = 0; k < NDIM; k++) {
		f4k[k] =
		    -2. * (Kt[0] *
			   (lconn[k][0][1] * Kt[1] +
			    lconn[k][0][2] * Kt[2] +
			    lconn[k][0][3] * Kt[3]) +
			   Kt[1] * (lconn[k][1][2] * Kt[2] +
				    lconn[k][1][3] * Kt[3]) +
			   lconn[k][2][3] * Kt[2] * Kt[3]
		    );

		f4k[k] -=
		    (lconn[k][0][0] * Kt[0] * Kt[0] +
		     lconn[k][1][1] * Kt[1] * Kt[1] +
		     lconn[k][2][2] * Kt[2] * Kt[2] +
		     lconn[k][3][3] * Kt[3] * Kt[3]
		    );
	}

	for (k = 0; k < NDIM; k++) {
		X[k] +=
		    0.166666666666667 * dl * (f1x[k] +
					      2. * (f2x[k] + f3x[k]) +
					      f4x[k]);
		K[k] +=
		    0.166666666666667 * dl * (f1k[k] +
					      2. * (f2k[k] + f3k[k]) +
					      f4k[k]);
	}

	init_dKdlam(X, K, dK);

	/* done */
}

void init_dKdlam(double X[], double Kcon[], double dK[])
{
	int k;
	double lconn[NDIM][NDIM][NDIM];

	get_connection(X, lconn);
	//fprintf(stderr, "Inside INIT_DKDLAM after X[0] = %lf, X[1] = %lf, X[2] = %lf, X[3] = %lf\n", X[0], X[1], X[2], X[3]);

	// for(int i = 0; i< 4; i++)
	// for(int j = 0; j < 4; j++)
	// for(int z = 0; z< 4; z++){
	// 	fprintf(stderr, "We are inside init_dKdlam function: lconn[%d][%d][%d] = %le\n", i, j, z, lconn[i][j][z]);
	// }
	for (k = 0; k < 4; k++) {

		dK[k] =
		    -2. * (Kcon[0] *
			   (lconn[k][0][1] * Kcon[1] +
			    lconn[k][0][2] * Kcon[2] +
			    lconn[k][0][3] * Kcon[3])
			   + Kcon[1] * (lconn[k][1][2] * Kcon[2] +
					lconn[k][1][3] * Kcon[3])
			   + lconn[k][2][3] * Kcon[2] * Kcon[3]
		    );

		dK[k] -=
		    (lconn[k][0][0] * Kcon[0] * Kcon[0] +
		     lconn[k][1][1] * Kcon[1] * Kcon[1] +
		     lconn[k][2][2] * Kcon[2] * Kcon[2] +
		     lconn[k][3][3] * Kcon[3] * Kcon[3]
		    );
	}


	return;
}
