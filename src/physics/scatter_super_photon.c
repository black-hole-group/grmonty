
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
	main scattering subroutine 

*/

#include <decs.h>

/* 
	scatter photon ph into photon php at same position 
*/

void scatter_super_photon(struct of_photon *ph, struct of_photon *php,
			  double Ne, double Thetae, double B,
			  double Ucon[NDIM], double Bcon[NDIM],
			  double Gcov[NDIM][NDIM])
{
	double P[NDIM], Econ[NDIM][NDIM], Ecov[NDIM][NDIM],
	    K_tetrad[NDIM], K_tetrad_p[NDIM], Bhatcon[NDIM], tmpK[NDIM];
	int k;

	/* quality control */

	if (isnan(ph->K[1])) {
		fprintf(stderr, "scatter: bad input photon\n");
		exit(0);
	}

	/* quality control */
	if (ph->K[0] > 1.e5 || ph->K[0] < 0. || isnan(ph->K[1])
	    || isnan(ph->K[0]) || isnan(ph->K[3])) {
		fprintf(stderr,
			"normalization problem, killing superphoton: %g \n",
			ph->K[0]);
		ph->K[0] = fabs(ph->K[0]);
		fprintf(stderr, "X1,X2: %g %g\n", ph->X[1], ph->X[2]);
		ph->w = 0.;
		return;
	}

	/* make trial vector for Gram-Schmidt orthogonalization in make_tetrad */
	/* note that B is in cgs but Bcon is in code units */
	if (B > 0.) {
		for (k = 0; k < NDIM; k++)
			Bhatcon[k] = Bcon[k] / (B / B_unit);
	} else {
		for (k = 0; k < NDIM; k++)
			Bhatcon[k] = 0.;
		Bhatcon[1] = 1.;
	}

	/* make local tetrad */
	make_tetrad(Ucon, Bhatcon, Gcov, Econ, Ecov);

	/* transform to tetrad frame */
	coordinate_to_tetrad(Ecov, ph->K, K_tetrad);

	/* quality control */
	if (K_tetrad[0] > 1.e5 || K_tetrad[0] < 0. || isnan(K_tetrad[1])) {
		fprintf(stderr,
			"conversion to tetrad frame problem: %g %g\n",
			ph->K[0], K_tetrad[0]);
/*		fprintf(stderr,"%g %g %g\n",ph->K[1], ph->K[2], ph->K[3]);
		fprintf(stderr,"%g %g %g\n",K_tetrad[1], K_tetrad[2], K_tetrad[3]);
		fprintf(stderr,"%g %g %g %g\n",Ucon[0], Ucon[1], Ucon[2], Ucon[3]);
		fprintf(stderr,"%g %g %g %g\n",Bhatcon[0], Bhatcon[1], Bhatcon[2], Bhatcon[3]);
		fprintf(stderr,"%g %g %g %g\n", Gcov[0][0], Gcov[0][1], Gcov[0][2], Gcov[0][3]) ;
		fprintf(stderr,"%g %g %g %g\n", Gcov[1][0], Gcov[1][1], Gcov[1][2], Gcov[1][3]) ;
		fprintf(stderr,"%g %g %g %g\n", Gcov[2][0], Gcov[2][1], Gcov[2][2], Gcov[2][3]) ;
		fprintf(stderr,"%g %g %g %g\n", Gcov[3][0], Gcov[3][1], Gcov[3][2], Gcov[3][3]) ;
		fprintf(stderr,"%g %g %g %g\n", Ecov[0][0], Ecov[0][1], Ecov[0][2], Ecov[0][3]) ;
		fprintf(stderr,"%g %g %g %g\n", Ecov[1][0], Ecov[1][1], Ecov[1][2], Ecov[1][3]) ;
		fprintf(stderr,"%g %g %g %g\n", Ecov[2][0], Ecov[2][1], Ecov[2][2], Ecov[2][3]) ;
		fprintf(stderr,"%g %g %g %g\n", Ecov[3][0], Ecov[3][1], Ecov[3][2], Ecov[3][3]) ;
		fprintf(stderr,"X1,X2: %g %g\n",ph->X[1],ph->X[2]) ;*/
		ph->w = 0.;
		return;
	}

	/* find the electron that we collided with */
	sample_electron_distr_p(K_tetrad, P, Thetae);

	/* given electron momentum P, find the new
	   photon momentum Kp */
	sample_scattered_photon(K_tetrad, P, K_tetrad_p);


	/* transform back to coordinate frame */
	tetrad_to_coordinate(Econ, K_tetrad_p, php->K);

	/* quality control */
	if (isnan(php->K[1])) {
		fprintf(stderr,
			"problem with conversion to coordinate frame\n");
		fprintf(stderr, "%g %g %g %g\n", Econ[0][0], Econ[0][1],
			Econ[0][2], Econ[0][3]);
		fprintf(stderr, "%g %g %g %g\n", Econ[1][0], Econ[1][1],
			Econ[1][2], Econ[1][3]);
		fprintf(stderr, "%g %g %g %g\n", Econ[2][0], Econ[2][1],
			Econ[2][2], Econ[2][3]);
		fprintf(stderr, "%g %g %g %g\n", Econ[3][0], Econ[3][1],
			Econ[3][2], Econ[3][3]);
		fprintf(stderr, "%g %g %g %g\n", K_tetrad_p[0],
			K_tetrad_p[1], K_tetrad_p[2], K_tetrad_p[3]);
		php->w = 0;
		return;
	}

	if (php->K[0] < 0) {
		fprintf(stderr, "K0, K0p, Kp, P[0]: %g %g %g %g\n",
			K_tetrad[0], K_tetrad_p[0], php->K[0], P[0]);
		php->w = 0.;
		return;
	}

	/* bookkeeping */
	K_tetrad_p[0] *= -1.;
	tetrad_to_coordinate(Ecov, K_tetrad_p, tmpK);

	php->E = php->E0s = -tmpK[0];
	php->L = tmpK[3];
	php->tau_abs = 0.;
	php->tau_scatt = 0.;
	php->b0 = B;

	php->X1i = ph->X[1];
	php->X2i = ph->X[2];
	php->X[0] = ph->X[0];
	php->X[1] = ph->X[1];
	php->X[2] = ph->X[2];
	php->X[3] = ph->X[3];
	php->ne0 = Ne;
	php->thetae0 = Thetae;
	php->E0 = ph->E;
	php->nscatt = ph->nscatt + 1;

	return;
}
