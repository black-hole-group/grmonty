
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

all functions related to creation and manipulation of tetrads

*/

#include <decs.h>


/* input and vectors are contravariant (index up) */
void coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM],
			  double K_tetrad[NDIM])
{
	int k;

	for (k = 0; k < 4; k++) {
		K_tetrad[k] =
		    Ecov[k][0] * K[0] +
		    Ecov[k][1] * K[1] +
		    Ecov[k][2] * K[2] + Ecov[k][3] * K[3];
	}
}

/* input and vectors are contravariant (index up) */
void tetrad_to_coordinate(double Econ[NDIM][NDIM], double K_tetrad[NDIM],
			  double K[NDIM])
{
	int l;

	for (l = 0; l < 4; l++) {
		K[l] = Econ[0][l] * K_tetrad[0] +
		    Econ[1][l] * K_tetrad[1] +
		    Econ[2][l] * K_tetrad[2] + Econ[3][l] * K_tetrad[3];
	}

	return;
}

#define SMALL_VECTOR	1.e-30

/* make orthonormal basis 
   first basis vector || U
   second basis vector || B
*/
void make_tetrad(double Ucon[NDIM], double trial[NDIM],
		 double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
		 double Ecov[NDIM][NDIM])
{
	int k, l;
	double norm;
	void normalize(double *vcon, double Gcov[4][4]);
	void project_out(double *vcona, double *vconb, double Gcov[4][4]);

	/* econ/ecov index explanation:
	   Econ[k][l]
	   k: index attached to tetrad basis
	   index down
	   l: index attached to coordinate basis 
	   index up
	   Ecov[k][l]
	   k: index attached to tetrad basis
	   index up
	   l: index attached to coordinate basis 
	   index down
	 */

	/* start w/ time component parallel to U */
	for (k = 0; k < 4; k++)
		Econ[0][k] = Ucon[k];
	normalize(Econ[0], Gcov);

	/*** done w/ basis vector 0 ***/

	/* now use the trial vector in basis vector 1 */
	/* cast a suspicious eye on the trial vector... */
	norm = 0.;
	for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
			norm += trial[k] * trial[l] * Gcov[k][l];
	if (norm <= SMALL_VECTOR) {	/* bad trial vector; default to radial direction */
		for (k = 0; k < 4; k++)	/* trial vector */
			trial[k] = delta(k, 1);
	}

	for (k = 0; k < 4; k++)	/* trial vector */
		Econ[1][k] = trial[k];

	/* project out econ0 */
	project_out(Econ[1], Econ[0], Gcov);
	normalize(Econ[1], Gcov);

	/*** done w/ basis vector 1 ***/

	/* repeat for x2 unit basis vector */
	for (k = 0; k < 4; k++)	/* trial vector */
		Econ[2][k] = delta(k, 2);
	/* project out econ[0-1] */
	project_out(Econ[2], Econ[0], Gcov);
	project_out(Econ[2], Econ[1], Gcov);
	normalize(Econ[2], Gcov);

	/*** done w/ basis vector 2 ***/

	/* and repeat for x3 unit basis vector */
	for (k = 0; k < 4; k++)	/* trial vector */
		Econ[3][k] = delta(k, 3);
	/* project out econ[0-2] */
	project_out(Econ[3], Econ[0], Gcov);
	project_out(Econ[3], Econ[1], Gcov);
	project_out(Econ[3], Econ[2], Gcov);
	normalize(Econ[3], Gcov);

	/*** done w/ basis vector 3 ***/

	/* now make covariant version */
	for (k = 0; k < 4; k++) {

		/* lower coordinate basis index */
		lower(Econ[k], Gcov, Ecov[k]);
	}

	/* then raise tetrad basis index */
	for (l = 0; l < 4; l++) {
		Ecov[0][l] *= -1.;
	}

	/* paranoia: check orthonormality */
	/*
	   double sum ;
	   int m ;
	   fprintf(stderr,"ortho check:\n") ;
	   for(k=0;k<NDIM;k++)
	   for(l=0;l<NDIM;l++) {
	   sum = 0. ;
	   for(m=0;m<NDIM;m++) {
	   sum += Econ[k][m]*Ecov[l][m] ;
	   }
	   fprintf(stderr,"%d %d %g\n",k,l,sum) ;
	   }
	   fprintf(stderr,"\n") ;
	   for(k=0;k<NDIM;k++)
	   for(l=0;l<NDIM;l++) {
	   fprintf(stderr,"%d %d %g\n",k,l,Econ[k][l]) ;
	   }
	   fprintf(stderr,"\n") ;
	 */


	/* done */

}

double delta(int i, int j)
{
	if (i == j)
		return (1.);
	else
		return (0.);
}

void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov)
{

	ucov[0] = Gcov[0][0] * ucon[0]
	    + Gcov[0][1] * ucon[1]
	    + Gcov[0][2] * ucon[2]
	    + Gcov[0][3] * ucon[3];
	ucov[1] = Gcov[1][0] * ucon[0]
	    + Gcov[1][1] * ucon[1]
	    + Gcov[1][2] * ucon[2]
	    + Gcov[1][3] * ucon[3];
	ucov[2] = Gcov[2][0] * ucon[0]
	    + Gcov[2][1] * ucon[1]
	    + Gcov[2][2] * ucon[2]
	    + Gcov[2][3] * ucon[3];
	ucov[3] = Gcov[3][0] * ucon[0]
	    + Gcov[3][1] * ucon[1]
	    + Gcov[3][2] * ucon[2]
	    + Gcov[3][3] * ucon[3];

	return;
}

void normalize(double *vcon, double Gcov[NDIM][NDIM])
{
	int k, l;
	double norm;

	norm = 0.;
	for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
			norm += vcon[k] * vcon[l] * Gcov[k][l];

	norm = sqrt(fabs(norm));
	for (k = 0; k < 4; k++)
		vcon[k] /= norm;

	return;
}

void project_out(double *vcona, double *vconb, double Gcov[NDIM][NDIM])
{

	double adotb, vconb_sq;
	int k, l;

	vconb_sq = 0.;
	for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
			vconb_sq += vconb[k] * vconb[l] * Gcov[k][l];

	adotb = 0.;
	for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
			adotb += vcona[k] * vconb[l] * Gcov[k][l];

	for (k = 0; k < 4; k++)
		vcona[k] -= vconb[k] * adotb / vconb_sq;

	return;
}

void normalize_null(double Gcov[NDIM][NDIM], double K[])
{
	int k, l;
	double A, B, C;

	/* pop K back onto the light cone */
	A = Gcov[0][0];
	B = 0.;
	for (k = 1; k < 4; k++)
		B += 2. * Gcov[k][0] * K[k];
	C = 0.;
	for (k = 1; k < 4; k++)
		for (l = 1; l < 4; l++)
			C += Gcov[k][l] * K[k] * K[l];

	K[0] = (-B - sqrt(fabs(B * B - 4. * A * C))) / (2. * A);

	return;
}
