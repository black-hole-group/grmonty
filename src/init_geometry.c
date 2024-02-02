
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

/* 

   set up the metric (indicies both up and down), and the
   connection coefficients on a grid, based
   on a HARM simulation grid.  
   
   In principle the geometry could be sampled on a finer or
   coarser grid than the simulation grid. 

   These routines are taken directly out of HARM.

   They require the code be compiled against the 
   Gnu scientific library (GSL).

   CFG 21 July 06
   
*/

gsl_matrix *gsl_gcov, *gsl_gcon;
gsl_permutation *perm;
#pragma omp threadprivate (gsl_gcov, gsl_gcon, perm)

/* assumes gcov has been set first; returns determinant */
double gdet_func(double gcov[][NDIM])
{
	double d;
	int k, l, signum;

	if (gsl_gcov == NULL) {
		gsl_gcov = gsl_matrix_alloc(NDIM, NDIM);
		gsl_gcon = gsl_matrix_alloc(NDIM, NDIM);
		perm = gsl_permutation_alloc(NDIM);
	}

	DLOOP gsl_matrix_set(gsl_gcov, k, l, gcov[k][l]);

	gsl_linalg_LU_decomp(gsl_gcov, perm, &signum);

	d = gsl_linalg_LU_det(gsl_gcov, signum);

	return (sqrt(fabs(d)));
}


/* invert gcov to get gcon */

/*
void gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
	int k, l, signum;
	
	if (gsl_gcov  == NULL) {
		gsl_gcov = gsl_matrix_alloc(NDIM, NDIM);
		gsl_gcon = gsl_matrix_alloc(NDIM, NDIM);
		perm = gsl_permutation_alloc(NDIM);
	}
	

	DLOOP gsl_matrix_set(gsl_gcov, k, l, gcov[k][l]);

	gsl_linalg_LU_decomp(gsl_gcov, perm, &signum);

	gsl_linalg_LU_invert(gsl_gcov, perm, gsl_gcon);

	DLOOP gcon[k][l] = gsl_matrix_get(gsl_gcon, k, l);

	// done!
}
*/

#undef DELTA
