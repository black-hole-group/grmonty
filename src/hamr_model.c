
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
#include "harm_model.h"
#include <math.h>


/*******************Functions used to calculate the grid*******************/
/* Matthews bl_coord*/
void vofx_matthewcoords(double *X, double *V){
	V[0] = X[0];
	double RTRANS =5000000.;
	double RB = 0.;
	double RADEXP = 1.0;
	double Xtrans = pow(log(RTRANS - RB), 1. / RADEXP);
	double BRAVO = 0.0;
	double TANGO = 1.0;
	double CHARLIE = 0.0;
	double DELTA = 3.0;
	if (X[1] < Xtrans){
		V[1] = exp(pow(X[1], RADEXP)) + RB;
	}
	else if (X[1] >= Xtrans && X[1]<1.01*Xtrans){
		V[1] = 10.*(X[1] / Xtrans - 1.)*((X[1] - Xtrans)*RADEXP*exp(pow(Xtrans, RADEXP))*pow(Xtrans, -1. + RADEXP) + RTRANS) +
			(1. - 10.*(X[1] / Xtrans - 1.))*(exp(pow(X[1], RADEXP)) + RB);
	}
	else{
		V[1] = (X[1] - Xtrans)*RADEXP*exp(pow(Xtrans, RADEXP))*pow(Xtrans, -1. + RADEXP) + RTRANS;
	}
	double A1 = 1. / (1. + pow(CHARLIE*(log(V[1]) / log(10.)), DELTA));
	double A2 = BRAVO*(log(V[1]) / log(10.)) + TANGO;
	double A3 = pow(0.5, 1. - A2);
	double sign = 1.;
	double X_2 =(X[2]+1.0)/2.0;
	double Xc = sqrt(pow(X_2, 2.));

	if (X_2 < 0.0){
		sign = -1.;
	}
	if (X_2 > 1.0){
		sign = -1.;
		Xc = 2. - Xc;
	}
	if (X_2 >= 0.5){
		Xc = 1. - Xc;
		V[2] = M_PI - sign*(A1* M_PI*Xc + M_PI*(1. - A1)*(A3*pow(Xc, A2) + 0.50 / M_PI*sin(M_PI + 2.*M_PI*(A3*pow(Xc, A2)))));
	}
	else{
		V[2] = sign*(A1* M_PI*Xc + M_PI*(1. - A1)*(A3*pow(Xc, A2) + 0.50 / M_PI*sin(M_PI + 2.*M_PI*(A3*pow(Xc, A2)))));
	}
	V[3] = X[3];
}

void bl_coord_hamr(double * X, double * r, double *th, double *phi)
{
	double V[4];
	double SINGSMALL = 1.e-20;
  	void (*vofx_function_pointer)(double*, double*);
    vofx_function_pointer = vofx_matthewcoords;
	vofx_function_pointer(X,V);
	// avoid singularity at polar axis
	if (fabs(V[2])<SINGSMALL){
		if (V[2] >= 0.0) V[2] = SINGSMALL;
		if (V[2]<0.0)  V[2] = -SINGSMALL;
	}
	if (fabs(M_PI - V[2]) <SINGSMALL){
		if (V[2] >= M_PI) V[2] = M_PI + SINGSMALL;
		if (V[2]<M_PI)  V[2] = M_PI -  SINGSMALL;
	}
	*r = V[1];
	*th = V[2];
	*phi = V[3];
	//fprintf(stderr, "bl_coord: V[1] = %le, V[2] = %le, V[3] = %le\n", V[1], V[2], V[3]);
	return ;
}

int LU_decompose( double A[][NDIM], int permute[] )
{
  double row_norm[NDIM];

  double absmin = 1.e-30; /* Value used instead of 0 for singular matrices */

  double  absmax, maxtemp, mintemp;

  int i, j, k, max_row;
  int n = NDIM;


  max_row = 0;

  /* Find the maximum elements per row so that we can pretend later
     we have unit-normalized each equation: */

  for( i = 0; i < n; i++ ) { 
    absmax = 0.;
    
    for( j = 0; j < n ; j++ ) { 
      
      maxtemp = fabs( A[i][j] ); 

      if( maxtemp > absmax ) { 
	absmax = maxtemp; 
      }
    }

    /* Make sure that there is at least one non-zero element in this row: */
    if( absmax == 0. ) { 
     //fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
      return(1);
    }

    row_norm[i] = 1. / absmax ;   /* Set the row's normalization factor. */
  }


  /* The following the calculates the matrix composed of the sum 
     of the lower (L) tridagonal matrix and the upper (U) tridagonal
     matrix that, when multiplied, form the original maxtrix.  
     This is what we call the LU decomposition of the maxtrix. 
     It does this by a recursive procedure, starting from the 
     upper-left, proceding down the column, and then to the next
     column to the right.  The decomposition can be done in place 
     since element {i,j} require only those elements with {<=i,<=j} 
     which have already been computed.
     See pg. 43-46 of "Num. Rec." for a more thorough description. 
  */

  /* For each of the columns, starting from the left ... */
  for( j = 0; j < n; j++ ) {

    /* For each of the rows starting from the top.... */

    /* Calculate the Upper part of the matrix:  i < j :   */
    for( i = 0; i < j; i++ ) {
      for( k = 0; k < i; k++ ) { 
	A[i][j] -= A[i][k] * A[k][j];
      }
    }

    absmax = 0.0;

    /* Calculate the Lower part of the matrix:  i <= j :   */

    for( i = j; i < n; i++ ) {

      for (k = 0; k < j; k++) { 
	A[i][j] -= A[i][k] * A[k][j];
      }

      /* Find the maximum element in the column given the implicit 
	 unit-normalization (represented by row_norm[i]) of each row: 
      */
      maxtemp = fabs(A[i][j]) * row_norm[i] ;

      if( maxtemp >= absmax ) {
	absmax = maxtemp;
	max_row = i;
      }

    }

    /* Swap the row with the largest element (of column j) with row_j.  absmax
       This is the partial pivoting procedure that ensures we don't divide
       by 0 (or a small number) when we solve the linear system.  
       Also, since the procedure starts from left-right/top-bottom, 
       the pivot values are chosen from a pool involving all the elements 
       of column_j  in rows beneath row_j.  This ensures that 
       a row  is not permuted twice, which would mess things up. 
    */
    if( max_row != j ) {

      /* Don't swap if it will send a 0 to the last diagonal position. 
	 Note that the last column cannot pivot with any other row, 
	 so this is the last chance to ensure that the last two 
	 columns have non-zero diagonal elements.
       */

      if( (j == (n-2)) && (A[j][j+1] == 0.) ) {
	max_row = j;
      }
      else { 
	for( k = 0; k < n; k++ ) { 

	  maxtemp       = A[   j   ][k] ; 
	  A[   j   ][k] = A[max_row][k] ;
	  A[max_row][k] = maxtemp; 

	}

	/* Don't forget to swap the normalization factors, too... 
	   but we don't need the jth element any longer since we 
	   only look at rows beneath j from here on out. 
	*/
	row_norm[max_row] = row_norm[j] ; 
      }
    }

    /* Set the permutation record s.t. the j^th element equals the 
       index of the row swapped with the j^th row.  Note that since 
       this is being done in successive columns, the permutation
       vector records the successive permutations and therefore
       index of permute[] also indexes the chronology of the 
       permutations.  E.g. permute[2] = {2,1} is an identity 
       permutation, which cannot happen here though. 
    */

    permute[j] = max_row;

    if( A[j][j] == 0. ) { 
      A[j][j] = absmin;
    }


  /* Normalize the columns of the Lower tridiagonal part by their respective 
     diagonal element.  This is not done in the Upper part because the 
     Lower part's diagonal elements were set to 1, which can be done w/o 
     any loss of generality.
  */
    if( j != (n-1) ) { 
      maxtemp = 1. / A[j][j]  ;
      
      for( i = (j+1) ; i < n; i++ ) {
	A[i][j] *= maxtemp;
      }
    }

  }

  return(0);

  /* End of LU_decompose() */

}

void LU_substitution( double A[][NDIM], double B[], int permute[] )
{
  int i, j ;
  int n = NDIM;
  double tmpvar,tmpvar2;

  
  /* Perform the forward substitution using the LU matrix. 
   */
  for(i = 0; i < n; i++) {

    /* Before doing the substitution, we must first permute the 
       B vector to match the permutation of the LU matrix. 
       Since only the rows above the currrent one matter for 
       this row, we can permute one at a time. 
    */
    tmpvar        = B[permute[i]];
    B[permute[i]] = B[    i     ];
    for( j = (i-1); j >= 0 ; j-- ) { 
      tmpvar -=  A[i][j] * B[j];
    }
    B[i] = tmpvar; 
  }
	   

  /* Perform the backward substitution using the LU matrix. 
   */
  for( i = (n-1); i >= 0; i-- ) { 
    for( j = (i+1); j < n ; j++ ) { 
      B[i] -=  A[i][j] * B[j];
    }
    B[i] /= A[i][i] ; 
  }

  /* End of LU_substitution() */

}

int invert_matrix( double Am[][NDIM], double Aminv[][NDIM] )  
{ 

  int i,j;
  int n = NDIM;
  int permute[NDIM]; 
  double dxm[NDIM], Amtmp[NDIM][NDIM];

  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  Amtmp[0][i] = Am[0][i]; }

  // Get the LU matrix:
  if( LU_decompose( Amtmp,  permute ) != 0  ) { 
    fprintf(stderr, "invert_matrix(): singular matrix encountered! \n");
    fprintf(stderr, "This is probably due to a nan value somewhere rather than determinant = 0. Investigate!\n");
	return(1);
  }

  for( i = 0; i < n; i++ ) { 
    for( j = 0 ; j < n ; j++ ) { dxm[j] = 0. ; }
    dxm[i] = 1.; 
    
    /* Solve the linear system for the i^th column of the inverse matrix: :  */
    LU_substitution( Amtmp,  dxm, permute );

    for( j = 0 ; j < n ; j++ ) {  Aminv[j][i] = dxm[j]; }

  }

  return(0);
}

/*This function has been tested against dxdxp in HAMR's python notebook. Seems to be working good*/
void dxdxp_func(double *X, double dxdxp[][NDIM])
{
	int i, j, k, l;
	double Xh[NDIM], Xl[NDIM];
	double Vh[NDIM], Vl[NDIM];
	//fprintf(stderr, "X variables entered function X[0] = %le, X[1] = %le, X[2] = %le, X[3] = %le \n",X[0], X[1], X[2], X[3]);
	for (k = 0; k<NDIM; k++) {
		for (l = 0; l<NDIM; l++) Xh[l] = X[l];
		for (l = 0; l<NDIM; l++) Xl[l] = X[l];
		Xh[k] += 0.00001;
		Xl[k] -= 0.00001;
		Vh[0] = Xh[0];
		Vl[0] = Xl[0];
		//fprintf(stderr, "Vh[0] = %le, Vl[0] = %le, Xh[0] = %le, Xl[0] = %le\n", Vh[0], Vl[0], Xh[0], Xl[0]);
		bl_coord_hamr(Xh, &Vh[1], &Vh[2], &Vh[3]);
		//fprintf(stderr, "Vh[1] = %le, Vh[2] = %le, Vh[3] = %le\n", Vh[1], Vh[2], Vh[3]);
		bl_coord_hamr(Xl, &Vl[1], &Vl[2], &Vl[3]);
		//fprintf(stderr, "Vl[1] = %le, Vl[2] = %le, Vl[3] = %le\n", Vl[1], Vl[2], Vl[3]);

		for (j = 0; j<NDIM; j++)
		{
			dxdxp[j][k] = (Vh[j] - Vl[j]) / (Xh[k] - Xl[k]);
			//fprintf(stderr, "dxdxp[%d][%d] = %le\n", j, k, dxdxp[j][k]);
		}
	}
}

/*This function has been tested against gcov in HAMR's python notebook. Seems to be working good*/
void gcov_func_hamr(double *X, double gcovp[][NDIM])
{
	int i, j, k, l;
	double sth, cth, s2, rho2;
	double del[NDIM];
	double r, th, phi;
	double a = 0.9375;
	double tilt = TILT_ANGLE / 180.*M_PI;
	//fprintf(stderr, "X = %lf, %lf, %lf, %lf\n", X[0], X[1], X[2], X[3]);

	for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++) gcovp[j][k] = 0.;
	bl_coord_hamr(X, &r, &th, &phi);
	cth = cos(th);
	sth = sin(th);

	s2 = sth*sth;
	rho2 = r*r + a*a*cth*cth;

	//compute Jacobian x1,x2,x3 -> r,th,phi (dr/dx1)
	gcovp[0][0] = (-1. + 2.*r / rho2);
	gcovp[0][1] = (2.*r / rho2)*r;
	gcovp[0][3] = (-2.*a*r*s2 / rho2);

	gcovp[1][0] = gcovp[0][1];
	gcovp[1][1] = (1. + 2.*r / rho2)*r*r;
	gcovp[1][3] = (-a*s2*(1. + 2.*r / rho2)) * r;

	gcovp[2][2] = rho2 * (M_PI/2) * (M_PI/2);

	gcovp[3][0] = gcovp[0][3];
	gcovp[3][1] = gcovp[1][3];
	gcovp[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r / rho2));
	//for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++) fprintf(stderr, "gcovp[%d][%d], %lf\n",j,k,gcovp[j][k]);


}

/* invert gcov to get gcon */
/*This function has been tested against gcon in HAMR's python notebook. Seems to be working good*/
void gcon_func_hamr(double gcov[][NDIM], double gcon[][NDIM])
{
  invert_matrix( gcov, gcon );
}

// void gcon_func_hamr(double *X, double gcon[][NDIM])
// {

// 	int k, l;
// 	double sth, cth, irho2;
// 	double r, th, phi;
// 	double hfac;
// 	/* required by broken math.h */
// 	void sincos(double in, double *sth, double *cth);

// 	DLOOP gcon[k][l] = 0.;
// 	#if(HAMR)
// 	bl_coord_hamr(X, &r, &th, &phi);
// 	#else
// 	bl_coord(X, &r, &th);
// 	#endif


// 	sincos(th, &sth, &cth);
// 	sth = fabs(sth) + SMALL;

// 	irho2 = 1. / (r * r + a * a * cth * cth);

// 	// transformation for Kerr-Schild -> modified Kerr-Schild 
// 	//hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
// 	hfac = 1;

// 	gcon[0][0] = -1. - 2. * r * irho2;
// 	gcon[0][1] = 2. * irho2;

// 	gcon[1][0] = gcon[0][1];
// 	gcon[1][1] = irho2 * (r * (r - 2.) + a * a) / (r * r);
// 	gcon[1][3] = a * irho2 / r;

// 	gcon[2][2] = irho2 / (hfac * hfac);

// 	gcon[3][1] = gcon[1][3];
// 	gcon[3][3] = irho2 / (sth * sth);
// }


void coord_hamr(int i, int j, int z, int loc, double * X)
{
	X[0] = 0.0;
	int j_local = j;
	if (j < 0) j_local = -j - 1;
	//if (j >= N2*pow(1 + REF_2, block[n][AMR_LEVEL2])) j_local = 2 * N2*pow(1 + REF_2, block[n][AMR_LEVEL2]) - 1 - j;
	//if (j == N2*pow(1 + REF_2, block[n][AMR_LEVEL2]) && loc == FACE2) j_local = j;
	if (j >= N2*pow(1 + REF_2, 0)) j_local = 2 * N2*pow(1 + REF_2, 0) - 1 - j;
	if (j == N2*pow(1 + REF_2, 0) && loc == FACE2) j_local = j;
	if (loc == FACE1) {
		X[1] = startx[1] + i*dx[1];
		X[2] = startx[2] + (j_local + 0.5)*dx[2];
		X[3] = startx[3] + (z + 0.5)*dx[3];
	}
	else if (loc == FACE2) {
		X[1] = startx[1] + (i + 0.5)*dx[1];
		X[2] = startx[2] + j_local*dx[2];
		X[3] = startx[3] + (z + 0.5)*dx[3];
	}
	else if (loc == FACE3) {
		X[1] = startx[1] + (i + 0.5)*dx[1];
		X[2] = startx[2] + (j_local + 0.5)*dx[2];
		X[3] = startx[3] + z*dx[3];
	}
	else if (loc == CENT) {
		X[1] = startx[1] + (i + 0.5)*dx[1];
		X[2] = startx[2] + (j_local + 0.5)*dx[2];
		X[3] = startx[3] + (z + 0.5)*dx[3]; 
		//X[3] = 0;
	}
	else {
		X[1] = startx[1] + i*dx[1];
		X[2] = startx[2] + j_local*dx[2];
		X[3] = startx[3] + z*dx[3];
	}

	#if(!CARTESIAN)
	if (j < 0){
		X[2] = X[2] + 1;
		X[2] = -X[2];
		X[2] = X[2] - 1;
	}
	if (j == N2*pow(1 + REF_2, 0) && loc == FACE2){
	}
	else if (j >= N2*pow(1 + REF_2, 0)){
		X[2] = X[2] + 1;
		X[2] = 4. - X[2];
		X[2] = X[2] - 1;
	}
	#endif
	//if (j == N2 * pow(1 + REF_2, block[n][AMR_LEVEL2]))fprintf(stderr, "test: %f %d \n" ,X[2], loc==FACE2);
	return;
}


