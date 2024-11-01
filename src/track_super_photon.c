
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

	main transport subroutine for tracking, absorbing,
	and scattering superphotons

	assumes superphotons do not step out of simulation then back in

*/

#include "decs.h"
#include "harm_model.h"


#define MAXNSTEP	1280000
int photon_scattered = 0;

void track_super_photon(struct of_photon *ph, int round, int photon_index)
{
	int bound_flag;
	double dtau_scatt, dtau_abs, dtau;
	double bi, bf;
	double alpha_scatti, alpha_scattf;
	double alpha_absi, alpha_absf;
	double dl, x1;
	double nu, Thetae, Ne, B, theta;
	struct of_photon php;
	double dtauK, frac;
	double bias = 0.;
	double Xi[NDIM], Ki[NDIM], dKi[NDIM], E0;
	double Gcov[NDIM][NDIM], Ucon[NDIM], Ucov[NDIM], Bcon[NDIM],
	    Bcov[NDIM];
	int nstep = 0;
	double var_id = 0;
	/* quality control */
	if (isnan(ph->X[0]) ||
	    isnan(ph->X[1]) ||
	    isnan(ph->X[2]) ||
	    isnan(ph->X[3]) ||
	    isnan(ph->K[0]) ||
	    isnan(ph->K[1]) ||
	    isnan(ph->K[2]) || isnan(ph->K[3]) || ph->w == 0.) {
		fprintf(stderr, "track_super_photon: bad input photon.\n");
		fprintf(stderr,
			"X0,X1,X2,X3,K0,K1,K2,K3,w,nscatt: %g %g %g %g %g %g %g %g %g %d\n",
			ph->X[0], ph->X[1], ph->X[2], ph->X[3], ph->K[0],
			ph->K[1], ph->K[2], ph->K[3], ph->w, ph->nscatt);
		return;
	}

	dtauK = 2. * M_PI * L_unit / (ME * CL * CL / HBAR);

	/* Initialize opacities */
	#if(HAMR)
	gcov_func_hamr(ph->X, Gcov);
	//gcov_func(ph->X, Gcov);
	#else
	gcov_func(ph->X, Gcov);
	#endif

	get_fluid_params(ph->X, Gcov, &Ne, &Thetae, &B, Ucon, Ucov, Bcon,
			 Bcov);

	theta = get_bk_angle(ph->X, ph->K, Ucov, Bcov, B);
	nu = get_fluid_nu(ph->X, ph->K, Ucov);
	alpha_scatti = alpha_inv_scatt(nu, Thetae, Ne);
	alpha_absi = alpha_inv_abs(nu, Thetae, Ne, B, theta);
	bi = bias_func(Thetae, ph->w);

	/* Initialize dK/dlam */
	init_dKdlam(ph->X, ph->K, ph->dKdlam);

	if(photon_index == 4821 && round == 0)
	{
		printf("alpha_Scatti = %le, ph->E0 = %le, nu = %le\n", alpha_scatti, ph->E0s, nu);
	}

	while (!stop_criterion(ph)) {
		/* Save initial position/wave vector */
		Xi[0] = ph->X[0];
		Xi[1] = ph->X[1];
		Xi[2] = ph->X[2];
		Xi[3] = ph->X[3];
		Ki[0] = ph->K[0];
		Ki[1] = ph->K[1];
		Ki[2] = ph->K[2];
		Ki[3] = ph->K[3];
		dKi[0] = ph->dKdlam[0];
		dKi[1] = ph->dKdlam[1];
		dKi[2] = ph->dKdlam[2];
		dKi[3] = ph->dKdlam[3];
		E0 = ph->E0s;

		/* evaluate stepsize */
		dl = stepsize(ph->X, ph->K);

		/* step the geodesic */
		// fprintf(stderr, "before pushing photon: X[0] = %lf, X[1] = %lf, X[2] = %lf, X[3] = %lf\n", ph->X[0], ph->X[1], ph->X[2], ph->X[3]);
		// fprintf(stderr, "Outside Push Photon function dl = %lf, dKcon[0] = %lf, dKcon[1] = %lf, dKcon[2] = %lf, dKcon[3] = %lf\n", dl, ph->dKdlam[0], ph->dKdlam[1], ph->dKdlam[2], ph->dKdlam[3]);

		push_photon(ph->X, ph->K, ph->dKdlam, dl, &(ph->E0s), 0);

		if (stop_criterion(ph))
			break;

		/* allow photon to interact with matter, */
		#if(HAMR)
		gcov_func_hamr(ph->X, Gcov);
		#else
		gcov_func(ph->X, Gcov);
		#endif
		get_fluid_params(ph->X, Gcov, &Ne, &Thetae, &B, Ucon, Ucov,
				 Bcon, Bcov);
		if(!(alpha_absi > 0. || alpha_scatti > 0. || Ne > 0.) && photon_index == 4821 && round == 0){
			printf("one of them is negative, alpha_absi = %le, alpha_scatti = %le, Ne = %le\n", alpha_absi, alpha_scatti, Ne);
			exit(1);
		}
		if (alpha_absi > 0. || alpha_scatti > 0. || Ne > 0.) {

			bound_flag = 0;
			if (Ne == 0.){
				bound_flag = 1;
				if(photon_index == 4821 && round == 0){
					printf("Ne = %le, boundflag reached!, (%le, %le, %le)\n", Ne, ph->X[1], ph->X[2], ph->X[3]);
				}
			}
			if (!bound_flag) {
				theta =
				    get_bk_angle(ph->X, ph->K, Ucov, Bcov,
						 B);
				nu = get_fluid_nu(ph->X, ph->K, Ucov);
				if (isnan(nu)) {
					fprintf(stderr,
						"isnan nu: track_super_photon dl,E0 %g %g\n",
						dl, E0);
					fprintf(stderr,
						"Xi, %g %g %g %g\n", Xi[0],
						Xi[1], Xi[2], Xi[3]);
					fprintf(stderr,
						"Ki, %g %g %g %g\n", Ki[0],
						Ki[1], Ki[2], Ki[3]);
					fprintf(stderr,
						"dKi, %g %g %g %g\n",
						dKi[0], dKi[1], dKi[2],
						dKi[3]);
					exit(1);
				}
			}

			/* scattering optical depth along step */
			if (bound_flag || nu < 0.) {
				dtau_scatt =
				    0.5 * alpha_scatti * dtauK * dl;
				dtau_abs = 0.5 * alpha_absi * dtauK * dl;
				alpha_scatti = alpha_absi = 0.;
				bias = 0.;
				bi = 0.;
			} else {
				alpha_scattf =
				    alpha_inv_scatt(nu, Thetae, Ne);
				dtau_scatt =
				    0.5 * (alpha_scatti +
					   alpha_scattf) * dtauK * dl;
				alpha_scatti = alpha_scattf;

				/* absorption optical depth along step */
				alpha_absf =
				    alpha_inv_abs(nu, Thetae, Ne, B,
						  theta);
				dtau_abs =
				    0.5 * (alpha_absi +
					   alpha_absf) * dtauK * dl;
				alpha_absi = alpha_absf;

				bf = bias_func(Thetae, ph->w);
				bias = 0.5 * (bi + bf);
				bi = bf;
			}

			
			x1 = -log(monty_rand());
			x1 = -log(0.99);
			php.w = ph->w / bias;
			if(round == 0 && photon_index == 4821){
				printf("dtau_scatt = %le, photon_index = %d, x1 = %le, alpha_absi = %le, alpha_scatti = %le, Ne = %le\n", dtau_scatt, photon_index, x1, alpha_absi, alpha_scatti, Ne);
				
			}
			if (bias * dtau_scatt > x1 && php.w > WEIGHT_MIN) {
				if (isnan(php.w) || isinf(php.w)) {
					fprintf(stderr,
						"w isnan in track_super_photon: Ne, bias, ph->w, php.w  %g, %g, %g, %g\n",
						Ne, bias, ph->w, php.w);
				}

				frac = x1 / (bias * dtau_scatt);

				/* Apply absorption until scattering event */
				dtau_abs *= frac;
				if (dtau_abs > 100)
					return;	/* This photon has been absorbed before scattering */

				dtau_scatt *= frac;
				dtau = dtau_abs + dtau_scatt;
				if (dtau_abs < 1.e-3)
					ph->w *=
					    (1. -
					     dtau / 24. * (24. -
							   dtau * (12. -
								   dtau *
								   (4. -
								    dtau))));
				else
					ph->w *= exp(-dtau);

				/* Interpolate position and wave vector to scattering event */
				push_photon(Xi, Ki, dKi, dl * frac, &E0,
					    0);
				ph->X[0] = Xi[0];
				ph->X[1] = Xi[1];
				ph->X[2] = Xi[2];
				ph->X[3] = Xi[3];
				ph->K[0] = Ki[0];
				ph->K[1] = Ki[1];
				ph->K[2] = Ki[2];
				ph->K[3] = Ki[3];
				ph->dKdlam[0] = dKi[0];
				ph->dKdlam[1] = dKi[1];
				ph->dKdlam[2] = dKi[2];
				ph->dKdlam[3] = dKi[3];
				ph->E0s = E0;

				/* Get plasma parameters at new position */
				#if(HAMR)
				gcov_func_hamr(ph->X, Gcov);
				//gcov_func(ph->X, Gcov);
				#else
				gcov_func(ph->X, Gcov);
				#endif
				get_fluid_params(ph->X, Gcov, &Ne, &Thetae,
						 &B, Ucon, Ucov, Bcon,
						 Bcov);

				if (Ne > 0.) {
					scatter_super_photon(ph, &php, Ne,
							     Thetae, B,
							     Ucon, Bcon,
							     Gcov);
					if (ph->w < 1.e-100) {	/* must have been a problem popping k back onto light cone */
						return;
					}

				if(round == 0){
				#pragma omp atomic
					photon_scattered +=1;
					if(var_id == 0 && photon_index == 4821){
						var_id = 123123;
						printf("photon_index scattered = %d\n", photon_index);
					}
				}
					track_super_photon(&php, round +1, photon_index);
				}

				theta =
				    get_bk_angle(ph->X, ph->K, Ucov, Bcov,
						 B);
				nu = get_fluid_nu(ph->X, ph->K, Ucov);
				if (nu < 0.) {
					alpha_scatti = alpha_absi = 0.;
				} else {
					alpha_scatti =
					    alpha_inv_scatt(nu, Thetae,
							    Ne);
					alpha_absi =
					    alpha_inv_abs(nu, Thetae, Ne,
							  B, theta);
				}
				bi = bias_func(Thetae, ph->w);

				ph->tau_abs += dtau_abs;
				ph->tau_scatt += dtau_scatt;

			} else {
				if (dtau_abs > 100)
					return;	/* This photon has been absorbed */
				ph->tau_abs += dtau_abs;
				ph->tau_scatt += dtau_scatt;
				dtau = dtau_abs + dtau_scatt;
				if (dtau < 1.e-3)
					ph->w *=
					    (1. -
					     dtau / 24. * (24. -
							   dtau * (12. -
								   dtau *
								   (4. -
								    dtau))));
				else
					ph->w *= exp(-dtau);
			
			}
		}

		nstep++;

		/* signs that something's wrong w/ the integration */
		if (nstep > MAXNSTEP) {
			fprintf(stderr,
				"X1,X2,K1,K2,bias: %g %g %g %g %g\n",
				ph->X[1], ph->X[2], ph->K[1], ph->K[2],
				bias);
			break;
		}

	}

	/* accumulate result in spectrum on escape */
	if(1){
	//if (record_criterion(ph) && nstep < MAXNSTEP)
		record_super_photon(ph);
	}
	/* done! */
	return;
}
