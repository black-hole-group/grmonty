
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

   grmonty Nph

   Using monte carlo method, estimate spectrum of an appropriately
   scaled axisymmetric GRMHD simulation as a function of 
   latitudinal viewing angle.

   Input simulation data is assumed to be in dump format provided by 
   HARM code.  Location of input file is, at present, hard coded
   (see init_sim_data.c).  

   Nph super-photons are generated in total and then allowed
   to propagate.  They are weighted according to the emissivity.
   The photons are pushed by the geodesic equation.
   Their weight decays according to the local absorption coefficient.
   The photons also scatter with probability related to the local
   scattering opacity.  

   The electrons are assumed to have a thermal distribution 
   function, and to be at the same temperature as the protons.

   CFG 8-17-06

   Implemented synchrotron sampling, 22 Jan 07

   fixed bugs in tetrad/compton scattering routines, 31 Jan 07

   Implemented normalization for output, 6 Feb 07

   Separated out different synchrotron sampling routines
   into separate files, 8 Mar 07

   fixed bug in energy recording; bug used Kcon[0] rather than 
   Kcov[0] as energy, 18 Mar 07

   major reorganization to encapsulate problem-dependent parts 5-6 Nov 07

 */

#include "decs.h"

/* defining declarations for global variables */
struct of_geom *geom;
int N1, N2, N3, n_within_horizon;
double F[N_ESAMP + 1], wgt[N_ESAMP + 1];
int Ns, N_superph_recorded, N_scatt;

/* some coordinate parameters */
double a;
double R0, Rin, Rh, Rout, Rms;
double hslope;
double startx[NDIM], stopx[NDIM], dx[NDIM];

double dlE, lE0;
double gam;
double dMsim;
double M_unit, L_unit, T_unit;
double RHO_unit, U_unit, B_unit, Ne_unit, Thetae_unit;
double max_tau_scatt, Ladv, dMact, bias_norm;

gsl_rng *r;
gsl_integration_workspace *w;

#pragma omp threadprivate(r)
#include <time.h>

int main(int argc, char *argv[])
{
	double Ntot, N_superph_made;
	int quit_flag, myid;
	struct of_photon ph;
	time_t currtime, starttime;
	#pragma omp parallel
	{
		#pragma omp single
		printf("num_threads in use = %d\n", omp_get_num_threads());
	}
	if (argc < 3) {
		fprintf(stderr, "usage: grmonty Ns infilename M_unit\n");
		exit(0);
	}
	sscanf(argv[1], "%lf", &Ntot);
	Ns = (int) Ntot;

	/* initialize random number generator */
#pragma omp parallel private(myid)
	{
		myid = omp_get_thread_num();
		init_monty_rand(139 * myid + time(NULL));	/* Arbitrarily picked initial seed */
	}

	/* spectral bin parameters */
	dlE = 0.25;		/* bin width */
	lE0 = log(1.e-12);	/* location of first bin, in electron rest-mass units */

	/* initialize model data, auxiliary variables */
	init_model(argv);

	/** main loop **/
	N_superph_made = 0;
	N_superph_recorded = 0;
	N_scatt = 0;
	starttime = time(NULL);
	quit_flag = 0;

	fprintf(stderr, "Entering main loop...\n");
	fflush(stderr);

#pragma omp parallel private(ph)
	{

		while (1) {

			/* get pseudo-quanta */
#pragma omp critical (MAKE_SPHOT)
			{
				if (!quit_flag)
					make_super_photon(&ph, &quit_flag);
			}
			if (quit_flag)
				break;

			/* push them around */
			track_super_photon(&ph);

			/* step */
#pragma omp atomic
			N_superph_made += 1;

			/* give interim reports on rates */
			if (((int) (N_superph_made)) % 100000 == 0
			    && N_superph_made > 0) {
				currtime = time(NULL);
				fprintf(stderr, "time %g, rate %g ph/s\n",
					(double) (currtime - starttime),
					N_superph_made / (currtime -
							  starttime));
			}
		}
	}
	currtime = time(NULL);
	fprintf(stderr, "Final time %g, rate %g ph/s\n",
		(double) (currtime - starttime),
		N_superph_made / (currtime - starttime));

#ifdef _OPENMP
#pragma omp parallel
	{
		omp_reduce_spect();
	}
#endif
	report_spectrum((int) N_superph_made);

	/* done! */
	return (0);

}
