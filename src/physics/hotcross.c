
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


#include <decs.h>

/* 

   given energy of photon in fluid rest frame w, in units of electron rest mass
   energy, and temperature of plasma, again in electron rest-mass units, return hot
   cross section in cgs.
   
   This has been checked against Wienke's Table 1, with some disagreement at
   the one part in 10^{-3} level, see wienke_table_1 in the subdirectory hotcross.
   It is not clear what this is due to, but Table 1 does appear to have been evaluated
   using Monte Carlo integration (!).

   A better way to do this would be to make a table in w*thetae and w/thetae; most
   	of the variation is accounted for by w*thetae.
   
*/

#define MINW	1.e-12
#define MAXW	1.e6
#define MINT	0.0001
#define MAXT	1.e4
#define NW	220
#define NT	80

#define HOTCROSS	"data/hotcross.dat"

double table[NW + 1][NT + 1];
double dlw, dlT, lminw, lmint;

void init_hotcross(void)
{
	int i, j, nread;
	double lw, lT;
	double total_compton_cross_num(double w, double thetae);
	FILE *fp;

	dlw = log10(MAXW / MINW) / NW;
	dlT = log10(MAXT / MINT) / NT;
	lminw = log10(MINW);
	lmint = log10(MINT);

	fp = fopen(HOTCROSS, "r");
	if (fp == NULL) {
		fprintf(stderr, "file %s not found.\n", HOTCROSS);
		fprintf(stderr,
			"making lookup table for compton cross section...\n");
#pragma omp parallel for private(i,j,lw,lT)
		for (i = 0; i <= NW; i++)
			for (j = 0; j <= NT; j++) {
				lw = lminw + i * dlw;
				lT = lmint + j * dlT;
				table[i][j] =
				    log10(total_compton_cross_num
					  (pow(10., lw), pow(10., lT)));
				if (isnan(table[i][j])) {
					fprintf(stderr, "%d %d %g %g\n", i,
						j, lw, lT);
					exit(0);
				}
			}
		fprintf(stderr, "done.\n\n");
		fprintf(stderr, "writing to file...\n");
		fp = fopen(HOTCROSS, "w");
		if (fp == NULL) {
			fprintf(stderr, "couldn't write to file\n");
			exit(0);
		}
		for (i = 0; i <= NW; i++)
			for (j = 0; j <= NT; j++) {
				lw = lminw + i * dlw;
				lT = lmint + j * dlT;
				fprintf(fp, "%d %d %g %g %15.10g\n", i, j,
					lw, lT, table[i][j]);
			}
		fprintf(stderr, "done.\n\n");
	} else {
		fprintf(stderr,
			"reading hot cross section data from %s...\n",
			HOTCROSS);
		for (i = 0; i <= NW; i++)
			for (j = 0; j <= NT; j++) {
				nread =
				    fscanf(fp, "%*d %*d %*lf %*lf %lf\n",
					   &table[i][j]);
				if (isnan(table[i][j]) || nread != 1) {
					fprintf(stderr,
						"error on table read: %d %d\n",
						i, j);
					exit(0);
				}
			}
		fprintf(stderr, "done.\n\n");
	}

	fclose(fp);

	return;
}



double total_compton_cross_lkup(double w, double thetae)
{
	int i, j;
	double lw, lT, di, dj, lcross;
	double total_compton_cross_num(double w, double thetae);
	double hc_klein_nishina(double we);

	/* cold/low-energy: just use thomson cross section */
	if (w * thetae < 1.e-6)
		return (SIGMA_THOMSON);

	/* cold, but possible high energy photon: use klein-nishina */
	if (thetae < MINT)
		return (hc_klein_nishina(w) * SIGMA_THOMSON);

	/* in-bounds for table */
	if ((w > MINW && w < MAXW) && (thetae > MINT && thetae < MAXT)) {

		lw = log10(w);
		lT = log10(thetae);
		i = (int) ((lw - lminw) / dlw);
		j = (int) ((lT - lmint) / dlT);
		di = (lw - lminw) / dlw - i;
		dj = (lT - lmint) / dlT - j;

		lcross =
		    (1. - di) * (1. - dj) * table[i][j] + di * (1. -
								dj) *
		    table[i + 1][j] + (1. - di) * dj * table[i][j + 1] +
		    di * dj * table[i + 1][j + 1];

		if (isnan(lcross)) {
			fprintf(stderr, "%g %g %d %d %g %g\n", lw, lT, i,
				j, di, dj);
		}

		return (pow(10., lcross));
	}

	fprintf(stderr, "out of bounds: %g %g\n", w, thetae);
	return (total_compton_cross_num(w, thetae));

}

#define MAXGAMMA	12.
#define DMUE		0.05
#define DGAMMAE		0.05

double total_compton_cross_num(double w, double thetae)
{
	double dmue, dgammae, mue, gammae, f, cross;
	double dNdgammae(double thetae, double gammae);
	double boostcross(double w, double mue, double gammae);
	double hc_klein_nishina(double we);

	if (isnan(w)) {
		fprintf(stderr, "compton cross isnan: %g %g\n", w, thetae);
		return (0.);
	}

	/* check for easy-to-do limits */
	if (thetae < MINT && w < MINW)
		return (SIGMA_THOMSON);
	if (thetae < MINT)
		return (hc_klein_nishina(w) * SIGMA_THOMSON);

	dmue = DMUE;
	dgammae = thetae * DGAMMAE;

	/* integrate over mu_e, gamma_e, where mu_e is the cosine of the
	   angle between k and u_e, and the angle k is assumed to lie,
	   wlog, along the z axis */
	cross = 0.;
	for (mue = -1. + 0.5 * dmue; mue < 1.; mue += dmue)
		for (gammae = 1. + 0.5 * dgammae;
		     gammae < 1. + MAXGAMMA * thetae; gammae += dgammae) {

			f = 0.5 * dNdgammae(thetae, gammae);

			cross +=
			    dmue * dgammae * boostcross(w, mue,
							gammae) * f;

			if (isnan(cross)) {
				fprintf(stderr, "%g %g %g %g %g %g\n", w,
					thetae, mue, gammae,
					dNdgammae(thetae, gammae),
					boostcross(w, mue, gammae));
			}
		}


	return (cross * SIGMA_THOMSON);
}

/* normalized (per unit proper electron number density)
   electron distribution */
double dNdgammae(double thetae, double gammae)
{
	double K2f;

	if (thetae > 1.e-2) {
		K2f = gsl_sf_bessel_Kn(2, 1. / thetae) * exp(1. / thetae);
	} else {
		K2f = sqrt(M_PI * thetae / 2.);
	}

	return ((gammae * sqrt(gammae * gammae - 1.) / (thetae * K2f)) *
		exp(-(gammae - 1.) / thetae));
}

double boostcross(double w, double mue, double gammae)
{
	double we, boostcross, v;
	double hc_klein_nishina(double we);

	/* energy in electron rest frame */
	v = sqrt(gammae * gammae - 1.) / gammae;
	we = w * gammae * (1. - mue * v);

	boostcross = hc_klein_nishina(we) * (1. - mue * v);

	if (boostcross > 2) {
		fprintf(stderr, "w,mue,gammae: %g %g %g\n", w, mue,
			gammae);
		fprintf(stderr, "v,we, boostcross: %g %g %g\n", v, we,
			boostcross);
		fprintf(stderr, "kn: %g %g %g\n", v, we, boostcross);
	}

	if (isnan(boostcross)) {
		fprintf(stderr, "isnan: %g %g %g\n", w, mue, gammae);
		exit(0);
	}

	return (boostcross);
}

double hc_klein_nishina(double we)
{
	double sigma;

	if (we < 1.e-3)
		return (1. - 2. * we);

	sigma = (3. / 4.) * (2. / (we * we) +
			     (1. / (2. * we) -
			      (1. + we) / (we * we * we)) * log(1. +
								2. * we) +
			     (1. + we) / ((1. + 2. * we) * (1. + 2. * we))
	    );

	return (sigma);

}
