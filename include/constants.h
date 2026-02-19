
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



/* == Definition of constants in CGS units == */

#define EE				(4.80320680e-10		)	/* electron charge */
#define CL				(2.99792458e10		)	/* speed of light */
#define ME				(9.1093826e-28		)	/* electron mass */
#define MP				(1.67262171e-24		)	/* proton mass */
#define MN				(1.67492728e-24		)	/* neutron mass */
#define AMU				(1.66053886e-24		)	/* atomic mass unit */
#define HPL				(6.6260693e-27		)	/* Planck constant */
#define HBAR			(HPL/(2.*M_PI)		)	/* Planck's consant / 2pi */
#define KBOL			(1.3806505e-16		)	/* Boltzmann constant */
#define GNEWT			(6.6742e-8			)	/* Gravitational constant */
#define SIG				(5.670400e-5		)	/* Stefan-Boltzmann constant */
#define RGAS			(8.3143e7			)	/* erg K^-1 mole^-1: ideal gas const */
#define EV				(1.60217653e-12		)	/* electron volt in erg */
#define SIGMA_THOMSON	(0.665245873e-24	)	/* Thomson cross section in cm^2 */
#define JY				(1.e-23				)	/* Jansky (flux/freq. unit) in cgs */

#define PC				(3.085678e18		)	/* parsec */
#define AU				(1.49597870691e13	)	/* Astronomical Unit */

#define YEAR			(31536000.			)	/* No. of seconds in year */
#define DAY				(86400.				)	/* No. of seconds in day  */
#define HOUR			(3600.				)	/* No. of seconds in hour */

#define MSUN			(1.989e33			)	/* solar mass */
#define RSUN			(6.96e10			)	/* Radius of Sun */
#define LSUN			(3.827e33				)	/* Luminousity of Sun */
#define TSUN			(5.78e3				)	/* Temperature of Sun's photosphere */

#define MEARTH			(5.976e27			)	/* Earth's mass */
#define REARTH			(6.378e8			)	/* Earth's radius */

#define DSGRA			(8.4e3 * PC			)	/* Distance from Earth to Sgr A*  */

#define TCBR			(2.726				)	/* CBR temperature, from COBE */

/* 
   abundances, from M & B, p. 99
*/

#define SOLX			(0.70				)	/* H */
#define SOLY			(0.28				)	/* He */
#define SOLZ			(0.02				)	/* Metals */

#define FOUR_PI	(1.2732395447351626862)	/* 4/pi */
#define FOUR_PISQ (0.40528473456935108578)	/* 4/pi^2 */
#define THREEPI_TWO (4.7123889803846898577)	/* 3pi/2 */

