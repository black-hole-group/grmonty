

#/***********************************************************************************
#    Copyright 2013 Joshua C. Dolence, Charles F. Gammie, Monika Mo\'scibrodzka,
#                   and Po Kin Leung
#
#                        GRMONTY  version 1.0   (released February 1, 2013)
#
#    This file is part of GRMONTY.  GRMONTY v1.0 is a program that calculates the
#    emergent spectrum from a model using a Monte Carlo technique.
#
#    This version of GRMONTY is configured to use input files from the HARM code
#    available on the same site.   It assumes that the source is a plasma near a
#    black hole described by Kerr-Schild coordinates that radiates via thermal 
#    synchrotron and inverse compton scattering.
#    
#    You are morally obligated to cite the following paper in any
#    scientific literature that results from use of any part of GRMONTY:
#
#    Dolence, J.C., Gammie, C.F., Mo\'scibrodzka, M., \& Leung, P.-K. 2009,
#        Astrophysical Journal Supplement, 184, 387
#
#    Further, we strongly encourage you to obtain the latest version of 
#    GRMONTY directly from our distribution website:
#    http://rainman.astro.illinois.edu/codelib/
#
#    GRMONTY is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    GRMONTY is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GRMONTY; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
#***********************************************************************************/
#
# requires an openmp-enabled version of gcc
#

CC = gcc
CFLAGS = -Wall -O2 -fopenmp -I$(INCDIR)
LDFLAGS = -lm -lgsl -lgslcblas -fopenmp

SRCDIR = src
INCDIR = include
BUILDDIR = build

SRCS = $(SRCDIR)/main.c \
       $(SRCDIR)/physics/compton.c \
       $(SRCDIR)/physics/radiation.c \
       $(SRCDIR)/physics/jnu_mixed.c \
       $(SRCDIR)/physics/hotcross.c \
       $(SRCDIR)/physics/scatter_super_photon.c \
       $(SRCDIR)/geometry/init_geometry.c \
       $(SRCDIR)/geometry/tetrads.c \
       $(SRCDIR)/geometry/geodesics.c \
       $(SRCDIR)/model/harm_model.c \
       $(SRCDIR)/model/harm_utils.c \
       $(SRCDIR)/model/init_harm_data.c \
       $(SRCDIR)/model/track_super_photon.c

OBJS = $(BUILDDIR)/main.o \
       $(BUILDDIR)/compton.o \
       $(BUILDDIR)/radiation.o \
       $(BUILDDIR)/jnu_mixed.o \
       $(BUILDDIR)/hotcross.o \
       $(BUILDDIR)/scatter_super_photon.o \
       $(BUILDDIR)/init_geometry.o \
       $(BUILDDIR)/tetrads.o \
       $(BUILDDIR)/geodesics.o \
       $(BUILDDIR)/harm_model.o \
       $(BUILDDIR)/harm_utils.o \
       $(BUILDDIR)/init_harm_data.o \
       $(BUILDDIR)/track_super_photon.o

INCS = $(INCDIR)/decs.h $(INCDIR)/constants.h $(INCDIR)/harm_model.h

grmonty : $(BUILDDIR) $(OBJS) $(INCS) makefile 
	$(CC) $(CFLAGS) -o grmonty $(OBJS) $(LDFLAGS)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/main.o: $(SRCDIR)/main.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/compton.o: $(SRCDIR)/physics/compton.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/radiation.o: $(SRCDIR)/physics/radiation.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/jnu_mixed.o: $(SRCDIR)/physics/jnu_mixed.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/hotcross.o: $(SRCDIR)/physics/hotcross.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/scatter_super_photon.o: $(SRCDIR)/physics/scatter_super_photon.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/init_geometry.o: $(SRCDIR)/geometry/init_geometry.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/tetrads.o: $(SRCDIR)/geometry/tetrads.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/geodesics.o: $(SRCDIR)/geometry/geodesics.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/harm_model.o: $(SRCDIR)/model/harm_model.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/harm_utils.o: $(SRCDIR)/model/harm_utils.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/init_harm_data.o: $(SRCDIR)/model/init_harm_data.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILDDIR)/track_super_photon.o: $(SRCDIR)/model/track_super_photon.c $(INCS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	/bin/rm -rf $(BUILDDIR) *.o grmonty.spec

.PHONY: clean
