# Makefile.am for the levmar Levenberg-Marquardt fitting library
# Copyright (C) 2007-2010 Emmanuel Bertin.
# levmar code by Manolis Lourakis <lourakis@ics.forth.gr> 
# http://www.ics.forth.gr/~lourakis/levmar/
noinst_LIBRARIES	= liblevmar.a
liblevmar_a_SOURCES	= Axb.c lmbc.c lm.c lmblec.c lmbleic.c lmlec.c misc.c \
			  compiler.h levmar.h misc.h
EXTRA_liblevmar_a_SOURCES = Axb_core.c lmbc_core.c lm_core.c \
			  lmblec_core.c lmbleic_core.c lmlec_core.c \
			  misc_core.c \
			  LICENSE README README.txt \
			  Makefile.icc Makefile.vc levmar.vcproj lmdemo.vcproj \
			  expfit.c lmdemo.c lm.h matlab
