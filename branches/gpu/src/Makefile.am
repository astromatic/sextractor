# Program Makefile for SEx
# Copyright (C) 2004-2010 Emmanuel Bertin.
if USE_MODEL
  FFTSOURCE		= fft.c
  PATTERNSOURCE		= pattern.c
  PROFITSOURCE		= profit.c
  LEVLIB		= $(top_builddir)/src/levmar/liblevmar.a
  LEVDIR		= levmar
endif
SUBDIRS			= fits $(LEVDIR) wcs
bin_PROGRAMS		= sex ldactoasc
check_PROGRAMS		= sex
sex_SOURCES		= analyse.c assoc.c astrom.c back.c bpro.c catout.c \
			  check.c clean.c extract.c $(FFTSOURCE) field.c \
			  filter.c fitswcs.c flag.c graph.c growth.c header.c \
			  image.c interpolate.c main.c makeit.c manobjlist.c \
			  misc.c neurro.c $(PATTERNSOURCE) pc.c photom.c \
			  plist.c prefs.c $(PROFITSOURCE) psf.c readimage.c \
			  refine.c retina.c scan.c som.c weight.c winpos.c \
			  xml.c \
			  assoc.h astrom.h back.h bpro.h check.h clean.h \
			  define.h extract.h fft.h field.h filter.h fitswcs.h \
			  flag.h globals.h growth.h header.h image.h \
			  interpolate.h key.h neurro.h param.h paramprofit.h \
			  pattern.h photom.h plist.h prefs.h preflist.h \
			  profit.h psf.h retina.h sexhead1.h sexhead.h \
			  sexheadsc.h som.h threads.h types.h wcscelsys.h \
			  weight.h winpos.h xml.h
ldactoasc_SOURCES 	= ldactoasc.c ldactoasc.h
sex_LDADD		= $(top_builddir)/src/fits/libfits.a \
			  $(top_builddir)/src/wcs/libwcs_c.a \
			  $(LEVLIB) 
ldactoasc_LDADD		= $(top_builddir)/src/fits/libfits.a
DATE=`date +"%Y-%m-%d"`

