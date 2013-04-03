#
#				Makefile.am
#
# src Makefile.am. Process this file with automake to generate a Makefile
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SExtractor
#
#	Copyright:		(C) 2002-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
#
#	License:		GNU General Public License
#
#	SExtractor is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SExtractor is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		03/04/2013
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if USE_MODEL
  FFTSOURCE		= fft.c
  PATTERNSOURCE		= pattern.c
  PROFITSOURCE		= profit.c
  LEVLIB		= $(srcdir)/levmar/liblevmar.a
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
sex_LDADD		= $(srcdir)/fits/libfits.a \
			  $(srcdir)/wcs/libwcs_c.a \
			  $(LEVLIB) 
ldactoasc_LDADD		= $(srcdir)/fits/libfits.a
DATE=`date +"%Y-%m-%d"`

