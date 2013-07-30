/*
*				check.h
*
* Include file for check.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SExtractor is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SExtractor is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		26/02/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

/*----------------------------- Internal constants --------------------------*/

#define CHECKINTERPW		6	/* Interpolation function range */
#define	CHECKINTERPFAC		3.0	/* Interpolation envelope factor */

#define	CHECKINTERPF(x)	(x<1e-5 && x>-1e-5? 1.0 \
			:(x>CHECKINTERPFAC?0.0:(x<-CHECKINTERPFAC?0.0 \
			:sinf(PI*x)*sinf(PI/CHECKINTERPFAC*x) \
				/(PI*PI/CHECKINTERPFAC*x*x))))
				/* Lanczos approximation */

/*--------------------------------- structures ------------------------------*/
/* Check-image parameters */
typedef struct structcheck
  {
  int		next;			/* Number of extensions */
  catstruct	*cat;			/* FITS file structure */
  void		*pix;			/* ptr to check-image pixmap */
  int		width, height, depth;	/* size of check-image */
  size_t	npix;			/* number of pixels in check-image */
  int		y;			/* current line in check-image */
  PIXTYPE	overlay;		/* intensity of the overlayed plots */
  void		*line;			/* buffered image line */
  checkenum	type;			/* CHECKIMAGE_TYPE */
  }	checkstruct;

/*------------------------------- functions ---------------------------------*/

checkstruct	*initcheck(char *, checkenum, int next);

void		addcheck(checkstruct *, float *, int,int, int,int, float),
		addcheck_resample(checkstruct *, float *, int,int, int,int,
			float, float),
		blankcheck(checkstruct *, PIXTYPE *, int,int,int,int,PIXTYPE),
		endcheck(checkstruct *),
		reendcheck(picstruct *field, checkstruct *),
		reinitcheck(picstruct *, checkstruct *),
		writecheck(checkstruct *, PIXTYPE *, int);
