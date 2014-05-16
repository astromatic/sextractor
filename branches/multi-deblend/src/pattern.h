/*
*				pattern.h
*
* Include file for pattern.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2007-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		11/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _PROFIT_H_
#include "profit.h"
#endif

#ifndef _PATTERN_H_
#define _PATTERN_H_

/*-------------------------------- flags ------------------------------------*/
/*-------------------------------- macros -----------------------------------*/
/*----------------------------- Internal constants --------------------------*/

#define	PATTERN_FMAX	4	/* Maximum pattern angular frequency */
#define	PATTERN_NCOMP	16	/* Default number of components (radii) */
#define	PATTERN_SCALE	5.0	/* Pattern scale in units of r_eff */
#define	PATTERN_MARGIN	0.2	/* Pattern margin in fractions of radius */
#define	PATTERN_BTMAX	0.6	/* Maximum B/T for pure disk scaling */

/* NOTES:
One must have:	PATTERN_SIZE > 1
		PATTERN_SCALE > 0.0
		PATTERN_BTMAX < 1.0
*/

/*--------------------------------- typedefs --------------------------------*/

typedef enum		{PATTERN_QUADRUPOLE, PATTERN_OCTOPOLE,
			PATTERN_POLARFOURIER, PATTERN_POLARSHAPELETS,
			PATTERN_NPATTERNS}
				pattypenum; /* Pattern code */

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  pattypenum	type;			/* Pattern code */
  int		ncomp;			/* Number of independent components */
  int		nmodes;			/* Number of modes per component */
  int		nfreq;			/* Number of waves per component */
  float		x[2];			/* Coordinate vector */
  float		rmax;			/* Largest radius in units of scale */
  float		*r;			/* Reduced radius */
  float		*norm;			/* Pattern vector norm */
  float		*coeff;			/* Fitted pattern coefficients */
  float		*mcoeff;		/* Modulus from pattern coefficients */
  float		*acoeff;		/* Argument from pattern coefficients */
  float		*modpix;		/* Pattern pixmaps */
  PIXTYPE	*lmodpix;		/* Low resolution pattern pixmaps */
  int		size[3];		/* Pixmap size for each axis */
  }	patternstruct;


/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/

patternstruct	*pattern_init(profitstruct *profit, pattypenum ptype, int nvec);

float		pattern_spiral(patternstruct *pattern);
void		pattern_compmodarg(patternstruct *pattern,profitstruct *profit),
		pattern_create(patternstruct *pattern, profitstruct *profit),
		pattern_end(patternstruct *pattern),
		pattern_fit(patternstruct *pattern, profitstruct *profit);

#endif

