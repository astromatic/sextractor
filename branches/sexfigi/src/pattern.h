 /*
 				pattern.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Authors:	E.BERTIN (IAP)
*
*	Contents:	Include file for pattern.c.
*
*	Last modify:	03/10/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
#define	PATTERN_RADIUS	1.0	/* Pattern radius over half the PATTERN_SIZE */
#define	PATTERN_SCALE	(0.2*PATTERN_RADIUS)	/* Pattern scale */

/* NOTES:
One must have:	PATTERN_SIZE > 1
		PATTERN_RADIUS > 0.0
		PATTERN_SCALE > 0.0
*/

/*--------------------------------- typedefs --------------------------------*/

typedef enum		{PATTERN_QUADRUPOLE, PATTERN_OCTOPOLE,
			PATTERN_POLARFOURIER,
			PATTERN_NPATTERNS}
				pattypenum; /* Pattern code */

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  pattypenum	type;			/* Pattern code */
  int		ncomp;			/* Number of independent components */
  int		nmodes;			/* Number of modes per component */
  int		nfreq;			/* Number of waves per component */
  double	x[2];			/* Coordinate vector */
  double	rmax;			/* Largest radius in units of scale */
  double	*r;			/* Reduced radius */
  double	*norm;			/* Pattern vector norm */
  double	*coeff;			/* Fitted pattern coefficients */
  double	*mcoeff;		/* Modulus from pattern coefficients */
  double	*acoeff;		/* Argument from pattern coefficients */
  double	*modpix;		/* Pattern pixmaps */
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

