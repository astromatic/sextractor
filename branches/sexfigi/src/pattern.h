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
*	Last modify:	09/09/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PATTERN_H_
#define _PATTERN_H_

/*-------------------------------- flags ------------------------------------*/
/*-------------------------------- macros -----------------------------------*/
/*----------------------------- Internal constants --------------------------*/

#define	PATTERN_SIZE	255     /* Pattern image size */
#define	PATTERN_RADIUS	1.0	/* Pattern radius over half the PATTERN_SIZE */
#define	PATTERN_SCALE	(0.2*PATTERN_RADIUS)	/* Pattern scale */

/* NOTES:
One must have:	PATTERN_SIZE > 1
		PATTERN_RADIUS > 0.0
		PATTERN_SCALE > 0.0
*/

/*--------------------------------- typedefs --------------------------------*/

typedef enum		{PATTERN_QUADRIPOLE, PATTERN_NPATTERNS}
				pattypenum; /* Pattern code */

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  pattypenum	type;			/* Pattern code */
  float		*pix;			/* Full pixmap of the model */
  int		size[3];		/* Pixmap size for each axis */
  }	patternstruct;


/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/

patternstruct	*pattern_init(pattypenum ptype, int nvec);

void		pattern_end(patternstruct *pattern);

#endif

