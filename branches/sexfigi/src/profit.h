 /*
 				profit.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Authors:	E.BERTIN (IAP)
*
*	Contents:	Include file for profit.c.
*
*	Last modify:	14/11/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

#define	PROFIT_NITER	1	/* Max. nb of iterations in profile fitting */
#define	PROFIT_MAXDIM	4	/* Max. nb of dimensions of profiles */

/* NOTES:
One must have:	PROFIT_NITER > 0
		PROFIT_MAXDIM > 2
*/

/*--------------------------------- typedefs --------------------------------*/

typedef enum		{SERSIC, DEVAUCOULEURS, EXPONENTIAL, DIRAC}
		codestruct; /* Profile code */

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  profcode	code;			/* Model code */
  PIXTYPE	*pix;			/* Pixmap of the model */
  int		naxis;			/* Number of pixmap dimensions */
  int		naxisn[PROFIT_MAXDIM];	/* Pixmap size for each axis */
/* Generic presentation parameters */
  double	*amp;			/* Amplitude */
  double	*x[PROFIT_NAXIS];	/* Pointer to coordinate vector */
  double	*rho[PROFIT_NAXIS];	/* Pointer to scaling vector */
  double	*theta;			/* Pointer to pos. angle (CCW/NAXIS1)*/
  double	*extra[PROFIT_MAXDIM-2]	/* Parameters along extra-dimension
  }	profstruct;

typedef struct
  {
  double	*param;		/* Vector of parameters to be fitted */
  int		nparam;		/* Number of parameters to be fitted */
  double	**parammin;	/* Parameter lower limit */
  double	**parammax;	/* Parameter upper limit */
  int		niter;		/* Number of iterations */
  profstruct	**prof;		/* Array of pointers to profiles */
  int		nprof;		/* Number of profiles to consider */
  }	profitstruct;

/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/
