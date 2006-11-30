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
*	Last modify:	30/11/2006
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
typedef enum	{INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}       interpenum;

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  profcode	code;			/* Model code */
  double	*comppix;		/* Composited pixmap of the model */
  double	*pix;			/* Pixmap of the model */
  int		naxis;			/* Number of pixmap dimensions */
  int		naxisn[PROFIT_MAXDIM];	/* Pixmap size for each axis */
  double	sizemax;		/* Maximum size in pixels */
/* Generic presentation parameters */
  double	*amp;			/* Amplitude */
  double	*x[2];			/* Pointer to coordinate vector */
  double	*scale[2];		/* Pointer to scaling vector */
  double	*posangle;		/* Pointer to pos. angle (CCW/NAXIS1)*/
  double	*extra[PROFIT_MAXDIM-2]	/* Parameters along extra-dimension */
  interpenum	interptype;		/* Interpolation type */
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
  psfstruct	*psf;		/* PSF */
  double	*psfdft;	/* Compressed Fourier Transform of the PSF */
  double	*fullpix;	/* Full resolution pixmap of the model */
  int		fullnaxisn[2];	/* Dimensions along each axis */
  double	objpix;		/* Copy of object pixmap */
  double	objweight;	/* Copy of object weight-map */
  int		objnaxisn[2];	/* Dimensions along each axis */
  double	*resi;		/* Vector of residuals */
  int		nresi;		/* Number of residual elements */
  }	profitstruct;

/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/
