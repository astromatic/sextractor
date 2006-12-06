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
*	Last modify:	06/12/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

#define	PROFIT_NITER	1	/* Max. nb of iterations in profile fitting */
#define	PROFIT_MAXEXTRA	2	/* Max. nb of extra free params of profiles */
#define PROFIT_RESOL	1024	/* Pixmap size of models */
/* NOTES:
One must have:	PROFIT_NITER > 0
		PROFIT_MAXEXTRA > 0
*/

/*--------------------------------- typedefs --------------------------------*/

typedef enum		{PROF_SERSIC, PROF_DEVAUCOULEURS, PROF_EXPONENTIAL,
			PROF_DIRAC}	proftypestruct; /* Profile code */
typedef enum	{INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}       interpenum;

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  proftypestruct	code;		/* Model code */
  double	*comppix;		/* Composited pixmap of the model */
  double	*pix;			/* Pixmap of the model */
  int		naxis;			/* Number of pixmap dimensions */
  int		naxisn[2+PROFIT_MAXEXTRA];	/* Pixmap size for each axis */
  double	sizemax;		/* Maximum size in pixels */
/* Generic presentation parameters */
  double	*amp;			/* Amplitude */
  double	*x[2];			/* Pointer to coordinate vector */
  double	*scale[2];		/* Pointer to scaling vector */
  double	*posangle;		/* Pointer to pos. angle (CCW/NAXIS1)*/
  double	*extra[PROFIT_MAXEXTRA]	/* Parameters along extra-dimension */
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
  double	*modpix;	/* Full resolution pixmap of the model */
  int		modnaxisn[2];	/* Dimensions along each axis */
  double	*objpix;	/* Copy of object pixmap */
  double	*objweight;	/* Copy of object weight-map */
  int		objnaxisn[2];	/* Dimensions along each axis */
  double	*resi;		/* Vector of residuals */
  int		nresi;		/* Number of residual elements */
  }	profitstruct;

/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/
