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
*	Last modify:	10/12/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

#define	PROFIT_MAXITER	200	/* Max. nb of iterations in profile fitting */
#define	PROFIT_MAXPROF	8	/* Max. nb of profile components */
#define	PROFIT_MAXEXTRA	2	/* Max. nb of extra free params of profiles */
#define PROFIT_PROFRES	1024	/* Pixmap size of model components */
#define PROFIT_PROFSRES	64	/* Number of model subcomponents */
#define INTERP_MAXKERNELWIDTH	8	/* Max. range of kernel (pixels) */
/* NOTES:
One must have:	PROFIT_NITER > 0
		PROFIT_MAXEXTRA > 0
*/

/*--------------------------------- typedefs --------------------------------*/

typedef enum		{PROF_SERSIC, PROF_DEVAUCOULEURS, PROF_EXPONENTIAL,
			PROF_DIRAC}	proftypenum; /* Profile code */
typedef enum	{INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}       interpenum;

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  proftypenum	code;			/* Model code */
  double	*pix;			/* Full pixmap of the model */
  int		naxis;			/* Number of pixmap dimensions */
  int		naxisn[2+PROFIT_MAXEXTRA];	/* Pixmap size for each axis */
  double	typscale;		/* Typical scale in prof pixels */
/* Generic presentation parameters */
  double	*amp;			/* Amplitude */
  double	*x[2];			/* Pointer to coordinate vector */
  double	*scale[2];		/* Pointer to scaling vector */
  double	*posangle;		/* Pointer to pos. angle (CCW/NAXIS1)*/
  double	*extra[PROFIT_MAXEXTRA];/* Parameters along extra-dimension */
  double	extrazero[PROFIT_MAXEXTRA]; /* Zero-point along extra-dim. */
  double	extrascale[PROFIT_MAXEXTRA]; /* Scaling along extra-dim. */
  int		extracycleflag[PROFIT_MAXEXTRA]; /* !=0 for cycling dim. */
  interpenum	interptype[2+PROFIT_MAXEXTRA];	/* Interpolation type */
  int		kernelwidth[2+PROFIT_MAXEXTRA];	/* Kernel size */
  double	*kernelbuf;		/* Kernel buffer */
  int		kernelnlines;		/* Number of interp kernel lines */
  }	profstruct;

typedef struct
  {
  double	*initparam;	/* Vector of parameter guesses */
  double	*param;		/* Vector of parameters to be fitted */
  int		nparam;		/* Number of parameters to be fitted */
  double	*paraminit;	/* Vector of parameter guesses */
  double	*parammin;	/* Parameter lower limit */
  double	*parammax;	/* Parameter upper limit */
  int		niter;		/* Number of iterations */
  profstruct	**prof;		/* Array of pointers to profiles */
  int		nprof;		/* Number of profiles to consider */
  psfstruct	*psf;		/* PSF */
  double	*psfdft;	/* Compressed Fourier Transform of the PSF */
  double	*modpix;	/* Full resolution pixmap of the model */
  int		modnaxisn[2];	/* Dimensions along each axis */
  PIXTYPE	*lmodpix;	/* Low resolution pixmap of the model */
  PIXTYPE	*objpix;	/* Copy of object pixmap */
  PIXTYPE	*objweight;	/* Copy of object weight-map */
  int		objnaxisn[2];	/* Dimensions along each axis */
  double	*resi;		/* Vector of residuals */
  int		nresi;		/* Number of residual elements */
  double	sigma;		/* Standard deviation of the pixel values */
  }	profitstruct;

/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/

profitstruct	*profit_init(psfstruct *psf, proftypenum *profcode,
			int nprof);

profstruct	*prof_init(profitstruct *profit, proftypenum profcode);

double		*profit_compresi(profitstruct *profit, picstruct *field,
			picstruct *wfield, objstruct *obj, double *resi),
		*profit_residuals(profitstruct *profit, picstruct *field,
			picstruct *wfield, objstruct *obj, double *param,
			double *resi);

PIXTYPE		*profit_resample(profitstruct *profit);

int		profit_copyobjpix(profitstruct *profit, picstruct *field,
				int ix, int iy);

void		prof_add(profstruct *prof, profitstruct *profit),
		prof_end(profstruct *prof),
		prof_fit(profitstruct *profit,
			picstruct *field, picstruct *wfield,
			objstruct *obj, obj2struct *obj2),
		profit_convolve(profitstruct *profit),
		profit_end(profitstruct *profit),
		profit_evaluate(double *par, int m_dat, double *fvec,
			void *data, int *info),
		profit_evaluate2(double *par, double *fvec, int m, int n,
			void *adata),
		profit_makedft(profitstruct *profit),
		profit_printout(int n_par, double* par, int m_dat, double* fvec,
			void *data, int iflag, int iter, int nfev );
