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
*	Last modify:	09/10/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PROFIT_H_
#define _PROFIT_H_

/*-------------------------------- flags ------------------------------------*/

#define		PROFLAG_MODSUB		0x0001
#define		PROFLAG_OBJSUB		0x0002
#define		PROFLAG_NOTCONST	0x0004

/*-------------------------------- macros -----------------------------------*/

#define		PROFIT_POW(x,a)		(x>0.01? exp(a*log(x)) : pow(x,a))
#define		PROFIT_POWF(x,a)	(x>0.01? expf(a*logf(x)) : powf(x,a))

/*----------------------------- Internal constants --------------------------*/

#define	PROFIT_MAXITER	1000	/* Max. nb of iterations in profile fitting */
#define	PROFIT_MAXPROF	8	/* Max. nb of profile components */
#define	PROFIT_OVERSAMP	5	/* Max. profile oversamp. factor on each axis */
#define	PROFIT_HIDEFRES	201	/* Resolution of the high def. model raster */
#define	PROFIT_REFFFAC	6.0	/* Factor in r_eff for measurement radius*/
#define	PROFIT_DYNPARAM	10.0	/* Dynamic compression param. in sigma units */
#define	PROFIT_MAXMODSIZE  512	/* Maximum size allowed for the model raster */
#define	PROFIT_MAXOBJSIZE  512	/* Maximum size allowed for the object raster */
#define	PROFIT_BARXFADE	0.1	/* Fract. of bar length crossfaded with arms */
#define	PROFIT_MAXEXTRA	2	/* Max. nb of extra free params of profiles */
#define PROFIT_PROFRES	256	/* Pixmap size of model components */
#define PROFIT_PROFSRES	64	/* Number of model subcomponents */
#define INTERP_MAXKERNELWIDTH	8	/* Max. range of kernel (pixels) */
/* NOTES:
One must have:	PROFIT_NITER > 0
		PROFIT_MAXEXTRA > 0
*/

/*--------------------------------- typedefs --------------------------------*/

typedef enum		{PROF_BACK, PROF_SERSIC, PROF_DEVAUCOULEURS,
			PROF_EXPONENTIAL, PROF_ARMS, PROF_BAR, PROF_INRING,
			PROF_OUTRING, PROF_SERSIC_TABEX, PROF_DIRAC, PROF_NPROF}
				proftypenum; /* Profile code */

typedef enum	{INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}       interpenum;

typedef enum	{PARAM_BACK, PARAM_X, PARAM_Y,
		PARAM_SPHEROID_FLUX, PARAM_SPHEROID_REFF, PARAM_SPHEROID_ASPECT,
		PARAM_SPHEROID_POSANG, PARAM_SPHEROID_SERSICN,
		PARAM_DISK_FLUX, PARAM_DISK_SCALE, PARAM_DISK_ASPECT,
		PARAM_DISK_POSANG,
		PARAM_ARMS_FLUX, PARAM_ARMS_QUADFRAC, PARAM_ARMS_SCALE,
		PARAM_ARMS_START, PARAM_ARMS_POSANG, PARAM_ARMS_PITCH,
		PARAM_ARMS_PITCHVAR, PARAM_ARMS_WIDTH,
		PARAM_BAR_FLUX, PARAM_BAR_ASPECT, PARAM_BAR_POSANG,
		PARAM_INRING_FLUX, PARAM_INRING_WIDTH, PARAM_INRING_ASPECT,
		PARAM_OUTRING_FLUX, PARAM_OUTRING_START, PARAM_OUTRING_WIDTH,
		PARAM_NPARAM}	paramenum;

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  proftypenum	code;			/* Model code */
  float		*pix;			/* Full pixmap of the model */
  int		naxis;			/* Number of pixmap dimensions */
  int		naxisn[3];		/* Pixmap size for each axis */
  float		typscale;		/* Typical scale in prof pixels */
  float		fluxfac;		/* Flux normalisation factor */
  float		lostfluxfrac;		/* Lost flux fraction */
/* Generic presentation parameters */
  float		*flux;			/* Integrated flux */
  float		*x[2];			/* Coordinate vector */
  float		*scale;			/* Scaling vector */
  float		*aspect;		/* Aspect ratio */
  float		*posangle;		/* Position angle (CCW/NAXIS1)*/
  float		*featfrac;		/* Feature flux fraction */
  float		*featscale;		/* Feature relative scalelength */
  float		*featstart;		/* Feature relative starting radius */
  float		*featposang;		/* Feature position angle */
  float		*featpitch;		/* Feature pitch */
  float		*featpitchvar;		/* Feature pitch variation */
  float		*featwidth;		/* Feature width */
  float		*feataspect;		/* Feature aspect ratio */
  float		*extra[PROFIT_MAXEXTRA];/* Parameters along extra-dimension */
  float		extrazero[PROFIT_MAXEXTRA]; /* Zero-point along extra-dim. */
  float		extrascale[PROFIT_MAXEXTRA]; /* Scaling along extra-dim. */
  int		extracycleflag[PROFIT_MAXEXTRA]; /* !=0 for cycling dim. */
  interpenum	interptype[2+PROFIT_MAXEXTRA];	/* Interpolation type */
  int		kernelwidth[2+PROFIT_MAXEXTRA];	/* Kernel size */
  float		*kernelbuf;		/* Kernel buffer */
  int		kernelnlines;		/* Number of interp kernel lines */
  }	profstruct;

typedef struct
  {
  objstruct	*obj;		/* Current object */
  obj2struct	*obj2;		/* Current object */
  int		nparam;		/* Number of parameters to be fitted */
  float		*paramlist[PARAM_NPARAM];	/* flat parameter list */
  int		paramindex[PARAM_NPARAM];/* Vector of parameter indices */
  float		param[PARAM_NPARAM];	/* Vector of parameters to be fitted */
  float		paraminit[PARAM_NPARAM];/* Parameter initial guesses */
  float		parammin[PARAM_NPARAM];	/* Parameter lower limits */
  float		parammax[PARAM_NPARAM];	/* Parameter upper limits */
  float		*covar;		/* Covariance matrix */
  float		paramerr[PARAM_NPARAM];	/* Std deviations of parameters */
  int		niter;		/* Number of iterations */
  profstruct	**prof;		/* Array of pointers to profiles */
  int		nprof;		/* Number of profiles to consider */
  struct psf	*psf;		/* PSF */
  float		pixstep;	/* Model/PSF sampling step */
  float		fluxfac;	/* Model flux scaling factor */
  float		subsamp;	/* Subsampling factor */
  float		*psfdft;	/* Compressed Fourier Transform of the PSF */
  float		*psfpix;	/* Full res. pixmap of the PSF */
  float		*modpix;	/* Full res. pixmap of the complete model */
  float		*pmodpix;	/* Full res. pixmap of the partial model */
  int		modnaxisn[3];	/* Dimensions along each axis */
  PIXTYPE	*lmodpix;	/* Low resolution pixmap of the model */
  PIXTYPE	*objpix;	/* Copy of object pixmap */
  PIXTYPE	*objweight;	/* Copy of object weight-map */
  int		objnaxisn[2];	/* Dimensions along each axis */
  int		ix, iy;		/* Integer coordinates of object pixmap */
  float		*resi;		/* Vector of residuals */
  int		nresi;		/* Number of residual elements */
  float		chi2;		/* Std error per residual element */
  float		sigma;		/* Standard deviation of the pixel values */
  float		flux;		/* Total flux in final convolved model */
  float		spirindex;	/* Spiral index (>0 for CCW) */
  }	profitstruct;

/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/

profitstruct	*profit_init(struct psf *psf);

profstruct	*prof_init(profitstruct *profit, proftypenum profcode);

float		*profit_compresi(profitstruct *profit, float dynparam,
				float *resi),
		*profit_residuals(profitstruct *profit, picstruct *field,
			picstruct *wfield, float dynparam,
			float *param, float *resi),
		prof_add(profstruct *prof, profitstruct *profit),
		profit_minradius(profitstruct *profit, float refffac),
		profit_spiralindex(profitstruct *profit);

int		profit_copyobjpix(profitstruct *profit, picstruct *field,
			picstruct *wfield),
		profit_minimize(profitstruct *profit, int niter),
		profit_setparam(profitstruct *profit, paramenum paramtype,
			float param, float parammin, float parammax);

void		prof_end(profstruct *prof),
		profit_addparam(profitstruct *profit, paramenum paramindex,
			float **param),
		profit_boundtounbound(profitstruct *profit, float *param),
		profit_fit(profitstruct *profit,
			picstruct *field, picstruct *wfield,
			objstruct *obj, obj2struct *obj2),
		profit_convolve(profitstruct *profit, float *modpix),
		profit_covarunboundtobound(profitstruct *profit),
		profit_end(profitstruct *profit),
		profit_evaluate(float *par, float *fvec, int m, int n,
			void *adata),
		profit_makedft(profitstruct *profit),
		profit_moments(profitstruct *profit, obj2struct *obj2),
		profit_printout(int n_par, float* par, int m_dat, float* fvec,
			void *data, int iflag, int iter, int nfev ),
		profit_psf(profitstruct *profit),
		profit_resample(profitstruct *profit, float *inpix,
			PIXTYPE *outpix, float factor),
		profit_resetparam(profitstruct *profit, paramenum paramtype),
		profit_resetparams(profitstruct *profit),
		profit_surface(profitstruct *profit, obj2struct *obj2,
			double lostfluxfrac),
		profit_unboundtobound(profitstruct *profit, float *param);

#endif
