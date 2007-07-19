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
*	Last modify:	19/07/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _PROFIT_H_
#define _PROFIT_H_

/*-------------------------------- flags ------------------------------------*/

#define		PROFIT_FLIPPED	0x0001

/*----------------------------- Internal constants --------------------------*/

#define	PROFIT_MAXITER	1000	/* Max. nb of iterations in profile fitting */
#define	PROFIT_OVERSAMP	1	/* Profile oversampling factor on each axis */
#define	PROFIT_MAXPROF	8	/* Max. nb of profile components */
#define	PROFIT_DYNPARAM	10.0	/* Dynamic compression param. in sigma units */
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
			PROF_EXPONENTIAL, PROF_ARMS, PROF_BAR,
			PROF_SERSIC_TABEX, PROF_NPROF}
				proftypenum; /* Profile code */
typedef enum	{INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}       interpenum;

typedef enum	{PARAM_BACK, PARAM_X, PARAM_Y,
		PARAM_DEVAUC_FLUX, PARAM_DEVAUC_MAJ, PARAM_DEVAUC_ASPECT,
		PARAM_DEVAUC_POSANG,
		PARAM_EXPO_FLUX, PARAM_EXPO_MAJ, PARAM_EXPO_ASPECT,
		PARAM_EXPO_POSANG,
		PARAM_SERSIC_FLUX, PARAM_SERSIC_MAJ, PARAM_SERSIC_ASPECT,
		PARAM_SERSIC_POSANG, PARAM_SERSIC_N,
		PARAM_ARMS_FLUX, PARAM_ARMS_QUADFRAC, PARAM_ARMS_SCALE,
		PARAM_ARMS_START, PARAM_ARMS_POSANG, PARAM_ARMS_PITCH,
		PARAM_ARMS_WIDTH,
		PARAM_BAR_FLUX, PARAM_BAR_ASPECT, PARAM_BAR_POSANG,
		PARAM_NPARAM}	paramenum;

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  proftypenum	code;			/* Model code */
  double	*pix;			/* Full pixmap of the model */
  int		naxis;			/* Number of pixmap dimensions */
  int		naxisn[2+PROFIT_MAXEXTRA];	/* Pixmap size for each axis */
  double	typscale;		/* Typical scale in prof pixels */
  double	scaling;		/* Scaling factor for lengths */
/* Generic presentation parameters */
  double	*flux;			/* Pointer to integrated flux */
  double	*x[2];			/* Pointer to coordinate vector */
  double	*scale;			/* Pointer to scaling vector */
  double	*aspect;		/* Pointer to aspect ratio */
  double	*posangle;		/* Pointer to pos. angle (CCW/NAXIS1)*/
  double	*armsquadfrac;		/* Pointer to the quad. arm fraction */
  double	*armscale;		/* Arm scalelength relative to disk */
  double	*armstart;		/* Pointer to arm starting radius */
  double	*armposang;		/* Pointer to arm position angle */
  double	*armpitch;		/* Pointer to arm pitch */
  double	*armwidth;		/* Pointer to arm width */
  double	*baramp;		/* Pointer to bar amplitude */
  double	*baraspect;		/* Pointer to bar aspect ratio */
  double	*barposang;		/* Pointer to bar position angle */
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
  int		nparam;		/* Number of parameters to be fitted */
  double	*paramlist[PARAM_NPARAM];	/* flat parameter list */
  double	param[PARAM_NPARAM];	/* Vector of parameters to be fitted */
  double	paraminit[PARAM_NPARAM];/* Parameter initial guesses */
  double	parammin[PARAM_NPARAM];	/* Parameter lower limits */
  double	parammax[PARAM_NPARAM];	/* Parameter upper limits */
  int		niter;		/* Number of iterations */
  profstruct	**prof;		/* Array of pointers to profiles */
  int		nprof;		/* Number of profiles to consider */
  struct psf	*psf;		/* PSF */
  double	*psfdft;	/* Compressed Fourier Transform of the PSF */
  double	*modpix;	/* Full res. pixmap of the complete model */
  double	*pmodpix;	/* Full res. pixmap of the partial model */
  int		modnaxisn[2];	/* Dimensions along each axis */
  PIXTYPE	*lmodpix;	/* Low resolution pixmap of the model */
  PIXTYPE	*objpix;	/* Copy of object pixmap */
  PIXTYPE	*objweight;	/* Copy of object weight-map */
  int		objnaxisn[2];	/* Dimensions along each axis */
  double	*resi;		/* Vector of residuals */
  int		nresi;		/* Number of residual elements */
  double	chi2;		/* Std error per residual element */
  double	sigma;		/* Standard deviation of the pixel values */
  double	flux;		/* Total flux in final convolved model */
  double	spirindex;	/* Spiral index (>0 for CCW) */
  }	profitstruct;

/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/

profitstruct	*profit_init(struct psf *psf, proftypenum *profcode,
			int nprof);

profstruct	*prof_init(profitstruct *profit, proftypenum profcode);

double		*profit_compresi(profitstruct *profit, picstruct *field,
			picstruct *wfield, objstruct *obj, double *resi),
		*profit_residuals(profitstruct *profit, picstruct *field,
			picstruct *wfield, objstruct *obj, double *param,
			double *resi),
		profit_spiralindex(profitstruct *profit, objstruct *obj,
			obj2struct *obj2);

PIXTYPE		*profit_resample(profitstruct *profit);

int		profit_copyobjpix(profitstruct *profit, picstruct *field,
				int ix, int iy),
		profit_minimize(profitstruct *profit, int niter),
		profit_resetparam(profitstruct *profit, proftypenum proftype,
			double param, double parammin, double parammax);

void		prof_add(profstruct *prof, profitstruct *profit),
		prof_end(profstruct *prof),
		profit_addparam(profitstruct *profit, paramenum paramindex,
			double **param),
		profit_boundtounbound(profitstruct *profit, double *param),
		profit_fit(profitstruct *profit,
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
			void *data, int iflag, int iter, int nfev ),
		profit_resetparams(profitstruct *profit, objstruct *obj,
			obj2struct *obj2),
		profit_unboundtobound(profitstruct *profit, double *param);

#endif
