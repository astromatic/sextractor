/*
*				profit.h
*
* Include file for profit.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2006-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		13/02/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _PROFIT_H_
#define _PROFIT_H_

/*-------------------------------- models -----------------------------------*/

#define		MODEL_NONE		0x0000
#define		MODEL_BACK		0x0001
#define		MODEL_DIRAC		0x0002
#define		MODEL_SERSIC		0x0004
#define		MODEL_DEVAUCOULEURS	0x0008
#define		MODEL_EXPONENTIAL	0x0010
#define		MODEL_ARMS		0x0020
#define		MODEL_BAR		0x0040
#define		MODEL_INRING		0x0080
#define		MODEL_OUTRING		0x0100
#define		MODEL_TABULATED		0x0200
#define		MODEL_NMAX		11

/*--------------------------- fitting flags ---------------------------------*/

#define		PROFLAG_MODSUB		0x0001
#define		PROFLAG_OBJSUB		0x0002
#define		PROFLAG_NOTCONST	0x0004
#define		PROFLAG_MINLIM		0x0008
#define		PROFLAG_MAXLIM		0x0010

/*------------------------- parameter type flags ----------------------------*/

#define		PROFPARAM_UNBOUNDED	1
#define		PROFPARAM_LINBOUNDED	2
#define		PROFPARAM_LOGBOUNDED	3

/*-------------------------------- macros -----------------------------------*/

#define		PROFIT_POW(x,a)		(x>0.01? exp(a*log(x)) : pow(x,a))
#define		PROFIT_POWF(x,a)	(x>0.01? expf(a*logf(x)) : powf(x,a))

/*----------------------------- Internal constants --------------------------*/

#define	PARAM_ALLPARAMS	(-1)	/* All parameters */
#define	PROFIT_MAXITER	1000	/* Max. nb of iterations in profile fitting */
#define	PROFIT_MAXPROF	8	/* Max. nb of profile components */
#define	PROFIT_HIDEFRES	201	/* Hi. def. model resol. (must be <MAXMODSIZE)*/
#define	PROFIT_REFFFAC	3.0	/* Factor in r_eff for measurement radius*/
#define	PROFIT_MAXR2MAX	1e6	/* Maximum r2_max for truncating profiles */
#define	PROFIT_DYNPARAM	10.0	/* Dynamic compression param. in sigma units */
#define	PROFIT_SMOOTHR	4.0	/* Profile smoothing radius (pixels) */
#define	PROFIT_MAXMODSIZE  512	/* Maximum size allowed for the model raster */
#define PROFIT_MAXSMODSIZE 64	/* Number of model planes */
#define	PROFIT_MAXOBJSIZE  512	/* Maximum size allowed for the object raster */
#define	PROFIT_BARXFADE	0.1	/* Fract. of bar length crossfaded with arms */
#define	PROFIT_MAXEXTRA	2	/* Max. nb of extra free params of profiles */
#define INTERP_MAXKERNELWIDTH	8	/* Max. range of kernel (pixels) */
/* NOTES:
One must have:	PROFIT_NITER > 0
		PROFIT_MAXEXTRA > 0
*/

/*--------------------------------- typedefs --------------------------------*/

typedef enum	{INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}       interpenum;

typedef enum	{PARAM_BACK,
		PARAM_DIRAC_FLUX, PARAM_X, PARAM_Y,
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

typedef enum 	{PARFIT_FIXED, PARFIT_UNBOUND, PARFIT_LINBOUND,
		PARFIT_LOGBOUND}	parfitenum;

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  unsigned int	code;			/* Model code */
  char		*name;			/* Model name */
  float		*pix;			/* Full pixmap of the model */
  int		naxis;			/* Number of pixmap dimensions */
  int		naxisn[3];		/* Pixmap size for each axis */
  int		npix;			/* Total number of prof pixels */
  float		typscale;		/* Typical scale in prof pixels */
  float		fluxfac;		/* Flux normalisation factor */
  float		lostfluxfrac;		/* Lost flux fraction */
  float		m0,mx2,my2,mxy;		/* 2nd order moments */
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
  int		paramrevindex[PARAM_NPARAM];/* Vector of reversed indices */
  parfitenum	parfittype[PARAM_NPARAM];/* Parameter fitting: fixed,bounded,.*/
  float		param[PARAM_NPARAM];	/* Vector of parameters to be fitted */
  float		paraminit[PARAM_NPARAM];/* Parameter initial guesses */
  float		parammin[PARAM_NPARAM];	/* Parameter lower limits */
  float		parammax[PARAM_NPARAM];	/* Parameter upper limits */
  double	dparampcen[PARAM_NPARAM];/* Parameter prior center */
  double	dparampsig[PARAM_NPARAM];/* Parameter prior sigma */
  int		nlimmin;	/* # of parameters that hit the min limit */
  int		nlimmax;	/* # of parameters that hit the max limit */
  float		paramerr[PARAM_NPARAM];	/* Std deviations of parameters */
  float		*covar;		/* Covariance matrix */
  int		iter;		/* Iteration counter */
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
  float		*modpix2;	/* 2nd full res. pixmap of the complete model */
  float		*cmodpix;	/* Full res. pixmap of the convolved model */
  int		modnaxisn[3];	/* Dimensions along each axis */
  int		nmodpix;	/* Total number of model pixels */
  PIXTYPE	*lmodpix;	/* Low resolution pixmap of the model */
  PIXTYPE	*lmodpix2;	/* 2nd low resolution pixmap of the model */
  PIXTYPE	*objpix;	/* Copy of object pixmap */
  PIXTYPE	*objweight;	/* Copy of object weight-map */
  int		objnaxisn[2];	/* Dimensions along each axis */
  int		nobjpix;	/* Total number of "final" pixels */
  int		ix, iy;		/* Integer coordinates of object pixmap */
  float		*resi;		/* Vector of residuals */
  int		nresi;		/* Number of residual elements */
  float		*presi;		/* Vector of prior residuals */
  int		npresi;		/* Number of prior residual elements */
  float		chi2;		/* Std error per residual element */
  float		sigma;		/* Standard deviation of the pixel values */
  float		flux;		/* Total flux in final convolved model */
  float		spirindex;	/* Spiral index (>0 for CCW) */
  float		guesssigbkg;	/* Best guess for background noise sigma */
  float		guessdx,guessdy;/* Best guess for relative source coordinates */
  float		guessflux;	/* Best guess for typical source flux (>0) */
  float		guessfluxmax;	/* Best guess for source flux upper limit (>0)*/
  float		guessradius;	/* Best guess for source half-light radius */
  float		guessaspect;	/* Best guess for source aspect ratio */
  float		guessposang;	/* Best guess for source position angle */
/* Buffers */
  double	dparam[PARAM_NPARAM];
  }	profitstruct;

/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/

profitstruct	*profit_init(struct psf *psf, unsigned int modeltype);

profstruct	*prof_init(profitstruct *profit, unsigned int modeltype);

float		*profit_compresi(profitstruct *profit, float dynparam,
				float *resi),
		*profit_presiduals(profitstruct *profit, double *dparam,
			float *presi),
		*profit_residuals(profitstruct *profit, picstruct *field,
			picstruct *wfield, float dynparam,
			float *param, float *resi),
		prof_add(profitstruct *profit, profstruct *prof,
			int extfluxfac_flag),
		profit_minradius(profitstruct *profit, float refffac),
		profit_noisearea(profitstruct *profit),
		profit_spiralindex(profitstruct *profit);

int		profit_boundtounbound(profitstruct *profit,
			float *param, double *dparam, int index),
		profit_copyobjpix(profitstruct *profit, picstruct *field,
			picstruct *wfield),
		profit_covarunboundtobound(profitstruct *profit,
			double *dparam, float *param),
		profit_minimize(profitstruct *profit, int niter),
		prof_moments(profitstruct *profit, profstruct *prof,
				double *jac),
		profit_resample(profitstruct *profit, float *inpix,
			PIXTYPE *outpix, float factor),
		profit_setparam(profitstruct *profit, paramenum paramtype,
			float param, float parammin, float parammax,
			parfitenum parfittype,
			float priorcen, float priorsig),
		profit_unboundtobound(profitstruct *profit,
			double *dparam, float *param, int index);

void		profit_dfit(profitstruct *profit, profitstruct *dprofit,
			picstruct *field, picstruct *dfield,
			picstruct *wfield, picstruct *dwfield,
			objstruct *obj, obj2struct *obj2),
		prof_end(profstruct *prof),
		profit_addparam(profitstruct *profit, paramenum paramindex,
			float **param),
		profit_fit(profitstruct *profit,
			picstruct *field, picstruct *wfield,
			objstruct *obj, obj2struct *obj2),
		profit_convmoments(profitstruct *profit, obj2struct *obj2),
		profit_convolve(profitstruct *profit, float *modpix),
		profit_end(profitstruct *profit),
		profit_evaluate(double *par, double *fvec, int m, int n,
			void *adata),
		profit_fluxcor(profitstruct *profit, objstruct *obj,
				obj2struct *obj2),
		profit_makedft(profitstruct *profit),
		profit_moments(profitstruct *profit, obj2struct *obj2),
		profit_printout(int n_par, float* par, int m_dat, float* fvec,
			void *data, int iflag, int iter, int nfev ),
		profit_psf(profitstruct *profit),
		profit_resetparam(profitstruct *profit, paramenum paramtype),
		profit_resetparams(profitstruct *profit),
		profit_surface(profitstruct *profit, obj2struct *obj2);

#endif
