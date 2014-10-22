/*
*				profit.h
*
* Include file for profit.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2006-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		14/10/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _PROFIT_H_
#define _PROFIT_H_

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _SUBIMAGE_H_
#include "subimage.h"
#endif

#ifndef _PSF_H_
#include "psf.h"
#endif

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
#define		MODEL_MOFFAT		0x0200
#define		MODEL_TABULATED		0x0400
#define		MODEL_NMAX		12

/*--------------------------- convolution flag ------------------------------*/

#define		PROFIT_NOCONV		0
#define		PROFIT_CONV		1

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
#define	PROFIT_DYNPARAM	3.0	/* Dynamic compression param. in sigma units */
#define	PROFIT_SMOOTHR	4.0	/* Profile smoothing radius (pixels) */
#define	PROFIT_MAXMODSIZE  5120	/* Maximum size allowed for the model raster */
#define PROFIT_MAXSMODSIZE 64	/* Number of model planes */
#define	PROFIT_MAXOBJSIZE  5120	/* Maximum size allowed for the object raster */
#define	PROFIT_BARXFADE	0.1	/* Fract. of bar length crossfaded with arms */
#define	PROFIT_MAXEXTRA	3	/* Max. nb of extra free params of profiles */
#define INTERP_MAXKERNELWIDTH	8	/* Max. range of kernel (pixels) */
#define	PROFIT_MAXNBAND	5
/* NOTES:
One must have:	PROFIT_NITER > 0
		PROFIT_MAXEXTRA > 0
*/

/*--------------------------------- typedefs --------------------------------*/

typedef enum	{INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}       interpenum;

typedef enum	{PARAM_BACK,
		PARAM_DIRAC_FLUX, PARAM_DIRAC_FLUX2,PARAM_DIRAC_FLUX3,
		PARAM_DIRAC_FLUX4,PARAM_DIRAC_FLUX5,
		PARAM_X, PARAM_Y,
		PARAM_SPHEROID_FLUX, PARAM_SPHEROID_FLUX2,PARAM_SPHEROID_FLUX3,
		PARAM_SPHEROID_FLUX4,PARAM_SPHEROID_FLUX5,
		PARAM_SPHEROID_REFF, PARAM_SPHEROID_ASPECT,
		PARAM_SPHEROID_POSANG, PARAM_SPHEROID_SERSICN,
		PARAM_DISK_FLUX, PARAM_DISK_FLUX2,PARAM_DISK_FLUX3,
		PARAM_DISK_FLUX4,PARAM_DISK_FLUX5,
		PARAM_DISK_SCALE, PARAM_DISK_ASPECT, PARAM_DISK_POSANG,
		PARAM_ARMS_FLUX, PARAM_ARMS_QUADFRAC, PARAM_ARMS_SCALE,
		PARAM_ARMS_START, PARAM_ARMS_POSANG, PARAM_ARMS_PITCH,
		PARAM_ARMS_PITCHVAR, PARAM_ARMS_WIDTH,
		PARAM_BAR_FLUX, PARAM_BAR_ASPECT, PARAM_BAR_POSANG,
		PARAM_INRING_FLUX, PARAM_INRING_WIDTH, PARAM_INRING_ASPECT,
		PARAM_OUTRING_FLUX, PARAM_OUTRING_START, PARAM_OUTRING_WIDTH,
		PARAM_MOFFAT_FLUX, PARAM_MOFFAT_ALPHA, PARAM_MOFFAT_ASPECT,
		PARAM_MOFFAT_POSANG, PARAM_MOFFAT_BETA,
		PARAM_MOFFAT_MINKP, PARAM_MOFFAT_OFFSET,
		PARAM_NPARAM}	paramenum;

typedef enum 	{PARFIT_FIXED, PARFIT_UNBOUND, PARFIT_LINBOUND,
		PARFIT_LOGBOUND}	parfitenum;

/*--------------------------- structure definitions -------------------------*/

typedef struct
  {
  paramenum	param_code;		/* Model parameter code */
  float		*value;			/* Vector of values */
  int		*index;			/* Indices to vector elements */
  parfitenum	fit_code;		/* Fitting behaviour code */
  float		min,max;		/* Min and max values */
  }	profparamstruct;

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
  profparamstruct	*profparam;	/* Model parameters */
  float		*flux[PROFIT_MAXNBAND];	/* Integrated flux */
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

typedef struct profit
  {
  int		nparam;		/* Number of parameters to be fitted */
  float		*paramlist[PARAM_NPARAM];	/* flat parameter list */
  int		paramindex[PARAM_NPARAM];/* Vector of parameter indices */
  int		paramrevindex[PARAM_NPARAM];/* Vector of reversed indices */
  parfitenum	parfittype[PARAM_NPARAM];/* Parameter fitting: fixed,bounded,.*/
  float		param[PARAM_NPARAM];	/* Vector of parameters to be fitted */
  float		paraminit[PARAM_NPARAM];/* Parameter initial guesses */
  float		parammin[PARAM_NPARAM];	/* Parameter lower limits */
  float		parammax[PARAM_NPARAM];	/* Parameter upper limits */
  int		paramsize[PARAM_NPARAM]; /* Parameter vector size */
  int		nlimmin;	/* # of parameters that hit the min limit */
  int		nlimmax;	/* # of parameters that hit the max limit */
  float		paramerr[PARAM_NPARAM];	/* Std deviations of parameters */
  float		*covar;		/* Covariance matrix */
  int		iter;		/* Iteration counter */
  int		niter;		/* Number of iterations */
  profstruct	**prof;		/* Array of pointers to profiles */
  int		nprof;		/* Number of profiles to consider */
  struct subprofit *subprofit;	/* Array of subprofile structures */
  int		nsubprofit;	/* Number of subprofiles */
  float		*resi;		/* Vector of residuals */
  int		nresi;		/* Number of residual elements */
  float		chi2;		/* Std error per residual element */
/* Buffers */
  double	dparam[PARAM_NPARAM];
  int		flag;		/* Model fitting flag */
  int		conv_flag;	/* Convolution flag */
  }	profitstruct;

typedef struct subprofit
  {
  int		index;		/* sub-profile index */
  struct field	*field;		/* Image field */
  struct field	*wfield;	/* Weight field */
  struct psf	*psf;		/* PSF */
  struct fftscratch
		*fftscratch;	/* FFT scratch space/plans */
  float		pixstep;	/* Model/PSF sampling step */
  float		fluxfac;	/* Model flux scaling factor */
  float		subsamp;	/* Subsampling factor */
  float		*psfpix;	/* Full res. pixmap of the PSF */
  float		*psfdft;	/* Compressed Fourier Transform of the PSF */
  PIXTYPE	*lmodpix;	/* Low resolution pixmaps of the model */
  PIXTYPE	*objpix;	/* Copy of object pixmaps */
  PIXTYPE	*objweight;	/* Copy of object weight-maps */
  int		objnaxisn[2];	/* Dimensions along each axis */
  int		nobjpix;	/* Total number of "final" pixels */
  float		*modpix;	/* Full res. pixmap of the complete model */
  float		*cmodpix;	/* Full res. pixmap of the convolved model */
  int		modnaxisn[3];	/* Dimensions along each axis */
  int		nmodpix;	/* Total number of model pixels */
  int		ix, iy;		/* Integer coordinates of object pixmap */
  float		spirindex;	/* Spiral index (>0 for CCW) */
  float		sigma;		/* Standard deviation of the pixel values */
  float		guesssigbkg;	/* Best guess for background noise sigma */
  float		guessdx,guessdy;/* Best guess for relative source coordinates */
  float		guessflux;	/* Best guess for typical source flux (>0) */
  float		guessfluxmax;	/* Best guess for source flux upper limit (>0)*/
  float		guessradius;	/* Best guess for source half-light radius */
  float		guessaspect;	/* Best guess for source aspect ratio */
  float		guessposang;	/* Best guess for source position angle */
  float		flux;		/* Total flux in final convolved models */
  }	subprofitstruct;

/*----------------------------- Global variables ----------------------------*/
/*-------------------------------- functions --------------------------------*/

profitstruct	*profit_init(objstruct *obj, subimagestruct *subimage,
			int nsubimage, unsigned int modeltype, int conv_flag);

profstruct	*prof_init(profitstruct *profit, unsigned int modeltype);

float		*profit_compresi(profitstruct *profit, float dynparam,
				float *resi),
		*profit_residuals(profitstruct *profit, float dynparam,
			float *param, float *resi),
		prof_add(subprofitstruct *subprofit, profstruct *prof,
			int extfluxfac_flag),
		profit_minradius(profitstruct *profit, float refffac),
		subprofit_noisearea(subprofitstruct *profit),
		subprofit_spiralindex(subprofitstruct *subprofit,
			obj2struct *obj2);

int		profit_boundtounbound(profitstruct *profit,
			float *param, double *dparam, int index),
		subprofit_copyobjpix(subprofitstruct *subprofit,
			subimagestruct *subimage),
		profit_covarunboundtobound(profitstruct *profit,
			double *dparam, float *param),
		profit_fit(profitstruct *profit, objstruct *obj),
		profit_minimize(profitstruct *profit, int niter),
		prof_moments(profitstruct *profit, profstruct *prof,
				double *jac),
		profit_resample(profitstruct *profit,
			subprofitstruct *subprofit,
			float *inpix, PIXTYPE *outpix, float factor),
		profit_setparam(profitstruct *profit, paramenum paramtype,
			float param, float parammin, float parammax,
			parfitenum parfittype),
		profit_unboundtobound(profitstruct *profit,
			double *dparam, float *param, int index);

void		prof_end(profstruct *prof),
		profit_addparam(profitstruct *profit, paramenum paramindex,
			float **param),
		subprofit_convmoments(subprofitstruct *subprofit,
			obj2struct *obj2),
		subprofit_convolve(subprofitstruct *subprofit, float *modpix),
		profit_end(profitstruct *profit),
		profit_evaluate(double *par, double *fvec, int m, int n,
			void *adata),
		subprofit_fluxcor(subprofitstruct *subprofit, obj2struct *obj2),
		subprofit_makedft(subprofitstruct *subprofit),
		profit_measure(profitstruct *profit, obj2struct *obj2),
		profit_moments(profitstruct *profit, obj2struct *obj2),
		profit_printout(int n_par, float* par, int m_dat, float* fvec,
			void *data, int iflag, int iter, int nfev ),
		subprofit_psf(subprofitstruct *subprofit, obj2struct *obj2),
		profit_resetparam(profitstruct *profit, paramenum paramtype),
		profit_resetparams(profitstruct *profit),
		profit_spread(profitstruct *profit,  fieldstruct *field,
			fieldstruct *wfield, objstruct *obj),
		subprofit_addmodpix(subprofitstruct *subprofitmod,
			PIXTYPE *pixout, int ix, int iy, int width, int height,
			float oversamp, float fac),
		subprofit_surface(profitstruct *profit,
			subprofitstruct *subprofit, obj2struct *obj2),
		subprofit_end(subprofitstruct *subprofit);
#endif
