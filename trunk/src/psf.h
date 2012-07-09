/*
*				psf.h
*
* Include file for psf.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1998-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		13/06/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*----------------------------- Internal constants --------------------------*/

#define	PSF_MAXSHIFT	20.0	/* Max shift from initial guess (pixels)*/
#define	PSF_MINSHIFT	1e-3	/* Min shift from previous guess (pixels)*/
#define PSF_NITER	20	/* Maximum number of iterations in fit */
#define PSF_NA		3	/* Number of fitted parameters per component */
#define PSF_NTOT	(PSF_NA*PSF_NPSFMAX)	/* Number of fitted parameters*/
#define PSF_DOUBLETOT   ((PSF_NA+1)*PSF_NPSFMAX)/* Nb of fitted parameters */
#define	PC_NITER	1	/* Maximum number of iterations in PC fit */

/* NOTES:
One must have:	PSF_MAXSHIFT > 0.0
		PSF_NPSF >= 1
		PSF_NITER >= 1
*/

/*--------------------------- structure definitions -------------------------*/

typedef struct code
  {
  float		*pc;
  float		**param;
  int		*parammod;
  int		ncode;
  int		nparam;
  }		codestruct;

typedef struct pc
  {
  char		name[MAXCHAR];	/* PC filename */
  int		npc;		/* Number of Principal Components */
  int		maskdim;	/* Dimensionality of the tabulated data */
  int		*masksize;	/* PC mask dimensions */
  int		masknpix;	/* Total number of involved PC pixels */
  float		*maskcomp; 	/* Convolved pix data (principal components) */
  int		omaskdim;	/* Dimensionality of the tabulated data */
  int		*omasksize;	/* PC mask dimensions */
  int		omasknpix;	/* Total number of involved PC pixels */
  float		*omaskcomp; 	/* Original pix data (principal components) */
  float		*maskcurr;	/* Current model */
  float		*mx2,*my2,*mxy;	/* 2nd order moments for each component */
  float		*flux;		/* Flux of each component */
  float		*bt;		/* B/T for each component */
  codestruct	*code;
  }	pcstruct;

typedef struct psf
  {
  char		name[MAXCHAR];	/* Name of the file containing the PSF data */
  int		maskdim;	/* Dimensionality of the tabulated data */
  int		*masksize;	/* PSF mask dimensions */
  int		masknpix;	/* Total number of involved PSF pixels */
  float		*maskcomp;      /* Complete pix. data (PSF components) */
  float		*maskloc;	/* Local PSF */
  double	**context;	/* Contexts */
  t_type	*contexttyp;	/* Context types */
  char		**contextname;	/* Array of context key-names */
  double	*contextoffset;	/* Offset to apply to context data */
  double	*contextscale;	/* Scaling to apply to context data */
  struct poly	*poly;		/* Polynom describing the PSF variations */
  pcstruct	*pc;		/* PC components */
  double	fwhm;		/* Typical PSF FWHM */
  float		pixstep;	/* PSF sampling step */
  int		build_flag;	/* Set if the current PSF has been computed */
  }	psfstruct;

typedef struct
  {
  int		niter;		/* Number of iterations required */
  int		npsf;		/* Number of fitted stars for this detection */
  double	*x,*y;		/* Position derived from the PSF-fitting */
  float		*flux;		/* Flux derived from the PSF-fitting */
  float		*fluxerr;	/* Flux error estimated from the PSF-fitting */
  }	psfitstruct;

/*----------------------------- Global variables ----------------------------*/
psfstruct	*psf,*thedpsf,*thepsf;
psfitstruct	*thepsfit,*thedpsfit;
PIXTYPE		*checkmask;

/*-------------------------------- functions --------------------------------*/
extern void	compute_pos(int *pnpsf,int *pconvflag,int *pnpsfflag,
			double radmin2, double radmax2,double r2, double *sol,
			double *flux , double *deltax,double *deltay,
			double *pdx,double *pdy),
		compute_pos_phot(int *pnpsf,double *sol,double *flux),
		compute_poserr(int j,double *var,double *sol,obj2struct *obj2,
			double *x2, double *y2,double *xy, int npsf),
		psf_build(psfstruct *psf),
		psf_end(psfstruct *psf, psfitstruct *psfit),
		psf_init(void),
		svdfit(double *a, float *b, int m, int n, double *sol,
			double *vmat, double *wmat),
		svdvar(double *vmat, double *wmat, int n, double *covmat);

extern double	*compute_gradient (float *weight,int width, int height,
			float *masks, float *maskx, float *masky,
			double *mat),
		*compute_gradient_phot(float *weight,int width, int height,
			float *masks, double *pm),
		psf_fwhm(psfstruct *psf);

extern psfstruct	*psf_load(char *filename, int ext);

extern void	pc_end(pcstruct *pc),
		pc_fit(psfstruct *psf, float *data, float *weight,
		int width, int height, int ix, int iy, float dx, float dy,
		int npc, float backrms),
		double_psf_fit(psfstruct *psf, picstruct *field,
			picstruct *wfield, objstruct *obj,
			psfstruct *dpsf, picstruct *dfield, picstruct *dwfield),
		psf_fit(psfstruct *psf, picstruct *field, picstruct *wfield,
		objstruct *obj),
		psf_readcontext(psfstruct *psf, picstruct *field);

extern pcstruct	*pc_load(catstruct *cat);

