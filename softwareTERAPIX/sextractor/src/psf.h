 /*
 				psf.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for psffit.c.
*
*	Last modify:	11/11/99
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

#define	PSF_MAXSHIFT	20.0	/* Max shift from initial guess (pixels)*/
#define	PSF_MINSHIFT	1e-3	/* Min shift from previous guess (pixels)*/
#define PSF_NITER	20	/* Maximum number of iterations in fit */
#define PSF_NA		3	/* Number of fitted parameters per component */
#define PSF_NTOT	(PSF_NA*PSF_NPSFMAX) /* Number of fitted parameters */
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
  double	*maskcurr;	/* Current model */
  double	*mx2,*my2,*mxy;	/* 2nd order moments for each component */
  double	*flux;		/* Flux of each component */
  double	*bt;		/* B/T for each component */
  codestruct	*code;
  }	pcstruct;

typedef struct
  {
  char		name[MAXCHAR];	/* Name of the file containing the PSF data */
  int		maskdim;	/* Dimensionality of the tabulated data */
  int		*masksize;	/* PSF mask dimensions */
  int		masknpix;	/* Total number of involved PSF pixels */
  float		*maskcomp;      /* Complete pix. data (PSF components) */
  double	*maskloc;	/* Local PSF */
  double	**context;	/* Contexts */
  t_type	*contexttyp;	/* Context types */
  char		**contextname;	/* Array of context key-names */
  double	*contextoffset;	/* Offset to apply to context data */
  double	*contextscale;	/* Scaling to apply to context data */
  struct poly	*poly;		/* Polynom describing the PSF variations */
  pcstruct	*pc;		/* PC components */
  double	fwhm;		/* Typical PSF FWHM */
  float		pixstep;	/* PSF sampling step */
  }	psfstruct;

typedef struct
  {
  int		niter;		/* Number of iterations required */
  int		npsf;		/* Number of fitted stars for this detection */
  float		*x,*y;		/* Position derived from the PSF-fitting */
  float		*flux;		/* Flux derived from the PSF-fitting */
  }	psfitstruct;

/*----------------------------- Global variables ----------------------------*/
psfstruct	*thepsf;
psfitstruct	*thepsfit;
PIXTYPE		*checkmask;

/*-------------------------------- functions --------------------------------*/
extern void	psf_build(psfstruct *psf),
		psf_end(psfstruct *psf),
		psf_init(psfstruct *psf),
		svdfit(double *a, double *b, int m, int n, double *sol,
			double *vmat, double *wmat),
		svdvar(double *vmat, double *wmat, int n, double *covmat);

extern psfstruct	*psf_load(char *filename);

extern void	pc_end(pcstruct *pc),
		pc_fit(psfstruct *psf, double *data, double *weight,
		int width, int height, int ix, int iy, double dx, double dy,
		int npc, float backrms),
		psf_fit(psfstruct *psf, picstruct *field, picstruct *wfield,
		objstruct *obj),
		psf_readcontext(psfstruct *psf, picstruct *field);

extern pcstruct	*pc_load(catstruct *cat);
