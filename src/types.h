/*
*				types.h
*
* Global type definitions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		05/07/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include <stdio.h>

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif
#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

/*-------------------------------- flags ------------------------------------*/

#define		OBJ_CROWDED	0x0001
#define		OBJ_MERGED	0x0002
#define		OBJ_SATUR	0x0004
#define		OBJ_TRUNC	0x0008
#define		OBJ_APERT_PB	0x0010
#define		OBJ_ISO_PB	0x0020
#define		OBJ_DOVERFLOW	0x0040
#define		OBJ_OVERFLOW	0x0080

/*----------------------------- weight flags --------------------------------*/

#define		OBJ_LOWWEIGHT	0x0001
#define		OBJ_LOWDWEIGHT	0x0002

/*---------------------------- preanalyse flags -----------------------------*/

#define		ANALYSE_FAST		0
#define		ANALYSE_FULL		1
#define		ANALYSE_ROBUST		2

/*--------------------------------- typedefs --------------------------------*/
typedef	unsigned char	BYTE;			/* a byte */
typedef	unsigned short	USHORT;			/* 0 to 65535 integers */
typedef	char		pliststruct;		/* Dummy type for plist */

typedef	int		LONG;
typedef	unsigned int	ULONG;

typedef  enum {BACK_RELATIVE, BACK_ABSOLUTE}
		backenum;				/* BACK_TYPE */

typedef  enum {CHECK_NONE, CHECK_IDENTICAL, CHECK_BACKGROUND,
        CHECK_BACKRMS, CHECK_MINIBACKGROUND, CHECK_MINIBACKRMS,
        CHECK_SUBTRACTED, CHECK_FILTERED, CHECK_OBJECTS, CHECK_APERTURES,
	CHECK_SEGMENTATION, CHECK_MASK, CHECK_SUBMASK, CHECK_ASSOC,
	CHECK_SUBOBJECTS, CHECK_PSFPROTOS, CHECK_SUBPSFPROTOS,
	CHECK_PCPROTOS, CHECK_SUBPCPROTOS, CHECK_PCOPROTOS,
	CHECK_MAPSOM, CHECK_PROFILES, CHECK_SUBPROFILES,
	CHECK_SPHEROIDS, CHECK_SUBSPHEROIDS, CHECK_DISKS, CHECK_SUBDISKS,
	CHECK_PATTERNS,CHECK_OTHER,
	MAXCHECK}
		checkenum;
	/* CHECK_IMAGE type */

typedef  enum {WEIGHT_NONE, WEIGHT_FROMBACK, WEIGHT_FROMRMSMAP,
		WEIGHT_FROMVARMAP, WEIGHT_FROMWEIGHTMAP, WEIGHT_FROMINTERP}
		weightenum;				/* WEIGHT_IMAGE type */

/*--------------------------------- objects ---------------------------------*/
/* I: "PIXEL" parameters */

typedef struct
  {
/* ---- basic parameters */
  int		number;				/* ID */
  int		fdnpix;				/* nb of extracted pix */
  int		dnpix;				/* nb of pix above thresh  */
  int		npix;				/* "" in measured frame */
  int		nzdwpix;			/* nb of zero-dweights around */
  int		nzwpix;				/* nb of zero-weights inside */
  float		fdflux;				/* integrated ext. flux */
  float		dflux;				/* integrated det. flux */
  float		flux;				/* integrated mes. flux */
  float		fluxerr;			/* integrated variance */
  PIXTYPE	fdpeak;				/* peak intensity (ADU) */
  PIXTYPE	dpeak;				/* peak intensity (ADU) */
  PIXTYPE	peak;				/* peak intensity (ADU) */
/* ---- astrometric data */
  int		peakx,peaky;			/* pos of brightest pix */
  double       	mx, my;				/* barycenter */
  double	poserr_mx2, poserr_my2,
		poserr_mxy;			/* Error ellips moments */
/* ---- morphological data */			
  int		xmin,xmax,ymin,ymax,ycmin,ycmax;/* x,y limits */
  PIXTYPE	*blank, *dblank; 	       	/* BLANKing sub-images  */
  int		*submap;			/* Pixel-index sub-map */
  int		subx,suby, subw,subh;		/* sub-image pos. and size */
  short		flag;				/* extraction flags */
  BYTE		wflag;				/* weighted extraction flags */
  FLAGTYPE	imaflag[MAXFLAG];		/* flags from FLAG-images */
  BYTE		singuflag;			/* flags for singularities */
  int		imanflag[MAXFLAG];     		/* number of MOST flags */
  double	mx2,my2,mxy;			/* variances and covariance */
  float		a, b, theta, abcor;		/* moments and angle */
  float		cxx,cyy,cxy;			/* ellipse parameters */
  int		firstpix;			/* ptr to first pixel */
  int		lastpix;			/* ptr to last pixel */
  float		bkg, dbkg, sigbkg, dsigbkg;	/* Background stats (ADU) */
  float		thresh;				/* measur. threshold (ADU) */
  float		dthresh;		       	/* detect. threshold (ADU) */
  float		mthresh;		       	/* max. threshold (ADU) */
  int		iso[NISO];			/* isophotal areas */
  float		fwhm;				/* IMAGE FWHM */
  
  }	objstruct;

/* II: "BLIND" parameters */
typedef struct
  {
/* ---- photometric data */
  float		flux_iso;			/* ISO integrated flux */
  float		fluxerr_iso;			/* RMS error on ISO flux */
  float		mag_iso;			/* ISO mag */
  float		magerr_iso;			/* ISO mag uncertainty */
  float		flux_isocor;			/* ISOCOR integrated flux */
  float		fluxerr_isocor;			/* RMS error on ISOCOR flux */
  float		mag_isocor;			/* ISOCOR mag */
  float		magerr_isocor;			/* ISOCOR mag uncertainty */
  float		kronfactor;			/* kron parameter */
  float		flux_auto;			/* AUTO integrated flux */
  float		fluxerr_auto;			/* RMS error on AUTO flux */
  float		mag_auto;			/* AUTO mag */
  float		magerr_auto;			/* AUTO mag uncertainty */
  float		petrofactor;			/* kron parameter */
  float		flux_petro;			/* AUTO integrated flux */
  float		fluxerr_petro;			/* RMS error on AUTO flux */
  float		mag_petro;			/* AUTO mag */
  float		magerr_petro;			/* AUTO mag uncertainty */
  float		flux_best;			/* BEST integrated flux */
  float		fluxerr_best;			/* RMS error on BEST flux */
  float		mag_best;			/* BEST mag */
  float		magerr_best;			/* BEST mag uncertainty */
  float		*flux_aper;			/* APER flux vector */
  float		*fluxerr_aper;			/* APER flux error vector  */
  float		*mag_aper;			/* APER magnitude vector */
  float		*magerr_aper;			/* APER mag error vector */
  float		flux_win;			/* WINdowed flux*/
  float		fluxerr_win;			/* WINdowed flux error */
  float		snr_win;			/* WINdowed SNR */
  float		mag_win;			/* WINdowed magnitude */
  float		magerr_win;			/* WINdowed magnitude error */
/* ---- astrometric data */
  BYTE		area_flagw;			/* World area or SB flag */
  double	posx,posy;			/* "FITS" pos. in pixels */
  double	jacob[NAXIS*NAXIS];		/* Local deproject. Jacobian */
  double	pixscale2;			/* Local pixel area */
  double	mamaposx,mamaposy;		/* "MAMA" pos. in pixels */
  float		sposx,sposy;			/* single precision pos. */
  float		poserr_a, poserr_b,
		poserr_theta;			/* Error ellips parameters */
  float		poserr_cxx, poserr_cyy,
		poserr_cxy;			/* pos. error ellipse */
  double	poserr_mx2w, poserr_my2w,
		poserr_mxyw;			/* WORLD error moments */
  float		poserr_aw, poserr_bw,
		poserr_thetaw;			/* WORLD error parameters */
  float		poserr_thetas;			/* native error pos. angle */
  float		poserr_theta2000;		/* J2000 error pos. angle */
  float		poserr_theta1950;		/* B1950 error pos. angle */
  float		poserr_cxxw, poserr_cyyw,
		poserr_cxyw;			/* WORLD error ellipse */
  double	mx2w,my2w,mxyw;			/* WORLD var. and covar. */
  double	peakxf, peakyf;			/* FOCAL of brightest pix */
  double	peakxw, peakyw;			/* WORLD of brightest pix */
  double	mxf, myf;			/* FOCAL barycenters */
  double	mxw, myw;			/* WORLD barycenters */
  double	alphas, deltas;			/* native alpha, delta */
  float		thetas;				/* native position angle E/N*/
  double	peakalphas, peakdeltas;		/* native for brightest pix */
  double	peakalpha2000, peakdelta2000;	/* J2000 for brightest pix */
  double	peakalpha1950, peakdelta1950;	/* B1950 for brightest pix */
  double	alpha2000, delta2000;		/* J2000 alpha, delta */
  float		theta2000;			/* J2000 position angle E/N */
  double	dtheta2000;			/* North J2000 - native angle*/
  double	alpha1950, delta1950;		/* B1950 alpha, delta */
  float		theta1950;			/* B1950 position angle E/N */
  double	dtheta1950;			/* North B1950 - native angle*/
  float		aw, bw;				/* WORLD ellipse size */
  float		thetaw;				/* WORLD position angle */
  float		cxxw,cyyw,cxyw;			/* WORLD ellipse parameters */
  float		npixw, fdnpixw;			/* WORLD isophotal areas */
  float		threshmu;			/* det. surface brightnees */
  float		maxmu;				/* max. surface brightnees */
  float		elong;				/* elongation */
  float		ellip;				/* ellipticity */
  float		polar;				/* Kaiser's "polarization" */
  float		polarw;				/* WORLD "polarization" */
  float		sprob;				/* Stellarity index */
  float		fwhmw;				/* WORLD FWHM */
  double	*assoc;				/* ASSOCiated data */
  int		assoc_number;			/* nb of ASSOCiated objects */
  float		*vignet;			/* Pixel data */
  float		*vigshift;			/* (Shifted) pixel data */

/* Windowed measurements */
  double	winpos_x,winpos_y;		/* Windowed barycenter */
  double	winposerr_mx2, winposerr_my2,
		winposerr_mxy;			/* Error ellips moments */
  float		winposerr_a, winposerr_b,
		winposerr_theta;		/* Error ellips parameters */
  float		winposerr_cxx, winposerr_cyy,
		winposerr_cxy;			/* pos. error ellipse */
  double	winposerr_mx2w, winposerr_my2w,
		winposerr_mxyw;			/* WORLD error moments */
  float		winposerr_aw, winposerr_bw,
		winposerr_thetaw;		/* WORLD error parameters */
  float		winposerr_thetas;		/* native error pos. angle */
  float		winposerr_theta2000;		/* J2000 error pos. angle */
  float		winposerr_theta1950;		/* B1950 error pos. angle */
  float		winposerr_cxxw, winposerr_cyyw,
		winposerr_cxyw;			/* WORLD error ellipse */
  double	win_mx2, win_my2,
		win_mxy;			/* Windowed moments */
  float		win_a, win_b,
		win_theta;			/* Windowed ellipse parameters*/
  float		win_polar;			/* Windowed "polarization" */
  float		win_cxx, win_cyy,
		win_cxy;			/* Windowed ellipse parameters*/
  double	win_mx2w, win_my2w,
		win_mxyw;			/* WORLD windowed moments */
  float		win_aw, win_bw,
		win_thetaw;			/* WORLD ellipse parameters */
  float		win_polarw;			/* WORLD WIN "polarization" */
  float		win_thetas;			/* native error pos. angle */
  float		win_theta2000;			/* J2000 error pos. angle */
  float		win_theta1950;			/* B1950 error pos. angle */
  float		win_cxxw, win_cyyw,
		win_cxyw;			/* WORLD ellipse parameters */
  double	winpos_xf, winpos_yf;		/* FOCAL coordinates */
  double	winpos_xw, winpos_yw;		/* WORLD coordinates */
  double	winpos_alphas, winpos_deltas;	/* native alpha, delta */
  double	winpos_alpha2000, winpos_delta2000;	/* J2000 alpha, delta */
  double	winpos_alpha1950, winpos_delta1950;	/* B1950 alpha, delta */
  short		winpos_niter;			/* Number of WIN iterations */
  short		win_flag;			/* 1:x2<0 2:xy=x2 4:flux<0 */

 /* ---- SOM fitting */
  float		flux_somfit;			/* Fitted amplitude */
  float		fluxerr_somfit;			/* RMS error on SOM flux */
  float		mag_somfit;			/* Magnitude from SOM fit */
  float		magerr_somfit;			/* Mag. err. from SOM fit */
  float		stderr_somfit;			/* Fitting reduced error */
  float		*vector_somfit;			/* SOM fit vector position */
/* ---- Growth curves and stuff */
  float		*flux_growth;			/* Cumulated growth_curve */
  float		flux_growthstep;		/* Growth-curve step */
  float		*mag_growth;			/* Cumulated growth_curve */
  float		mag_growthstep;			/* Growth-curve step */
  float		*flux_radius;			/* f-light-radii */
  float		hl_radius;			/* Scalar half-light radius */
/* ---- PSF-fitting */
  float		fwhm_psf;			/* PSF FWHM */
  float		fwhmw_psf;			/* WORLD PSF FWHM */
  float		flux_psf;			/* Flux from PSF-fitting */
  float		fluxerr_psf;			/* RMS error on PSF flux */
  float		mag_psf;			/* Mag from PSF-fitting */
  float		magerr_psf;			/* RMS mag from PSF-fitting */
  double	x_psf, y_psf;			/* Coords from PSF-fitting */
  double	xf_psf, yf_psf;			/* FOCAL coordinates */
  short		niter_psf;			/* # of PSF-fitting iterat. */
  short		npsf;				/* # of fitted PSFs */
  float		chi2_psf;			/* Red. chi2 of PSF-fitting */
  double	xw_psf, yw_psf;			/* WORLD coords */
  double	alphas_psf, deltas_psf;		/* native alpha, delta */
  double	alpha2000_psf, delta2000_psf;	/* J2000 alpha, delta */
  double	alpha1950_psf, delta1950_psf;	/* B1950 alpha, delta */
  double	poserrmx2_psf, poserrmy2_psf,
		poserrmxy_psf;			/* Error ellips moments */
  float		poserra_psf, poserrb_psf,
		poserrtheta_psf;		/* Error ellips parameters */
  float		poserrcxx_psf, poserrcyy_psf,
		poserrcxy_psf;			/* pos. error ellipse */
  double	poserrmx2w_psf, poserrmy2w_psf,
		poserrmxyw_psf;			/* WORLD error moments */
  float		poserraw_psf, poserrbw_psf,
		poserrthetaw_psf;		/* WORLD error parameters */
  float		poserrthetas_psf;		/* native error pos. angle */
  float		poserrtheta2000_psf;		/* J2000 error pos. angle */
  float		poserrtheta1950_psf;		/* B1950 error pos. angle */
  float		poserrcxxw_psf, poserrcyyw_psf,
		poserrcxyw_psf;			/* WORLD error ellipse */
/* ---- PC-fitting */
  double	mx2_pc,my2_pc,mxy_pc;		/* PC 2nd-order parameters */
  float		a_pc,b_pc,theta_pc;		/* PC shape parameters */
  float		*vector_pc;			/* Principal components */
  float		gdposang;			/* Gal. disk position angle */
  float		gdscale;			/* Gal. disk scalelength */
  float		gdaspect;			/* Gal. disk aspect-ratio */
  float		gde1,gde2;			/* Gal. disk ellipticities */
  float		gbratio;			/* Galaxy B/T */
  float		gbposang;			/* Gal. bulge position angle */
  float		gbscale;			/* Gal. bulge scalelength */
  float		gbaspect;			/* Gal. bulge aspect-ratio */
  float		gbe1,gbe2;			/* Gal. bulge ellipticities */
  float		flux_galfit;			/* Galaxy tot. flux from fit */
  float		fluxerr_galfit;			/* RMS error on galfit flux */
  float		mag_galfit;			/* Galaxy tot. mag from fit */
  float		magerr_galfit;			/* RMS error on galfit mag */
/* ---- Profile-fitting */
  float		*prof_vector;			/* Model parameters */
  float		*prof_errvector;		/* Model parameter errors */
  float		*prof_errmatrix;		/* Model parameter covariances*/
  float		prof_chi2;			/* Reduced chi2 */
  BYTE		prof_flag;			/* Model-fitting flags */
  BYTE		prof_flagw;			/* Model-fitting WORLD flag */
  short		prof_niter;			/* # of model-fitting iter. */
  float		flux_prof;			/* Flux from model-fitting */
  float		fluxerr_prof;			/* RMS error on model flux */
  float		mag_prof;			/* Mag from model-fitting */
  float		magerr_prof;			/* RMS mag from model-fitting */
  float		peak_prof;			/* Model peak flux */
  float		fluxeff_prof;			/* Effective model flux */
  float		fluxmean_prof;			/* Mean effective model flux */
  float		mumax_prof;			/* Model peak surf. bri. */
  float		mueff_prof;			/* Model effective surf. bri. */
  float		mumean_prof;			/* Mean model effective SB */
  float		fluxcor_prof;			/* Hybrid model-fitting flux */
  float		fluxcorerr_prof;		/* RMS error on hybrid flux */
  float		magcor_prof;			/* Hybrid model-fitting mag */
  float		magcorerr_prof;			/* RMS error on hybrid mag */
  double	x_prof, y_prof;			/* Coords from model-fitting*/
  double	xf_prof, yf_prof;		/* FOCAL coordinates */
  double	xw_prof, yw_prof;		/* WORLD coords */
  double	alphas_prof, deltas_prof;	/* native alpha, delta */
  double	alpha2000_prof, delta2000_prof;	/* J2000 alpha, delta */
  double	alpha1950_prof, delta1950_prof;	/* B1950 alpha, delta */
  double	poserrmx2_prof, poserrmy2_prof,
		poserrmxy_prof;			/* Error ellips moments */
  float		poserra_prof, poserrb_prof,
		poserrtheta_prof;		/* Error ellips parameters */
  float		poserrcxx_prof, poserrcyy_prof,
		poserrcxy_prof;			/* pos. error ellipse */
  double	poserrmx2w_prof, poserrmy2w_prof,
		poserrmxyw_prof;		/* WORLD error moments */
  float		poserraw_prof, poserrbw_prof,
		poserrthetaw_prof;		/* WORLD error parameters */
  float		poserrthetas_prof;		/* native error pos. angle */
  float		poserrtheta2000_prof;		/* J2000 error pos. angle */
  float		poserrtheta1950_prof;		/* B1950 error pos. angle */
  float		poserrcxxw_prof, poserrcyyw_prof,
		poserrcxyw_prof;		/* WORLD error ellipse */
  double	prof_mx2, prof_my2, prof_mxy;	/* Model-fitting moments */
  double	prof_mx2cov[9];			/* Mod-fit moment cov. matrix */
  float		prof_a, prof_b,
		prof_theta;			/* Model-fitting ellip. params*/
  float		prof_cxx, prof_cyy,
		prof_cxy;			/* Model-fitting ellip. params*/
  double	prof_convmx2, prof_convmy2,
		prof_convmxy;			/* Convolved model moments */
  float		prof_convcxx, prof_convcyy,
		prof_convcxy;			/* Conv. model ellipse params */
  float		prof_conva, prof_convb,
		prof_convtheta;			/* Conv. model ellipse params */
  float		prof_pol1, prof_pol2;		/* Model-fitting pol. vector*/
  float		prof_pol1err, prof_pol2err,
		prof_pol12corr;			/* Model-fitting pol. errors */
  float		prof_e1, prof_e2;		/* Model-fitting ellip.vector*/
  float		prof_e1err, prof_e2err,
		prof_e12corr;			/* Model-fitting ellip. errors */
  double	prof_mx2w, prof_my2w,
		prof_mxyw;			/* WORLD model-fitting moments*/
  float		prof_aw, prof_bw,
		prof_thetaw;			/* WORLD ellipse parameters */
  float		prof_thetas;			/* native error pos. angle */
  float		prof_theta2000;			/* J2000 error pos. angle */
  float		prof_theta1950;			/* B1950 error pos. angle */
  float		prof_cxxw, prof_cyyw,
		prof_cxyw;			/* WORLD ellipse parameters */
  float		prof_pol1w, prof_pol2w;		/* WORLD polarisation vector*/
  float		prof_pol1errw, prof_pol2errw,
		prof_pol12corrw;		/* WORLD polarisation errors */
  float		prof_e1w, prof_e2w;		/* WORLD ellipticity vector*/
  float		prof_e1errw, prof_e2errw,
		prof_e12corrw;			/* WORLD ellipticity errors */
  float		prof_class_star;		/* Mod.-fitting star/gal class*/
  float		prof_concentration;		/* Model-fitting concentration*/
  float		prof_concentrationerr;		/* RMS error */
  float		prof_noisearea;			/* Equivalent noise area */
  float		prof_offset_flux;		/* Background offset */
  float		prof_offset_fluxerr;		/* RMS error */
  float		prof_dirac_flux;		/* Point source total flux */
  float		prof_dirac_fluxerr;		/* RMS error */
  float		prof_dirac_fluxratio;		/* Point source flux ratio */
  float		prof_dirac_fluxratioerr;	/* RMS error */
  float		prof_dirac_mag;			/* Point source "total" mag */
  float		prof_dirac_magerr;		/* RMS error */
  float		prof_spheroid_flux;		/* Spheroid total flux */
  float		prof_spheroid_fluxerr;		/* RMS error */
  float		prof_spheroid_peak;		/* Spheroid peak flux */
  float		prof_spheroid_fluxeff;		/* Spheroid effective flux */
  float		prof_spheroid_fluxmean;		/* Spheroid mean effect. flux */
  float		prof_spheroid_mag;		/* Spheroid "total" mag */
  float		prof_spheroid_magerr;		/* RMS error */
  float		prof_spheroid_mumax;		/* Spehroid peak surf. brigh.*/
  float		prof_spheroid_mueff;		/* Spheroid effect. surf. bri.*/
  float		prof_spheroid_mumean;		/* Spheroid mean eff. su. bri.*/
  float		prof_spheroid_fluxratio;	/* Spheroid flux ratio */
  float		prof_spheroid_fluxratioerr;	/* RMS error */
  float		prof_spheroid_reff;		/* Spheroid effective radius */
  float		prof_spheroid_refferr;		/* RMS error */
  float		prof_spheroid_reffw;		/* WORLD spheroid eff. radius */
  float		prof_spheroid_refferrw;		/* RMS error */
  float		prof_spheroid_aspect;		/* Spheroid aspect ratio */
  float		prof_spheroid_aspecterr;	/* RMS error */
  float		prof_spheroid_aspectw;		/* WORLD spheroid aspect ratio*/
  float		prof_spheroid_aspecterrw;	/* RMS error */
  float		prof_spheroid_theta;		/* Spheroid position angle */
  float		prof_spheroid_thetaerr;		/* RMS error */
  float		prof_spheroid_thetaw;		/* WORLD spheroid pos. angle */
  float		prof_spheroid_thetaerrw;	/* RMS error */
  float		prof_spheroid_thetas;		/* Sky spheroid pos. angle */
  float		prof_spheroid_theta2000;	/* J2000 spheroid pos. angle */
  float		prof_spheroid_theta1950;	/* B1950 spheroid pos. angle */
  float		prof_spheroid_sersicn;		/* Spheroid Sersic index */
  float		prof_spheroid_sersicnerr;	/* RMS error */
  float		prof_disk_flux;			/* Disk total flux */
  float		prof_disk_fluxerr;		/* RMS error */
  float		prof_disk_peak;			/* Disk peak flux */
  float		prof_disk_fluxeff;		/* Disk effective flux */
  float		prof_disk_fluxmean;		/* Disk mean effective flux */
  float		prof_disk_mag;			/* Disk "total" mag */
  float		prof_disk_magerr;		/* RMS error */
  float		prof_disk_mumax;		/* Disk peak surf. brightness */
  float		prof_disk_mueff;		/* Disk effective surf. bri. */
  float		prof_disk_mumean;		/* Disk mean eff. surf. bri. */
  float		prof_disk_fluxratio;		/* Disk flux ratio */
  float		prof_disk_fluxratioerr;		/* RMS error */
  float		prof_disk_scale;		/* Disk scale length */
  float		prof_disk_scaleerr;		/* RMS error */
  float		prof_disk_scalew;		/* WORLD disk scale length */
  float		prof_disk_scaleerrw;		/* RMS error */
  float		prof_disk_aspect;		/* Disk aspect ratio */
  float		prof_disk_aspecterr;		/* RMS error */
  float		prof_disk_aspectw;		/* WORLD disk aspect ratio */
  float		prof_disk_aspecterrw;		/* RMS error */
  float		prof_disk_inclination;		/* Disk inclination */
  float		prof_disk_inclinationerr;	/* RMS error */
  float		prof_disk_theta;		/* Disk position angle */
  float		prof_disk_thetaerr;		/* RMS error */
  float		prof_disk_thetaw;		/* WORLD disk position angle */
  float		prof_disk_thetaerrw;		/* RMS error */
  float		prof_disk_thetas;		/* Sky disk position angle */
  float		prof_disk_theta2000;		/* J2000 disk position angle */
  float		prof_disk_theta1950;		/* B1950 disk position angle */
  float		*prof_disk_patternvector;	/* Disk pattern coefficients */
  float		*prof_disk_patternmodvector;	/* Disk pattern moduli */
  float		*prof_disk_patternargvector;	/* Disk pattern arguments */
  float		prof_disk_patternspiral;	/* Disk pattern spiral index */
  float		prof_bar_flux;			/* Galactic bar total flux */
  float		prof_bar_fluxerr;		/* RMS error */
  float		prof_bar_mag;			/* Bar "total" magnitude */
  float		prof_bar_magerr;		/* RMS error */
  float		prof_bar_fluxratio;		/* Bar flux ratio */
  float		prof_bar_fluxratioerr;		/* RMS error */
  float		prof_bar_length;		/* Bar length */
  float		prof_bar_lengtherr;		/* RMS error */
  float		prof_bar_lengthw;		/* WORLD bar length */
  float		prof_bar_lengtherrw;		/* RMS error */
  float		prof_bar_aspect;		/* Bar aspect ratio */
  float		prof_bar_aspecterr;		/* RMS error */
  float		prof_bar_aspectw;		/* WORLD bar aspect ratio */
  float		prof_bar_aspecterrw;		/* RMS error */
  float		prof_bar_posang;		/* Bar true prosition angle */
  float		prof_bar_posangerr;		/* RMS error */
  float		prof_bar_theta;			/* Bar projected angle */
  float		prof_bar_thetaerr;		/* RMS error */
  float		prof_bar_thetaw;		/* WORLD bar projected angle */
  float		prof_bar_thetaerrw;		/* RMS error */
  float		prof_bar_thetas;		/* Sky bar projected angle */
  float		prof_bar_theta2000;		/* J2000 bar projected angle */
  float		prof_bar_theta1950;		/* B1950 bar projected angle */
  float		prof_arms_flux;			/* Spiral arms total flux */
  float		prof_arms_fluxerr;		/* RMS error */
  float		prof_arms_mag;			/* Arms "total" magnitude */
  float		prof_arms_magerr;		/* RMS error */
  float		prof_arms_fluxratio;		/* Arms flux ratio */
  float		prof_arms_fluxratioerr;		/* RMS error */
  float		prof_arms_scale;		/* Arms scalelength */
  float		prof_arms_scaleerr;		/* RMS error */
  float		prof_arms_scalew;		/* WORLD arms scalelength */
  float		prof_arms_scaleerrw;		/* RMS error */
  float		prof_arms_posang;		/* Arms true position angle */
  float		prof_arms_posangerr;		/* RMS error */
  float		prof_arms_thetaw;		/* WORLD arms position angle */
  float		prof_arms_thetas;		/* Sky arms position angle */
  float		prof_arms_theta2000;		/* J2000 arms position angle */
  float		prof_arms_theta1950;		/* B1950 arms position angle */
  float		prof_arms_pitch;		/* Arms pitch angle */
  float		prof_arms_pitcherr;		/* RMS error */
  float		prof_arms_start;		/* Arms starting radius */
  float		prof_arms_starterr;		/* RMS error */
  float		prof_arms_startw;		/* WORLD arms starting radius */
  float		prof_arms_starterrw;		/* RMS error */
  float		prof_arms_quadfrac;		/* Arms quadrature fraction */
  float		prof_arms_quadfracerr;		/* RMS error */
  float		dprof_chi2;			/* Det. model fit reduced chi2*/
  BYTE		dprof_flag;			/* Detection model flags*/
  short		dprof_niter;			/* # of detection model iter. */
  float		flux_dprof;			/* Flux in detection model*/
  float		fluxerr_dprof;			/* Error on detect model flux */
  float		mag_dprof;			/* Mag from detection model */
  float		magerr_dprof;			/* RMS mag from detect. model */
/* ---- MEF ----*/
  short		ext_number;			/* FITS extension number */
/* ---- Time ---- */
  float		analtime;			/* Analysis time (s) */
  }	obj2struct;

/*----------------------------- lists of objects ----------------------------*/
typedef struct
  {
  int		nobj;			/* number of objects in list */
  objstruct	*obj;			/* pointer to the object array */
  int		npix;			/* number of pixels in pixel-list */
  pliststruct	*plist;			/* pointer to the pixel-list */
  PIXTYPE	dthresh;		/* detection threshold */
  PIXTYPE	thresh;			/* analysis threshold */
  }	objliststruct;


/*----------------------------- image parameters ----------------------------*/
typedef struct pic
  {
  char		filename[MAXCHAR];	/* pointer to the image filename */
  char		*rfilename;		/* pointer to the reduced image name */
  char		hfilename[MAXCHAR];	/* external header filename */
  int		headflag;		/* external header found? */
  char		ident[MAXCHAR];		/* field identifier (read from FITS)*/
  char		rident[MAXCHAR];	/* field identifier (relative) */
  catstruct	*cat;			/* FITS structure */
  tabstruct	*tab;			/* FITS extension structure */
  FILE		*file;			/* pointer the image file structure */
/* ---- main image parameters */
  int		bitpix, bytepix;	/* nb of bits and bytes per pixel */
  int		bitsgn;			/* non-zero if signed integer data */
  int		width, height;		/* x,y size of the field */
  KINGSIZE_T	npix;			/* total number of pixels */
  double	bscale, bzero;		/* FITS scale and offset */
  double	ngamma;			/* normalized photo gamma */
  int		nlevels;		/* nb of quantification levels */
  float		pixmin, pixmax;		/* min and max values in frame */
  int		y;			/* y current position in field */
  int		ymin;			/* y limit (lowest accessible) */
  int		ymax;			/* y limit (highest accessible+1) */
  int		yblank;			/* y blanking limit (highest+1) */
  PIXTYPE	*strip;			/* pointer to the image buffer */
  FLAGTYPE	*fstrip;		/* pointer to the FLAG buffer */
  int		stripheight;		/* height  of a strip (in lines) */
  int		stripmargin;		/* number of lines in margin */
  int		stripstep;		/* number of lines at each read */
  int		stripy;			/* y position in buffer */
  int		stripylim;		/* y limit in buffer */
  int		stripysclim;		/* y scroll limit in buffer */
/* ---- basic astrometric parameters */
   double	pixscale;		/* pixel size in arcsec.pix-1 */
   double	epoch;			/* epoch of coordinates */
/* ---- basic photometric parameters */
   double	gain;			/* conversion factor in e-/ADU */
   double	satur_level;		/* saturation level in ADUs */
/* ---- background parameters */
  float		*back;			/* ptr to the background map in mem */
  float		*dback;			/* ptr to the background deriv. map */
  float		*sigma;			/* ptr to the sigma map */
  float		*dsigma;		/* Ptr to the sigma deriv. map */
  int		backw, backh;		/* x,y size of a bkgnd mesh */
  int		nbackp;			/* total nb of pixels per bkgnd mesh */
  int		nbackx, nbacky;		/* x,y number of bkgnd meshes */
  int		nback;			/* total number of bkgnd meshes */
  int		nbackfx, nbackfy;	/* x,y size of bkgnd filtering mask */
  float		backmean;		/* median bkgnd value in image */
  float		backsig;		/* median bkgnd rms in image */
  float		sigfac;			/* scaling RMS factor (for WEIGHTs) */
  PIXTYPE	*backline;		/* current interpolated bkgnd line */
  PIXTYPE	dthresh;		/* detection threshold */
  PIXTYPE	thresh;			/* analysis threshold */
  backenum	back_type;		/* Background type */
/* ---- astrometric parameters */
  struct wcs	*wcs;			/* astrometric data */
  struct structassoc	*assoc;		/* ptr to the assoc-list */
  int		flags;			/* flags defining the field type */
/* ---- image interpolation */
  int		interp_flag;		/* interpolation for this field? */
  PIXTYPE	*interp_backup;		/* backup line for interpolation */
  PIXTYPE	weight_thresh;		/* interpolation threshold */
  int		*interp_ytimeoutbuf;	/* interpolation timeout line buffer */
  int		interp_xtimeout;	/* interpolation timeout value in x */
  int		interp_ytimeout;	/* interpolation timeout value in y */
  struct pic	*reffield;	       	/* pointer to a reference field */
  OFF_T		mefpos;			/* Position in a MEF file */
  }	picstruct;


/*-------------------------------- catalog  ---------------------------------*/

typedef struct
  {
  int		ndetect;				/* nb of detections */
  int		ntotal;					/* Total object nb */
  int		nparam;					/* Nb of parameters */
/*----- Misc. strings defining the extraction */
  char		prefs_name[MAXCHAR];			/* Prefs filename*/
  char		image_name[MAXCHAR];			/* image filename*/
  char		psf_name[MAXCHAR];			/* PSF filename*/
  char		nnw_name[MAXCHAR];			/* NNW name */
  char		filter_name[MAXCHAR];			/* Filter name */
  char		soft_name[MAXCHAR];			/* Sextractor version*/
/*----- time */
  char		ext_date[16],ext_time[16];		/* date and time */
  double	ext_elapsed;				/* processing time */
/*----- MEF */
  int		currext;				/* current extension */
  int		next;					/* Nb of extensions */
  }		sexcatstruct;

