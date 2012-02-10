/*
*				prefs.h
*
* Include file for prefs.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-20112 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		18/01/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _PROFIT_H_
#include        "profit.h"
#endif

#ifndef _PATTERN_H_
#include        "pattern.h"
#endif

#ifndef _PREFS_H_
#define _PREFS_H_

/*----------------------------- Internal constants --------------------------*/

#define	MAXCHARL	16384		/* max. nb of chars in a string list */
#define	MAXLIST		2048		/* max. nb of list members */
#define	MAXLISTSIZE	2000000		/* max size of list */
#define	MAXNPARAMAXIS	4		/* max parameter dimensionality */

/* NOTES:
One must have:	MAXLIST >= 1 (preferably >= 16!)
*/

/*------------------------------- preferences -------------------------------*/
typedef struct
  {
  char		**command_line;				/* Command line */
  int		ncommand_line;				/* nb of params */
  char		prefs_name[MAXCHAR];			/* prefs filename*/
  char		*(image_name[MAXIMAGE]);		/* image filenames */
  int		nimage;					/* nb of images */
  char		cat_name[MAXCHAR];			/* catalog filename*/
  char		head_suffix[MAXCHAR];			/* ext. header suffix */
/*----- catalog output */
  char		*(param[MAXLIST]);			/* catalog parameters*/
  int		nparam;					/* nb of params */
  enum	{CAT_NONE, ASCII, ASCII_HEAD, ASCII_SKYCAT, ASCII_VO,
	FITS_LDAC, FITS_TPX, FITS_10}
		cat_type;				/* type of catalog */
/*----- thresholding */
  double	dthresh[2];				/* detect. threshold */
  int		ndthresh;				/* (1 or 2 entries) */
  double	thresh[2];				/* analysis thresh. */
  int		nthresh;				/* (1 or 2 entries) */
  enum	{THRESH_RELATIVE, THRESH_ABSOLUTE}
					thresh_type[2];	/* bkgnd type */
  int		nthresh_type;				/* nb of params */
/*----- extraction */
  int		ext_minarea;				/* min area in pix. */
  int		ext_maxarea;				/* max area in pix. */
  int		deb_maxarea;				/* max deblend. area */
  int		filter_flag;				/* smoothing on/off */
  char		filter_name[MAXCHAR];			/* mask filename */
  double	filter_thresh[2];			/* Filter thresholds */
  int		nfilter_thresh;				/* nb of params */
  int		deblend_nthresh;			/* threshold number */
  double	deblend_mincont;			/* minimum contrast */
  char		satur_key[8];				/* saturation keyword */
  double	satur_level;				/* saturation level */
  enum	{CCD, PHOTO}			detect_type;	/* detection type */
/*----- Flagging */
  char		*(fimage_name[MAXFLAG]);		/* flagmap filenames */
  int		nfimage_name;				/* nb of params */
  enum	{FLAG_OR, FLAG_AND, FLAG_MIN, FLAG_MAX, FLAG_MOST}
				flag_type[MAXFLAG];	/* flag combination */
  int		imaflag_size;				/* requested iso nb1 */
  int		imanflag_size;				/* requested iso nb2 */
  int		nimaisoflag;				/* effective iso nb */
  int		nimaflag;				/* effective ima nb */
/*----- cleaning */
  int		clean_flag;				/* allow cleaning ? */
  double	clean_param;				/* cleaning effic. */
/*----- Weighting */
  char		*(wimage_name[2]);       		/* weight filenames */
  int		nwimage_name;				/* nb of params */
  weightenum	weight_type[MAXIMAGE];			/* weighting scheme */
  int		nweight_type;				/* nb of params */
  int		weight_flag[MAXIMAGE];			/* do we weight ? */
  int		weightgain_flag;			/* weight gain? */
  int		wscale_flag[2];		/* Weight rescaling */
  int		nwscale_flag;				/* nb of params */
/*----- photometry */
  int		nphotinstru;				/* nb of channels */
  char		*(photinstru_key[72]);		/* Photom instrument keywords*/
  int		nphotinstru_key;			/* nb of params */
  char		**photinstrustr;		/* Photom instrument labels */
  int		nphotinstrumax;				/* max nb of params */
  double	apert[MAXNAPER];			/* apert size (pix) */
  int		naper;					/* effective apert. */
  int		aper_size[2];				/* requested apert. */
  double	autoparam[2];				/* Kron parameters */
  int		nautoparam;				/* nb of Kron params */
  double	petroparam[2];				/* Kron parameters */
  int		npetroparam;				/* nb of Kron params */
  double	autoaper[2];				/* minimum apertures */
  int		nautoaper;				/* nb of min. aperts */
  double	mag_zeropoint;				/* magnitude offsets */
  double	mag_gamma;				/* for emulsions */
  double	gain;					/* only for CCD */
  char		gain_key[8];				/* gain keyword */
/*----- S/G separation */
  double	pixel_scale;				/* in arcsec */
  double	seeing_fwhm;				/* in arcsec */
  char		nnw_name[MAXCHAR];			/* nnw filename */
/*----- background */
  enum	{IMAGE, AFILE}			back_origin;	/* origin of bkgnd */
  char		back_name[MAXCHAR];			/* bkgnd filename */
  backenum	back_type[2];				/* bkgnd type */
  int		nback_type;				/* nb of params */
  double	back_val[2];				/* user-def. bkg */
  int		nback_val;				/* nb of params */
  int		backsize[2];				/* bkgnd mesh size */
  int		nbacksize;				/* nb of params */
  int		backfsize[2];				/* bkgnd filt. size */
  int		nbackfsize;				/* nb of params */
  double	backfthresh;				/* bkgnd fil. thresh */
  enum	{GLOBAL, LOCAL}			pback_type;	/* phot. bkgnd type */
  int		pback_size;				/* rect. ann. width */
/*----- memory */
  int		clean_stacksize;			/* obj buffer size*/
  int		obj2_stacksize;				/* obj2 buffer size*/
  int		mem_pixstack;				/* pixel stack size */
  int		mem_bufsize;				/* strip height */
/*----- miscellaneous */
  int		pipe_flag;				/* allow piping ? */
  enum	{QUIET, NORM, WARN, FULL}      	verbose_type;	/* display type */
  int		xml_flag;				/* Write XML file? */
  char		xml_name[MAXCHAR];			/* XML file name */
  char		xsl_name[MAXCHAR];			/* XSL file name (or URL) */
  char		sdate_start[12];			/* SCAMP start date */
  char		stime_start[12];			/* SCAMP start time */
  char		sdate_end[12];				/* SCAMP end date */
  char		stime_end[12];				/* SCAMP end time */
  double	time_diff;				/* Execution time */
/*----- CHECK-images */
  int		check_flag;				/* CHECK-image flag */
  checkenum    	check_type[MAXCHECK];		       	/* check-image types */
  int		ncheck_type;				/* nb of params */
  char		*(check_name[MAXCHECK]);       		/* check-image names */
  int		ncheck_name;				/* nb of params */
  struct structcheck	*(check[MAXCHECK]);		/* check-image ptrs */
/*----- ASSOCiation */
  int		assoc_flag;				/* ASSOCiation flag */
  char		assoc_name[MAXCHAR];			/* ASSOC-list name */
  int		assoc_param[3];				/* ASSOC param cols */
  int		nassoc_param;				/* nb of params */
  int		assoc_data[MAXNASSOC];		       	/* ASSOC data cols */
  int		nassoc_data;				/* nb of params */
  enum	{ASSOC_FIRST, ASSOC_NEAREST, ASSOC_MEAN, ASSOC_MAGMEAN,
	 ASSOC_SUM, ASSOC_MAGSUM, ASSOC_MIN, ASSOC_MAX}
		assoc_type;				/* type of assoc. */
  enum	{ASSOCCOORD_PIXEL, ASSOCCOORD_WORLD}
		assoccoord_type;		       	/* type of coords */
  enum	{ASSOCSELEC_ALL, ASSOCSELEC_MATCHED, ASSOCSELEC_NOMATCHED}
		assocselec_type;		       	/* type of assoc. */
  double	assoc_radius;				/* ASSOC range */
  int		assoc_size;				/* nb of parameters */
  char		retina_name[MAXCHAR];			/* retina filename */
  int		vignet_size[3];				/* vignet size */
  int		vigshift_size[3];			/* shift-vignet size */
  int		cleanmargin;				/* CLEANing margin */
/*----- SOM-fitting */
  char		som_name[MAXCHAR];			/* SOM filename */
  int		somfit_flag;				/* ASSOCiation flag */
  int		somfit_vector_size[2];			/* SOMfit vec. size */
/*----- masking */
  enum {MASK_NONE, MASK_BLANK, MASK_CORRECT}
		mask_type;				/* type of masking */
  int		blank_flag;				/* BLANKing flag */
  double	weight_thresh[MAXIMAGE];		/* weight thresholds */
  int		nweight_thresh;				/* nb of params */
/*----- interpolation */
  enum	{INTERP_NONE, INTERP_VARONLY, INTERP_ALL}
		interp_type[MAXIMAGE];			/* interpolat. type */
  int		ninterp_type;				/* nb of params */
  int		interp_xtimeout[MAXIMAGE]; 		/* interp. x timeout */
  int		ninterp_xtimeout;		        /* nb of params */
  int		interp_ytimeout[MAXIMAGE];		/* interp. y timeout */
  int		ninterp_ytimeout;       		/* nb of params */
/*----- astrometry */
  int		world_flag;				/* WORLD required */
  char		coosys[16];				/* VOTable coord.sys */
  double	epoch;					/* VOTable epoch */
/*----- growth curve */
  int		growth_flag;				/* gr. curve needed */
  int		flux_growth_size[2];			/* number of elem. */
  int		mag_growth_size[2];			/* number of elem. */
  int		flux_radius_size[2];			/* number of elem. */
  double	growth_step;				/* step size (pix) */
  double	flux_frac[MAXNAPER];			/* for FLUX_RADIUS */
  int		nflux_frac;       			/* number of elem. */
/*----- PSF-fitting */
  int		psf_flag;				/* PSF needed */
  int		dpsf_flag;				/* detectiob PSF */
  int		psffit_flag;				/* PSF-fit needed */
  int		dpsffit_flag;				/* dual image PSF-fit */
  char		*(psf_name[2]);				/* PSF filename */
  int		npsf_name;				/* nb of params */
  int		psf_npsfmax;				/* Max # of PSFs */
  int		psf_xsize,psf_ysize;			/* nb of params */
  int		psf_xwsize,psf_ywsize;			/* nb of params */
  int		psf_alphassize,psf_deltassize;		/* nb of params */
  int		psf_alpha2000size,psf_delta2000size;	/* nb of params */
  int		psf_alpha1950size,psf_delta1950size;	/* nb of params */
  int		psf_fluxsize;				/* nb of params */
  int		psf_fluxerrsize;			/* nb of params */
  int		psf_magsize;				/* nb of params */
  int		psf_magerrsize;				/* nb of params */
  int		prof_flag;				/* Profile-fitting */
  int		pattern_flag;				/* Pattern-fitting */
/*----- model-fitting */
  int		prof_modelflags;			/* model flags */
  int		prof_vectorsize;			/* nb of params */
  int		prof_errvectorsize;			/* nb of params */
  int		prof_errmatrixsize[2];			/* nb of params */
  int		prof_disk_patternvectorsize;		/* nb of params */
  int		prof_disk_patternncomp;			/* nb of params */
  int		prof_disk_patternmodvectorsize;		/* nb of params */
  int		prof_disk_patternmodncomp;		/* nb of params */
  int		prof_disk_patternargvectorsize;		/* nb of params */
  int		prof_disk_patternargncomp;		/* nb of params */
/*----- Pattern-fitting */
  pattypenum	pattern_type;				/* Disk pattern type */
/*----- customize */
  int		fitsunsigned_flag;			/* Force unsign FITS */
  int		next;			     /* Number of extensions in file */
/* Multithreading */
  int		nthreads;			/* Number of active threads */
  }	prefstruct;

prefstruct	prefs;

/*-------------------------------- protos -----------------------------------*/
extern char	*list_to_str(char *listname);

extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	prefs_dump(int state),
		prefs_end(void),
		prefs_tune(void),
		prefs_read(char *filename,char **argkey,char **argval,int narg),
		prefs_use(void);
#endif
