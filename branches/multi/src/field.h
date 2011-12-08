/*
*				field.h
*
* Include file for field.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		08/12/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FIELD_H_
#define _FIELD_H_

/*------------------------------ field flags -------------------------------*/
#define		DETECT_FIELD	0x01	/* Detection */
#define		MEASURE_FIELD	0x02	/* Measurement */
#define		FLAG_FIELD	0x04	/* Flagging */
#define		RMS_FIELD	0x08	/* Weighting with std deviations */
#define		VAR_FIELD	0x10	/* Weighting with variances */
#define		WEIGHT_FIELD	0x20	/* Weighting with weights */
#define		BACKRMS_FIELD	0x40	/* Weighting from a backrms matrix */
#define		INTERP_FIELD	0x80	/* Purely interpolated data */


/*------------------------------- structures --------------------------------*/
/* image parameters */
typedef struct field
  {
  struct field	*reffield;	       	/* pointer to a reference field */
  char		filename[MAXCHAR];	/* pointer to the image filename */
  char		*rfilename;		/* pointer to the reduced image name */
  char		hfilename[MAXCHAR];	/* external header filename */
  int		headflag;		/* external header found? */
  char		ident[MAXCHAR];		/* field identifier (read from FITS)*/
  char		rident[MAXCHAR];	/* field identifier (relative) */
  catstruct	*cat;			/* FITS structure */
  tabstruct	*tab;			/* FITS extension structure */
  int		bitpix;			/* FITS image bitpix (for metadata) */
/* ---- main image parameters */
  int		width, height;		/* x,y size of the field */
  KINGSIZE_T	npix;			/* total number of pixels */
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
/* ---- PSF */
  struct psf	*psf;			/* image PSF model */
  }	fieldstruct;


/*------------------------------- functions ---------------------------------*/

void		field_end(fieldstruct *field);

fieldstruct	*field_inherit(fieldstruct *infield, int flags),
		*field_init(char *filename, int flags, int ext);

#endif
