/*
*				back.c
*
* Functions dealing with the image background.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		30/06/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

// TODO: put this block in new file
#include <stdbool.h>
#include <sys/mman.h>
#include <time.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"

//#include        "objmask.h"

#include	"back.h"
#include	"field.h"
#include	"misc.h"
#include	"weight.h"

// REMARK: EITHER include 'readimage.h' OR define 'readimage_loadstrip()',
//         otherwise the return value is NOT a pointer value, and the code crashes;
//#include "readimage.h"
void *readimage_loadstrip(fieldstruct *field, fieldstruct *wfield, const int write_checkimgs);

/****************************back_map*****************************************/
/**
 * Function: back_map
 * Computes a map of a field (image extension) and stores it in the field
 * structure.
 *
 * @author E. Bertin (IAP), MK
 * @date   30th June 2014
 *
 * @param[in,out] field       the field
 * @param[in]     wfield      the weight field
 * @param[in]     wscale_flag weight scaling flag
 */
void    back_map(fieldstruct *field, fieldstruct *wfield, int wscale_flag)
{
  // set the object mask to NULL and call the masked function
  objmaskstruct *nullmask=NULL;
  back_map_mask(field, wfield, wscale_flag, nullmask);
}

/****************************back_map_mask************************************/
/**
 * Function: back_map_mask
 * Computes a map of a field (image extension) and stores it in the field
 * structure. Pixel which are as part of an object masked in the object masked
 * are *NOT* for the background and noise determination.
 *
 * @author E. Bertin (IAP), MK
 * @date   30th June 2014
 *
 * @param[in,out] field       the field
 * @param[in]     wfield      the weight field
 * @param[in]     wscale_flag weight scaling flag
 * @param[in]     omask       the mask for the field
 */
void    back_map_mask(fieldstruct *field, fieldstruct *wfield, int wscale_flag, objmaskstruct *omask)

    {
   backstruct	*backmesh,*wbackmesh, *bm,*wbm;
   PIXTYPE	*buf,*wbuf, *buft,*wbuft;
   OFF_T	fcurpos,wfcurpos, wfcurpos2,fcurpos2, bufshift, jumpsize;
   long int     ocurpos, ocurpos2;
   size_t	bufsize, bufsize2,
		size,meshsize;
   int		i,j,k,m,n, step, nlines,
		w,bw, bh, nx,ny,nb,
		lflag, nr;
   float	*ratio,*ratiop, *weight, *sigma,
		sratio, sigfac;

/* If the weight-map is not an external one, no stats are needed for it */
  if (wfield && wfield->flags&(INTERP_FIELD|BACKRMS_FIELD))
    wfield= NULL;

  w = field->width;
  bw = field->backw;
  bh = field->backh;
  nx = field->nbackx;
  ny = field->nbacky;
  nb = field->nback;

  NFPRINTF(OUTPUT, "Setting up background maps");

/* Decide if it is worth displaying progress each 16 lines */

  lflag = (field->width*field->backh >= (size_t)65536);

/* Save current positions in files */

  wfcurpos = wfcurpos2 = 0;		/* to avoid gcc -Wall warnings */
  QFTELL(field->cat->file, fcurpos, field->cat->filename);
  ocurpos = ocurpos2 = 0;
  if (omask)
    ocurpos = tell_objmask(omask);
  if (wfield)
    QFTELL(wfield->cat->file, wfcurpos, wfield->cat->filename);

/* Allocate a correct amount of memory to store pixels */

  bufsize = (OFF_T)w*bh;
  meshsize = (size_t)bufsize;
  nlines = 0;
  if (bufsize > (size_t)BACK_BUFSIZE)
    {
    nlines = BACK_BUFSIZE/w;
    step = (field->backh-1)/nlines+1;
    bufsize = (size_t)(nlines = field->backh/step)*w;
    bufshift = (step/2)*(OFF_T)w;
    jumpsize = (step-1)*(OFF_T)w;
    }
  else
    bufshift = jumpsize = 0;		/* to avoid gcc -Wall warnings */

  /* Allocate some memory */
  QMALLOC(backmesh, backstruct, nx);		/* background information */
  QMALLOC(buf, PIXTYPE, bufsize);		/* pixel buffer */
  free(field->back);
  QMALLOC(field->back, float, nb);		/* background map */
  free(field->backline);
  QMALLOC(field->backline, PIXTYPE, w);		/* current background line */
  free(field->sigma);
  QMALLOC(field->sigma, float, nb);		/* sigma map */
  if (wfield)
    {
    QMALLOC(wbackmesh, backstruct, nx);		/* background information */
    QMALLOC(wbuf, PIXTYPE, bufsize);		/* pixel buffer */
    free(wfield->back);
    QMALLOC(wfield->back, float, nb);		/* background map */
    free(wfield->backline);
    QMALLOC(wfield->backline, PIXTYPE, w);	/* current background line */
    free(wfield->sigma);
    QMALLOC(wfield->sigma, float, nb);		/* sigma map */
    wfield->sigfac = 1.0;
    }
  else
    {
    wbackmesh = NULL;
    wbuf = NULL;
    }

  // =1 prints to file,
  // =0 prints to screen
  int fileindex=1;

/* Loop over the data packets */

  for (j=0; j<ny; j++)
    {
    if (lflag && j)
      NPRINTF(OUTPUT, "\33[1M> Setting up background map at line:%5d\n\33[1A",
	      j*bh);
    if (!nlines)
      {
/*---- The image is small enough so that we can make exhaustive stats */
      if (j == ny-1 && field->npix%bufsize)
        bufsize = field->npix%bufsize;
      read_body(field->tab, buf, bufsize);
      if (omask)
        apply_objmask(buf, bufsize, omask);
      if (wfield)
        {
        read_body(wfield->tab, wbuf, bufsize);
        weight_to_var(wfield, wbuf, bufsize);
        }
/*---- Build the histograms */
      back_stat(backmesh, wbackmesh, buf, wbuf, bufsize,nx, w, bw,
	wfield?wfield->weight_thresh:0.0);
      bm = backmesh;
      for (m=nx; m--; bm++)
        if (bm->mean <= -BIG)
          bm->histo=NULL;
        else
          QCALLOC(bm->histo, LONG, bm->nlevels);
      if (wfield)
        {
        wbm = wbackmesh;
        for (m=nx; m--; wbm++)
          if (wbm->mean <= -BIG)
            wbm->histo=NULL;
          else
            QCALLOC(wbm->histo, LONG, wbm->nlevels);
        }
      back_histo(backmesh, wbackmesh, buf, wbuf, bufsize,nx, w, bw,
	wfield?wfield->weight_thresh:0.0);
      }
    else
      {
/*---- Image size too big, we have to skip a few data !*/
      QFTELL(field->cat->file, fcurpos2, field->cat->filename);
      if (omask)
        ocurpos2 = tell_objmask(omask);
      if (wfield)
        QFTELL(wfield->cat->file, wfcurpos2, wfield->cat->filename);
      if (j == ny-1 && (n=field->height%field->backh))
        {
        meshsize = n*(size_t)w;
        nlines = BACK_BUFSIZE/w;
        step = (n-1)/nlines+1;
        bufsize = (nlines = n/step)*(size_t)w;
        bufshift = (step/2)*(OFF_T)w;
        jumpsize = (step-1)*(OFF_T)w;
        free(buf);
        QMALLOC(buf, PIXTYPE, bufsize);		/* pixel buffer */
        if (wfield)
          {
          free(wbuf);
          QMALLOC(wbuf, PIXTYPE, bufsize);	/* pixel buffer */
          }
        }

/*---- Read and skip, read and skip, etc... */
      QFSEEK(field->cat->file, bufshift*(OFF_T)field->tab->bytepix,
		SEEK_CUR, field->cat->filename);
      if (omask)
        seek_objmask(bufshift, SEEK_CUR, omask);

      buft = buf;
      for (i=nlines; i--; buft += w)
        {
        read_body(field->tab, buft, w);
        if (omask)
          apply_objmask(buft, w, omask);
        if (i)
          {
            QFSEEK(field->cat->file, jumpsize*(OFF_T)field->tab->bytepix,
                SEEK_CUR, field->cat->filename);
            if (omask)
              seek_objmask(jumpsize, SEEK_CUR, omask);
          }
        }

      if (wfield)
        {
/*------ Read and skip, read and skip, etc... now on the weight-map */
        QFSEEK(wfield->cat->file, bufshift*(OFF_T)wfield->tab->bytepix,
		SEEK_CUR, wfield->cat->filename);
        wbuft = wbuf;
        for (i=nlines; i--; wbuft += w)
          {
          read_body(wfield->tab, wbuft, w);
          weight_to_var(wfield, wbuft, w);
          if (i)
            QFSEEK(wfield->cat->file, jumpsize*(OFF_T)wfield->tab->bytepix,
		SEEK_CUR, wfield->cat->filename);
          }
        }
      back_stat(backmesh, wbackmesh, buf, wbuf, bufsize, nx, w, bw,
	wfield?wfield->weight_thresh:0.0);
      QFSEEK(field->cat->file, fcurpos2, SEEK_SET, field->cat->filename);
      if (omask)
        seek_objmask(ocurpos2, SEEK_SET, omask);
      bm = backmesh;
      for (m=nx; m--; bm++)
        if (bm->mean <= -BIG)
          bm->histo=NULL;
        else
          QCALLOC(bm->histo, LONG, bm->nlevels);
      if (wfield)
        {
        QFSEEK(wfield->cat->file, wfcurpos2, SEEK_SET, wfield->cat->filename);
        wbm = wbackmesh;
        for (m=nx; m--; wbm++)
          if (wbm->mean <= -BIG)
            wbm->histo=NULL;
          else
            QCALLOC(wbm->histo, LONG, wbm->nlevels);
        }
/*---- Build (progressively this time) the histograms */
      for(size=meshsize, bufsize2=bufsize; size>0; size -= bufsize2)
        {
        if (bufsize2>size)
          bufsize2 = size;
        read_body(field->tab, buf, bufsize2);
        if (omask)
          apply_objmask(buf, bufsize2, omask);
        if (wfield)
          {
          read_body(wfield->tab, wbuf, bufsize2);
          weight_to_var(wfield, wbuf, bufsize2);
          }
        back_histo(backmesh, wbackmesh, buf, wbuf, bufsize2, nx, w, bw,
		wfield?wfield->weight_thresh:0.0);
        }
      }

    // print the background mesh
    //back_printmeshs(backmesh, nx,&fileindex);

    /*-- Compute background statistics from the histograms */
    bm = backmesh;
    for (m=0; m<nx; m++, bm++)
      {
      k = m+nx*j;
      back_guess(bm, field->back+k, field->sigma+k);
      free(bm->histo);
      }
    if (wfield)
      {
      wbm = wbackmesh;
      for (m=0; m<nx; m++, wbm++)
        {
        k = m+nx*j;
        back_guess(wbm, wfield->back+k, wfield->sigma+k);
        free(wbm->histo);
        }
      }
    }

/* Free memory */
  free(buf);
  free(backmesh);
  if (wfield)
    {
    free(wbackmesh);
    free(wbuf);
    }

/* Go back to the original position */
  QFSEEK(field->cat->file, fcurpos, SEEK_SET, field->cat->filename);
  if (omask){
    seek_objmask(ocurpos, SEEK_SET, omask);
  }
 if (wfield)
    QFSEEK(wfield->cat->file, wfcurpos, SEEK_SET, wfield->cat->filename);

/* Median-filter and check suitability of the background map */
  NFPRINTF(OUTPUT, "Filtering background map(s)");
  back_filter(field);
  if (wfield)
    back_filter(wfield);

/* Compute normalization for variance- or weight-maps*/
  if (wfield && wfield->flags&(VAR_FIELD|WEIGHT_FIELD))
    {
    nr = 0;
    QMALLOC(ratio, float, wfield->nback);
    ratiop = ratio;
    weight = wfield->back;
    sigma = field->sigma;
    for (i=wfield->nback; i--; sigma++)
      if ((sratio=*(weight++)) > 0.0
		&& (sratio = *sigma/sqrt(sratio)) > 0.0)
        {
        *(ratiop++) = sratio;
        nr++;
        }
    sigfac = fqmedian(ratio, nr);
    for (i=0; i<nr && ratio[i]<=0.0; i++);
    if (i<nr)
      sigfac = fqmedian(ratio+i, nr-i);
    else
      {
      warning("Null or negative global weighting factor:","defaulted to 1");
      sigfac = 1.0;
      } 
    free(ratio);

    if (wscale_flag)
      wfield->sigfac = sigfac;
    else
      {
      wfield->sigfac = 1.0;
      field->backsig /= sigfac;
      }
    }


/* Compute 2nd derivatives along the y-direction */
  NFPRINTF(OUTPUT, "Computing background d-map");
  free(field->dback);
  field->dback = back_makespline(field, field->back);
  NFPRINTF(OUTPUT, "Computing background-noise d-map");
  free(field->dsigma);
  field->dsigma = back_makespline(field, field->sigma);
/* If asked for, force the backmean parameter to the supplied value */
  if (field->back_type == BACK_ABSOLUTE)
    field->backmean = (float)prefs.back_val[(field->flags&DETECT_FIELD)?0:1];

/* Set detection/measurement threshold */
  if (prefs.ndthresh > 1)
    {
     double	dval;

    if (fabs(dval=prefs.dthresh[0] - prefs.dthresh[1])> 70.0)
      error(EXIT_FAILURE,
	"*Error*: I cannot deal with such extreme thresholds!", "");

    field->dthresh = field->pixscale*field->pixscale*pow(10.0, -0.4*dval);
    }
  else if (prefs.thresh_type[0]==THRESH_ABSOLUTE)
    field->dthresh = prefs.dthresh[0];
  else
    field->dthresh = prefs.dthresh[0]*field->backsig;
  if (prefs.nthresh > 1)
    {
     double	dval;

    if (fabs(dval=prefs.thresh[0] - prefs.thresh[1]) > 70.0)
      error(EXIT_FAILURE,
	"*Error*: I cannot deal with such extreme thresholds!", "");

    field->thresh = field->pixscale*field->pixscale*pow(10.0, -0.4*dval);
    }
  else if (prefs.thresh_type[1]==THRESH_ABSOLUTE)
    field->thresh = prefs.thresh[0];
  else
    field->thresh = prefs.thresh[0]*field->backsig;

#ifdef	QUALITY_CHECK
  printf("%-10g %-10g %-10g\n", field->backmean, field->backsig,
	(field->flags & DETECT_FIELD)? field->dthresh : field->thresh);
#endif
  if (field->dthresh<=0.0 || field->thresh<=0.0)
    error(EXIT_FAILURE,
	"*Error*: I cannot deal with zero or negative thresholds!", "");

  if (field->detector_type == DETECTOR_PHOTO
	&& field->backmean+3*field->backsig > 50*field->ngamma)
    error(EXIT_FAILURE,
	"*Error*: The density range of this image is too large for ",
	"PHOTO mode");

  return;
  }


/****** back_stat ***********************************************************
PROTO	void back_stat(backstruct *backmesh, backstruct *wbackmesh,
		PIXTYPE *buf, PIXTYPE *wbuf, size_t bufsize,
			int n, int w, int bw, PIXTYPE wthresh)
PURPOSE	Compute robust statistical estimators of the background in a row of
	meshes.
INPUT	Pointer to background image mesh structure,
	pointer to background variance mesh structure,
	pointer to image buffer,
	pointer to variance buffer,
	buffer size (number of pixels),
	number of meshes
	frame width (pixels),
	mesh width (pixels),
	weight threshold.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2011
 ***/
void	back_stat(backstruct *backmesh, backstruct *wbackmesh,
		PIXTYPE *buf, PIXTYPE *wbuf, size_t bufsize,
			int n, int w, int bw, PIXTYPE wthresh)

  {
   backstruct	*bm, *wbm;
   double	pix,wpix, sig, mean,wmean, sigma,wsigma, step;
   PIXTYPE	*buft,*wbuft,
		lcut,wlcut, hcut,whcut;
   int		m,h,x,y, npix,wnpix, offset, lastbite;

  h = bufsize/w;
  bm = backmesh;
  wbm = wbackmesh;
  offset = w - bw;
  step = sqrt(2/PI)*QUANTIF_NSIGMA/QUANTIF_AMIN;
  wmean = wsigma = wlcut = whcut = 0.0;	/* to avoid gcc -Wall warnings */
  for (m = n; m--; bm++,buf+=bw)
    {
    if (!m && (lastbite=w%bw))
      {
      bw = lastbite;
      offset = w-bw;
      }
    mean = sigma = 0.0;
    buft=buf;
/*-- We separate the weighted case at this level to avoid penalty in CPU */
    npix = 0;
    if (wbackmesh)
      {
      wmean = wsigma = 0.0;
      wbuft = wbuf;
      for (y=h; y--; buft+=offset,wbuft+=offset)
        for (x=bw; x--;)
          {
          pix = *(buft++);
          if ((wpix = *(wbuft++)) < wthresh && pix > -BIG)
            {
            wmean += wpix;
            wsigma += wpix*wpix;
            mean += pix;
            sigma += pix*pix;
            npix++;
            }
	  }
      }
    else
      for (y=h; y--; buft+=offset)
        for (x=bw; x--;)
          if ((pix = *(buft++)) > -BIG)
            {
            mean += pix;
            sigma += pix*pix;
            npix++;
            }
/*-- If not enough valid pixels, discard this mesh */
    if ((float)npix < (float)(bw*h*BACK_MINGOODFRAC))
      {
      bm->mean = bm->sigma = -BIG;
      if (wbackmesh)
        {
        wbm->mean = wbm->sigma = -BIG;
        wbm++;
        wbuf += bw;
        }
      continue;
      }
    if (wbackmesh)
      {
      wmean /= (double)npix;
      wsigma = (sig = wsigma/npix - wmean*wmean)>0.0? sqrt(sig):0.0;
      wlcut = wbm->lcut = (PIXTYPE)(wmean - 2.0*wsigma);
      whcut = wbm->hcut = (PIXTYPE)(wmean + 2.0*wsigma);
      }
    mean /= (double)npix;
    sigma = (sig = sigma/npix - mean*mean)>0.0? sqrt(sig):0.0;
    lcut = bm->lcut = (PIXTYPE)(mean - 2.0*sigma);
    hcut = bm->hcut = (PIXTYPE)(mean + 2.0*sigma);
    mean = sigma = 0.0;
    npix = wnpix = 0;
    buft = buf;
    if (wbackmesh)
      {
      wmean = wsigma = 0.0;
      wbuft=wbuf;
      for (y=h; y--; buft+=offset, wbuft+=offset)
        for (x=bw; x--;)
          {
          pix = *(buft++);
          if ((wpix = *(wbuft++))<wthresh && pix<=hcut && pix>=lcut)
            {
            mean += pix;
            sigma += pix*pix;
            npix++;
            if (wpix<=whcut && wpix>=wlcut)
              {
              wmean += wpix;
              wsigma += wpix*wpix;
              wnpix++;
              }
            }
          }
      }
    else
      for (y=h; y--; buft+=offset)
        for (x=bw; x--;)
          {
          pix = *(buft++);
          if (pix<=hcut && pix>=lcut)
            {
            mean += pix;
            sigma += pix*pix;
            npix++;
            }
          }

    bm->npix = npix;
    mean /= (double)npix;
    sig = sigma/npix - mean*mean;
    sigma = sig>0.0 ? sqrt(sig):0.0;
    bm->mean = mean;
    bm->sigma = sigma;
    if ((bm->nlevels = (int)(step*npix+1)) > QUANTIF_NMAXLEVELS)
      bm->nlevels = QUANTIF_NMAXLEVELS;
    bm->qscale = sigma>0.0? 2*QUANTIF_NSIGMA*sigma/bm->nlevels : 1.0;
    bm->qzero = mean - QUANTIF_NSIGMA*sigma;
    if (wbackmesh)
      {
      wbm->npix = wnpix;
      wmean /= (double)wnpix;
      sig = wsigma/wnpix - wmean*wmean;
      wsigma = sig>0.0 ? sqrt(sig):0.0;
      wbm->mean = wmean;
      wbm->sigma = wsigma;
      if ((wbm->nlevels = (int)(step*wnpix+1)) > QUANTIF_NMAXLEVELS)
        wbm->nlevels = QUANTIF_NMAXLEVELS;
      wbm->qscale = wsigma>0.0? 2*QUANTIF_NSIGMA*wsigma/wbm->nlevels : 1.0;
      wbm->qzero = wmean - QUANTIF_NSIGMA*wsigma;
      wbm++;
      wbuf += bw;
      }
    }

  return;
  }


/****** back_histo ***********************************************************
PROTO	void back_histo(fieldstruct *field, fieldstruct *wfield,
		PIXTYPE *buf, PIXTYPE *wbuf, size_t bufsize,
			int n, int w, int bw, PIXTYPE wthresh)
PURPOSE	Compute image and variance histograms in a row of meshes.
INPUT	Pointer to background image mesh structure,
	pointer to background variance mesh structure,
	pointer to image buffer,
	pointer to variance buffer,
	buffer size (number of pixels),
	frame width (pixels),
	mesh width (pixels),
	weight threshold.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2011
 ***/
void	back_histo(backstruct *backmesh, backstruct *wbackmesh,
		PIXTYPE *buf, PIXTYPE *wbuf, size_t bufsize,
			int n, int w, int bw, PIXTYPE wthresh)

  {
   backstruct	*bm,*wbm;
   PIXTYPE	*buft,*wbuft;
   float	qscale,wqscale, cste,wcste, wpix;
   LONG		*histo,*whisto;
   int		h,m,x,y, nlevels,wnlevels, lastbite, offset, bin;

  h = bufsize/w;
  bm = backmesh;
  wbm = wbackmesh;
  offset = w - bw;
  for (m=0; m++<n; bm++ , buf+=bw)
    {
    if (m==n && (lastbite=w%bw))
      {
      bw = lastbite;
      offset = w-bw;
      }
/*-- Skip bad meshes */
    if (bm->mean <= -BIG)
      {
      if (wbackmesh)
        {
        wbm++;
        wbuf += bw;
        }
      continue;
      }
    nlevels = bm->nlevels;
    histo = bm->histo;
    qscale = bm->qscale;
    cste = 0.499999 - bm->qzero/qscale;
    buft = buf;
    if (wbackmesh)
      {
      wnlevels = wbm->nlevels;
      whisto = wbm->histo;
      wqscale = wbm->qscale;
      wcste = 0.499999 - wbm->qzero/wqscale;
      wbuft = wbuf;
      for (y=h; y--; buft+=offset, wbuft+=offset)
        for (x=bw; x--;)
          {
          bin = (int)(*(buft++)/qscale + cste);
          if ((wpix = *(wbuft++))<wthresh && bin<nlevels && bin>=0)
            {
        	(*(histo+bin))++;
            bin = (int)(wpix/wqscale + wcste);
            if (bin>=0 && bin<wnlevels)
              (*(whisto+bin))++;
            }
          }
      wbm++;
      wbuf += bw;
      }
    else
      for (y=h; y--; buft += offset)
        for (x=bw; x--;)
          {//QPRINTF(OUTPUT, " %.5g", *buft);
          bin = (int)(*(buft++)/qscale + cste);
          if (bin>=0 && bin<nlevels)
        	(*(histo+bin))++;
          //QPRINTF(OUTPUT, "->%i", bin);
          }
    }

  return;
  }


/****** back_guess ***********************************************************
PROTO	float back_guess(backstruct *bkg, float *mean, float *sigma)
PURPOSE	Estimate the background value and its variance from a histogram.
INPUT	Pointer to background structure,
	pointer to the background level estimate (output),
	pointer to the background variance estimate (output).
OUTPUT	Estimate of the background level.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2011
 ***/
float	back_guess(backstruct *bkg, float *mean, float *sigma)

  {
   LONG		*histo, *hilow, *hihigh, *histot;
   unsigned long lowsum, highsum, sum;
   double	ftemp, mea, sig, sig1, med, dpix;
   int		i, n, lcut,hcut, nlevelsm1, pix;

/* Leave here if the mesh is already classified as `bad' */
  if (bkg->mean<=-BIG)
    {
    *mean = *sigma = -BIG;
    return -BIG;
    }

  histo = bkg->histo;
  hcut = nlevelsm1 = bkg->nlevels-1;
  lcut = 0;

  sig = 10.0*nlevelsm1;
  sig1 = 1.0;
  mea = med = bkg->mean;
  for (n=100; n-- && (sig>=0.1) && (fabs(sig/sig1-1.0)>BACK_EPS);)
    {
    sig1 = sig;
    sum = mea = sig = 0.0;
    lowsum = highsum = 0;
    histot = hilow = histo+lcut;
    hihigh = histo+hcut;
    for (i=lcut; i<=hcut; i++)
      {
      if (lowsum<highsum)
        lowsum += *(hilow++);
      else
        highsum +=  *(hihigh--);
      sum += (pix = *(histot++));
      mea += (dpix = (double)pix*i);
      sig += dpix*i;
      }

    med = hihigh>=histo?
	((hihigh-histo)+0.5+((double)highsum-lowsum)/(2.0*(*hilow>*hihigh?
                                                *hilow:*hihigh)))
       : 0.0;
    if (sum)
      {
      mea /= (double)sum;
      sig = sig/sum - mea*mea;
      }
    sig = sig>0.0?sqrt(sig):0.0;
    lcut = (ftemp=med-3.0*sig)>0.0 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5):0;
    hcut = (ftemp=med+3.0*sig)<nlevelsm1 ?(int)(ftemp>0.0?ftemp+0.5:ftemp-0.5)
								: nlevelsm1;
    }
  *mean = fabs(sig)>0.0? (fabs(bkg->sigma/(sig*bkg->qscale)-1) < 0.0 ?
			    bkg->qzero+mea*bkg->qscale
			    :(fabs((mea-med)/sig)< 0.3 ?
			      bkg->qzero+(2.5*med-1.5*mea)*bkg->qscale
			     :bkg->qzero+med*bkg->qscale))
                       :bkg->qzero+mea*bkg->qscale;

  *sigma = sig*bkg->qscale;

  return *mean;
  }

/****** back_filter **********************************************************
PROTO	void back_filter(fieldstruct *field)
PURPOSE	Median filter a background map to remove the contribution from bright
	sources.
INPUT	Pointer to image field structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2011
 ***/
void	back_filter(fieldstruct *field)

  {
   float	*back,*sigma, *back2,*sigma2, *bmask,*smask, *sigmat,
		d2,d2min, fthresh, med, val,sval;
   int		i,j,px,py, np, nx,ny, npx,npx2, npy,npy2, dpx,dpy, x,y, nmin;

  fthresh = prefs.backfthresh;
  nx = field->nbackx;
  ny = field->nbacky;
  np = field->nback;
  npx = field->nbackfx/2;
  npy = field->nbackfy/2;
  npy *= nx;

  QMALLOC(bmask, float, (2*npx+1)*(2*npy+1));
  QMALLOC(smask, float, (2*npx+1)*(2*npy+1));
  QMALLOC(back2, float, np);
  QMALLOC(sigma2, float, np);

  back = field->back;
  sigma = field->sigma;
  val = sval = 0.0;			/* to avoid gcc -Wall warnings */

/* Look for `bad' meshes and interpolate them if necessary */
  for (i=0,py=0; py<ny; py++)
    for (px=0; px<nx; px++,i++)
      if ((back2[i]=back[i])<=-BIG)
        {
/*------ Seek the closest valid mesh */
        d2min = BIG;
        nmin = 0.0;
        for (j=0,y=0; y<ny; y++)
          for (x=0; x<nx; x++,j++)
            if (back[j]>-BIG)
              {
              d2 = (float)(x-px)*(x-px)+(y-py)*(y-py);
              if (d2<d2min)
                {
                val = back[j];
                sval = sigma[j];
                nmin = 1;
                d2min = d2;
                }
              else if (d2==d2min)
                {
                val += back[j];
                sval += sigma[j];
                nmin++;
                }
              }
        back2[i] = nmin? val/nmin: 0.0;
        sigma[i] = nmin? sval/nmin: 1.0;
        }
  memcpy(back, back2, (size_t)np*sizeof(float));

/* Do the actual filtering */
  for (py=0; py<np; py+=nx)
    {
    npy2 = np - py - nx;
    if (npy2>npy)
      npy2 = npy;
    if (npy2>py)
      npy2 = py;
    for (px=0; px<nx; px++)
      {
      npx2 = nx - px - 1;
      if (npx2>npx)
        npx2 = npx;
      if (npx2>px)
        npx2 = px;
      i=0;
      for (dpy = -npy2; dpy<=npy2; dpy+=nx)
        {
        y = py+dpy;
        for (dpx = -npx2; dpx <= npx2; dpx++)
          {
          x = px+dpx;
          bmask[i] = back[x+y];
          smask[i++] = sigma[x+y];
          }
        }
      if (fabs((med=fqmedian(bmask, i))-back[px+py])>=fthresh)
        {
        back2[px+py] = med;
        sigma2[px+py] = fqmedian(smask, i);
        }
      else
        {
        back2[px+py] = back[px+py];
        sigma2[px+py] = sigma[px+py];
        }
      }
    }

  free(bmask);
  free(smask);
  memcpy(back, back2, np*sizeof(float));
  field->backmean = fqmedian(back2, np);
  free(back2);
  memcpy(sigma, sigma2, np*sizeof(float));
  field->backsig = fqmedian(sigma2, np);

  if (field->backsig<=0.0)
    {
    sigmat = sigma2+np;
    for (i=np; i-- && *(--sigmat)>0.0;);
    if (i>=0 && i<(np-1))
      field->backsig = fqmedian(sigmat+1, np-1-i);
    else
      {
      if (field->flags&(DETECT_FIELD|MEASURE_FIELD))
        warning("Image contains mainly constant data; ",
		"I'll try to cope with that...");
      field->backsig = 1.0;
      }
    }

  free(sigma2);

  return;
  }


/****** back_local ***********************************************************
PROTO	float back_local(fieldstruct *field, objstruct *obj, float *sigma)
PURPOSE	Compute local background if possible.
INPUT	Pointer to image field structure,
	pointer to the source obj structure,
	pointer to local background variance estimate (output).
OUTPUT	Local background level estimate.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	21/12/2011
 ***/
float	back_local(fieldstruct *field, objstruct *obj, float *sigma)

  {
   static backstruct	backmesh;
   int			bxmin,bxmax, bymin,bymax, ixmin,ixmax, iymin,iymax,
			bxnml,bynml, oxsize,oysize, npix,
			i, x,y, bin, w,sh, bmn, pbs;
   float		dbkg, bqs,cste;
   LONG			*bmh;
   PIXTYPE		*backpix, *bp, *strip, *st,
			pix;

  strip = field->strip;
  w = field->width;
  sh = field->stripheight;
  pbs = prefs.pback_size;

/* Estimate background in a 'rectangular annulus' around the object */
  oxsize = obj->xmax - obj->xmin + 1;
  oysize = obj->ymax - obj->ymin + 1;
  bxnml = oxsize<w/2? oxsize/4 : (w-oxsize)/4;
  bynml = oysize<field->height/2? oysize/4 : (field->height-oysize)/4;
  bxmin = (ixmin = obj->xmin - bxnml) - pbs;
  bxmax = (ixmax = obj->xmax+1 + bxnml) + pbs;
  bymin = (iymin = obj->ymin - bynml) - pbs;
  bymax = (iymax = obj->ymax+1 + bynml) + pbs;

  if (bymin>=field->ymin && bymax<field->ymax
	&& bxmin>=0 && bxmax<w)
    {
    npix = (bxmax-bxmin)*(bymax-bymin) - (ixmax-ixmin)*(iymax-iymin);

    QMALLOC(backpix, PIXTYPE, npix);
    bp = backpix;

/*--store all the pixels*/
    npix = 0;
    for (y=bymin; y<bymax; y++)
      {
      st = strip + (y%sh)*w + bxmin;
      for (x=pbs; x--;)
        if ((pix=*(st++))>-BIG)
          {
          *(bp++) = pix;
          npix++;
          }
      st += ixmax-ixmin;
      for (x=pbs; x--;)
        if ((pix=*(st++))>-BIG)
          {
          *(bp++) = pix;
          npix++;
          }
      }

    for (y=bymin; y<iymin; y++)
      {
      st = strip + (y%sh)*w + ixmin;
      for (x=ixmax-ixmin; x--;)
        if ((pix=*(st++))>-BIG)
          {
          *(bp++) = pix;
          npix++;
          }
      }
    for (y=iymax; y<bymax; y++)
      {
      st = strip + (y%sh)*w + ixmin;
      for (x=ixmax-ixmin; x--;)
        if ((pix=*(st++))>-BIG)
          {
          *(bp++) = pix;
          npix++;
          }
      }

    if (npix)
      {
      back_stat(&backmesh, NULL, backpix, NULL, npix, 1, 1, 1, 0.0);
      QCALLOC(backmesh.histo, LONG, backmesh.nlevels);
      bmh = backmesh.histo;
      bmn = backmesh.nlevels;
      cste = 0.499999 - backmesh.qzero/(bqs = backmesh.qscale);
      bp = backpix;
      for (i=npix; i--;)
        {
        bin = (int)(*(bp++)/bqs + cste);
        if (bin>=0 && bin<bmn)
          (*(bmh+bin))++;
        }
      back_guess(&backmesh, &dbkg, sigma);
      free(backmesh.histo);
      }
    else
      {
      dbkg = 0.0;
      *sigma = -1.0;
      }
    free(backpix);
    }
  else
    {
    dbkg = 0.0;
    *sigma = -1.0;
    }

  return dbkg;
  }


/****** back_interpolate *****************************************************
PROTO	PIXTYPE back_interpolate(fieldstruct *field, double x, double y)
PURPOSE	Return interpolated background value at pixel position x,y (linear
	interpolation between background map vertices).
INPUT	Pointer to image field structure,
	pixel x coordinate,
	pixel y coordinate.
OUTPUT	Local background level estimate.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/12/2011
 ***/
PIXTYPE	back_interpolate(fieldstruct *field, double x, double y)

  {
   int		nx,ny, xl,yl, pos;
   double	dx,dy, cdx;
   float	*b, b0,b1,b2,b3;

  b = field->back;
  nx = field->nbackx;
  ny = field->nbacky;

  dx = x/field->backw - 0.5;
  dy = y/field->backh - 0.5;
  dx -= (xl = (int)dx);
  dy -= (yl = (int)dy);

  if (xl<0)
    {
    xl = 0;
    dx -= 1.0;
    }
  else if (xl>=nx-1)
    {
    xl = nx<2 ? 0 : nx-2;
    dx += 1.0;
    }

  if (yl<0)
    {
    yl = 0;
    dy -= 1.0;
    }
  else if (yl>=ny-1)
    {
    yl = ny<2 ? 0 : ny-2;
    dy += 1.0;
    }

  pos = yl*nx + xl;
  cdx = 1 - dx;

  b0 = *(b+=pos);		/* consider when nbackx or nbacky = 1 */
  b1 = nx<2? b0:*(++b);
  b2 = ny<2? *b:*(b+=nx);
  b3 = nx<2? *b:*(--b);

  return (PIXTYPE)((1-dy)*(cdx*b0 + dx*b1) + dy*(dx*b2 + cdx*b3));
  }


/****** back_makespline ******************************************************
PROTO	float *back_makespline(fieldstruct *field, float *map)
PURPOSE	Pre-compute spline 2nd derivatives along the y direction at background
	nodes.
INPUT	Pointer to image field structure,
	pointer to the background map.
OUTPUT	Pointer to the array of 2nd derivatives.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2011
 ***/
float *back_makespline(fieldstruct *field, float *map)
  {
   int		x,y, nbx,nby,nbym1;
   float	*dmap,*dmapt,*mapt, *u, temp;

  nbx = field->nbackx;
  nby = field->nbacky;
  nbym1 = nby - 1;
  QMALLOC(dmap, float, field->nback);
  for (x=0; x<nbx; x++)
    {
    mapt = map+x;
    dmapt = dmap+x;
    if (nby>1)
      {
      QMALLOC(u, float, nbym1);	/* temporary array */
      *dmapt = *u = 0.0;	/* "natural" lower boundary condition */
      mapt += nbx;
      for (y=1; y<nbym1; y++, mapt+=nbx)
        {
        temp = -1/(*dmapt+4);
        *(dmapt += nbx) = temp;
        temp *= *(u++) - 6*(*(mapt+nbx)+*(mapt-nbx)-2**mapt);
        *u = temp;
        }
      *(dmapt+=nbx) = 0.0;	/* "natural" upper boundary condition */
      for (y=nby-2; y--;)
        {
        temp = *dmapt;
        dmapt -= nbx;
        *dmapt = (*dmapt*temp+*(u--))/6.0;
        }
      free(u);
      }
    else
      *dmapt = 0.0;
    }

  return dmap;
  }


/****** back_subline *********************************************************
PROTO	void back_subline(fieldstruct *field, int y, PIXTYPE *line)
PURPOSE	Interpolate background at line y (bicubic spline interpolation between
	background map vertices) and subtract it from the current line.
INPUT	Pointer to image field structure,
	y coordinate of image line,
	pointer to the image line.

OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/04/2012
 ***/
void	back_subline(fieldstruct *field, int y, int xmin, int width,
		PIXTYPE *line)

  {
   int		i,j,x,yl, nbx,nbxm1,nby, nx, ystep, changepoint;
   float	dx,dx0,dy,dy3, cdx,cdy,cdy3, temp, xstep,
		*node,*nodep,*dnode, *blo,*bhi,*dblo,*dbhi, *u;
   PIXTYPE	*backline, bval;

  backline = field->backline;

  if (field->back_type==BACK_ABSOLUTE)
    {
/*-- In absolute background mode, just subtract a cste */
    bval = field->backmean;
    line += xmin;
    for (i=width; i--;)
      *(line++) -= ((*backline++)=bval);
    return;
    }

  nbx = field->nbackx;
  nbxm1 = nbx - 1;
  nby = field->nbacky;
  if (nby > 1)
    {
    dy = (float)y/field->backh - 0.5;
    dy -= (yl = (int)dy);
    if (yl<0)
      {
      yl = 0;
      dy -= 1.0;
      }
    else if (yl>=nby-1)
      {
      yl = nby<2 ? 0 : nby-2;
      dy += 1.0;
      }
/*-- Interpolation along y for each node */
    cdy = 1 - dy;
    dy3 = (dy*dy*dy-dy);
    cdy3 = (cdy*cdy*cdy-cdy);
    ystep = nbx*yl;
    blo = field->back + ystep;
    bhi = blo + nbx;
    dblo = field->dback + ystep;
    dbhi = dblo + nbx;
    QMALLOC(node, float, nbx);	/* Interpolated background */
    nodep = node;
    for (x=nbx; x--;)
      *(nodep++) = cdy**(blo++) + dy**(bhi++) + cdy3**(dblo++) + dy3**(dbhi++);

/*-- Computation of 2nd derivatives along x */
    QMALLOC(dnode, float, nbx);	/* 2nd derivative along x */
    if (nbx>1)
      {
      QMALLOC(u, float, nbxm1);	/* temporary array */
      *dnode = *u = 0.0;	/* "natural" lower boundary condition */
      nodep = node+1;
      for (x=nbxm1; --x; nodep++)
        {
        temp = -1/(*(dnode++)+4);
        *dnode = temp;
        temp *= *(u++) - 6*(*(nodep+1)+*(nodep-1)-2**nodep);
        *u = temp;
        }
      *(++dnode) = 0.0;	/* "natural" upper boundary condition */
      for (x=nbx-2; x--;)
        {
        temp = *(dnode--);
        *dnode = (*dnode*temp+*(u--))/6.0;
        }
      free(u);
      dnode--;
      }
    }
  else
    {
/*-- No interpolation and no new 2nd derivatives needed along y */
    node = field->back;
    dnode = field->dback;
    }

/*-- Interpolation along x */
  if (nbx>1)
    {
    nx = field->backw;
    xstep = 1.0/nx;
    changepoint = nx/2;
    dx  = (xstep - 1)/2;	/* dx of the first pixel in the row */
    dx0 = ((nx+1)%2)*xstep/2;	/* dx of the 1st pixel right to a bkgnd node */
    blo = node;
    bhi = node + 1;
    dblo = dnode;
    dbhi = dnode + 1;
    for (x=i=0,j=xmin; j--; i++, dx += xstep)
      {
      if (i==changepoint && x>0 && x<nbxm1)
        {
        blo++;
        bhi++;
        dblo++;
        dbhi++;
        dx = dx0;
        }
      cdx = 1 - dx;
      if (i==nx)
        {
        x++;
        i = 0;
        }
      }

    for (j=width; j--; i++, dx += xstep)
      {
      if (i==changepoint && x>0 && x<nbxm1)
        {
        blo++;
        bhi++;
        dblo++;
        dbhi++;
        dx = dx0;
        }
      cdx = 1 - dx;
      
      *(line++) -= (*(backline++) = (PIXTYPE)(cdx*(*blo+(cdx*cdx-1)**dblo)
			+ dx*(*bhi+(dx*dx-1)**dbhi)));
      if (i==nx)
        {
        x++;
        i = 0;
        }
      }
    }
  else
    for (j=width; j--;)
      *(line++) -= (*(backline++) = (PIXTYPE)*node);

  if (nby>1)
    {
    free(node);
    free(dnode);
    }

  return;
  }


/****** back_rmsline *********************************************************
PROTO	void back_rmsline(fieldstruct *field, int y, PIXTYPE *line)
PURPOSE	Bicubic-spline interpolation of the background noise along the current
	scanline (y).
INPUT	Measurement or detection field pointer,
	current line position. 
	where to put the data. 
OUTPUT	-.
NOTES	Most of the code is a copy of back_subline(), for optimization reasons.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2011
 ***/
void	back_rmsline(fieldstruct *field, int y, PIXTYPE *line)

  {
   int		i,j,x,yl, nbx,nbxm1,nby, nx,width, ystep, changepoint;
   float	dx,dx0,dy,dy3, cdx,cdy,cdy3, temp, xstep,
		*node,*nodep,*dnode, *blo,*bhi,*dblo,*dbhi, *u;

  nbx = field->nbackx;
  nbxm1 = nbx - 1;
  nby = field->nbacky;
  if (nby > 1)
    {
    dy = (float)y/field->backh - 0.5;
    dy -= (yl = (int)dy);
    if (yl<0)
      {
      yl = 0;
      dy -= 1.0;
      }
    else if (yl>=nby-1)
      {
      yl = nby<2 ? 0 : nby-2;
      dy += 1.0;
      }
/*-- Interpolation along y for each node */
    cdy = 1 - dy;
    dy3 = (dy*dy*dy-dy);
    cdy3 = (cdy*cdy*cdy-cdy);
    ystep = nbx*yl;
    blo = field->sigma + ystep;
    bhi = blo + nbx;
    dblo = field->dsigma + ystep;
    dbhi = dblo + nbx;
    QMALLOC(node, float, nbx);	/* Interpolated background */
    nodep = node;
    for (x=nbx; x--;)
      *(nodep++) = cdy**(blo++) + dy**(bhi++) + cdy3**(dblo++) + dy3**(dbhi++);

/*-- Computation of 2nd derivatives along x */
    QMALLOC(dnode, float, nbx);	/* 2nd derivative along x */
    if (nbx>1)
      {
      QMALLOC(u, float, nbxm1);	/* temporary array */
      *dnode = *u = 0.0;	/* "natural" lower boundary condition */
      nodep = node+1;
      for (x=nbxm1; --x; nodep++)
        {
        temp = -1/(*(dnode++)+4);
        *dnode = temp;
        temp *= *(u++) - 6*(*(nodep+1)+*(nodep-1)-2**nodep);
        *u = temp;
        }
      *(++dnode) = 0.0;	/* "natural" upper boundary condition */
      for (x=nbx-2; x--;)
        {
        temp = *(dnode--);
        *dnode = (*dnode*temp+*(u--))/6.0;
        }
      free(u);
      dnode--;
      }
    }
  else
    {
/*-- No interpolation and no new 2nd derivatives needed along y */
    node = field->sigma;
    dnode = field->dsigma;
    }

/*-- Interpolation along x */
  width = field->width;
  if (nbx>1)
    {
    nx = field->backw;
    xstep = 1.0/nx;
    changepoint = nx/2;
    dx  = (xstep - 1)/2;	/* dx of the first pixel in the row */
    dx0 = ((nx+1)%2)*xstep/2;	/* dx of the 1st pixel right to a bkgnd node */
    blo = node;
    bhi = node + 1;
    dblo = dnode;
    dbhi = dnode + 1;
    for (x=i=0,j=width; j--; i++, dx += xstep)
      {
      if (i==changepoint && x>0 && x<nbxm1)
        {
        blo++;
        bhi++;
        dblo++;
        dbhi++;
        dx = dx0;
        }
      cdx = 1 - dx;
      *(line++) = (PIXTYPE)(cdx*(*blo+(cdx*cdx-1)**dblo)
			+ dx*(*bhi+(dx*dx-1)**dbhi));
      if (i==nx)
        {
        x++;
        i = 0;
        }
      }
    }
  else
    for (j=width; j--;)
      *(line++) = (PIXTYPE)*node;

  if (nby>1)
    {
    free(node);
    free(dnode);
    }

  return;
  }


/****** back_copy ************************************************************
PROTO	void back_copy(fieldstruct *infield, fieldstruct *outfield)
PURPOSE Copy sub-structures related to background procedures (mainly copying
	memory).
INPUT	Input field pointer to be copied,
	destination field pointer. 
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2011
 ***/
void	back_copy(fieldstruct *infield, fieldstruct *outfield)

  {
  QMEMCPY(infield->back, outfield->back, float, infield->nback);
  QMEMCPY(infield->dback, outfield->dback, float, infield->nback);
  QMEMCPY(infield->sigma, outfield->sigma, float, infield->nback);
  QMEMCPY(infield->dsigma, outfield->dsigma, float, infield->nback);
  QMEMCPY(infield->backline, outfield->backline, PIXTYPE, infield->width);

  return;
  }


/****** back_end *************************************************************
PROTO	void back_end(fieldstruct *field)
PURPOSE Terminate background procedures (freeing memory).
INPUT	Image field pointer. 
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2011
 ***/
void	back_end(fieldstruct *field)

  {
  free(field->back);
  free(field->dback);
  free(field->sigma);
  free(field->dsigma);
  free(field->backline);

  return;
  }

/****** back_printmeshs ******************************************************
PROTO	void back_printmeshs(const backstruct *backmesh, const int nmeshs)
PURPOSE Print info on a mesh structure.
INPUT	Pointer to a  mesh structure, number of meshes
OUTPUT	-.
NOTES	-.
AUTHOR	M. Kuemmel (LMU)
VERSION	12/02/2014
 ***/
void	back_printmeshs(backstruct *backmesh, const int nmeshs, int *tofile)
{
  char filename[MAXCHAR];
  int index;
  FILE *outstream;
  backstruct *bm;

  // go over the meshes
  bm=backmesh;
  for (index=0; index<nmeshs; index++, bm++){
      if (*tofile){
          sprintf(filename, "histoMsk%i.txt", *tofile);
          outstream = fopen(filename, "w+");
          back_printmesh(bm, outstream);
          fclose(outstream);
          (*tofile)++;
      }
      else{
          // print one mesh
          back_printmesh(bm, OUTPUT);
      }
  }
}

/****** back_printmesh ******************************************************
PROTO	void back_printmesh(const backstruct *backmesh)
PURPOSE Print info on a mesh structure.
INPUT	Pointer to a  mesh structure
OUTPUT	-.
NOTES	-.
AUTHOR	M. Kuemmel (LMU)
VERSION	12/02/2014
 ***/
void	back_printmesh(const backstruct *backmesh, FILE * outstream)
{
  int index;
  QPRINTF(outstream, "# Mode: %.5g, Mean: %.5g, Sigma: %.5g, Npixel: %i\n", backmesh->mode, backmesh->mean, backmesh->sigma, backmesh->npix);
  QPRINTF(outstream, "# Lcut: %.5g, Hcut: %.5g, Qzero: %.5g Qscale: %.5g\n", backmesh->lcut, backmesh->hcut, backmesh->qzero, backmesh->qscale);
  QPRINTF(outstream, "# Nlevels: %i, \n", backmesh->nlevels);
  for (index=0; index<backmesh->nlevels; index++)
    //if (backmesh->histo[index])
      QPRINTF(outstream, "%i %.5g %i\n", index, backmesh->qzero+(float)index*backmesh->qscale, backmesh->histo[index]);
}

// TODO: put everything from here into new file
/****************************create_objmask***********************************/
/**
 * Function: create_objmask
 *
 * Create an object mask of the specified size and in the desired structure.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] nx       number of pixels in x
 * @param[in] ny       number of pixels in y
 * @param[in] swapflag mapp to a file [1] or keep in memory [0]
 * @param[in] index    file index, which will be part of the name
 *
 * @return the generated object mask
 */
objmaskstruct *create_objmask(const size_t nx, const size_t ny, const int swapflag, const int index)
{
  objmaskstruct *omask;
  //size_t index;
  bool *abool;

  // basic allocation
  omask = (objmaskstruct *)malloc(sizeof(objmaskstruct));

  // copy the "constants"
  omask->width  = nx;
  omask->height = ny;
  omask->npix   = nx*ny;

  // store the array size
  omask->size = omask->npix*sizeof(bool);

  // check for the swap flag
  if (!swapflag)
    {
      // malloc the big array
      omask->maskdata = (bool *)malloc(omask->npix*sizeof(bool));
      if (!omask->maskdata)
        {
          error(EXIT_FAILURE, "Problems allocating memory", "");
        }

      // mark the other elements as not being used
      omask->swapfile=NULL;
      sprintf(omask->swapname, "");
    }
  else
    {
      // get the swap filename
      get_obmask_name(index, omask->swapname);

      // open the swap file
      omask->swapfile = fopen(omask->swapname, "w+");
      if (!omask->swapfile)
        {
          error(EXIT_FAILURE, "Can not open swap file: ", omask->swapname);
        }

      // stretch the file to the desired size
      if (fseek(omask->swapfile, omask->size-1, SEEK_SET) == -1) {
          fclose(omask->swapfile);
          error(EXIT_FAILURE, "Error calling fseek() to 'stretch' the file: ", omask->swapname);
      }

      //int result=0;
      // the fwrite() line should work as well, but somehow
      // it dumps away, probably someparams are wrong
      //if (fwrite("", 1, 1, omask->swapfile) != 1) {
      //result = write(fileno(omask->swapfile), "", 1);
      //if (result != 1) {
      //result = fwrite("o", 1, 1, omask->swapfile);
      //if (result !=1){
      //fprintf(stderr, "rrrr: %i\n", result);

      // write something at that point to get the file size
      if (write(fileno(omask->swapfile), "", 1) != 1) {
          fclose(omask->swapfile);
          error(EXIT_FAILURE, "Error writing last byte of the file: ", omask->swapname);
      }

      // make the memory mapping and ceck that it worked
      omask->maskdata = mmap(0, omask->size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(omask->swapfile), 0);
      if (omask->maskdata == MAP_FAILED) {
          fclose(omask->swapfile);
          error(EXIT_FAILURE, "Error mmapping the file: ", omask->swapname);
      }
    }

  // set swap flag
  omask->swapflag=swapflag;

  // set the stepper position
  omask->actindex=0;
  omask->actpos=omask->maskdata;

  // initialize all values
  if (!memset(omask->maskdata, 0, omask->size))
    {
      if (omask->swapfile)
        fclose(omask->swapfile);
      error(EXIT_FAILURE, "Error setting the memory to FALSE.", "");
    }
  // that's a classical variant of the initialization
  //for (index=0; index<omask->npix; index++)
  //  omask->maskdata[index]=true;

  return omask;
}

/****************************free_objmask***********************************/
/**
 * Function: free_objmask
 *
 * Free the memory in an object mask
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] omask  the object mask
 */
void free_objmask(objmaskstruct *omask)
{

  // check the swap flag
  if (!omask->swapflag)
    {
      // free the data in memory
      free(omask->maskdata);
      omask->maskdata=NULL;
    }
  else
    {
      // unmap the memory
      if (munmap(omask->maskdata, omask->size) == -1) {
          fclose(omask->swapfile);
          error(EXIT_FAILURE, "Problems unmapping the swap file: ", omask->swapname);
      }

      // close the file and set to NULL
      if (fclose(omask->swapfile))
        warning("Problems closing the swap file: ", omask->swapname);
      omask->swapfile=NULL;

      // delete  the swap file
      if (unlink(omask->swapname))
        warning("Problems deleting swap file: ", omask->swapname);

    }

  // free the struct
  free(omask);
  omask=NULL;
}

/****************************get_obmask_name**********************************/
/**
 * Function: get_obmask_name
 *
 * Creates an unique file name for the object mask. The file is located in the
 * tmp-directory and contains the creation date/time, process number and an
 * index, making the name rather unique.
 *
 * @author MK
 * @date   Jan 2015
 *
 * @param[in] index    file index
 * @param[in] filename file name   y-coordinate of element
 */
void get_obmask_name(const int index, char *filename)
{
  time_t    thetime=time(NULL);
  struct tm *tm= localtime(&thetime);

  // get the swap filename
  sprintf(filename, "./objmask_%02d-%02d_%02d:%02d:%02d_%05ld_%03i.tmp",
      tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min, tm->tm_sec,
      (long)getpid(), index);

  // check whether the object exists
  // already, which really should not be
  if(fopen(filename, "r"))
    {
      error(EXIT_FAILURE, "File does already exist: ", filename);
    }

return;
}

/****************************set_objmask_value*********************************/
/**
 * Function: set_objmask_value
 *
 * Set an entry in an object mask. The integer input value is transformed
 * to true/false.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] xpos   x-coordinate of element
 * @param[in] ypos   y-coordinate of element
 * @param[in] value  value to set the element to (0=false, !=0 = true)
 * @param[in] omask  the object mask of the element
 */
void set_objmask_value(const size_t xpos, const size_t ypos, const int value, objmaskstruct *omask)
{
  // check for a valid element
  if (xpos >= omask->width || ypos >= omask->height) {
      char err_msg[MAXCHAR];
      sprintf (err_msg, "Array size: (%zu,%zu), index (%zu,%zu) does not exist!", omask->width, omask->width, xpos, ypos);
      error(EXIT_FAILURE, err_msg, "");
    }
  omask->maskdata[ypos*omask->width+xpos] = value ? true : false;
}

/****************************get_objmask_value*********************************/
/**
 * Function: get_objmask_value
 *
 * Get an entry in an object mask. True/false is translated to 1/0.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] xpos   x-coordinate of element
 * @param[in] ypos   y-coordinate of element
 * @param[in] omask  the object mask of the element
 *
 * @return 1/0 for true/false
 */
int get_objmask_value(const size_t xpos, const size_t ypos, const objmaskstruct *omask)
{
  // check for a valid element
  if (xpos >= omask->width || ypos >= omask->height) {
      char err_msg[MAXCHAR];
      sprintf (err_msg, "Array size: (%zu,%zu), index (%zu,%zu) does not exist!", omask->width, omask->width, xpos, ypos);
      error(EXIT_FAILURE, err_msg, "");
    }
  return omask->maskdata[ypos*omask->width+xpos] ? 1 : 0;
}

/****************************reset_objmask************************************/
/**
 * Function: reset_objmask
 *
 * Initialize all positional values in an object mask.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] omask  the object mask
 */
void reset_objmask(objmaskstruct *omask){
  // put the index and position to the beginning
  omask->actindex=0;
  omask->actpos=omask->maskdata;
  return;
}

/****************************tell_objmask*************************************/
/**
 * Function: tell_objmask
 *
 * Kind of a parallel function to "ftell()" which returns the current index
 * in the object mask;
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] omask  the object mask
 *
 * @return the current index position in the object mask
 */
long int tell_objmask(const objmaskstruct *omask)
{
  return omask->actindex;
}

/****************************seek_objmask*************************************/
/**
 * Function: seek_objmask
 *
 * Kind of a parallel function to "fseek()" which changes the position in an
 * object mask. If the new position would be outside of the array,  an error
 * is given.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] offset the object mask
 * @param[in] origin the object mask
 * @param[in] omask  the object mask
 */
void seek_objmask(const long int offset, const int origin, objmaskstruct *omask)
{
  switch (origin)
  {
  case SEEK_SET:
    // make sure the position exists
    if (offset<0)
      error(EXIT_FAILURE, "Can not seek position before the start of mask: ", omask->swapname);
    else if (offset>omask->npix)
      error(EXIT_FAILURE, "Can not seek position after the end of mask: ", omask->swapname);

    // change to the desired position
    omask->actindex = offset;
    omask->actpos   = omask->maskdata+offset;
    break;
  case SEEK_CUR:
    // make sure the position exists
    if (offset<0 && (-offset>omask->actindex))
      error(EXIT_FAILURE, "Can not seek position before the start of mask: ", omask->swapname);
    else if (offset>=0 && omask->actindex+offset>omask->npix)
      error(EXIT_FAILURE, "Can not seek position after the end of mask: ", omask->swapname);

    // change to the desired position
    omask->actindex += offset;
    omask->actpos   += offset;
    break;
  case SEEK_END:
    error(EXIT_FAILURE, "The origin SEEK_END is not supported!%s", "");
    break;
  default:
    {
      // some unknown origin, nothing can be done
      char errint[MAXCHAR];
      sprintf(errint, "%d", origin);
      error(EXIT_FAILURE, "The chosen origin is not implemented: ", errint);
    break;
    }
  }
}

/****************************read_objmaskOld**********************************/
/**
 * Function: read_objmaskOld
 * Old version of "read_objmask". Could be removed at some point.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] nobj   number of objects to read
 * @param[in] omask  the object mask
 *
 * @return the current position to read nobject elements from
 */
bool *read_objmaskOld(const size_t nobj, objmaskstruct *omask)
{
  bool *actpos;

  // make sure there are enough elements to access
  if ((size_t)omask->actindex+nobj > omask->npix)
    error(EXIT_FAILURE, "Not enough elements in swapfile:", omask->swapname);

  // store the current position
  actpos = omask->actpos;

  // prepare the next position
  omask->actindex += (long int)nobj;
  omask->actpos   += nobj;

  // return the current position
  return actpos;
}

/****************************read_objmask**********************************/
/**
 * Function: read_objmask
 *
 * Kind of a fake function which gives back the actual stored position in an
 * object mask together with the guaranty that nobj elements can be read
 * from it. Naming is similar to "read_tab()" to be used in parallel to
 * this function.
 * @author MK
 * @date   June 2014
 *
 * @param[in] nobj   number of objects to read
 * @param[in] bbuff  address of actual position to read from
 * @param[in] omask  the object mask
 */
void read_objmask(const size_t nobj, bool **boolbuff, objmaskstruct *omask)
{
  // make sure that nobj elements can be read
  if ((size_t)omask->actindex+nobj > omask->npix)
    error(EXIT_FAILURE, "Not enough elements in swapfile:", omask->swapname);

  // copy the address of the position
  *boolbuff = omask->actpos;

  // prepare the next position
  omask->actindex += (long int)nobj;
  omask->actpos   += nobj;
}

/****************************apply_objmask**********************************/
/**
 * Function: apply_objmask
 * Applies the object mask to the values in the pixel buffer. The mask values
 * starting from the current position of the internal index of the mask is
 * applied to the buffer values. This means that the masks need to be in sync
 * with the buffer. The mask is applied by setting the values to "-BIG", such
 * that they are not used in the background detemrination.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in,out] buffer   number of objects to read
 * @param[in]     npix  address of actual position to read from
 * @param[in]     omask  the object mask
 */
void apply_objmask(PIXTYPE *buffer, size_t npix, objmaskstruct *omask)
{
  PIXTYPE *actbuff;
  bool    *actbool;
  size_t  index;

  // "read" the next mask pixels
  read_objmask(npix, &actbool, omask);

  // go over the buffer
  actbuff = buffer;
  for (index=0; index<npix; index++, actbuff++, actbool++)
    {
      if (*actbool)
        {
          //set masked pixels to -BIG
          *actbuff = -BIG;
        }
    }
}

/****************************print_objmask_elem*********************************/
/**
 * Function: print_objmask_elem
 *
 * Print an element in an object mask.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] xpos   x-coordinate of element
 * @param[in] ypos   y-coordinate of element
 * @param[in] omask  the object mask of the element
 */
void print_objmask_elem(const size_t xpos, const size_t ypos, const objmaskstruct *omask)
{
  fprintf(OUTPUT, "%s ", get_objmask_value(xpos, ypos, omask)? "true" : "false");
}

/****************************print_objmask*********************************/
/**
 * Function: print_objmask
 *
 * Print all element in an object mask.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] omask  the object mask of the element
 */
void print_objmask(const objmaskstruct *omask)
{
  size_t i, j;
  for (j=0; j<omask->height; j++)
    for (i=0; i<omask->width; i++)
      print_objmask_elem(i, j, omask);
  fprintf(OUTPUT, "\n");
}

/****************************objmask_info*********************************/
/**
 * Function: objmask_info
 *
 * Print information on an object mask.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in] omask  the object mask
 */
void objmask_info(const objmaskstruct *omask)
{
  size_t i;
  size_t n_true=0, n_false=0;
  bool   *act_bool;

  // go over the data and count the true's and false's
  act_bool=omask->maskdata;
  for (i=0, act_bool=omask->maskdata; i < omask->npix; i++)
    if (*(act_bool++))
      n_true++;
    else
      n_false++;

  // print the basic information
  // on the object mask
  if (omask->swapfile)
    QPRINTF(OUTPUT, "Object mask in file: %s\n", omask->swapname);
  QPRINTF(OUTPUT, "Object mask dimension: %zux%zu pix\n", omask->width, omask->height);
  QPRINTF(OUTPUT, "fraction of masked pixels: %.2f%%\n", 100.0*(double)n_true/(double)omask->npix);
  QPRINTF(OUTPUT, "# of masked pixels: %zu, unmasked: %zu\n", n_true, n_false);
}


/****************************populate_objmask*********************************/
/**
 * Function: populate_objmask
 * Populate an object mask by going through a field line by line, convolving
 * with the filter kernel and identifying pixels above a certain threshold.
 *
 * @author MK
 * @date   June 2014
 *
 * @param[in]     dfield  the detection field
 * @param[in]     dwfield the weight for the detection field
 * @param[in,out] omask   the object mask
 */
void populate_objmask(fieldstruct *dfield, fieldstruct *dwfield, objmaskstruct *omask)
{
  size_t width, height;
  size_t y_act, x_act;
  int varthreshflag;
  PIXTYPE wthresh, relthresh, thresh;
  PIXTYPE *scan=NULL, *wscan=NULL;
  PIXTYPE *cscan=NULL, *cwscan=NULL;
  OFF_T   fcurpos,wfcurpos;

  // If WEIGHTing and no absolute thresholding, activate threshold scaling */
  varthreshflag = (dwfield && prefs.thresh_type[0]!=THRESH_ABSOLUTE);
  relthresh = varthreshflag ? prefs.dthresh[0] : 0.0;// to avoid gcc warnings
  width  = (size_t)dfield->width;
  height = (size_t)dfield->height;

  // store the threshold value
  thresh = dfield->dthresh;

  // do some initialization on the fields,
  // to prepare the reading
  scan_initmarkers(dfield);
  scan_initmarkers(dwfield);

  // memorize the original position in the file
  QFTELL(dfield->cat->file, fcurpos, dfield->cat->filename);
  if (dwfield){
      QFTELL(dwfield->cat->file, wfcurpos, dwfield->cat->filename);
  }

  // allocate memory for buffers
  if (prefs.filter_flag)
    {
      QMALLOC(cscan, PIXTYPE, width);
      if (dwfield)
        {
          QCALLOC(cwscan, PIXTYPE, width);
        }
    }

  // go over all lines
  for (y_act=0; y_act<height; )
    {
      // read in the weight image
      if (dwfield){
        wscan = (dwfield->stripy==dwfield->stripysclim)
            ? (PIXTYPE *)readimage_loadstrip(dwfield, (fieldstruct *)NULL, 0)
            : &dwfield->strip[dwfield->stripy*dwfield->width];
      }

      // read in the detection image
      scan = (dfield->stripy==dfield->stripysclim)
          ? (PIXTYPE *)readimage_loadstrip(dfield, dwfield, 0)
          : &dfield->strip[dfield->stripy*dfield->width];

      if (prefs.filter_flag)
        {
          // filter the detection image row
          filter(dfield, cscan, dfield->y);
          if (dwfield)
            {
              // filter the weight image row
              filter(dwfield, cwscan, (int)y_act);
            }
        }
      else
        {
          // use the unfiltered rows
          cscan = scan;
          cwscan = wscan;
        }

      // go over all columns
      for (x_act=0; x_act<width; x_act++)
        {
          // if necessary, adjust the threshold
          if (varthreshflag)
            thresh = relthresh*sqrt(cwscan[x_act]);

          // mark a detected pixel in the object mask
          if (cscan[x_act]>thresh)
            set_objmask_value(x_act, y_act, 1, omask);

        }

      // prepare markers for the next line
      y_act++;
      scan_updatemarkers(dfield, y_act);
      scan_updatemarkers(dwfield, y_act);
    }

  // go back to the original position in the files, allowing them to be re-read
  QFSEEK(dfield->cat->file, fcurpos, SEEK_SET, dfield->cat->filename);
  if (dwfield)
    QFSEEK(dwfield->cat->file, wfcurpos, SEEK_SET, dwfield->cat->filename);

  // release the buffer memory
  if (prefs.filter_flag){
      free(cscan);
      cscan=NULL;
      if (dwfield){
          free(cwscan);
          cwscan=NULL;
      }
  }

  return;
}
