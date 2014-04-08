/*
*				filter.c
*
* 
* Filter image rasters (for detection).
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		11/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"bpro.h"
#include	"filter.h"
#include	"psf.h"
#include	"wcs/poly.h"
#include	"image.h"

/**
 * Function: convolve
 *
 * Convolves a scan line with an array. The row number to convole is given
 * as input. The convoluted line is returned, hence it does NOT work
 * in place.
 * The main iteration is NOT over the image pixels but over the filter
 * values.
 *
 * @param[in]  field  structure with the scan data
 * @param[out] mscan the convolved image line
 * @param[in]  y     image row to convolve
 *
 */
void convolve(fieldstruct *field, PIXTYPE *mscan, int y)

{
        int             mw,mw2,m0,me,m,mx,dmx, y0, dy, sw,sh;
        float   *mask;
        PIXTYPE *mscane, *s,*s0, *d,*de, mval;

        sw = field->width;        // width of the image
        sh = field->stripheight;  // height of the data available
        mw = thefilter->convw;    // certainly the width of the kernel
        mw2 = mw/2;               // half width of the kernel
        mscane = mscan+sw;        // the end of the result vector

        // fix the y-value to start the convolution
        y0 = y - (thefilter->convh/2);

        // check whether the starting
        // y-value is available
        if ((dy = field->ymin-y0) > 0)
        {
                // fix the smallest available y-value
                y0 = field->ymin;

                // restrict the convolution kernel
                // to the usable area.
                // NOTE: the lower part of y in the image corresponds
                //       to the upper end in the kernel, hence the
                //       upper end needs to be adjusted
                me = mw*(thefilter->convh-dy);
        }
        else
                // us the entire kernel as default
                me = mw*thefilter->convh;

        // check whether the upper
        // y-value is available
        if ((dy = field->ymax - y0) < thefilter->convh){
                // if the high y-values are not available,
                // restrict the convolution kernel to the usable area.
                // NOTE: the upper part of y in the image corresponds
                //       to the lower part in the kernel, hence the
                //       lower end needs to be adjusted
                m0 = mw*(thefilter->convh-dy);
        }
        else
                // us the entire kernel as default
                m0 = 0;

        // m0 and me mark the pointer indices in
        // the convolution kernel that can be used

        // initialize the resulting vector
        memset(mscan, 0, sw*sizeof(PIXTYPE));
        s0 = NULL;                              /* To avoid gcc -Wall warnings */

        // set mask to the end of the kernel
        mask = thefilter->conv+me;

        // iterate over all usable filter
        // values; for each value add its
        // contribution to the corresponding
        // value of the filtered array
        for (m = m0, mx = 0; m<me; m++, mx++)
        {
                // check for a new row
                // in the convolution kernel
                if (mx==mw)
                        mx = 0;

                // jump to a new row in the
                // image data
                if (!mx)
                        s0 = field->strip+sw*((y0++)%sh);

                // make reasonable start
                // and end values
                if ((dmx = mx-mw2)>=0)
                {
                        s = s0 + dmx;
                        d = mscan;
                        de = mscane - dmx;
                }
                else
                {
                        s = s0;
                        d = mscan - dmx;
                        de = mscane;
                }

                // get the kernel value
                mval = *(--mask);

                // go over the row
                // and add the contribution
                // to the result vector
                while (d<de)
                        *(d++) += mval**(s++);
        }

        return;
}


/**
 * Function: getconv
 *
 * Read the convolution filter from a file. Accepted file formats are an ASCII
 * filter or a PSFEx file.
 *
 * @param[in] filename the filename
 *
 * @return flag indicating success or error
 */
int	getconv(const char *filename)
{
	// check whether the filename ends with ".psf"
	if(strlen(filename) > 4 && !strcmp(filename + strlen(filename) - 4, ".psf"))
		// load the convolution filter from a psf file
		return getPSFExconv(filename);
	else
		// load the convolution filter from an ASCII
		return getASCIIconv(filename);

	// it should NEVER arrive here
	return RETURN_ERROR;
}

/**
 * Function: getASCIIconv
 *
 * Read the convolution mask from an ASCII file.
 *
 * @param[in] filename the filename
 *
 * @return flag indicating success or error
 */
int	getASCIIconv(const char *filename)

  {
  FILE		*file;
  char		str[MAXCHAR], *sstr, *null = NULL;
  double	sum, var, pix;
  int		i,j, n, normflag;


/* Open the file which may contain a convolution mask */

  if (!(file = fopen(filename, "r")))
    error(EXIT_FAILURE, "*Error*: cannot open ", filename);

/* Check it is a convolution mask. Otherwise, exit */
  fgets(str,MAXCHAR,file);
  if (strncmp(str,"CONV",4))
    {
    fclose(file);
    return RETURN_ERROR;
    }

  // set the flag for normalization
  if (strstr(str,"NORM"))
    normflag = strstr(str,"NONORM")?0:1;
  else
    {
    warning("No normalization info found in convolution file (old format?)\n",
	"> => I will assume you want the mask to be normalized...");
    normflag = 1;
    }

/* Allocate memory for storing mask elements */
  QMALLOC(thefilter->conv, float, MAXMASK);

  for (i=0, n=0; fgets(str,MAXCHAR,file);)
    {
    j=0;
    sstr = strtok(str, " \t\n");
    if (sstr && sstr[0]!=(char)'#')
      do
        {
        j++;
        thefilter->conv[i++] = (float)atof(sstr);
        if (i>MAXMASK)
          error(EXIT_FAILURE, "*Error*: Convolution mask too large in ",
		filename);
        }	while ((sstr = strtok(null, " \t\n")));
    if (j>n)
      n = j;
    }

  fclose(file);

  if ((thefilter->convw = n)<1)
    error(EXIT_FAILURE, "*Error*: unappropriate convolution mask width in ",
	filename);
  if (i%n)
    error(EXIT_FAILURE, "*Error*: unappropriate convolution mask line in ",
	filename);

  QREALLOC(thefilter->conv, float, i);

  thefilter->convh = i/n;
  var = 0.0;
  for (j=0, sum=0.0; j<i; j++)
    {
    sum += fabs(pix = thefilter->conv[j]);
    var += pix*pix;
    }

  thefilter->varnorm = (float)(var=sqrt(var));

  if (normflag)
    {
    if (sum == 0.0)
      {
      warning("Zero-sum filter: ", "Normalization switched to variance-mode");
      sum = var;
      }
    for (j=0; j<i; j++)
      thefilter->conv[j] /= sum;
    }

  thefilter->nconv = thefilter->convw*thefilter->convh;

  return RETURN_OK;
  }

/**
 * Function: getPSFExconv
 *
 * Read the convolution mask a PSFEx file
 *
 * @param[in] filename the filename
 *
 * @return flag indicating success or error
 */
int getPSFExconv(const char *filename){
	psfstruct *psf;
	int ndim, index, posx, posy;
	char *ctxname;
	double pos[POLY_MAXDIM];
	int w2, h2, ixmean, iymean;
	float pixstep, *pix2, stats[4];
	int xstart, xend, ii, ystart, yend, jj, icent;
	float *rpos, *actpos, *fpos, sum;
	float FRAC=0.95;

	// load the psf file
	psf=psf_load(filename);

	// make sure there are at most two dimensions
	ndim = psf->poly->ndim;
	if (ndim>2)
		error(EXIT_FAILURE, "*Error*: more than two dimensions in: ", filename);

	// make sure the context variables
	// mark positions by checking the names
	posx=-1;
	posy=-1;
	for (index=0; index<ndim; index++){
		ctxname = psf->contextname[index];

		// make sure the context name ends with "_IMAGE"
		if(strlen(ctxname) < 6 || strcmp(ctxname + strlen(ctxname) - 6, "_IMAGE"))
			error(EXIT_FAILURE, "*Error*: no image coordinate in context: ", ctxname);

		// check for the allowed x-values; make sure only one context points to x
		if (!strncmp(ctxname, "X", 1) || !strncmp(ctxname, "XMIN", 4) || !strncmp(ctxname, "XMAX", 4) || !strncmp(ctxname, "XPEAK", 5))
			if (posx>-1)
				error(EXIT_FAILURE, "*Error*: two context parameters point to X, the last being: ", ctxname);
			else
				posx=index;

		// check for the allowed y-values; make sure only one context points to y
		if (!strncmp(ctxname, "Y", 1) || !strncmp(ctxname, "YMIN", 4) || !strncmp(ctxname, "YMAX", 4) || !strncmp(ctxname, "YPEAK", 5))
			if (posy>-1)
				error(EXIT_FAILURE, "*Error*: two context parameters point to Y the last being: ", ctxname);
			else
				posy=index;
	}

	// build the psf at the offset position,
	// which is in the middle of the field
	for (index=0; index < ndim; index++)
		pos[index] = psf->contextoffset[index];
	psf_buildpos(psf, pos, ndim);

	// get something onto the screen
	//psf_print(psf, 1, 1);

	// compute the center of gravity
	getImageStats(psf->maskloc, psf->masksize[0], psf->masksize[1], stats);
	// print the stats
	//QPRINTF(OUTPUT, "\n sum: %.6g", stats[0]);
	//QPRINTF(OUTPUT, "\n cogx %.6g", stats[1]);
	//QPRINTF(OUTPUT, "\n cogy %.6g\n", stats[2]);

	// compute a reasonable size for the re-sampled image;
	// force odd axis lengths
	// allocate memory for the re-sampled image;
	// re-sample the image at the new center and the 'true' pixel size
	// REMARK: not sure whether the parameters are correct...
	//w2=psf->masksize[0];
	//h2=psf->masksize[1];
	w2=(int)ceil((double)psf->masksize[0]*psf->pixstep+1.0);
	h2=(int)ceil((double)psf->masksize[1]*psf->pixstep+1.0);
	w2 = !(w2 % 2) ? w2+1 : w2;
	h2 = !(h2 % 2) ? h2+1 : h2;
	pixstep=1./psf->pixstep;
	QMALLOC(pix2, float, w2*h2);
	vignet_resample(psf->maskloc, psf->masksize[0], psf->masksize[1], pix2, w2, h2, (stats[1]-(float)(psf->masksize[0]/2))*pixstep, (stats[2]-(float)(psf->masksize[1]/2))*pixstep, pixstep);
	//vignet_resample(
	// float *pix1, data for input image
	// int    w1,   width of input image
	// int    h1,   height of input image
	// float *pix2, data of  re-sampled image NEEDS TO BE ALLOCATED
	// int    w2,   width of re-sampled image;   the cutout starts from the middle!!
	// int    h2,   height of re-sampled image;  the cutout starts from the middle!!
	// float  dx,   x-shift, but what shift??
	// float  dy,   y-shift, but what shift??
	// float step2  old vs new steps?? used as "1.0/psf->pixstep"
	// )



	// print the stats again
	//QPRINTF(OUTPUT, "\n\nre-centered psf: %ix%i", w2,h2);
	getImageStats(pix2, w2, h2, stats);
	//QPRINTF(OUTPUT, "\n sum: %.6g", stats[0]);
	//QPRINTF(OUTPUT, "\n cogx %.6g", stats[1]);
	//QPRINTF(OUTPUT, "\n cogy %.6g\n", stats[2]);

	// compute the center pixel
	// take the target values from the re-sampling
	ixmean = w2/2;
	iymean = h2/2;

	// get the bigger of the x/y center index;
	// define a progressively larger area around
	// the center and sum up the pixel values
	icent = ixmean > iymean ? ixmean : iymean;
	index=0;
	sum=0.0;
	while (index<=icent && sum/stats[0]<FRAC){
		// set the start and end values for x/y
		xstart = ixmean-index < 0 ? 0 : ixmean-index;
		ystart = iymean-index < 0 ? 0 : iymean-index;
		xend = ixmean+index+1 < w2 ? ixmean+index+1 : w2;
		yend = iymean+index+1 < h2 ? iymean+index+1 : h2;

		// sum up the psf values
		sum=0.0;
		for (jj=ystart, rpos=pix2+(ystart*w2); jj<yend; jj++, rpos+=w2)
			for (ii=xstart, actpos=rpos+xstart; ii<xend; ii++, actpos++)
				sum += *actpos;
		//QPRINTF(OUTPUT, "index=%i, npix=%i, %i<=x<%i, %i<=y<%i, sum=%.6g\n", index, 2*index+1, xstart, xend, ystart, yend, sum/stats[0]);

		index+=1;
	}

	// set some filter parameters
	// (sum, xstart, xend, ystar, yend are still valid)
	thefilter->convw = xend-xstart;
	thefilter->convh = yend-ystart;
	thefilter->nconv = thefilter->convw*thefilter->convh;

	// allocate memory, copy the normalized values
	// from the psf to the filter
	QMALLOC(thefilter->conv, float, thefilter->nconv);
	fpos = thefilter->conv;
	for (jj=ystart, rpos=pix2+(ystart*w2); jj<yend; jj++, rpos+=w2)
		for (ii=xstart, actpos=rpos+xstart; ii<xend; ii++, actpos++)
			*(fpos++) = *actpos / sum;

	// compute the center of gravity;
	// print the stats
	//getImageStats(thefilter->conv, thefilter->convw, thefilter->convh, stats);
	//QPRINTF(OUTPUT, "\n sum: %.6g", stats[0]);
	//QPRINTF(OUTPUT, "\n cogx %.6g", stats[1]);
	//QPRINTF(OUTPUT, "\n cogy %.6g\n", stats[2]);

	// release memory
	psf_end(psf, NULL);
	free(pix2);

	return RETURN_OK;
}

/**
 * Function: getImageStats
 *
 * Read the convolution mask a PSFEx file
 *
 * @param[in]  pix    the pointer to the pixekl values
 * @param[in]  width  the width of the array
 * @param[in]  height the height of the array
 * @param[out] stats  vector for the stat values
 *
 */
void getImageStats(const float *pix, const int width, const int height, float stats[]){
	long nx, ny, index;
	float *actpix;

	// iterate over pixels
	actpix =  pix;
	stats[0]=0.0;
	stats[1]=0.0;
	stats[2]=0.0;
	for (index=0; index<(long)(width*height); index++){
		// compute nx and ny
		ny = index / width;
		nx = index - ny*width;

		// update the fields
		stats[0] += *actpix;
		stats[1] += (float)nx* *actpix;
		stats[2] += (float)ny* *(actpix++);
	}
	// normalize the COG
	stats[1] =  stats[1] / stats[0];
	stats[2] =  stats[2] / stats[0];
}


/**
 * Function: getfilter
 *
 * Reads in filter. This can be either a convolution mask or an ANN file
 *
 * @param[in] filename the filename of the filter
 */
void	getfilter(const char *filename)
  {
  QCALLOC(thefilter, filterstruct, 1);
  if (getconv(filename) != RETURN_OK && getneurfilter(filename) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: not a suitable filter in ", filename);
  return;
  }

/**
 * Function: getneurfilter
 *
 * Read an ANN RETINA-filter file.
 *
 * @param[in] filename the ANN RETINA-filter filter
 *
 * @return flag indicating success or error
 */
int	getneurfilter(const char *filename)
  {
#define	FILTEST(x) \
        if (x != RETURN_OK) \
	  error(EXIT_FAILURE, "*Error*: incorrect filter header in ", filename)

   catstruct	*fcat;
   tabstruct	*ftab;
   int		ival;

/* We first map the catalog */
  if (!(fcat = read_cat(filename)))
    return RETURN_ERROR;
/* Test if the requested table is present */
  if (!(ftab = name_to_tab(fcat, "BP-ANN", 0)))
    error(EXIT_FAILURE, "*Error*: no BP-ANN info found in ", filename);
  FILTEST(fitsread(ftab->headbuf, "BPTYPE  ", gstr,H_STRING,T_STRING));
  if (strcmp(gstr, "RETINA"))
    error(EXIT_FAILURE, "*Error*: not a retina-filter in ", filename);
  FILTEST(fitsread(ftab->headbuf, "RENAXIS ", &ival ,H_INT, T_LONG));
  if (ival != 2) 
    error(EXIT_FAILURE, "*Error*: not a 2D retina in ", filename);
  FILTEST(fitsread(ftab->headbuf, "RENAXIS1", &thefilter->convw ,H_INT,T_LONG));
  FILTEST(fitsread(ftab->headbuf, "RENAXIS2", &thefilter->convh ,H_INT,T_LONG));
  thefilter->nconv = thefilter->convw*thefilter->convh;
  QMALLOC(thefilter->conv, float, thefilter->nconv);
  thefilter->bpann = loadtab_bpann(ftab, filename);

  close_cat(fcat);
  free_cat(&fcat,1);

  return RETURN_OK;
  }


/**
 * Function: endfilter
 *
 * Terminate filtering or convolution procedures by releasing memory allocated
 * in the global variable 'thefilter'.
 *
 */
void	endfilter()
  {
  QFREE(thefilter->conv);
  if (thefilter->bpann)
    {
    free_bpann(thefilter->bpann);
    thefilter->bpann = NULL;
    }

  QFREE(thefilter);

  return;
  }


/********************************** filter ***********************************/
/*
Switch to the appropriate filtering routine.
*/
void	filter(fieldstruct *field, PIXTYPE *mscan, int y)

  {
  if (thefilter->bpann)
    neurfilter(field, mscan, y);
  else
    convolve(field, mscan, y);

  return;
  }


/******************************** neurfilter *********************************/
/*
Filter a scan line using an artificial retina.
*/
void	neurfilter(fieldstruct *field, PIXTYPE *mscan, int y)

  {
   PIXTYPE	cval;
   float	*pix, resp, sig, threshlow, threshhigh;
   int		i,x, tflag;

  sig = field->backsig;
  if ((tflag = (prefs.nfilter_thresh>0)))
    {
    threshlow = prefs.filter_thresh[0]*field->backsig;
    threshhigh = (prefs.nfilter_thresh<2)?
		BIG : prefs.filter_thresh[1]*field->backsig;
    }
  else
    threshlow = threshhigh = 0.0;	/* To avoid gcc -Wall warnings */

  for (x=0;x<field->width;x++)
    {
    if (tflag)
      {
      cval = PIX(field, x, y);
      if (cval<threshlow || cval>threshhigh)
        {
        *(mscan++) = cval;
        continue;
	}
      }
/*-- Copy the surrounding image area to the retina */
    copyimage(field, thefilter->conv, thefilter->convw, thefilter->convh,
		x, y);
    pix = thefilter->conv;
/*-- Apply a transform of the intensity scale */
    for (i=thefilter->nconv; i--; pix++)
      *pix = *pix>0.0?log(1+*pix/sig):-log(1-*pix/sig);
    play_bpann(thefilter->bpann, thefilter->conv, &resp);
    if (resp>70.0)
      resp = 70.0;
    else if (resp<-70.0)
      resp = -70.0;
/*-- Recover the linear intensity scale */
    *(mscan++) = resp>0.0?sig*(exp(resp)-1):sig*(exp(-resp)-1);
    }

  return;
  }


/******************************* convolve_image ******************************/
/*
Convolve a vignet with an array.
*/
void	convolve_image(fieldstruct *field, float *vig1,
		float *vig2, int width, int height)

  {
   int		mw,mw2,m0,me,m,mx,dmx, y, y0,dy;
   float	*mask, *mscane, *s,*s0, *vig3, *d,*de, mval, val;

/* If no filtering, just return a copy of the input data */
  if (!thefilter || thefilter->bpann)
    {
    memcpy(vig2, vig1, width*height*sizeof(float));
    return;
    }    
  s0 = NULL;				/* To avoid gcc -Wall warnings */
  mw = thefilter->convw;
  mw2 = mw/2;
  memset(vig2, 0, width*height*sizeof(float));
  vig3 = vig2;
  for (y=0; y<height; y++, vig3+=width)
    {
    mscane = vig3+width;
    y0 = y - (thefilter->convh/2);
    if (y0 < 0)
      {
      m0 = -mw*y0;
      y0 = 0;
      }
    else
      m0 = 0;
    if ((dy = height - y0) < thefilter->convh)
      me = mw*dy;
    else
      me = mw*thefilter->convh;
    mask = thefilter->conv+m0;
    for (m = m0, mx = 0; m<me; m++, mx++)
      {
      if (mx==mw)
        mx = 0;
      if (!mx)
        s0 = vig1+width*(y0++);
      if ((dmx = mx-mw2)>=0)
        {
        s = s0 + dmx;
        d = vig3;
        de = mscane - dmx;
        }
      else
        {
        s = s0;
        d = vig3 - dmx;
        de = mscane;
        }

      mval = *(mask++);
      while (d<de)
        *(d++) += ((val = *(s++))>-BIG? mval*val:0.0);
      }
    }

  return;
  }


