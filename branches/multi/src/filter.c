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
#include        "psf.h"
#include	"filter.h"
#include	"wcs/poly.h"
#include	"image.h"

/**************************convolve*******************************************/
/**
 *
 * Function: convolve
 *
 * Convolves a scan line with an array. The row number to convole is given
 * as input. The convoluted line is returned, hence it does NOT work
 * in place.
 * The main iteration is NOT over the image pixels but over the filter
 * values.
 *
 * @author EB, MK
 * @date   April 2014
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
        //QPRINTF(OUTPUT, "convolve\n");

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

/****************************convolve_var*************************************/
/**
 *
 * Function: convolve_var
 *
 * Convolves a scan line with an array using a 2D variable psf. The row number
 * to convole is given as input. The convoluted line is returned, hence it does NOT
 * work in place.
 * The main iteration is NOT over the image pixels but over the filter
 * values.
 *
 * @author MK
 * @date   April 2014
 *
 * @param[in]  field  structure with the scan data
 * @param[out] mscan the convolved image line
 * @param[in]  y     image row to convolve
 *
 */
void convolve_var(fieldstruct *field, PIXTYPE *mscan, int y)

{
  double        pos[POLY_MAXDIM];
  int     mw,mw2,m0,me,m,mx,dmx, y0, dy, sw,sh, index, ii;
  float   *mask;
  PIXTYPE *mscane, *s,*s0, *d,*de, mval;
  double  *dvec;
  varconv *varconvlov;

  float   *actcomp;
  double  *actbasis, actzero;

  varconvlov = thefilter->varpsf; // the variable convolution filter
  sw = field->width;              // width of the image
  sh = field->stripheight;        // height of the data available
  mw = thefilter->convw;          // certainly the width of the kernel
  mw2 = mw/2;                     // half width of the kernel
  mscane = mscan+sw;              // the end of the result vector

  // allocate space for the x/y-dependencies
  if (!varconvlov->basis){
      QMALLOC(varconvlov->basis, double, sw*varconvlov->psf->poly->ncoeff);
  }

  // go over all x values
  for (index=0; index<sw; index++)
    {
      // set the current x/y positions in IRAF notations
      *varconvlov->xpos = (double)index+1.0;
      *varconvlov->ypos = (double)y+1.0;

      // transform to relative position
      for (ii=0; ii<varconvlov->psf->poly->ndim; ii++)
        {
          varconvlov->pos[ii] = (varconvlov->pos[ii] - varconvlov->psf->contextoffset[ii]) / varconvlov->psf->contextscale[ii];
        }

      // evaluate the polynomial
      poly_func(varconvlov->psf->poly, varconvlov->pos);

      // fill into the memory
      for (ii=0, dvec=varconvlov->basis+index; ii<varconvlov->psf->poly->ncoeff; ii++, dvec+=sw)
        {
          *dvec = varconvlov->psf->poly->basis[ii];
         }
   }

  /*for (index=10; index<20; index++)
    {
      //dvec=varconvlov->basis;
      for (ii=0, dvec=varconvlov->basis+index; ii<varconvlov->psf->poly->ncoeff; ii++, dvec+=sw)
        QPRINTF(OUTPUT, " %i: %.5g", ii, *dvec);
        //QPRINTF(OUTPUT, " %i: %.5g", ii, varconvlov->basis[ii*sw+index]);
      QPRINTF(OUTPUT, "\n");
    }
  */

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

  // set mask in the first layer to
  // the end of the kernel
  mask = varconvlov->psf->maskcomp+me;

  // iterate over all usable kernel
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

      // proceed to the
      // next kernel value
      --mask;

      // go over all polynomial components of the kernel
      for (index=0, actcomp=mask; index<varconvlov->psf->poly->ncoeff; index++, actcomp+=thefilter->nconv)
        {

          // make reasonable start
          // and end values
          if ((dmx = mx-mw2)>=0)
            {
              s = s0 + dmx;
              d = mscan;
              de = mscane - dmx;

              // go to the correct position in the polynomials
              actbasis = varconvlov->basis + index*sw;
              //          ^^^^^^^^^^^^^       ^^^^^^
              //       x-position as in "d"   polynomial offset
            }
          else
            {
              s = s0;
              d = mscan - dmx;
              de = mscane;

              // go to the correct position in the polynomials
              actbasis = varconvlov->basis - dmx + index*sw;
              //          ^^^^^^^^^^^^^       ^^^^^^
              //       x-position as in "d"   polynomial offset
            }

          // go over the row
          // and add the contribution
          // to the result vector
          //if (index==0){
          while (d<de)
            {
              *(d++) += *actcomp*(float)*(actbasis++) * *(s++);
              //           ^^^^^^^^^^^^^^^^^^^^^^       ^^^^^^
              //        current kernel layer factor  value of neighbour pixel
            }
          //}
        }
    }

  return;
}

/******************************getconv****************************************/
/**
 *
 * Function: getconv
 *
 * Read the convolution filter from a file. Accepted file formats are an ASCII
 * filter or a PSFEx file.
 *
 * @author EB, MK
 * @date   April 2014
 *
 * @param[in] filename the filename
 *
 * @return flag indicating success or error
 */
int	getconv(const char *filename)
{
  // check whether the filename ends with ".psf"
  if(strlen(filename) > 4 && !strcmp(filename + strlen(filename) - 4, ".psf")){
      // load a psf file and set up the convolution
      return getPSFExConv(filename);
  }
  else
    // load the convolution filter from an ASCII file
    return getASCIIconv(filename);

  // it should NEVER arrive here
  return RETURN_ERROR;
}

/****************************getASCIIconv*************************************/
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

  // assign the filter function
  thefilter->filterFunc = convolve;

  return RETURN_OK;
  }

/*****************************getPSFExConv************************************/
/**
 *
 * Function: getPSFExConv
 *
 * Read the convolution mask a PSFEx file. The psf structure is then re-sampled
 * to the native image sixe size for the convolution. The optimal size is
 * determined on the variable parameters, and the structure is then trimmed
 * to the optimal size and trimmed to that size. The final result is then
 * transferred to the filter structure.
 *
 * @author MK
 * @date   April 2014
 *
 * @param[in] filename the filename
 *
 * @return flag indicating success or error
 *
 * PSFSCALE is the full hub of the values covered,
 * meaning the PSF is valid in each direction for:
 * psf->contextoffset[index] - 0.5 * psf->contextscale[index] < x < psf->contextoffset[index] + 0.5 * psf->contextscale[index]
 * psf->contextscale[index]
 */
int getPSFExConv(const char *filename){
  psfstruct *psf;
  int ndim, index, posx, posy;
  char *ctxname;
  double maxScale=1.0, tmp;
  double pos[POLY_MAXDIM];
  double startPos[POLY_MAXDIM];
  double stepSize;
  psfstruct *res_psf=NULL;
  //float FRAC=0.95;
  int optsize, optsizeold;
  int NTESTGRID=10;

  // load the psf file
  psf=psf_load(filename);

  // TODO:
  // - testing
  // - "psf_resample()" uses "interpolate_pix()" and "make_kernel()" for the sinc
  //   interpolation. Those programs only give results only in the inner area
  //   (no zero-extrapolation), which is sub-optimal. EB wants to provide functions
  //   that make this extrapolation;
  // - when trimming the psf, taking into account the shape (trimming for surface
  //   brightness rather than fractional flux
  // - damping the PSF to zero, as EB suggests

  // make sure there are at most two dimensions
  ndim = psf->poly->ndim;
  if (ndim>2)
    error(EXIT_FAILURE, "*Error*: more than two dimensions in: ", filename);

  // make sure the context variables
  // mark positions by checking the names
  posx=POLY_MAXDIM-1;
  posy=POLY_MAXDIM-1;
  for (index=0; index<ndim; index++){
      ctxname = psf->contextname[index];

      // make sure the context name ends with "_IMAGE"
      if(strlen(ctxname) < 6 || strcmp(ctxname + strlen(ctxname) - 6, "_IMAGE")){
        error(EXIT_FAILURE, "*Error*: no image coordinate in context: ", ctxname);
      }

      // check for the allowed x-values; make sure only one context points to x
      if (!strncmp(ctxname, "X", 1) || !strncmp(ctxname, "XMIN", 4) || !strncmp(ctxname, "XMAX", 4) || !strncmp(ctxname, "XPEAK", 5)){
        if (posx<2) {
          error(EXIT_FAILURE, "*Error*: two context parameters point to X, the last being: ", ctxname);
        }
        else {
          posx=index;
        }
      }

      // check for the allowed y-values; make sure only one context points to y
      if (!strncmp(ctxname, "Y", 1) || !strncmp(ctxname, "YMIN", 4) || !strncmp(ctxname, "YMAX", 4) || !strncmp(ctxname, "YPEAK", 5)){
        if (posy<2) {
          error(EXIT_FAILURE, "*Error*: two context parameters point to Y the last being: ", ctxname);
        }
        else {
          posy=index;
        }
      }
  }

  // re-sample the psf to the image pixels size
  res_psf = psf_resample(psf, 1./psf->pixstep);

  // get the maximum scale at all dimensions
  for (index=0; index < ndim; index++)
    maxScale = maxScale < psf->contextscale[index] ? psf->contextscale[index] : maxScale ;

  // get the set size in the longest direction
  stepSize = maxScale/(double)NTESTGRID;

  // get the maximum scale at all dimensions
  for (index=0; index < ndim; index++) {
      // get the starting position in each dimensions
      tmp = floor(psf->contextscale[index]/stepSize);
      startPos[index] = psf->contextoffset[index] - stepSize*tmp/2.0;
      //QPRINTF(OUTPUT, "\n dim: %i, startPos: %.5g, stepsize: %5g\n", index, startPos[index], stepSize);

  }
  //QPRINTF(OUTPUT, "\n cogy %.6g\n", stats[2]);
  //exit(0);

  optsize=0;
  optsizeold=0;
  if (ndim==0){
      //QPRINTF(OUTPUT, "\nNo need to evaluate something....\n\n\n");
      //psf->build_flag =0;
      //psf_buildpos(psf, pos, ndim);
      //optPSFExSize(psf, &optsizeold);

      res_psf->build_flag =0;
      psf_buildpos(res_psf, pos, ndim);
      optPSFExSize(res_psf, &optsize);
  }
  else if (ndim==1){
      for (pos[0]=startPos[index]; pos[0] <= psf->contextoffset[0]+0.5*psf->contextscale[0]; pos[0]+=stepSize){
          //QPRINTF(OUTPUT, "1D evaluation at (%.5g)\n", pos[0]);
          //psf->build_flag =0;
          //psf_buildpos(psf, pos, ndim);
          //optPSFExSize(psf, &optsizeold);

          res_psf->build_flag =0;
          psf_buildpos(res_psf, pos, ndim);
          optPSFExSize(res_psf, &optsize);
      }
  }
  else if (ndim==2){
      for (pos[0]=startPos[0]; pos[0] <= psf->contextoffset[0]+0.5*psf->contextscale[0]; pos[0]+=stepSize){
          for (pos[1]=startPos[1]; pos[1] <= psf->contextoffset[1]+0.5*psf->contextscale[1]; pos[1]+=stepSize){
             // QPRINTF(OUTPUT, "2D evaluation at (%.5g,%.5g)\n", pos[0], pos[1]);
              //psf->build_flag =0;
              //psf_buildpos(psf, pos, ndim);
              //optPSFExSize(psf, &optsizeold);

              //QPRINTF(OUTPUT, "start size %i", optsize);
              res_psf->build_flag =0;
              psf_buildpos(res_psf, pos, ndim);
              optPSFExSize(res_psf, &optsize);
              //QPRINTF(OUTPUT, " end size %i\n", optsize);
              //psf_print(new_psf, 0, 1);
          }
      }
  }

  // trim the psf to the useful size
  psf_trim(res_psf, optsize);
  psf_normalize(res_psf);

  // release memory
  psf_end(psf, NULL);
  psf=NULL;

  // transfer the size values
  thefilter->convw = res_psf->masksize[0];
  thefilter->convh = res_psf->masksize[1];
  thefilter->nconv = thefilter->convw*thefilter->convh;

  if (ndim)
    {
      float *fpos1, *fpos2;
      //float sum=0.;
      //float stats[4];

      // save the psf in the filter structure
      QMALLOC(thefilter->varpsf, varconv, 1);
      thefilter->varpsf->psf = res_psf;
      thefilter->varpsf->basis = NULL;

      // connect x/y positions to the context indices in the psf
      thefilter->varpsf->xpos = &thefilter->varpsf->pos[posx];
      thefilter->varpsf->ypos = &thefilter->varpsf->pos[posy];

      // just to avoid a core dump
      QMALLOC(thefilter->conv, float, thefilter->nconv);
      fpos1 = thefilter->conv;
      fpos2 = res_psf->maskcomp;
      for (index=0; index<thefilter->nconv; index++)
        *(fpos1++) = *(fpos2++);

      // assign the filter function
      thefilter->filterFunc=convolve_var;

      //fpos1 = thefilter->conv;
      //for (index=0; index<thefilter->nconv; index++)
      //  sum +=  *(fpos1++);
      //QPRINTF(OUTPUT, "\nDimensions %ix%i, sum: %g, \n", thefilter->convw, thefilter->convh, sum);
      // compute the center of gravity;
      // print the stats
      //getImageStats(thefilter->conv, thefilter->convw, thefilter->convh, stats);
      //QPRINTF(OUTPUT, "\n summmm: %.6g", stats[0]);
      //QPRINTF(OUTPUT, "\n cogxxx %.6g %i ", stats[1], thefilter->convw);
      //QPRINTF(OUTPUT, "\n cogyyy %.6g %i\n", stats[2], thefilter->convh);
    }
  else
    {
      // the psf is constant, copy its content
      // to the filter structure

      float *fpos1, *fpos2;
      //float sum=0.;
      //float stats[4];

      // allocate memory for the convolution mask
      // and transfer the values from the psf
      QMALLOC(thefilter->conv, float, thefilter->nconv);
      fpos1 = thefilter->conv;
      fpos2 = res_psf->maskcomp;
      for (index=0; index<thefilter->nconv; index++)
        *(fpos1++) = *(fpos2++);

      // assign the filter function
      thefilter->filterFunc=convolve;

      //fpos = thefilter->conv;
      //for (index=0; index<thefilter->nconv; index++)
      //  sum +=  *(fpos++);
      //QPRINTF(OUTPUT, "\nDimensions %ix%i, sum: %g, \n", thefilter->convw, thefilter->convh, sum);
      // compute the center of gravity;
      // print the stats
      //getImageStats(thefilter->conv, thefilter->convw, thefilter->convh, stats);
      //QPRINTF(OUTPUT, "\n summm: %.6g", stats[0]);
      //QPRINTF(OUTPUT, "\n cogxx %.6g %i ", stats[1], thefilter->convw);
      //QPRINTF(OUTPUT, "\n cogyy %.6g %i\n", stats[2], thefilter->convh);

      // release memory
      psf_end(res_psf, NULL);
      res_psf=NULL;
    }

  return RETURN_OK;
}

/****************************optPSFExSize*************************************/
/**
 *
 * Function: optPSFExSize
 *
 * Determine the optimal size for a psf-struct, meaning the size than encloses
 * more than a certain fraction of the flux.
 *
 * @author MK
 * @date   April 2014
 *
 * @param[in]     psf      the psf-structure
 * @param[in,out] optsize  the optimal (half-)size
 *
 * @return flag indicating success or error
 */
int optPSFExSize(const psfstruct *psf, int *optsize){
  int ixmean, iymean, icent, inloop=0;
  int xstart, xend, ii, ystart, yend, jj;
  float sum;
  float *rpos, *actpos;
  float stats[4];
  float FRAC=0.95;

  // compute the center of gravity and total sum and so on
  getImageStats(psf->maskloc, psf->masksize[0], psf->masksize[1], stats);
  //QPRINTF(OUTPUT, "name: %s, sum: %.9g, cog: (%.6g,%.6g) size (%i,%i)\n", psf->name, stats[0], stats[1], stats[2], psf->masksize[0], psf->masksize[1]);

  // compute the center pixel
  // take the target values from the re-sampling
  ixmean = psf->masksize[0]/2;
  iymean = psf->masksize[1]/2;

  // get the bigger of the x/y center index;
  // define a progressively larger area around
  // the center and sum up the pixel values
  // Remark: A quadratic PSF is assumed. One could improve
  //         upon that by assuming a rectangular PSF progressing
  //         asymetrical according to the side ratios.
  icent = ixmean > iymean ? ixmean : iymean;
  sum=0.0;
  while (*optsize<=icent && sum/stats[0]<FRAC){
      // set the start and end values for x/y
      xstart = ixmean-*optsize < 0 ? 0 : ixmean-*optsize;
      ystart = iymean-*optsize < 0 ? 0 : iymean-*optsize;
      xend = ixmean+*optsize+1 < psf->masksize[0] ? ixmean+*optsize+1 : psf->masksize[0];
      yend = iymean+*optsize+1 < psf->masksize[1] ? iymean+*optsize+1 : psf->masksize[1];

      // sum up the psf values
      sum=0.0;
      for (jj=ystart, rpos=psf->maskloc+(ystart*psf->masksize[1]); jj<yend; jj++, rpos+=psf->masksize[1]){
        for (ii=xstart, actpos=rpos+xstart; ii<xend; ii++, actpos++){
          sum += *actpos;
        }
      }

      // increment the optsize
      *optsize+=1;

      // mark the loop passing
      inloop=1;
  }

  // correct for the last
  // increment in the loop
  if (inloop)
    *optsize-=1;

  return RETURN_OK;
}

void optPSFExSizeOld(const psfstruct *psf){
	int ixmean, iymean, icent, index;
	int xstart, xend, ii, ystart, yend, jj;
	float sum;
	float *rpos, *actpos;
	float stats[4];
	float FRAC=0.95;

	// compute the center of gravity and total sum and so on
	getImageStats(psf->maskloc, psf->masksize[0], psf->masksize[1], stats);
	QPRINTF(OUTPUT, "name: %s, sum: %.9g, cog: (%.6g,%.6g) size (%i,%i)\n", psf->name, stats[0], stats[1], stats[2], psf->masksize[0], psf->masksize[1]);

	// compute the center pixel
	// take the target values from the re-sampling
	ixmean = psf->masksize[0]/2;
	iymean = psf->masksize[1]/2;

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
		xend = ixmean+index+1 < psf->masksize[0] ? ixmean+index+1 : psf->masksize[0];
		yend = iymean+index+1 < psf->masksize[1] ? iymean+index+1 : psf->masksize[1];

		// sum up the psf values
		sum=0.0;
		for (jj=ystart, rpos=psf->maskloc+(ystart*psf->masksize[1]); jj<yend; jj++, rpos+=psf->masksize[1])
			for (ii=xstart, actpos=rpos+xstart; ii<xend; ii++, actpos++)
				sum += *actpos;

		// enhance the index
		index+=1;
	}
}

/*************************getConstPSFExConv***********************************/
/**
 * Function: getConstPSFExConv
 *
 * Determines a constant PSF from a psf-file. For variable PSF's the psf at
 * the center is evaluated and transferred to the filter structure.
 *
 * @author MK
 * @date   April 2014
 *
 * @param[in] filename the filename
 *
 * @return flag indicating success or error
 */
int getConstPSFExConv(const char *filename){
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
  if (ndim>2){
    error(EXIT_FAILURE, "*Error*: more than two dimensions in: ", filename);
  }

  // make sure the context variables
  // mark positions by checking the names
  posx=-1;
  posy=-1;
  for (index=0; index<ndim; index++){
      ctxname = psf->contextname[index];

      // make sure the context name ends with "_IMAGE"
      if(strlen(ctxname) < 6 || strcmp(ctxname + strlen(ctxname) - 6, "_IMAGE")){
        error(EXIT_FAILURE, "*Error*: no image coordinate in context: ", ctxname);
      }

      // check for the allowed x-values; make sure only one context points to x
      if (!strncmp(ctxname, "X", 1) || !strncmp(ctxname, "XMIN", 4) || !strncmp(ctxname, "XMAX", 4) || !strncmp(ctxname, "XPEAK", 5)){
        if (posx>-1) {
          error(EXIT_FAILURE, "*Error*: two context parameters point to X, the last being: ", ctxname);
        }
        else {
          posx=index;
        }
      }

      // check for the allowed y-values; make sure only one context points to y
      if (!strncmp(ctxname, "Y", 1) || !strncmp(ctxname, "YMIN", 4) || !strncmp(ctxname, "YMAX", 4) || !strncmp(ctxname, "YPEAK", 5)){
        if (posy>-1) {
          error(EXIT_FAILURE, "*Error*: two context parameters point to Y the last being: ", ctxname);
        }
        else {
          posy=index;
        }
      }
  }

  // build the psf at the offset position,
  // which is in the middle of the field
  for (index=0; index < ndim; index++){
    pos[index] = psf->contextoffset[index];
  }
  psf_buildpos(psf, pos, ndim);

  // get something onto the screen
  //psf_print(psf, 1, 1);

  // compute the center of gravity
  getImageStats(psf->maskloc, psf->masksize[0], psf->masksize[1], stats);
  //QPRINTF(OUTPUT, "sum: %.9g, cog: (%.6g,%.6g) size (%i,%i) at pos: (%g,%g)\n", stats[0], stats[1], stats[2], psf->masksize[0], psf->masksize[1], pos[0], pos[1]);
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
      for (jj=ystart, rpos=pix2+(ystart*w2); jj<yend; jj++, rpos+=w2){
        for (ii=xstart, actpos=rpos+xstart; ii<xend; ii++, actpos++){
          sum += *actpos;
        }
      }
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
  for (jj=ystart, rpos=pix2+(ystart*w2); jj<yend; jj++, rpos+=w2){
    for (ii=xstart, actpos=rpos+xstart; ii<xend; ii++, actpos++){
      *(fpos++) = *actpos / sum;
    }
  }

  // compute the center of gravity;
  // print the stats
  getImageStats(thefilter->conv, thefilter->convw, thefilter->convh, stats);
  QPRINTF(OUTPUT, "\n sum: %.6g", stats[0]);
  QPRINTF(OUTPUT, "\n cogx %.6g %i ", stats[1], thefilter->convw);
  QPRINTF(OUTPUT, "\n cogy %.6g %i\n", stats[2], thefilter->convh);

  // release memory
  psf_end(psf, NULL);
  free(pix2);

  // assign the filter function
  thefilter->filterFunc=convolve_var;

  return RETURN_OK;
}

/************************getImageStats****************************************/
/**
 *
 * Function: getImageStats
 *
 * Compute some simple image statistics: the center of gravity in x and y,
 * and the sum of all pixel values.
 *
 * @author MK
 * @date   April 2014
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


/***************************getfilter*****************************************/
/**
 *
 * Function: getfilter
 *
 * Reads in filter. This can be either a convolution mask or a psf file
 * or an ANN file
 *
 * @author EB, MK
 * @date   April 2014
 *
 * @param[in] filename the filename of the filter
 */
void	getfilter(const char *filename)
  {
  // allocate the filter structure
  // and explicitly initialize the pointers
  QCALLOC(thefilter, filterstruct, 1);
  thefilter->varpsf     = NULL;
  thefilter->bpann      = NULL;
  thefilter->filterFunc = NULL;
  thefilter->conv       = NULL;

  // load in a convolution filter or a neural network filter
  if (getconv(filename) != RETURN_OK && getneurfilter(filename) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: not a suitable filter in ", filename);

  return;
  }

/***************************getneurfilter*************************************/
/**
 *
 * Function: getneurfilter
 *
 * Read an ANN RETINA-filter file.
 *
 * @author EB, MK
 * @date   April 2014
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

  // assign the filter function
  thefilter->filterFunc=neurfilter;

  return RETURN_OK;
  }


/************************endfilter********************************************/
/**
 *
 * Function: endfilter
 *
 * Terminate filtering or convolution procedures by releasing memory allocated
 * in the global variable 'thefilter'.
 *
 * @author EB, MK
 * @date   April 2014
 */
void	endfilter()
  {
  if (thefilter->conv)
    {
      // release convolution data memory
      QFREE(thefilter->conv);
      thefilter->conv = NULL;
    }

  if (thefilter->bpann)
    {
      // release ANN memory
      free_bpann(thefilter->bpann);
      thefilter->bpann = NULL;
    }

  if (thefilter->varpsf)
    {
      // release variable psf memory
      psf_end(thefilter->varpsf->psf, NULL);
      thefilter->varpsf->psf = NULL;
      if (thefilter->varpsf->basis)
        {
          QFREE(thefilter->varpsf->basis);
          thefilter->varpsf->basis=NULL;
        }
      QFREE(thefilter->varpsf);
      thefilter->varpsf=NULL;
    }

  // release filter memory
  QFREE(thefilter);
  thefilter=NULL;

  return;
  }


/**************************filter*********************************************/
/**
 *
 * Function: filter
 *
 * General filter function. Switches to the appropriate routine,
 * which is either convolution or neural newtork filtering.
 * Remark: The function is rather superfluous now. The calls to the functions
 *         could be replaced with the line below.
 *
 * @author EB, MK
 * @date   April 2014
 *
 * @param[in]  field  structure with the scan data
 * @param[out] mscan  the convolved image line
 * @param[in]  y      image row to convolve
 *
 */
void filter(fieldstruct *field, PIXTYPE *mscan, int y)
{
  /* the old version
  if (thefilter->bpann)
    neurfilter(field, mscan, y);
  else
    convolve(field, mscan, y);
 */

  // TODO: needs to be ckecked for
  // ANN filter (neural network filters)

  // new version with the function assigned
  (*thefilter->filterFunc)(field, mscan, y);

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


