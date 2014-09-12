/*
*				psf.c
*
* Fit a PSF model to an image.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1998-2014 IAP/CNRS/UPMC
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
*	Last modified:		11/09/2014
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
#include	"check.h"
#include	"field.h"
#include	"filter.h"
#include	"image.h"
#include	"wcs/poly.h"
#include	"subimage.h"
#include	"psf.h"

static float	interpolate_pixZERO(float *posin, float *pix, int *naxisn,
		interpenum interptype);
static float    interpolate_pix(float *posin, float *pix, int *naxisn,
                interpenum interptype);

static void	make_kernel(float pos, float *kernel, interpenum interptype);

/*------------------------------- variables ---------------------------------*/
const int	interp_kernwidth_psf[5]={1,2,4,6,8};


extern keystruct	obj2key[];
extern objstruct	outobj;

/********************************* psf_init **********************************/
/*
Allocate memory and stuff for the PSF-fitting.
*/
void	psf_init(psfstruct *psf)
  {
  QMALLOC(thepsfit, psfitstruct, 1);
  QMALLOC(thepsfit->x, double, prefs.psf_npsfmax);
  QMALLOC(thepsfit->y, double, prefs.psf_npsfmax);
  QMALLOC(thepsfit->flux, float, prefs.psf_npsfmax);
  QMALLOC(thepsfit->fluxerr, float, prefs.psf_npsfmax);
  QMALLOC(ppsfit, psfitstruct, 1); /*?*/
  QMALLOC(ppsfit->x, double, prefs.psf_npsfmax);
  QMALLOC(ppsfit->y, double, prefs.psf_npsfmax);
  QMALLOC(ppsfit->flux, float, prefs.psf_npsfmax);
  QMALLOC(ppsfit->fluxerr, float, prefs.psf_npsfmax);

  return;
  }  


/********************************* psf_end ***********************************/
/*
Free memory occupied by the PSF-fitting stuff.
*/
void	psf_end(psfstruct *psf, psfitstruct *psfit)
  {
   int	d, ndim;

  ndim = psf->poly->ndim;
  for (d=0; d<ndim; d++)
    free(psf->contextname[d]);
  free(psf->context);
  free(psf->contextname);
  free(psf->contextoffset);
  free(psf->contextscale);
  free(psf->contexttyp);
  poly_end(psf->poly);
  free(psf->maskcomp);
  free(psf->maskloc);
  free(psf->masksize);
  free(psf);

  if (psfit)
    {
    free(psfit->x);
    free(psfit->y);
    free(psfit->flux);
    free(psfit->fluxerr);
    free(psfit);
    }

  return;
  }

/****************************psf_load*****************************************/
/**
 * Function: psf_load
 *
 * Read the PSF data from a FITS file in the PSFEx format.
 *
 * @param[in] filename name of the PSFEx file
 *
 * @return the generated PSF structure
 */
psfstruct	*psf_load(char *filename)
  {
   static obj2struct	saveobj2;
   psfstruct		*psf;
   catstruct		*cat;
   tabstruct		*tab;
   keystruct		*key;
   char			*head, *ci,*co;
   int			deg[POLY_MAXDIM], group[POLY_MAXDIM], ndim, ngroup,
			i,k;

/* Open the cat (well it is not a "cat", but simply a FITS file */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: PSF file not found: ", filename);

/* OK, we now allocate memory for the PSF structure itself */
  QCALLOC(psf, psfstruct, 1);

/* Store a short copy of the PSF filename */
  if ((ci=strrchr(filename, '/')))
    strcpy(psf->name, ci+1);
  else
    strcpy(psf->name, filename);

  if (!(tab = name_to_tab(cat, "PSF_DATA", 0)))
    error(EXIT_FAILURE, "*Error*: PSF_DATA table not found in catalog ",
	filename);

  head = tab->headbuf;

/*-- Dimension of the polynomial */
  if (fitsread(head, "POLNAXIS", &ndim, H_INT,T_LONG) == RETURN_OK
	&& ndim)
    {
/*-- So we have a polynomial description of the PSF variations */
    if (ndim > POLY_MAXDIM)
        {
        sprintf(gstr, "*Error*: The POLNAXIS parameter in %s exceeds %d",
		psf->name, POLY_MAXDIM);
        error(EXIT_FAILURE, gstr, "");
        }

    QMALLOC(psf->contextname, char *, ndim);
    QMALLOC(psf->context, double *, ndim);
    QMALLOC(psf->contexttyp, t_type, ndim);
    QMALLOC(psf->contextoffset, double, ndim);
    QMALLOC(psf->contextscale, double, ndim);

/*-- We will have to use the flagobj2 struct, so we first save its content */
    saveobj2 = flagobj2;
/*-- flagobj2 is used as a FLAG array, so we initialize it to 0 */
    memset(&flagobj2, 0, sizeof(flagobj2));
    for (i=0; i<ndim; i++)
      {
/*---- Polynomial groups */
      sprintf(gstr, "POLGRP%1d", i+1);
      if (fitsread(head, gstr, &group[i], H_INT,T_LONG) != RETURN_OK)
        goto headerror;

/*---- Contexts */
      QMALLOC(psf->contextname[i], char, 80);
      sprintf(gstr, "POLNAME%1d", i+1);
      if (fitsread(head,gstr,psf->contextname[i],H_STRING,T_STRING)!=RETURN_OK)
        goto headerror;
      if (*psf->contextname[i]==(char)':')
/*------ It seems we're facing a FITS header parameter */
        psf->context[i] = NULL;	/* This is to tell we'll have to load */
				/* a FITS header context later on */
      else
/*------ The context element is a dynamic object parameter */
        {
        if ((k = findkey(psf->contextname[i], (char *)obj2key,
		sizeof(keystruct)))==RETURN_ERROR)
          {
          sprintf(gstr, "*Error*: %s CONTEXT parameter in %s unknown",
		psf->contextname[i], psf->name);
          error(EXIT_FAILURE, gstr, "");
          }
        key = obj2key+k;
        psf->context[i] = key->ptr;	/* !CHECK should depend on channel */
        psf->contexttyp[i] = key->ttype;
/*------ Declare the parameter "active" to trigger computation by SExtractor */
        *((char *)key->ptr) = (char)'\1';
        }
/*---- Scaling of the context parameter */
      sprintf(gstr, "POLZERO%1d", i+1);
      if (fitsread(head, gstr, &psf->contextoffset[i], H_EXPO, T_DOUBLE)
		!=RETURN_OK)
        goto headerror;
      sprintf(gstr, "POLSCAL%1d", i+1);
      if (fitsread(head, gstr, &psf->contextscale[i], H_EXPO, T_DOUBLE)
		!=RETURN_OK)
        goto headerror;
      }

/*-- Number of groups */
    if (fitsread(head, "POLNGRP ", &ngroup, H_INT, T_LONG) != RETURN_OK)
      goto headerror;

    for (i=0; i<ngroup; i++)
      {
/*---- Polynomial degree for each group */
      sprintf(gstr, "POLDEG%1d", i+1);
      if (fitsread(head, gstr, &deg[i], H_INT,T_LONG) != RETURN_OK)
        goto headerror;
      }

    psf->poly = poly_init(group, ndim, deg, ngroup);

/*-- Restore previous flagobj2 content */
    flagobj2 = saveobj2;
    }
  else
    {
/*-- This is a simple, constant PSF */
    psf->poly = poly_init(group, 0, deg, 0);
    psf->context = NULL;
    }

/* Dimensionality of the PSF mask */
  if (fitsread(head, "PSFNAXIS", &psf->maskdim, H_INT, T_LONG) != RETURN_OK)
    goto headerror;
  if (psf->maskdim<2 || psf->maskdim>3)
    error(EXIT_FAILURE, "*Error*: wrong dimensionality for the PSF "
	"mask in ", filename);
  QMALLOC(psf->masksize, int, psf->maskdim);
  for (i=0; i<psf->maskdim; i++)
    psf->masksize[i] = 1;
  psf->masknpix = 1;
  for (i=0; i<psf->maskdim; i++)
    {
    sprintf(gstr, "PSFAXIS%1d", i+1);
    if (fitsread(head, gstr, &psf->masksize[i], H_INT,T_LONG) != RETURN_OK)
      goto headerror;
    psf->masknpix *= psf->masksize[i];
    }

/* PSF FWHM: defaulted to 3 pixels */
 if (fitsread(head, "PSF_FWHM", &psf->fwhm, H_FLOAT,T_DOUBLE) != RETURN_OK)
    psf->fwhm = 3.0;

/* PSF oversampling: defaulted to 1 */
  if (fitsread(head, "PSF_SAMP", &psf->pixstep,H_FLOAT,T_FLOAT) != RETURN_OK
	|| psf->pixstep <= 0.0)
    psf->pixstep = 1.0;

/* Load the PSF mask data */
  key = read_key(tab, "PSF_MASK");
  psf->maskcomp = key->ptr;

  QMALLOC(psf->maskloc, float, psf->masksize[0]*psf->masksize[1]);

/* But don't touch my arrays!! */
  blank_keys(tab);

  free_cat(&cat, 1);

  return psf;

headerror:
  error(EXIT_FAILURE, "*Error*: Incorrect or obsolete PSF data in ", filename);
  return NULL;
  }


/***************************** psf_readcontext *******************************/
/*
Read the PSF context parameters in the FITS header.
*/
void	psf_readcontext(psfstruct *psf, fieldstruct *field)
  {
   static double	contextval[POLY_MAXDIM];
   int			i, ndim;

  ndim = psf->poly->ndim;
  for (i=0; i<ndim; i++)
    if (!psf->context[i])
      {
      psf->context[i] = &contextval[i];
      psf->contexttyp[i] = T_DOUBLE;
      if (fitsread(field->tab->headbuf, psf->contextname[i]+1, &contextval[i],
		H_FLOAT,T_DOUBLE) == RETURN_ERROR)
        {
        sprintf(gstr, "*Error*: %s parameter not found in the header of ",
		psf->contextname[i]+1);
        error(EXIT_FAILURE, gstr, field->rfilename);
        }
      }

  return;
  }

/*****************************psf_print***************************************/
/**
 *
 * Function: psf_print
 *
 * Prints information on a PSFEx structure. Depending on the input only
 * metadata or also the polynomial data and, if available, the evaluated
 * PSF in the storage is printed.
 *
 * @author MK
 * @date   April 2014
 *
 * @param[in] psf        the PSF struct
 * @param[in] printdata  print the polynomial data?
 * @param[in] printloc   print the evaluated psf?
 *
 */
void psf_print(const psfstruct *psf, const int printdata, const int printloc){
  int index, ndim, ngroup;
  QPRINTF(OUTPUT, "PSF name: %s\n", psf->name);
  QPRINTF(OUTPUT, "maskdim: %i\n", psf->maskdim);
  for (index=0; index < psf->maskdim; index++)
    QPRINTF(OUTPUT, "masksize[%i]: %i\n", index, psf->masksize[index]);
  ndim   = psf->poly->ndim;
  ngroup = psf->poly->ngroup;
  QPRINTF(OUTPUT, "number of dimensions: %i\n", ndim);
  for (index=0; index < ndim; index++){
      if (!psf->context[index]){
          QPRINTF(OUTPUT, "context[%i] not (yet) defined!\n", index);
      }
      else{
          QPRINTF(OUTPUT, "context[i].name:   %s\n", psf->contextname[index]);
          QPRINTF(OUTPUT, "context[i].offset: %.3g\n", psf->contextoffset[index]);
          QPRINTF(OUTPUT, "context[i].scale:  %.3g\n", psf->contextscale[index]);
      }
  }

  for (index=0; index < ndim; index++){
      QPRINTF(OUTPUT, "group index %i:   %i\n", index, psf->poly->group[index]);
  }

  for (index=0; index < ngroup; index++){
      QPRINTF(OUTPUT, "degree index %i:   %i\n", index, psf->poly->degree[index]);
  }

  if (printdata) {
      QPRINTF(OUTPUT, "masknpix: %i\n", psf->masknpix);
      for (index=0; index < psf->masknpix; index++)
        QPRINTF(OUTPUT, " %.4g", psf->maskcomp[index]);
  }
  QPRINTF(OUTPUT, "\nFWHM: %.3g\n", psf->fwhm);
  QPRINTF(OUTPUT, "pixstep: %.3g\n", psf->pixstep);
  QPRINTF(OUTPUT, "build_flag: %i\n", psf->build_flag);
  if (printloc && psf->build_flag){
      QPRINTF(OUTPUT, "local psf:\n");
      for (index=0; index<psf->masksize[0]*psf->masksize[1]; index++)
        QPRINTF(OUTPUT, " %.4g", psf->maskloc[index]);
  }
}


/******************************psf_resample***********************************/
/**
 *
 * Function: psf_resample
 *
 * "factor" is the ratio of the new pixel size over the old pixel size, meaning
 * "factor>1.0" means the first two dimension of the new PSF struct get smaller,
 * "factor<1.0" means the first two dimension of the new PSF struct get bigger!
 * As a consequence, "factor=1.0/psf->pixstep" re-samples the PSF to the pixel
 * size of the image it was derived from.
 *
 * @author MK
 * @date   April 2014
 *
 * @param[in] psf     the PSF structure
 * @param[in] factor  the re-sampling factor
 *
 * @return re-sampled PSF structure
 */
psfstruct *psf_resample(const psfstruct *psf, const float factor){
  psfstruct *res_psf;
  int w2, h2;
  float xcin, ycin, xcout, ycout;
  float *pix1, *pix2, flux, factorsquare;
  float	posin[2], posout[2], dnaxisn[2];
  int d, i, index, ii;
  int	deg[POLY_MAXDIM], group[POLY_MAXDIM], ndim, ngroup;

  // square of the resample factor
  factorsquare = factor*factor;

  // determine the new size
  w2=(int)ceil((double)psf->masksize[0]/factor);
  h2=(int)ceil((double)psf->masksize[1]/factor);
  w2 = !(w2 % 2) ? w2+1 : w2;
  h2 = !(h2 % 2) ? h2+1 : h2;
  //QPRINTF(OUTPUT, "size1: (%i,%i) --> size2: (%i,%i)\n\n", psf->masksize[0], psf->masksize[1], w2, h2);

  // compute the image centers using FITS convention
  xcout = (float)(w2/2) + 1.0;	/* FITS convention */
  ycout = (float)(h2/2) + 1.0;	/* FITS convention */
  xcin = (float)(psf->masksize[0]/2) + 1.0;			/* FITS convention */
  ycin = (float)(psf->masksize[1]/2) + 1.0;			/* FITS convention */

  // allocate memory for the PSF structure itself
  QCALLOC(res_psf, psfstruct, 1);

  // define and set the dimensions
  // in the new PSF
  res_psf->maskdim = psf->maskdim;
  QMALLOC(res_psf->masksize, int, res_psf->maskdim);
  res_psf->masksize[0]=w2;
  res_psf->masksize[1]=h2;
  res_psf->masksize[2]=psf->masksize[2];
  res_psf->masknpix = res_psf->masksize[0]*res_psf->masksize[1]*res_psf->masksize[2];

  // allocate memory for the data arrays
  QCALLOC(res_psf->maskcomp, float, res_psf->masknpix);
  QCALLOC(res_psf->maskloc, float, res_psf->masksize[0]*res_psf->masksize[1]);

  // copy some static metadata over
  strcpy(res_psf->name, "<REBIN>");
  strcat(res_psf->name, psf->name);
  res_psf->pixstep    = factor*psf->pixstep;
  res_psf->fwhm       = psf->fwhm/factor;
  res_psf->build_flag = 0;

  // allocate memory for dynamic metadata
  ndim   = psf->poly->ndim;
  ngroup = psf->poly->ngroup;
  QMALLOC(res_psf->contextname, char *, ndim);
  QMALLOC(res_psf->context, double *, ndim);
  QMALLOC(res_psf->contexttyp, t_type, ndim);
  QMALLOC(res_psf->contextoffset, double, ndim);
  QMALLOC(res_psf->contextscale, double, ndim);

  // copy over the dynamic metadata
  for (index=0; index<ndim; index++){
      // copy context names
      QMALLOC(res_psf->contextname[index], char, 80);
      strcpy(res_psf->contextname[index], psf->contextname[index]);

      // not sure whether it is so easy
      res_psf->context[index]       = psf->context[index];
      res_psf->contexttyp[index]    = psf->contexttyp[index];

      // copy offsets and scales
      res_psf->contextoffset[index] = psf->contextoffset[index];
      res_psf->contextscale[index]  = psf->contextscale[index];

      // copy group numbers
      // the "+1" is necessary since there i "-1"
      // in "poly_init()" down
      group[index] = psf->poly->group[index]+1;
  }

  // copy degree numbers
  for (index=0; index<ngroup; index++)
    deg[index] = psf->poly->degree[index];

  // initialize the polynomials
  res_psf->poly = poly_init(group, ndim, deg, ngroup);

  pix1 = psf->maskcomp;
  pix2 = res_psf->maskcomp;
  for (ii=0; ii<res_psf->masksize[2]; ii++, pix1+=psf->masksize[0]*psf->masksize[1]){
      for (posout[1]=1.0; posout[1] < (float)res_psf->masksize[1]+0.5; posout[1]+=1.0){
          for (posout[0]=1.0; posout[0] < (float)res_psf->masksize[0]+0.5; posout[0]+=1.0){
              //for (i=w2*h2; i--;){
              //QPRINTF(OUTPUT, " (%.5g,%.5g)", posout[0], posout[1]);

              posin[0] = (posout[0] - xcout)*factor + xcin;
              posin[1] = (posout[1] - ycout)*factor + ycin;
              //QPRINTF(OUTPUT, " (%.5g,%.5g)<->(%.5g,%.5g %.2g)", posout[0], posout[1], posin[0], posin[1], factor);
              //flux += ((*(pixout++) = interpolate_pix(posin, psf->maskloc, psf->masksize, INTERP_LANCZOS3)));
              //flux += ((*(pix2++) = interpolate_pix(posin, psf->maskcomp, psf->masksize, INTERP_LANCZOS3)));
              //*(pix2++) = interpolate_pix(posin, psf->maskcomp, psf->masksize, INTERP_LANCZOS3);

              // fill in with the Lanczos-resampled value times factor^2 to preserve normalization
              //*(pix2++) = interpolate_pix(posin, pix1, psf->masksize, INTERP_LANCZOS3)*factorsquare;
              *(pix2++) = interpolate_pixZERO(posin, pix1, psf->masksize, INTERP_LANCZOS3)*factorsquare;

              // for (posout[0],posout[1])=(x,y) iterates over all x,y values,
              // varying over x more rapidly:
              // (posout[0],posout[1]) = (1,1), (2,1), (3,1),....., (2,1), (2,2), ...
              //    for (d=0; d<2; d++)
              //      if ((posout[d]+=1.0) < dnaxisn[d])
              //        break;
              //      else
              //        posout[d] = 1.0;
              //}
          }
      }
  }

  return res_psf;
}


/******************************psf_trim***************************************/
/**
 *
 * Function: psf_trim
 *
 * Trim the part of a psf-struct, leaving the center where
 * most of the signal is. Adjust the other metadata of the
 * structure as well.
 * The definitions on how to compute the trimmed dimensions
 * from "optsize" follow "filter.c:optPSFExSize()"
 *
 * @author MK
 * @date   April 2014
 *
 * @param[in,out] psf     the PSF structure
 * @param[in]     optsize  the desired size
 */
void psf_trim(psfstruct *psf, const int optsize){
	int xstart, ystart, xend, yend;
	int xnewsize, ynewsize;
	int ii, jj, kk;
	float *pix1, *pix2, *pix2save, *rowpos, *actpos;

	// set the start and end values for x/y
	xstart = psf->masksize[0]/2-optsize < 0 ? 0 : psf->masksize[0]/2-optsize;
	ystart = psf->masksize[1]/2-optsize < 0 ? 0 : psf->masksize[1]/2-optsize;
	xend   = psf->masksize[0]/2+optsize+1 < psf->masksize[0] ? psf->masksize[0]/2+optsize+1 : psf->masksize[0];
	yend   = psf->masksize[1]/2+optsize+1 < psf->masksize[1] ? psf->masksize[1]/2+optsize+1 : psf->masksize[1];

	// return if there is nothing to do
	if (xstart==0 && ystart==0 && xend==psf->masksize[0] && yend==psf->masksize[1]){
		return;
	}

	// determine the old and new sizes
	xnewsize = xend-xstart;
	ynewsize = yend-ystart;

	// allocate memory for the new data vector
	QMALLOC(pix2, float, xnewsize*ynewsize*psf->masksize[2]);

	// transfer the pixel values
	pix1 = psf->maskcomp;
	pix2save = pix2;
	for (kk=0, pix1=psf->maskcomp; kk<psf->masksize[2]; kk++, pix1+=psf->masksize[0]*psf->masksize[1]){
		for (jj=ystart, rowpos=pix1+(ystart*psf->masksize[1]); jj<yend; jj++, rowpos+=psf->masksize[1]){
			for (ii=xstart, actpos=rowpos+xstart; ii<xend; ii++){
				*(pix2++) = *(actpos++);
			}
		}
	}

	// transfer the memory to
	// the structure
	free(psf->maskcomp);
	psf->maskcomp = pix2save;

	// set the new size
	psf->masksize[0] = xnewsize;
	psf->masksize[1] = ynewsize;
	psf->masknpix = psf->masksize[0]*psf->masksize[1]*psf->masksize[2];

	// re-set the local storage
	free(psf->maskloc);
	QMALLOC(psf->maskloc, float, psf->masksize[0]*psf->masksize[1]);
	memset(psf->maskloc, psf->masksize[0]*psf->masksize[1], sizeof(float));

	// mark as empty
	psf->build_flag=0;
}

/******************************psf_normalize**********************************/
/**
 *
 * Function: psf_normalize
 *
 * Normalize a psf structure. This is done by evaluating the psf at
 * the offset position, summing over the local psf and then dividing
 * all components by that sum.
 *
 * @author MK
 * @date   April 2014
 *
 * @param[in,out] psf     the PSF structure
 */
void psf_normalize(psfstruct *psf){
  double pos[POLY_MAXDIM];
  float *actpix;
  float sum=0.0;
  int index=0;

  // build the psf at the offset position,
  // which is in the middle of the field
  for (index=0; index < psf->poly->ndim; index++){
    pos[index] = psf->contextoffset[index];
  }
  psf->build_flag=0;
  psf_buildpos(psf, pos, psf->poly->ndim);

  // sum up the evaluated psf
  actpix = psf->maskloc;
  for (index=0; index<psf->masksize[0]*psf->masksize[1]; index++){
      sum += *(actpix++);
  }

  // normalize the data values
  actpix = psf->maskcomp;
  for (index=0; index<psf->masknpix; index++){
      *(actpix++) /= sum;
  }
}

/******************************** psf_fit ***********************************/
/*                   standard PSF fit for one component                     */
/****************************************************************************/

void	psf_fit(psfstruct *psf, fieldstruct *field, fieldstruct *wfield,
		obj2struct *obj2)
{
  checkstruct		*check;
  static double	x2[PSF_NPSFMAX],y2[PSF_NPSFMAX],xy[PSF_NPSFMAX],
			deltax[PSF_NPSFMAX],
			deltay[PSF_NPSFMAX],
			flux[PSF_NPSFMAX],fluxerr[PSF_NPSFMAX],
			deltaxb[PSF_NPSFMAX],deltayb[PSF_NPSFMAX],
			fluxb[PSF_NPSFMAX],fluxerrb[PSF_NPSFMAX],
			sol[PSF_NTOT], covmat[PSF_NTOT*PSF_NTOT], 
			vmat[PSF_NTOT*PSF_NTOT], wmat[PSF_NTOT];
  float			*data, *data2, *data3, *weight, *d, *w;
  double		*mat,
			*m, *var,
			dx,dy,
			pix,pix2, wthresh,val,
			backnoise2, gain, radmin2,radmax2,satlevel,chi2,
			r2, valmax, psf_fwhm;
  float			**psfmasks, **psfmaskx,**psfmasky,
			*ps, *dh, *wh, pixstep;
  PIXTYPE		*datah, *weighth;
  int			i,j,p, npsf,npsfmax, npix, nppix, ix,iy,niter,
			width, height, pwidth,pheight, x,y,
			xmax,ymax, wbad, gainflag, convflag, npsfflag,
			ival,kill=0;
    
  dx = dy = 0.0;
  niter = 0;
  npsfmax = prefs.psf_npsfmax;
  pixstep = 1.0/psf->pixstep;
  gain = (field->gain >0.0? field->gain: 1e30);
  backnoise2 = field->backsig*field->backsig;
  satlevel = field->satur_level - obj2->bkg[0];
  wthresh = wfield?wfield->weight_thresh:BIG;
  gainflag = prefs.weightgain_flag[0];
  psf_fwhm = psf->fwhm*psf->pixstep;

 
  /* Initialize outputs */
  thepsfit->niter = 0;
  thepsfit->npsf = 0;
  for (j=0; j<npsfmax; j++) 
    {
      thepsfit->x[j] = obj2->posx[0];
      thepsfit->y[j] = obj2->posy[0];
      thepsfit->flux[j] = 0.0;
      thepsfit->fluxerr[j] = 0.0;
    }

  /* Scale data area with object "size" */
  ix = (obj2->xmax+obj2->xmin+1)/2;
  iy = (obj2->ymax+obj2->ymin+1)/2;
  width = obj2->xmax-obj2->xmin+1+psf_fwhm;
  if (width < (ival=(int)(psf_fwhm*2)))
    width = ival;
  height = obj2->ymax-obj2->ymin+1+psf_fwhm;
  if (height < (ival=(int)(psf_fwhm*2)))
    height = ival;
  npix = width*height;
  radmin2 = PSF_MINSHIFT*PSF_MINSHIFT;
  radmax2 = npix/2.0;

  /* Scale total area with PSF FWHM */
  pwidth = (int)(psf->masksize[0]*psf->pixstep)+width;;
  pheight = (int)(psf->masksize[1]*psf->pixstep)+height;
  nppix = pwidth*pheight;

  QMALLOC(weighth, PIXTYPE, npix);
  QMALLOC(weight, float, npix);
  QMALLOC(datah, PIXTYPE, npix);
  QMALLOC(data, float, npix);
  QMALLOC(data2, float, npix);
  QMALLOC(data3, float, npix);
  QMALLOC(mat, double, npix*PSF_NTOT);
  if (prefs.check[CHECK_SUBPSFPROTOS] || prefs.check[CHECK_PSFPROTOS])
    {
      QMALLOC(checkmask, PIXTYPE, nppix);
    }

  QMALLOC(psfmasks, float *, npsfmax);
  QMALLOC(psfmaskx, float *, npsfmax);
  QMALLOC(psfmasky, float *, npsfmax);
  for (i=0; i<npsfmax; i++)
    {
      QMALLOC(psfmasks[i], float, npix);
      QMALLOC(psfmaskx[i], float, npix);
      QMALLOC(psfmasky[i], float, npix);
    }

  /* Compute weights */
  wbad = 0;
  if (wfield)
    {
    psf_copyobjpix(datah, weighth, width, height, ix,iy, obj2,
		!(field->flags&MEASURE_FIELD));
    for (wh=weighth, w=weight, dh=datah,p=npix; p--;)
      if ((pix=*(wh++)) < wthresh && pix>0
	&& (pix2=*(dh++))>-BIG
	&& pix2<satlevel)
        *(w++) = 1/sqrt(pix+(pix2>0.0?
			(gainflag? pix2*pix/backnoise2:pix2)/gain
			:0.0));
        else
          {
          *(w++) = 0.0;
          wbad++;
          }
    }
  else
    {
    psf_copyobjpix(datah, NULL, width, height, ix,iy, obj2,
		!(field->flags&MEASURE_FIELD));
    for (w=weight, dh=datah, p=npix; p--;)
      if ((pix=*(dh++))>-BIG && pix<satlevel)
        *(w++) = 1.0/sqrt(backnoise2+(pix>0.0?pix/gain:0.0));
      else
        {
          *(w++) = 0.0;
          wbad++;
        }
    }

  /* Special action if most of the weights are zero!! */
  if (wbad>=npix-3)
    return;

  /* Weight the data */
  dh = datah;
  val = obj2->dbkg[0];      /* Take into account a local background change */
  d = data;
  w = weight;
  for (p=npix; p--;)
    *(d++) = (*(dh++)-val)**(w++);

  /* Get the local PSF */
  psf_build(psf, obj2);

  npsfflag = 1;
  r2 = psf_fwhm*psf_fwhm/2.0;
  fluxb[0] = fluxerrb[0] = deltaxb[0] = deltayb[0] = 0.0;

  for (npsf=1; npsf<=npsfmax && npsfflag; npsf++)
    {
      kill=0;
/*-- First compute an optimum initial guess for the positions of components */
      if (npsf>1)
        {
/*---- Subtract previously fitted components */
          d = data2;
          dh = datah;
          for (p=npix; p--;)
            *(d++) = (double)*(dh++);
          for (j=0; j<npsf-1; j++)
            {
              d = data2;
              ps = psfmasks[j];
              for (p=npix; p--;)
                *(d++) -= flux[j]**(ps++);
            }
          convolve_image(field, data2, data3, width,height);
/*---- Ignore regions too close to stellar cores */
          for (j=0; j<npsf-1; j++)
            {
              d = data3;
              dy = -((double)(height/2)+deltay[j]);
              for (y=height; y--; dy += 1.0)
                {
                  dx = -((double)(width/2)+deltax[j]);
                  for (x=width; x--; dx+= 1.0, d++)
                    if (dx*dx+dy*dy<r2)
                      *d = -BIG;
                }
            }
/*---- Now find the brightest pixel (poor man's guess, to be refined later) */
          d = data3;
          valmax = -BIG;
          xmax = width/2;
          ymax = height/2;
          for (y=0; y<height; y++)
            for (x=0; x<width; x++)
              {
                if ((val = *(d++))>valmax)
                  {
                    valmax = val;
                    xmax = x;
                    ymax = y;
                  }
              }
          deltax[npsf-1] = (double)(xmax - width/2);
          deltay[npsf-1] = (double)(ymax - height/2);
        }
      else
        {
/*---- Only one component to fit: simply use the barycenter as a guess */
          deltax[npsf-1] = obj2->mx - ix;
          deltay[npsf-1] = obj2->my - iy;
        }

      niter = 0;
      convflag = 1;
      for (i=0; i<PSF_NITER && convflag; i++)
        {
          convflag = 0,niter++,m=mat;
          for (j=0; j<npsf; j++)
            {
/*------ Resample the PSFs here for the 1st iteration */
              vignet_resample(psf->maskloc, psf->masksize[0], psf->masksize[1],
                              psfmasks[j], width, height,
                              -deltax[j]*pixstep, -deltay[j]*pixstep,
                              pixstep);       
              m=compute_gradient(weight,width,height,
                                 psfmasks[j],psfmaskx[j],psfmasky[j],m);
            }
          
          
          svdfit(mat, data, npix, npsf*PSF_NA, sol, vmat, wmat);
          
          compute_pos( &npsf, &convflag, &npsfflag,radmin2,radmax2,
                       r2, sol,flux, deltax, deltay,&dx,&dy);
        }

/*-- Compute variances and covariances */
      svdvar(vmat, wmat, npsf*PSF_NA, covmat);
      var = covmat;
      for (j=0; j<npsf; j++, var += (npsf*PSF_NA+1)*PSF_NA)
        {
/*---- First, the error on the flux estimate */      
          fluxerr[j] = sqrt(*var)>0.0?  sqrt(*var):999999.0;
          //if (flux[j]<12*fluxerr && j>0)
          //  npsfmax--,flux[j]=0;
          if (flux[j]<12*fluxerr[j] && j>0)
                 {
                   flux[j]=0,kill++,npsfmax--;
                   //if(j==npsfmax-1)
                   //  kill++;             
                 } 
        }
      if (npsfflag)
        {
/*--- If we reach this point we know the data are worth backuping */
          for (j=0; j<npsf; j++)
            {
              deltaxb[j] = deltax[j];
              deltayb[j] = deltay[j];
              fluxb[j] = flux[j];
              fluxerrb[j]=fluxerr[j];
            }
        }
    }
  npsf=npsf-1-kill;

/* Now keep only fitted stars that fall within the current detection area */
  i = 0;
  for (j=0; j<npsf; j++)
    {      
      x = (int)(deltaxb[j]+0.4999)+width/2;
      y = (int)(deltayb[j]+0.4999)+height/2;
      if (x<0 || x>=width || y<0 || y>=height)
        continue;
      if (weight[y*width+x] < 1/BIG)
        continue;
      if (10*fluxb[j]<fluxb[0] )
        continue;
      if (fluxb[j]<=0 )
        continue; 

      if (FLAG(obj2.poserrmx2_psf))
        compute_poserr(j,covmat,sol,obj2,x2,y2,xy, npsf);
      
      deltax[i] = deltaxb[j];
      deltay[i] = deltayb[j];
      flux[i] = fluxb[j];
      fluxerr[i++] = fluxerrb[j];
    }

  npsf = i;

  /* Compute chi2 if asked to 
  if (FLAG(obj2.chi2_psf))
    {
      for (j=0; j<npsf; j++)
        {
          chi2 = 0.0;
          for (d=data,w=weight,p=0; p<npix; w++,p++)
            {
              pix = *(d++);
              pix -=  psfmasks[j][p]*flux[j]**w;
              chi2 += pix*pix;
              if (chi2>1E29) chi2=1E28;
            }
          obj2->chi2_psf = obj2->sigbkg>0.?
            chi2/((npix - 3*npsf)*obj2->sigbkg*obj2->sigbkg):999999;

        }
      
    }*/
 /* Compute relative chi2 if asked to */
    if (FLAG(obj2.chi2_psf))
    {
      for (j=0; j<npsf; j++)
        {
          chi2 = 0.0;
          for (d=data,w=weight,p=0; p<npix; w++,p++)
            {
              pix = *(d++)/flux[j];
              pix -=  psfmasks[j][p]**w;
              chi2 += pix*pix;
              if (chi2>1E29) chi2=1E28;
            }
          obj2->chi2_psf = flux[j]>0?
		chi2/((npix - 3*npsf)*obj2->sigbkg[0]*obj2->sigbkg[0]):999999;

        }
      
    }
  /* CHECK images */
  if (prefs.check[CHECK_SUBPSFPROTOS] || prefs.check[CHECK_PSFPROTOS])
    for (j=0; j<npsf; j++)
      {
        vignet_resample(psf->maskloc, psf->masksize[0], psf->masksize[1],
                        checkmask, pwidth, pheight,
                        -deltax[j]*pixstep, -deltay[j]*pixstep, pixstep);
        if ((check = prefs.check[CHECK_SUBPSFPROTOS]))
          check_add(check, checkmask, pwidth,pheight, ix,iy,-flux[j]);
        if ((check = prefs.check[CHECK_PSFPROTOS]))
          check_add(check, checkmask, pwidth,pheight, ix,iy,flux[j]);
      }

  thepsfit->niter = niter;
  thepsfit->npsf = npsf;
  for (j=0; j<npsf; j++)
    {
      thepsfit->x[j] = ix+deltax[j]+1.0;
      thepsfit->y[j] = iy+deltay[j]+1.0;
      thepsfit->flux[j] = flux[j];
      thepsfit->fluxerr[j] = fluxerr[j];
    }  
  
  for (i=0; i<prefs.psf_npsfmax; i++)
    {
      QFREE(psfmasks[i]);
      QFREE(psfmaskx[i]);
      QFREE(psfmasky[i]);
    }

  QFREE(psfmasks);
  QFREE(psfmaskx);
  QFREE(psfmasky);
  QFREE(datah);
  QFREE(data);
  QFREE(data2);
  QFREE(data3);
  QFREE(weighth);
  QFREE(weight);
  QFREE(data);
  QFREE(mat);

  if (prefs.check[CHECK_SUBPSFPROTOS] || prefs.check[CHECK_PSFPROTOS])
    QFREE(checkmask);

  return;
}


/****************************** double_psf_fit ******************************/
/* double fit to make psf detection on one image and photometry on another  */
/****************************************************************************/

void    double_psf_fit(psfstruct *ppsf, fieldstruct *pfield, fieldstruct *pwfield,
                       obj2struct *obj2,
			psfstruct *psf, fieldstruct *field, fieldstruct *wfield)
{
  static double      /* sum[PSF_NPSFMAX]*/ pdeltax[PSF_NPSFMAX],
    pdeltay[PSF_NPSFMAX],psol[PSF_NPSFMAX], pcovmat[PSF_NPSFMAX*PSF_NPSFMAX], 
    pvmat[PSF_NPSFMAX*PSF_NPSFMAX], pwmat[PSF_NPSFMAX],pflux[PSF_NPSFMAX],
    pfluxerr[PSF_NPSFMAX];

    double *pmat,
     *pm, /* *pps,  *px, *py,*/
    dx,dy,pdx,pdy, /* x1,y1, mx,my,mflux, */
    val, ppix,ppix2, /* dflux, */
    gain, radmin2,radmax2,satlevel
    ,chi2,pwthresh,pbacknoise2, /* mr, */
    r2=0, psf_fwhm,ppsf_fwhm ;
  float         **ppsfmasks, **ppsfmaskx,**ppsfmasky, *pps;
  float         *pdata, *pdata2, *pdata3, *pweight, *pd, *pw, 
		*pdh, *pwh, pixstep,ppixstep;
  PIXTYPE       *pdatah, *pweighth;
  int                   i,j,k,p, npsf, npix,ix,iy,
    width, height, /* hw,hh, */
    x,y, /* yb, */
    wbad, gainflag,
    ival,npsfmax;
  double *pvar;
  
  pdx = pdy =dx = dy = 0.0;
  ppixstep = 1.0/ppsf->pixstep;
  pixstep = 1.0/psf->pixstep;
  gain = (field->gain >0.0? field->gain: 1e30);
  npsfmax=prefs.psf_npsfmax;
  pbacknoise2 = pfield->backsig*pfield->backsig;
  satlevel = field->satur_level - obj2->bkg[0];
  gainflag = prefs.weightgain_flag[0];
  psf_fwhm = psf->fwhm*psf->pixstep;
  ppsf_fwhm = ppsf->fwhm*ppsf->pixstep;
  pwthresh = pwfield?pwfield->weight_thresh:BIG;

  /* Initialize outputs */
  ppsfit->niter = 0;
  ppsfit->npsf = 0;
  for (j=0; j<npsfmax; j++) 
    {
      ppsfit->x[j] = 999999.0;
      ppsfit->y[j] = 999999.0;
      ppsfit->flux[j] = 0.0;
      ppsfit->fluxerr[j] = 0.0;
      pdeltax[j]= pdeltay[j]=psol[j]= pwmat[j]=pflux[j]=pfluxerr[j]=0.0;
   
    }

  ix = (obj2->xmax+obj2->xmin+1)/2;
  iy = (obj2->ymax+obj2->ymin+1)/2;
  width = obj2->xmax-obj2->xmin+1+psf_fwhm;
  if (width < (ival=(int)(psf_fwhm*2)))
    width = ival;
  height = obj2->ymax-obj2->ymin+1+psf_fwhm;
  if (height < (ival=(int)(psf_fwhm*2)))
    height = ival;
  npix = width*height;
  radmin2 = PSF_MINSHIFT*PSF_MINSHIFT;
  radmax2 = npix/2.0;
  psf_fit(psf, field, wfield, obj2);
  npsf=thepsfit->npsf;
  
  QMALLOC(ppsfmasks,float *,npsfmax);
  QMALLOC(ppsfmaskx,float *,npsfmax);
  QMALLOC(ppsfmasky,float *,npsfmax);

  for (i=0; i<npsfmax; i++)
    {
      QMALLOC(ppsfmasks[i],float,npix);
      QMALLOC(ppsfmaskx[i],float,npix);
      QMALLOC(ppsfmasky[i],float,npix);
    }

  QMALLOC(pweighth, PIXTYPE, npix);
  QMALLOC(pweight, float, npix);
  QMALLOC(pdatah, PIXTYPE, npix);
  QMALLOC(pdata, float, npix);
  QMALLOC(pdata2, float, npix);
  QMALLOC(pdata3, float, npix);
  QMALLOC(pmat, double, npix*npsfmax);
  
   for (j=0; j<npsf; j++)
    {
      pdeltax[j] =thepsfit->x[j]-ix-1 ;
      pdeltay[j] =thepsfit->y[j]-iy-1 ;
      ppsfit->flux[j] = 0;
      ppsfit->fluxerr[j] = 0;
    }

/*-------------------  Now the photometry fit ---------------------*/
   /* Compute photometry weights */
  wbad = 0;
  if (pwfield)
    {
    psf_copyobjpix(pdatah, pweighth, width, height, ix,iy, obj2, 0);
    for (pwh=pweighth, pw=pweight, pdh=pdatah,p=npix; p--;)
      {
      if ((ppix=*(pwh++)) < pwthresh && ppix>0
		&& (ppix2=*(pdh++))>-BIG  && ppix2<satlevel)
        *(pw++) = 1/sqrt(ppix+(ppix2>0.0?
			(gainflag? ppix2*ppix/pbacknoise2:ppix2)/gain : 0.0));
      else
        {
        *(pw++) = 0.0;          
        wbad++;
        }
      }
    }
  else
    {
    psf_copyobjpix(pdatah, NULL, width, height, ix,iy, obj2, 0);
    for (pw=pweight, pdh=pdatah, p=npix; p--;)
      if ((ppix=*(pdh++))>-BIG && ppix<satlevel)
        *(pw++) = 1.0/sqrt(pbacknoise2+(ppix>0.0?ppix/gain:0.0));
      else
        {
        *(pw++) = 0.0;
        wbad++;
        }
    }

  /* Special action if most of the weights are zero!! */
  if (wbad>=npix-3)
    return;

  /* Weight the data */
  pdh = pdatah;
  pd = pdata;
  pw = pweight;
  val = obj2->dbkg[0];
  for (p=npix; p--;)
    *(pd++) = (*(pdh++)-val)**(pw++);

 
  /* Get the photmetry PSF */
   psf_build(ppsf, obj2);
  for (j=1; j<=npsf; j++)
    {
      if (j>1)
        {
          /*---- Subtract //previously fitted components in photometry image */
          pd = pdata2;
          pdh = pdatah;
          for (p=npix; p--;)
            *(pd++) = (double)*(pdh++);
          for (k=0; k<j-1; k++)
            {
              pd = pdata2;
              pps = ppsfmasks[k];
              for (p=npix; p--;)
                *(pd++) -= pflux[k]**(pps++);
            }
          convolve_image(pfield, pdata2, pdata3, width,height);
         /*---- Ignore regions too close to stellar cores */
          for (k=0; k<j-1; k++)
            {
              pd = pdata3;
              dy = -((double)(height/2)+pdeltay[k]);
              for (y=height; y--; dy += 1.0)
                {
                  dx = -((double)(width/2)+pdeltax[k]);
                  for (x=width; x--; dx+= 1.0, pd++)
                    if (dx*dx+dy*dy<r2) /*?*/
                      *pd = -BIG;
                }
            } 
        }
   
      pm=pmat;
      for (k=0; k<j; k++)
            {
              /*------ Resample the PSFs here for the 1st iteration */
              vignet_resample(ppsf->maskloc,
			ppsf->masksize[0], ppsf->masksize[1],
			ppsfmasks[k], width, height,
			-pdeltax[k]*ppixstep, -pdeltay[k]*ppixstep,
			ppixstep);              
              pm=compute_gradient_phot(pweight,width,height, ppsfmasks[k],pm);
            }
      
      svdfit(pmat, pdata, npix, j, psol, pvmat, pwmat);  
      compute_pos_phot( &j, psol,pflux);
   
  for (k=0; k<j; k++)
        {
          svdvar(pvmat, pwmat, j, pcovmat);
          pvar = pcovmat;
          pfluxerr[k]= sqrt(*pvar)>0.0 && sqrt(*pvar)<99? 
            sqrt(*pvar):99;
        }
    }
  /* Compute chi2 if asked to 
  if (FLAG(obj2.chi2_psf))
    {
      for (j=0; j<npsf; j++)
        {
          chi2 = 0.0;
          for (pd=pdata,pw=pweight,p=0; p<npix; pw++,p++)
            {
              ppix = *(pd++);
              ppix -=  ppsfmasks[j][p]*pflux[j]**pw;
              chi2 += ppix*ppix;
              if (chi2>1E29) chi2=1E28;
            }
          obj2->chi2_psf = obj2->sigbkg>0.?
            chi2/((npix - 3*npsf)*obj2->sigbkg*obj2->sigbkg):999999;

        }
      
    }
 */
 /* Compute relative error if asked to */
  if (FLAG(obj2.chi2_psf))
  {
      for (j=0; j<npsf; j++)
        {
          chi2 = 0.0;
          for (pd=pdata,pw=pweight,p=0; p<npix; pw++,p++)
            {
              ppix = *(pd++)/pflux[j];
              ppix -=  ppsfmasks[j][p]**pw;
              chi2 += ppix*ppix;
              if (chi2>1E29) chi2=1E28;
            }
          obj2->chi2_psf = pflux[j]>0?
		chi2/((npix - 3*npsf)*obj2->sigbkg[0]*obj2->sigbkg[0]):999999;

        }
      
    }
  ppsfit->niter = thepsfit->niter;
  ppsfit->npsf = npsf;

  for (j=0; j<npsf; j++)
    {
      thepsfit->x[j] = ix+pdeltax[j]+1.0;
      thepsfit->y[j] = iy+pdeltay[j]+1.0;
      thepsfit->flux[j] = pflux[j];
      thepsfit->fluxerr[j] = pfluxerr[j];
      ppsfit->x[j] = ix+pdeltax[j]+1.0;
      ppsfit->y[j] = iy+pdeltay[j]+1.0;
      ppsfit->flux[j] = pflux[j];
      ppsfit->fluxerr[j] = pfluxerr[j];
    }
    
  
  for (i=0; i<npsfmax; i++)
    {
      QFREE(ppsfmasks[i]);
      QFREE(ppsfmaskx[i]);
      QFREE(ppsfmasky[i]);
    }

  
  QFREE(ppsfmasks);
  QFREE(ppsfmaskx);
  QFREE(ppsfmasky);
  QFREE(pdatah);
  QFREE(pdata);
  QFREE(pdata2);
  QFREE(pdata3);
  QFREE(pweighth);
  QFREE(pweight);
  QFREE(pdata);
  QFREE(pmat);
   
  if (prefs.check[CHECK_SUBPSFPROTOS] || prefs.check[CHECK_PSFPROTOS])
    {
      QFREE(checkmask);
    }
  return;
  }


/****** psf_copyobjpix ******************************************************
PROTO	int psf_copyobjpix(PIXTYPE *image, PIXTYPE *imvar,
			int wout, int hout, int ix, int iy,
			obj2struct *obj2, int detect_flag);
PURPOSE	Copy a piece of the input object image/weights to local arrays.
INPUT	Pointer to the output image array,
	pointer to the output image variance array,
	output frame width,
	output frame height,
	integer x coordinate,
	integer y coordinate,
	pointer to the obj2 structure,
	detection field flag (non 0 only for pure detection images).
OUTPUT	RETURN_ERROR if the coordinates are outside object image,
	RETURN_OK otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/09/2014
 ***/
int	psf_copyobjpix(PIXTYPE *image, PIXTYPE *imvar,
			int wout, int hout, int ix, int iy,
			obj2struct *obj2, int detect_flag)
  {
   subimagestruct	*subimage;
   PIXTYPE		*imaget,*imvart, *imagein,*imvarin;
   int			i,y, win,hin, w2, xmin,xmax,ymin,ymax;

/* Set output to -BIG */
  imaget = image;
  for (i=wout*hout; i--;)
    *(imaget++) = -BIG;
  if (imvar)
    memset(imvar, 0, wout*hout*sizeof(PIXTYPE));

  subimage = obj2->subimage;	/* !CHECK */
  ix -= subimage->xmin[0];
  iy -= subimage->xmin[1];
  win = subimage->size[0];
  hin = subimage->size[1];

/* Don't go further if out of frame!! */
  if (ix<0 || ix>=win || iy<0 || iy>=hin)
    return RETURN_ERROR;

/* Set the image boundaries */
  w2 = wout;
  ymin = iy-hout/2;
  ymax = ymin + hout;
  if (ymin<0)
    {
    image -= ymin*wout;
    if (imvar)
      imvar -= ymin*wout;
    ymin = 0;
    }
  if (ymax>hin)
    ymax = hin;

  xmin = ix-wout/2;
  xmax = xmin + wout;
  if (xmax>win)
    {
    w2 -= xmax-win;
    xmax = win;
    }
  if (xmin<0)
    {
    image -= xmin;
    if (imvar)
      imvar -= xmin;
    w2 += xmin;
    xmin = 0;
    }

/* Copy the right pixels to the destination */
  imagein = subimage->image;
  for (y=ymin; y<ymax; y++, image += wout)
    memcpy(image, imagein + xmin+y*win, w2*sizeof(PIXTYPE));
  if (imvar)
    {
    imvarin = subimage->imvar;
    for (y=ymin; y<ymax; y++, imvar += wout)
      memcpy(imvar, imvarin + xmin+y*win, w2*sizeof(PIXTYPE));
    }


  return RETURN_OK;
  }



/******************************* psf_build **********************************/
/*
Build the local PSF (function of "context").
*/
void	psf_build(psfstruct *psf, obj2struct *obj2)
  {
   double	pos[POLY_MAXDIM];
   double	*basis, fac;
   float	*ppc, *pl;
   int		i,n,p, ndim, npix;

  if (psf->build_flag)
    return;

  npix = psf->masksize[0]*psf->masksize[1];

/* Reset the Local PSF mask */
  memset(psf->maskloc, 0, npix*sizeof(float));

/* Grab the context vector */
  ndim = psf->poly->ndim;
  for (i=0; i<ndim; i++)
    {
    ttypeconv((char *)obj2 + ((void *)psf->context[i]-(void *)&flagobj2),
		&pos[i], psf->contexttyp[i],T_DOUBLE);
    pos[i] = (pos[i] - psf->contextoffset[i]) / psf->contextscale[i];
    }
  poly_func(psf->poly, pos);

  basis = psf->poly->basis;

  ppc = psf->maskcomp;
/* Sum each component */
  for (n = (psf->maskdim>2?psf->masksize[2]:1); n--;)
    {
    pl = psf->maskloc;
    fac = *(basis++);
    for (p=npix; p--;)
      *(pl++) +=  fac**(ppc++);
    }

  psf->build_flag = 1;

  return;
  }

/**************************psf_buildpos***************************************/
/**
 * Function: psf_buildpos
 *
 * Evaluates the PSF struct at a given position. The PSF values
 * are written to the local storage of the PSF struct.
 *
 * @author MK
 * @date   April 2014
 *
 * @param[in,out] psf    the PSF struct
 * @param[in]     pos    position to evaluate the PSF struct
 * @param[in]     inndim number of dimensions/values in 'pos'
 *
 */
void psf_buildpos(psfstruct *psf, const double *posin, const int inndim)
{
  double	pos[POLY_MAXDIM];
  double	*basis, fac;
  float	*ppc, *pl;
  int		n, p, npix;

  // looks like this can only
  // be un-set from the 'outside'
  if (psf->build_flag){
      QPRINTF(OUTPUT, "nothing to do, build-flag is: %i\n", psf->build_flag);
      return;
  }
  // assure the dimensions match
  if (inndim != psf->poly->ndim){
      error(EXIT_FAILURE, "*Error*: the dimensions differ ", "for the PSF!");
  }

  // compute the number of pixels and reset the Local PSF mask
  npix = psf->masksize[0]*psf->masksize[1];
  memset(psf->maskloc, 0, npix*sizeof(float));

  // normalize the positional values
  for (n=0; n<psf->poly->ndim; n++){
      pos[n] = (posin[n] - psf->contextoffset[n]) / psf->contextscale[n];
  }

  // evaluate the polynomial
  poly_func(psf->poly, pos);
  basis = psf->poly->basis;

  // Sum each component
  ppc = psf->maskcomp;
  for (n = (psf->maskdim>2?psf->masksize[2]:1); n--;){
      pl = psf->maskloc;
      fac = *(basis++);
      for (p=npix; p--;)
        *(pl++) +=  fac**(ppc++);
  }

  // set the build flag
  psf->build_flag = 1;

  return;
}


/******************************** psf_fwhm **********************************/
/*
Return the local PSF FWHM.
*/
double	psf_fwhm(psfstruct *psf, obj2struct *obj2)
  {
   float	*pl,
		val, max;
   int		n,p, npix;

  if (!psf->build_flag)
    psf_build(psf, obj2);

  npix = psf->masksize[0]*psf->masksize[1];
  max = -BIG;
  pl = psf->maskloc;
  for (p=npix; p--;)
    if ((val=*(pl++)) > max)
      max = val;
  pl = psf->maskloc;
  max /= 2.0;
  n = 0;
  for (p=npix; p--;)
    if (*(pl++) >= max)
      n++;

  return 2.0*sqrt(n/PI)*psf->pixstep;
  }


/*****************************compute_gradient*********************************/

double *compute_gradient(float *weight,int width, int height,
                         float *masks,float *maskx,float *masky
                        ,double *m)
{
  int x,y;
  float	*w, *ps,*px,*py;
    
  /*------ copy of the (weighted) PSF, with outer ring set to zero */
      ps = masks;
      w = weight;
      for (y=0; y<height; y++)
        for (x=0; x<width; x++, ps++, w++)
          *(m++) = y?(y>=(height-1)?0:(x?(x>=(width-1)?0:*ps**w):0)):0;
      /*------ (weighted) PSF gradient in x (kernel for first moment in x) */
      ps = masks;
      px = maskx;
      w = weight;
      for (y=0; y<height; y++)
        for (x=0; x<width; x++, ps++, w++)
          *(m++) = ((*px++) = (x?(x>=(width-1)?0:*(ps+1)-*(ps-1)):0))**w/2;
      /*------ (weighted) PSF gradient in y (kernel for first moment in y) */
      ps = masks; 
      py = masky;
      w = weight;
      for (y=0; y<height; y++)
        for (x=0; x<width; x++, ps++, w++)
          *(m++) = (*(py++)=(y?(y>=(height-1)?0:*(ps+width)-*(ps-width)):0))
            **w/2;
  return m;
}

/*****************************compute_gradient_phot*****************************
****/

double *compute_gradient_phot(float  *pweight,int width, int height,
                         float *pmasks,double *pm)

{
  int x,y;
  float  *pw, *pps;
    
  /*------ copy of the (weighted) PSF, with outer ring set to zero */
      pps = pmasks;
      pw = pweight;
      for (y=0; y<height; y++)
        for (x=0; x<width; x++, pps++, pw++)
          *(pm++) = y?(y>=(height-1)?0:(x?(x>=(width-1)?0:*pps**pw):0)):0;

  return pm;
}

/**************************compute_pos********************************/

void compute_pos(int *pnpsf,int *pconvflag,int *pnpsfflag,double radmin2,
                         double radmax2,double r2,double *sol,double *flux 
                        ,double *deltax,double *deltay,double *pdx,double *pdy)
{
  int j,k,convflag,npsfflag,npsf; 
  double dx,dy;

  dx=*pdx;
  dy=*pdy;
  convflag=*pconvflag;
  npsfflag=*pnpsfflag;
  npsf=*pnpsf;
  for (j=0; j<npsf; j++)
    {
      flux[j] = sol[j*PSF_NA];
      /*------ Update the PSF shifts */
      if (fabs(flux[j])>0.0)
        {
          dx = -sol[j*PSF_NA+1]/((npsf>1?2:1)*flux[j]);
          dy = -sol[j*PSF_NA+2]/((npsf>1?2:1)*flux[j]);
        }
      
      deltax[j] += dx;
      deltay[j] += dy;
      /*------ Continue until all PSFs have come to a complete stop */
      if ((dx*dx+dy*dy) > radmin2)
        convflag = 1;
    }
  for (j=0; j<npsf; j++)
    {
      /*------ Exit if too much decentering or negative flux */
      for (k=j+1; k<npsf; k++)
        {
          dx = deltax[j]-deltax[k];
          dy = deltay[j]-deltay[k];
          if (dx*dx+dy*dy<r2/4.0)
            {
              flux[j] = -BIG;
              break;
            }
        }
      if (flux[j]<0.0
          || (deltax[j]*deltax[j] + deltay[j]*deltay[j]) > radmax2)
        {
          npsfflag = 0;
          convflag = 0;
          npsf--;
          break;
        }
    }
  *pdx=dx;
  *pdy=dy;
  *pconvflag=convflag;
  *pnpsfflag= npsfflag;
  *pnpsf=npsf;
  return;
}

/**************************compute_pos_phot********************************/

void compute_pos_phot(int *pnpsf,double *sol,double *flux)
{
  int j,npsf;   
  npsf=*pnpsf;
  for (j=0; j<npsf; j++)
    {
      flux[j] = sol[j];     
    }
  *pnpsf=npsf;
  return;
}


/************************************compute_poserr*****************************
*********/

void compute_poserr( int j,double *var,double *sol,obj2struct *obj2,double *x2,
                    double *y2,double *xy, int npsf)
{
  double vara,covab,varb, f2;

  /*------ Variances and covariance along x and y */
  vara = *(var += (PSF_NA*npsf+1)*(j*PSF_NA+1));
  covab = *(++var);
  varb = *(var += PSF_NA*npsf);
  f2 = sol[PSF_NA*j];
  f2 *= f2;
  obj2->poserrmx2_psf = vara/f2;
  obj2->poserrmy2_psf = varb/f2;
  obj2->poserrmxy_psf = covab/f2;

  /*------ If requested, translate variances to major and minor error axes... */
  if (FLAG(obj2.poserra_psf))
    {
      double    pmx2,pmy2,temp,theta;
      
      if (fabs(temp=obj2->poserrmx2_psf-obj2->poserrmy2_psf) > 0.0)
        theta = atan2(2.0 * obj2->poserrmxy_psf,temp) / 2.0;
      else
        theta = PI/4.0;
      
      temp = sqrt(0.25*temp*temp+obj2->poserrmxy_psf*obj2->poserrmxy_psf);
      pmy2 = pmx2 = 0.5*(obj2->poserrmx2_psf+obj2->poserrmy2_psf);
      pmx2+=temp;
      pmy2-=temp;

      obj2->poserra_psf = (float)sqrt(pmx2);
      obj2->poserrb_psf = (float)sqrt(pmy2);
      obj2->poserrtheta_psf = theta*180.0/PI;
    }
  
  /*------ ...Or ellipse parameters */
  if (FLAG(obj2.poserr_cxx))
    {
      double    xm2,ym2, xym, temp;
      
      xm2 = obj2->poserrmx2_psf;
      ym2 = obj2->poserrmy2_psf;
      xym = obj2->poserrmxy_psf;
      obj2->poserrcxx_psf = (float)(ym2/(temp=xm2*ym2-xym*xym));
      obj2->poserrcyy_psf = (float)(xm2/temp);
      obj2->poserrcxy_psf = (float)(-2*xym/temp);
    }
  return;
}


/******************************** svdfit ************************************/
/*
General least-square fit A.x = b, based on Singular Value Decomposition (SVD).
Loosely adapted from Numerical Recipes in C, 2nd Ed. (p. 671).
Note: the a and v matrices are transposed with respect to the N.R. convention.
*/
void svdfit(double *a, float *b, int m, int n, double *sol,
	double *vmat, double *wmat)
  {
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define	PYTHAG(a,b)	((at=fabs(a)) > (bt=fabs(b)) ? \
				  (ct=bt/at,at*sqrt(1.0+ct*ct)) \
				: (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define	TOL		1.0e-11

   int			flag,i,its,j,jj,k,l,nm,mmi,nml;
   double		c,f,h,s,x,y,z,
			anorm, g, scale,
			at,bt,ct,maxarg1,maxarg2,
			thresh, wmax,
			*w,*ap,*ap0,*ap1,*ap10,*rv1p,*vp,*vp0,*vp1,*vp10,
			*tmpp, *rv1,*tmp;
   float		*bp;

  anorm = g = scale = 0.0;
  if (m < n)
    error(EXIT_FAILURE, "*Error*: Not enough rows for solving the system ",
	"in svdfit()");
  
  QMALLOC(rv1, double, n);
  QMALLOC(tmp, double, n);
  l = nm = nml = 0;			/* To avoid gcc -Wall warnings */
  for (i=0;i<n;i++)
    {
    l = i+1;
    nml = n-l;
    rv1[i] = scale*g;
    g = s = scale = 0.0;
    if ((mmi = m - i) > 0)
      {
      ap = ap0 = a+i*(m+1);
      for (k=mmi;k--;)
        scale += fabs(*(ap++));
      if (scale)
        {
        for (ap=ap0,k=mmi; k--; ap++)
          {
          *ap /= scale;
          s += *ap**ap;
          }
        f = *ap0;
        g = -SIGN(sqrt(s),f);
        h = f*g-s;
        *ap0 = f-g;
        ap10 = a+l*m+i;
        for (j=nml;j--; ap10+=m)
          {
          for (s=0.0,ap=ap0,ap1=ap10,k=mmi; k--;)
            s += *(ap1++)**(ap++);
          f = s/h;
          for (ap=ap0,ap1=ap10,k=mmi; k--;)
            *(ap1++) += f**(ap++);
          }
        for (ap=ap0,k=mmi; k--;)
          *(ap++) *= scale;
        }
      }
    wmat[i] = scale*g;
    g = s = scale = 0.0;
    if (i < m && i+1 != n)
      {
      ap = ap0 = a+i+m*l;
      for (k=nml;k--; ap+=m)
        scale += fabs(*ap);
      if (scale)
        {
        for (ap=ap0,k=nml;k--; ap+=m)
          {
          *ap /= scale;
          s += *ap**ap;
          }
        f=*ap0;
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        *ap0=f-g;
        rv1p = rv1+l;
        for (ap=ap0,k=nml;k--; ap+=m)
          *(rv1p++) = *ap/h;
        ap10 = a+l+m*l;
        for (j=m-l; j--; ap10++)
          {
          for (s=0.0,ap=ap0,ap1=ap10,k=nml; k--; ap+=m,ap1+=m)
            s += *ap1**ap;
          rv1p = rv1+l;
          for (ap1=ap10,k=nml;k--; ap1+=m)
            *ap1 += s**(rv1p++);
          }
        for (ap=ap0,k=nml;k--; ap+=m)
          *ap *= scale;
        }
      }
    anorm=MAX(anorm,(fabs(wmat[i])+fabs(rv1[i])));
    }

  for (i=n-1;i>=0;i--)
    {
    if (i < n-1)
      {
      if (g)
        {
        ap0 = a+l*m+i;
        vp0 = vmat+i*n+l;
        vp10 = vmat+l*n+l;
        g *= *ap0;
        for (ap=ap0,vp=vp0,j=nml; j--; ap+=m)
          *(vp++) = *ap/g;
        for (j=nml; j--; vp10+=n)
          {
          for (s=0.0,ap=ap0,vp1=vp10,k=nml; k--; ap+=m)
            s += *ap**(vp1++);
          for (vp=vp0,vp1=vp10,k=nml; k--;)
            *(vp1++) += s**(vp++);
          }
        }
      vp = vmat+l*n+i;
      vp1 = vmat+i*n+l;
      for (j=nml; j--; vp+=n)
        *vp = *(vp1++) = 0.0;
      }
    vmat[i*n+i]=1.0;
    g=rv1[i];
    l=i;
    nml = n-l;
    }

  for (i=(m<n?m:n); --i>=0;)
    {
    l=i+1;
    nml = n-l;
    mmi=m-i;
    g=wmat[i];
    ap0 = a+i*m+i;
    ap10 = ap0 + m;
    for (ap=ap10,j=nml;j--;ap+=m)
      *ap=0.0;
    if (g)
      {
      g=1.0/g;
      for (j=nml;j--; ap10+=m)
        {
        for (s=0.0,ap=ap0,ap1=ap10,k=mmi; --k;)
              s += *(++ap)**(++ap1);
        f = (s/(*ap0))*g;
        for (ap=ap0,ap1=ap10,k=mmi;k--;)
          *(ap1++) += f**(ap++);
        }
      for (ap=ap0,j=mmi;j--;)
        *(ap++) *= g;
      }
    else
      for (ap=ap0,j=mmi;j--;)
        *(ap++)=0.0;
    ++(*ap0);
    }

  for (k=n; --k>=0;)
      {
      for (its=0;its<100;its++)
        {
        flag=1;
        for (l=k;l>=0;l--)
          {
          nm=l-1;
          if (fabs(rv1[l])+anorm == anorm)
            {
            flag=0;
            break;
            }
          if (fabs(wmat[nm])+anorm == anorm)
            break;
          }
        if (flag)
          {
          c=0.0;
          s=1.0;
          ap0 = a+nm*m;
          ap10 = a+l*m;
          for (i=l; i<=k; i++,ap10+=m)
            {
            f=s*rv1[i];
            if (fabs(f)+anorm == anorm)
              break;
            g=wmat[i];
            h=PYTHAG(f,g);
            wmat[i]=h;
            h=1.0/h;
            c=g*h;
            s=(-f*h);
            for (ap=ap0,ap1=ap10,j=m; j--;)
              {
              z = *ap1;
              y = *ap;
              *(ap++) = y*c+z*s;
              *(ap1++) = z*c-y*s;
              }
            }
          }
        z=wmat[k];
        if (l == k)
          {
          if (z < 0.0)
            {
            wmat[k] = -z;
            vp = vmat+k*n;
            for (j=n; j--; vp++)
              *vp = (-*vp);
            }
          break;
          }
        if (its == 99)
          error(EXIT_FAILURE, "*Error*: No convergence in 100 SVD iterations ",
		"in svdfit()");
        x=wmat[l];
        nm=k-1;
        y=wmat[nm];
        g=rv1[nm];
        h=rv1[k];
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
        g=PYTHAG(f,1.0);
        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
        c=s=1.0;
        ap10 = a+l*m;
        vp10 = vmat+l*n;
        for (j=l;j<=nm;j++,ap10+=m,vp10+=n)
          {
          i=j+1;
          g=rv1[i];
          y=wmat[i];
          h=s*g;
          g=c*g;
          z=PYTHAG(f,h);
          rv1[j]=z;
          c=f/z;
          s=h/z;
          f=x*c+g*s;
          g=g*c-x*s;
          h=y*s;
          y=y*c;
          for (vp=(vp1=vp10)+n,jj=n; jj--;)
            {
            z = *vp;
            x = *vp1;
            *(vp1++) = x*c+z*s;
            *(vp++) = z*c-x*s;
            }
          z=PYTHAG(f,h);
          wmat[j]=z;
          if (z)
            {
            z=1.0/z;
            c=f*z;
            s=h*z;
            }
          f=c*g+s*y;
          x=c*y-s*g;
          for (ap=(ap1=ap10)+m,jj=m; jj--;)
            {
            z = *ap;
            y = *ap1;
            *(ap1++) = y*c+z*s;
            *(ap++) = z*c-y*s;
            }
          }
        rv1[l]=0.0;
        rv1[k]=f;
        wmat[k]=x;
        }
      }

  wmax=0.0;
  w = wmat;
  for (j=n;j--; w++)
    if (*w > wmax)
      wmax=*w;
  thresh=TOL*wmax;
  w = wmat;
  for (j=n;j--; w++)
    if (*w < thresh)
      *w = 0.0;

  w = wmat;
  ap = a;
  tmpp = tmp;
  for (j=n; j--; w++)
    {
    s=0.0;
    if (*w)
      {
      bp = b;
      for (i=m; i--;)
        s += *(ap++)**(bp++);
      s /= *w;
      }
    else
      ap += m;
    *(tmpp++) = s;
    }

  vp0 = vmat;
  for (j=0; j<n; j++,vp0++)
    {
    s=0.0;
    tmpp = tmp;
    for (vp=vp0,jj=n; jj--; vp+=n)
      s += *vp**(tmpp++);
    sol[j]=s;
    }

/* Free temporary arrays */
  free(tmp);
  free(rv1);

  return;
  }

#undef SIGN
#undef MAX
#undef PYTHAG
#undef TOL

/******************************** svdvar ************************************/
/*
Computation of the covariance matrix from the SVD vmat and wmat matrices.A
dapted from Numerical Recipes in C, 2nd Ed. (p. 679).
*/
void svdvar(double *v, double *w, int n, double *cov)
  {
   static double	wti[PSF_NTOT];
   double		sum;
   int			i,j,k;

  for (i=0; i<n; i++)
    wti[i] = w[i]? 1.0/(w[i]*w[i]) : 0.0;

  for (i=0; i<n; i++)
    for (j=0; j<=i; j++)
      {
      for (sum=0.0,k=0; k<n; k++)
        sum += v[k*n+i]*v[k*n+j]*wti[k];
      cov[j*n+i] = cov[i*n+j] = sum;
      }

  return;
  }


/******************************interpolate_pixZERO****************************/
/**
 *
 * Function: interpolate_pixZERO
 *
 * Basically a copy of 'interpolate_pix()', but with a different boundary
 * condition: here the grid given in 'pix' is extrapolated with values 0.0.
 * 'interpolate_pix()' only interpolates positions which don't have outside
 * positions within the kernel area.
 *
 * @author EB, MK
 * @date   April 2014
 *
 * @param[in] posin      - the position to get the value for
 * @param[in] pix        - the interpolation matrix
 * @param[in] naxisn     - the axes length of the interp. matrix
 * @param[in] interptype - the type of interpolation
 *
 * @return the interpolated value
 */
static float    interpolate_pixZERO(float *posin, float *pix, int *naxisn,
    interpenum interptype)
{
  float        buffer[INTERP_MAXKERNELWIDTH],
  kernel[INTERP_MAXKERNELWIDTH], dpos[2],
  *kvector, *pixin, *pixout,
  val;
  int          kwidth, i,j, n, ijstart[2], iact, jact;

  // get the characteristic kernel length
  kwidth = interp_kernwidth_psf[interptype];

  // go over all dimension
  for (n=0; n<2; n++)
    {
      // get the in-value
      val = *(posin++);

      //get the integer part of the current coordinate or nearest neighbour
      ijstart[n] = (interptype==INTERP_NEARESTNEIGHBOUR)? (int)(val-0.50001):(int)val;

      // store the fractional part of the current coordinate
      dpos[n] = val - ijstart[n];

      // store the starting index
      ijstart[n]-=kwidth/2;
    }

  // compute the interpolation kernel for x
  make_kernel(dpos[0], kernel, interptype);

  // first step: interpolate along NAXIS1 from the data themselves
  // go over all y-values
  pixout = buffer;
  for (j=kwidth, jact=ijstart[1]; j--; jact++)
    {
      // outside pixels do not contribute (have value=0.0)
      if (jact<0 || jact>=naxisn[1]){
          // enhance the pointer
          *(pixout++) = 0.0;
      }
      else {
          // go over all x-values
          val = 0.0;
          kvector = kernel;
          for (i=kwidth, iact=ijstart[0]; i--; iact++){

              // outside pixels do not contribute (have value=0.0)
              if (iact<0 || iact>=naxisn[0]){
                  // enhance the pointer
                  kvector++;
              }
              else {
                  // add the contribution
                  val += *(kvector++)*pix[jact*naxisn[0]+iact];
              }
          }
          // transfer the final value
          *(pixout++) = val;
      }
    }

  // as of now, buffer contains the values interpolated at
  // the correct x-value for all relevant y-values

  // compute the interpolation kernel for y
  make_kernel(dpos[1], kernel, interptype);

  // apply the kernel y-kernel values
  // to the values stored in buffer
  pixin = buffer;
  val = 0.0;
  kvector = kernel;
  for (i=kwidth; i--;){
      val += *(kvector++)**(pixin++);
  }

  // return the interpolated value
  return val;
}

/****** interpolate_pix ******************************************************
PROTO   void interpolate_pix(float *posin, float *pix, int naxisn,
                interpenum interptype)
PURPOSE Interpolate a model profile at a given position.
INPUT   Profile structure,
        input position vector,
        input pixmap dimension vector,
        interpolation type.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 07/12/2006
 ***/
static float    interpolate_pix(float *posin, float *pix, int *naxisn,
                        interpenum interptype)
  {
   float        buffer[INTERP_MAXKERNELWIDTH],
                kernel[INTERP_MAXKERNELWIDTH], dpos[2],
                *kvector, *pixin, *pixout,
                val;
   int          fac, ival, kwidth, start, width, step,
                i,j, n;

  kwidth = interp_kernwidth_psf[interptype];
  start = 0;
  fac = 1;
  for (n=0; n<2; n++)
    {
    val = *(posin++);
    width = naxisn[n];
/*-- Get the integer part of the current coordinate or nearest neighbour */
    ival = (interptype==INTERP_NEARESTNEIGHBOUR)? (int)(val-0.50001):(int)val;
/*-- Store the fractional part of the current coordinate */
    dpos[n] = val - ival;
/*-- Check if interpolation start/end exceed image boundary... */
    ival-=kwidth/2;
    if (ival<0 || ival+kwidth<=0 || ival+kwidth>width)
      return 0.0;
/*-- Update starting pointer */
    start += ival*fac;
/*-- Update step between interpolated regions */
    fac *= width;
    }

/* First step: interpolate along NAXIS1 from the data themselves */
  make_kernel(dpos[0], kernel, interptype);
  step = naxisn[0]-kwidth;
  pixin = pix+start;
  pixout = buffer;
  for (j=kwidth; j--;)
    {
    val = 0.0;
    kvector = kernel;
    for (i=kwidth; i--;)
      val += *(kvector++)**(pixin++);
    *(pixout++) = val;
    pixin += step;
    }

/* Second step: interpolate along NAXIS2 from the interpolation buffer */
  make_kernel(dpos[1], kernel, interptype);
  pixin = buffer;
  val = 0.0;
  kvector = kernel;
  for (i=kwidth; i--;)
    val += *(kvector++)**(pixin++);

  return val;
  }


/****** make_kernel **********************************************************
PROTO	void make_kernel(float pos, float *kernel, interpenum interptype)
PURPOSE	Conpute interpolation-kernel data
INPUT	Position,
	Pointer to the output kernel data,
	Interpolation method.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/07/2011
 ***/
void	make_kernel(float pos, float *kernel, interpenum interptype)
  {
   float	x, val, sinx1,sinx2,sinx3,cosx1;

  if (interptype == INTERP_NEARESTNEIGHBOUR)
    *kernel = 1;
  else if (interptype == INTERP_BILINEAR)
    {
    *(kernel++) = 1.0-pos;
    *kernel = pos;
    }
  else if (interptype == INTERP_LANCZOS2)
    {
    if (pos<1e-5 && pos>-1e-5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/2.0*(pos+1.0);
#ifdef HAVE_SINCOSF
      sincosf(x, &sinx1, &cosx1);
#else
      sinx1 = sinf(x);
      cosx1 = cosf(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/2.0;
      val += (*kernel = cosx1/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else if (interptype == INTERP_LANCZOS3)
    {
    if (pos<1e-5 && pos>-1e-5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/3.0*(pos+2.0);
#ifdef HAVE_SINCOSF
      sincosf(x, &sinx1, &cosx1);
#else
      sinx1 = sinf(x);
      cosx1 = cosf(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx2=-0.5*sinx1-0.866025403785*cosx1)
				/ (x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx3=-0.5*sinx1+0.866025403785*cosx1)
				/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = sinx1/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/3.0;
      val += (*kernel = sinx3/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else if (interptype == INTERP_LANCZOS4)
    {
    if (pos<1e-5 && pos>-1e-5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/4.0*(pos+3.0);
#ifdef HAVE_SINCOSF
      sincosf(x, &sinx1, &cosx1);
#else
      sinx1 = sinf(x);
      cosx1 = cosf(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/4.0;
      val +=(*(kernel++) = -(sinx2=0.707106781186*(sinx1+cosx1))
				/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = cosx1/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -(sinx3=0.707106781186*(cosx1-sinx1))/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/4.0;
      val += (*kernel = sinx3/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else
    error(EXIT_FAILURE, "*Internal Error*: Unknown interpolation type in ",
		"make_kernel()");

  return;
  }
