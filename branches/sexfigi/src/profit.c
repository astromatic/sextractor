 /*
 				profit.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Authors:	E.BERTIN (IAP)
*
*	Contents:	Fit an arbitrary profile combination to a detection.
*
*	Last modify:	15/11/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
#include	"image.h"
#include	"profit.h"
#include	"psf.h"

/*------------------------------- variables ---------------------------------*/

/****** prof_init *************************************************************
PROTO	profstruct prof_init(void)
PURPOSE	Allocate and initialize a new profile structure.
INPUT	-.
OUTPUT	A pointer to an allocated prof structure.
AUTHOR	E. Bertin (IAP)
VERSION	09/11/2006
 ***/
profstruct	*prof_init(void)
  {
   profstruct		*profit;

  QMALLOC(prof, profstruct, 1);

  return prof;
  }  


/****** prof_end **************************************************************
PROTO	void prof_end(profstruct *prof)
PURPOSE	End (deallocate) a profile structure.
INPUT	Prof structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/11/2006
 ***/
void	prof_end(profstruct *profit)
  {
  free(prof);

  return;
  }


/****** prof_residuals ********************************************************
PROTO	double *prof_residuals(profstruct *prof, int nprof, psfstruct *psf,
			picstruct *field, picstruct *wfield, objstruct *obj)
PURPOSE	Compute the vector of residuals between the data and the galaxy
	profile model.
INPUT	Array of profile structures,
	Number of profiles,
	Pointer to a PSF structure,
	Pointer to the field,
	Pointer to the field weight,
	Pointer to the obj.
OUTPUT	Vector of residuals.
AUTHOR	E. Bertin (IAP)
VERSION	15/11/2006
 ***/
double	*prof_residuals(profitstruct *profit, picstruct *field,
		picstruct *wfield, objstruct *obj)
  {

  nprof = profit->nprof;
  pixout = profit->pix;
  memset(pixout, 0, profit->naxisn[0]*profit->naxisn[1]*sizeof(double));
  for (p=0; p<nprof; p++)
    {
    prof_add(prof[p], profit);
    }
  convolve_profit(profit);
  compresi_profit(profit, field, wfield, obj);

  return profit->resi;
  }


/****** prof_add **************************************************************
PROTO	void prof_add(profstruct *prof, profitstruct *profit)
PURPOSE	Add a model profile to an image.
INPUT	Profile structure,
	profile-fitting structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/11/2006
 ***/
void	prof_add(profstruct *prof, profitstruct *profit)
  {
  naxis = 2;
  npix = prof->naxisn[0]*prof->naxisn[1];
  nextra = prof->naxis-2;
  if (nextra>0)
/*-- Create an interpolated image of model */
    prof_interpextra(prof);
  else
/*-- No extra axis: one single pixmap */
    memcpy(prof->comppix, prof->pix, npix*sizeof(double));

/* Compute Profile CD matrix */
  ctheta = cos(*prof->posangle*DEG);
  stheta = sin(*prof->posangle*DEG);
  cd11 = *prof->scale[0]*ctheta;
  cd12 =-*prof->scale[1]*stheta;
  cd21 = *prof->scale[0]*stheta;
  cd22 = *prof->scale[1]*ctheta;
  dx1 = prof->x[0];
  dx2 = prof->x[1];

/* Initialize multi-dimensional counters */
  for (d=0; d<naxis; d++)
    {
    posout[d] = 0.0;
    dnaxisn[d] = prof->naxisn[d] - 0.00001;
    }

/* Remap each pixel */
  for (i=npix; i--;)
    {
    x1 = posout[0] - dx1;
    x2 = posout[1] - dx2;
    posin[0] = cd11*x1 + cd12*x2;
    posin[1] = cd21*x1 + cd22*x2;
    for (d=0; d<naxis; d++)
      if ((posout[d]+=1.0) < dnaxisn[d])
        break;
      else
        posout[d] = 0.0;
    }

  return;
  }


/****** prof_mapcoords ********************************************************
PROTO	void prof_mapcoords(profstruct *prof, double *posout, double *posin)
PURPOSE	Add a model profile to an image.
INPUT	Profile structure,
	profile-fitting structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/11/2006
 ***/
void	prof_mapcoords(profstruct *prof, profitstruct *profit)
  {
  naxis = prof->naxis;
  npix = prof->naxisn[0]*prof->naxisn[1];
  nextra = prof->naxis-2;
  if (nextra>0)
/*-- Create an interpolated image of model */
    prof_interpextra(prof);
  else
/*-- No extra axis: one single pixmap */
    memcpy(prof->comppix, prof->pix, npix*sizeof(double));

/* Initialize multi-dimensional counters */
  for (d=0; d<naxis; d++)
    {
    posout[d] = 0.0;
    dnaxisn[d] = prof->naxisn[d] - 0.00001;
    }

/* Remap each pixel */
  for (i=npix; i--;)
    {
    prof_coordmap(prof, posout, posin);
    for (d=0; d<naxis; d++)
      if ((posout[d]+=1.0) < dnaxisn[d])
        break;
      else
        posout[d] = 0.0;
    }

  return;
  }


/****** prof_interpextra ******************************************************
PROTO	void prof_interpextra(profstruct *prof)
PURPOSE	Interpolate the "extra" components of a model profile.
INPUT	Profile structure.
OUTPUT	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/11/2006
 ***/
void	prof_interpextra(profstruct *prof)
  {

  npix = prof->naxisn[0]*prof->naxisn[1];
  nextra = prof->naxis-2;
/* Interpolate pixmap over all extra axes */
  inc = npix;
  pixpoint = prof->pix;
  for (e=0; e<nextra; e++)
    {
    extralen = prof->naxisn[e+2];
/*-- Compute position along axis */
    pos = prof->extrapos[e]*prof->extrascale[e]+prof->extraoffset[e];
    proint = (int)pos;
/*-- Keep position within boundaries and let interpolation do the rest */
    if (prof->extracycleflag)
      {
      proint = proint % extralen;
      if (proint < 0)
        proint += extralen;
      }
    else
      {
      if (proint < 0)
        proint = 0;
      else if (proint >= extralen)
        proint = extralen - 1;
      }
    profrac[e] = pos - proint;
    pixpoint += proint*inc;
    pixstep[e] = inc;
    inc *= extralen;
    count[e] = 0;
    }
  nprof = 2 << nextra;	/* Nb of pixmaps involved in interpolation */
  for (n=0; n<nprof; n++)
    {
    fac = 0.0;
    pixint = pixpoint;
    for (e=0; e<nextra; e++)
      {
      fac *= count[e]? 1.0 - profac[e]: profac[e];
      pixint += count[e]*pixstep[e];
      }
/*-- The actual multi-linear interpolation */
    pixout = prof->comppix;
    for (i=npix; i--;)
      *(pixout++) += fac**(pixint++);
    for (e=0; e<nextra; e++)
      if ((++count[e])<2)
        break;
      else
        count[e] = 0;
    }

  return;
  }


/****** prof_fit ************************************************************
PROTO	fitstruct *prof_fit(profstruct *prof, int nprof, psfstruct *psf,
		picstruct *field, picstruct *wfield, objstruct *obj)
PURPOSE	Fit profile(s) convolved with the PSF to a detected object.
INPUT	Array of profile structures,
	Number of profiles,
	Pointer to a PSF structure,
	Pointer to the field,
	Pointer to the field weight,
	Pointer to the obj.
OUTPUT	Pointer to an allocated fit structure (containing details about the
	fit).
AUTHOR	E. Bertin (IAP)
VERSION	10/11/2006
 ***/
fitstruct	*prof_fit(profstruct *prof, int nprof, psfstruct *psf,
		picstruct *field, picstruct *wfield, objstruct *obj)
  {

  return;
  }
