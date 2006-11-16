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
VERSION	16/11/2006
 ***/
void	prof_add(profstruct *prof, profitstruct *profit)
  {
  npix = prof->naxisn[0]*prof->naxisn[1];
  nextra = prof->naxis-2;
  if (nextra>0)
    {
    inc = npix;
    for (i=0; i<nextra; i++)
      {
      profrac1[i] = 1.0 - prof->extra[i] + (proint=(int)prof->extra[i]));
      point1[i] = prof->pix[proint*inc];
      point2[i] = point1[i] + inc;
      for (i=npix; i--;)
        {
      inc *= prof->naxisn[i+2];
      }
      
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
