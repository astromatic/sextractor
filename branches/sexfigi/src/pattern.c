 /*
 				pattern.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Authors:	E.BERTIN (IAP)
*
*	Contents:	Generate and handle image patterns for image fitting.
*
*	Last modify:	09/09/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#define _GNU_SOURCE
#include <math.h>
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"pattern.h"
#include	"profit.h"

/*------------------------------- variables ---------------------------------*/


/****** pattern_init ***********************************************************
PROTO	patternstruct pattern_init(pattern_type ptype, int nvec)
PURPOSE	Allocate and initialize a new pattern structure.
INPUT	Pattern type,
	Number of vectors.
OUTPUT	Pointer to the new pattern structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	12/09/2008
 ***/
patternstruct	*pattern_init(pattypenum ptype, int nvec)
  {
   patternstruct	*pattern;

  QCALLOC(pattern, patternstruct, 1);
  pattern->type = ptype;
  pattern->size[2] = nvec*2;

  return pattern;
  }  


/****** pattern_end ***********************************************************
PROTO	void pattern_end(patternstruct *pattern)
PURPOSE	End (deallocate) a pattern structure.
INPUT	Pattern structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/09/2008
 ***/
void	pattern_end(patternstruct *pattern)
  {
  free(pattern->pix);
  free(pattern);

  return;
  }


/****** pattern_fit ******************************************************
PROTO	void pattern_resample(patternstruct *pattern)
PURPOSE	Resample a pattern structure.
INPUT	Pointer to pattern structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	12/09/2008
 ***/
void	pattern_fit(patternstruct *pattern, profitstruct *profit)
  {
   int	npix;

  npix = profit->modnaxisn[0]*profit->modnaxisn[1] * pattern->size[2];
  pattern->aspect = *profit->paramlist[PARAM_DISK_ASPECT];
  pattern->posangle = fmod_m90_p90(*profit->paramlist[PARAM_DISK_POSANG]);
  pattern->scale = *profit->paramlist[PARAM_DISK_SCALE]/profit->pixstep;
  QMALLOC(profit->patpix, PIXTYPE, npix);
  QMALLOC(profit->lpatpix, PIXTYPE, npix);

  pattern_add(pattern, profit);

  free(profit->patpix);
  free(profit->lpatpix);

  return;
  }


/****** pattern_add ******************************************************
PROTO	void pattern_add(patternstruct *pattern, profstruct *profit)
PURPOSE	Resample a pattern structure.
INPUT	Pointer to pattern structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	12/09/2008
 ***/
void	pattern_add(patternstruct *pattern, profitstruct *profit)
  {
catstruct *cat;
   double		x1,x2, x1t,x2t, r2,r2min,r2max, lr, lr0, 
			mod,ang, cosang,sinang, angcoeff,
			ctheta,stheta, saspect,xscale,yscale,
			cd11,cd12,cd21,cd22, x1cout,x2cout;
   PIXTYPE		*cpix,*spix;
   int			p, ix1,ix2, nvec, npix;

/* Compute Profile CD matrix */
  ctheta = cos(pattern->posangle*DEG);
  stheta = sin(pattern->posangle*DEG);
  saspect = sqrt(fabs(pattern->aspect));
  xscale = (pattern->scale==0.0)? 0.0 : 1.0/fabs(pattern->scale);
  yscale = (pattern->scale*saspect == 0.0)?
			0.0 : 1.0/fabs(pattern->scale*saspect);
  cd11 = xscale*ctheta;
  cd12 = xscale*stheta;
  cd21 = -yscale*stheta;
  cd22 = yscale*ctheta;
 
  x1cout = (double)(profit->modnaxisn[0]/2);
  x2cout = (double)(profit->modnaxisn[1]/2);

  switch(pattern->type)
    {
    case PATTERN_QUADRUPOLE:
    case PATTERN_OCTOPOLE:
      nvec = pattern->size[2]/2;
      npix = profit->modnaxisn[0]*profit->modnaxisn[1];
      r2min = fabs(cd11*cd22-cd12*cd21)/10.0;
      r2max = BIG;
      cpix = profit->patpix;
      spix = profit->patpix+npix;
      angcoeff = (pattern->type==PATTERN_OCTOPOLE)? 4.0 : 2.0;
      for (p=0; p<nvec; p++, cpix+=npix, spix+=npix)
        {
        x1 = -x1cout;
        x2 = -x2cout;
        lr0 = log(5.0*(p+1)*(p+1)/nvec/nvec);
        for (ix2=profit->modnaxisn[1]; ix2--; x2+=1.0)
          {
          x1t = cd12*x2 + cd11*x1;
          x2t = cd22*x2 + cd21*x1;
          for (ix1=profit->modnaxisn[0]; ix1--;)
            {
            r2 = x1t*x1t+x2t*x2t;
            if (r2<r2max)
              {
              lr = 8.0*(0.5*log(r2 > r2min ? r2 : r2min)-lr0);
              mod = exp(-0.5*lr*lr);
              ang = angcoeff*atan2(x2t,x1t);
#ifdef HAVE_SINCOS
              sincos(ang, &sinang, &cosang);
#else
              sinang = sin(ang);
              cosang = cos(ang);
#endif
              *(cpix++) = mod*cosang;
              *(spix++) = mod*sinang;
              }
            else
              *(cpix++) = *(spix++) = 0.0;
            x1t += cd11;
            x2t += cd21;
            }
          }
        }
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown Pattern type","");
    }

bswapflag =1;
cat=new_cat(1);
init_cat(cat);
cat->tab->naxis=3;
QMALLOC(cat->tab->naxisn, int, 3);
cat->tab->naxisn[0]=profit->modnaxisn[0];
cat->tab->naxisn[1]=profit->modnaxisn[1];
cat->tab->naxisn[2]=pattern->size[2];
cat->tab->bitpix=BP_FLOAT;
cat->tab->bytepix=4;
cat->tab->bodybuf=profit->patpix;
cat->tab->tabsize=cat->tab->naxisn[0]*cat->tab->naxisn[1]*cat->tab->naxisn[2]*sizeof(PIXTYPE);
save_cat(cat, "toto2.fits");
cat->tab->bodybuf=NULL;
free_cat(&cat, 1);

  return;
  }

