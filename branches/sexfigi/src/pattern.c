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
*	Last modify:	16/09/2008
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
#include	ATLAS_LAPACK_H

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"pattern.h"
#include	"profit.h"

/*------------------------------- variables ---------------------------------*/

/****** pattern_init ***********************************************************
PROTO	patternstruct pattern_init(profitstruct *profit, pattern_type ptype,
				int nvec)
PURPOSE	Allocate and initialize a new pattern structure.
INPUT	Pointer to a profit structure,
	Pattern type,
	Number of vectors.
OUTPUT	Pointer to the new pattern structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/09/2008
 ***/
patternstruct	*pattern_init(profitstruct *profit, pattypenum ptype, int nvec)
  {
   patternstruct	*pattern;
   int	npix;

  QCALLOC(pattern, patternstruct, 1);
  pattern->type = ptype;
  pattern->size[0] = profit->modnaxisn[0];
  pattern->size[1] = profit->modnaxisn[1];
  pattern->size[2] = nvec*2;
  npix = pattern->size[0]*pattern->size[1] * pattern->size[2];
  pattern->aspect = *profit->paramlist[PARAM_DISK_ASPECT];
  pattern->posangle = fmod_m90_p90(*profit->paramlist[PARAM_DISK_POSANG]);
  pattern->scale = *profit->paramlist[PARAM_DISK_SCALE]/profit->pixstep;
  QMALLOC(pattern->modpix, double, npix);
  QMALLOC(pattern->lmodpix, PIXTYPE, npix);
  QMALLOC(pattern->coeff, double, pattern->size[2]);

  return pattern;
  }  


/****** pattern_end ***********************************************************
PROTO	void pattern_end(patternstruct *pattern)
PURPOSE	End (deallocate) a pattern structure.
INPUT	Pattern structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/09/2008
 ***/
void	pattern_end(patternstruct *pattern)
  {
  free(pattern->modpix);
  free(pattern->lmodpix);
  free(pattern->coeff);
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
VERSION	15/09/2008
 ***/
void	pattern_fit(patternstruct *pattern, profitstruct *profit)
  {
catstruct *cat;
char	name[MAXCHAR];
static int number;
   double	*inpix, *doutpix1, *alpha,*beta,
		dval, dprod;
   PIXTYPE	*outpix,*outpix1,*outpix2;
   PIXTYPE	*weightpix;
   int		n,p,p2, nvec, ninpix, noutpix;

  nvec = pattern->size[2];
  pattern_create(pattern);
  QMALLOC(alpha, double, nvec*nvec);
  beta = pattern->coeff;
  inpix = pattern->modpix;
  ninpix = pattern->size[0]*pattern->size[1];
  outpix = pattern->lmodpix;
  noutpix = profit->objnaxisn[0]*profit->objnaxisn[1];
  for (p=0; p<nvec; p++)
    {
    profit_convolve(profit, inpix);
    profit_resample(profit, inpix, outpix);
    outpix1 = pattern->lmodpix;
    for (p2=0; p2<=p; p2++)
      {
      weightpix = profit->objweight;
      outpix2 = outpix;
      dval = 0.0;
      for (n=noutpix; n--;)
        {
        dprod = *(outpix1++)**(outpix2++);
        if (*(weightpix++)>0.0)
          dval += dprod;
        }
      alpha[p*nvec+p2] = alpha[p2*nvec+p] = dval;
      }
    weightpix = profit->objweight;
    doutpix1 = profit->resi;
    outpix2 = outpix;
    dval = 0.0;
    for (n=noutpix; n--;)
      {
      dprod = *doutpix1**(outpix2++);
      if (*(weightpix++)>0.0)
        {
        dval += dprod;
        doutpix1++;
        }
      }
    beta[p] = dval;
    inpix += ninpix;
    outpix += noutpix;
    }

/* Solve the system */
  clapack_dpotrf(CblasRowMajor,CblasUpper,nvec,alpha,nvec);
  clapack_dpotrs(CblasRowMajor,CblasUpper,nvec,1,alpha,nvec,beta,nvec);

  free(alpha);

QCALLOC(outpix, PIXTYPE, noutpix);

outpix2 = pattern->lmodpix;
for (p=0; p<nvec; p++)
{
dval = pattern->coeff[p];
outpix1 = outpix;
for (n=noutpix; n--;)
*(outpix1++) += dval**(outpix2++);
}
cat=new_cat(1);
init_cat(cat);
cat->tab->naxis=2;
QMALLOC(cat->tab->naxisn, int, 3);
cat->tab->naxisn[0]=profit->objnaxisn[0];
cat->tab->naxisn[1]=profit->objnaxisn[1];
cat->tab->naxisn[2]=1;
cat->tab->bitpix=BP_FLOAT;
cat->tab->bytepix=4;
cat->tab->bodybuf=outpix;
cat->tab->tabsize=cat->tab->naxisn[0]*cat->tab->naxisn[1]*cat->tab->naxisn[2]*sizeof(PIXTYPE);
sprintf(name, "toto_%02d_c.fits", ++number);
save_cat(cat, name);
cat->tab->bodybuf=NULL;
free_cat(&cat, 1);

weightpix = profit->objweight;
doutpix1 = profit->resi;
outpix1 = outpix;
for (n=noutpix; n--;)
if ((*weightpix++)>0.0)
*(outpix1++) = *(doutpix1++);
else
*(outpix1++)=0.0;
cat=new_cat(1);
init_cat(cat);
cat->tab->naxis=2;
QMALLOC(cat->tab->naxisn, int, 3);
cat->tab->naxisn[0]=profit->objnaxisn[0];
cat->tab->naxisn[1]=profit->objnaxisn[1];
cat->tab->naxisn[2]=1;
cat->tab->bitpix=BP_FLOAT;
cat->tab->bytepix=4;
cat->tab->bodybuf=outpix;
cat->tab->tabsize=cat->tab->naxisn[0]*cat->tab->naxisn[1]*cat->tab->naxisn[2]*sizeof(PIXTYPE);
sprintf(name, "toto_%02d_b.fits", number);
save_cat(cat, name);
cat->tab->bodybuf=NULL;
free_cat(&cat, 1);
cat=new_cat(1);
init_cat(cat);
cat->tab->naxis=2;
QMALLOC(cat->tab->naxisn, int, 3);
cat->tab->naxisn[0]=profit->objnaxisn[0];
cat->tab->naxisn[1]=profit->objnaxisn[1];
cat->tab->naxisn[2]=1;
cat->tab->bitpix=BP_FLOAT;
cat->tab->bytepix=4;
cat->tab->bodybuf=profit->objpix;
cat->tab->tabsize=cat->tab->naxisn[0]*cat->tab->naxisn[1]*cat->tab->naxisn[2]*sizeof(PIXTYPE);
sprintf(name, "toto_%02d_a.fits", number);
save_cat(cat, name);
cat->tab->bodybuf=NULL;
free_cat(&cat, 1);

free(outpix);
QCALLOC(outpix, PIXTYPE, noutpix*nvec/2);
outpix2 = pattern->lmodpix;
for (p=0; p<nvec/2; p++)
{
dval = pattern->coeff[2*p];
outpix1 = outpix+noutpix*p;
for (n=noutpix; n--; )
*(outpix1++) += dval**(outpix2++);
dval = pattern->coeff[2*p+1];
outpix1 = outpix+noutpix*p;
for (n=noutpix; n--;)
*(outpix1++) += dval**(outpix2++);
}
//outpix1 = outpix;
//for (n=noutpix*nvec/2; n--; outpix1++)
//*outpix1 *= *outpix1;

cat=new_cat(1);
init_cat(cat);
cat->tab->naxis=3;
QMALLOC(cat->tab->naxisn, int, 3);
cat->tab->naxisn[0]=profit->objnaxisn[0];
cat->tab->naxisn[1]=profit->objnaxisn[1];
cat->tab->naxisn[2]=nvec/2;
cat->tab->bitpix=BP_FLOAT;
cat->tab->bytepix=4;
cat->tab->bodybuf=outpix;
cat->tab->tabsize=cat->tab->naxisn[0]*cat->tab->naxisn[1]*cat->tab->naxisn[2]*sizeof(PIXTYPE);
sprintf(name, "tata_%02d.fits", number);
save_cat(cat, name);
cat->tab->bodybuf=NULL;
free_cat(&cat, 1);

free(outpix);

  return;
  }


/****** pattern_create ******************************************************
PROTO	void pattern_create(patternstruct *pattern)
PURPOSE create a pattern basis.
INPUT	Pointer to pattern structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	15/09/2008
 ***/
void	pattern_create(patternstruct *pattern)
  {
   double		x1,x2, x1t,x2t, r2,r2min,r2max, lr, lr0, 
			mod,ang, cosang,sinang, angcoeff,
			ctheta,stheta, saspect,xscale,yscale,
			cd11,cd12,cd21,cd22, x1cout,x2cout;
   double		*cpix,*spix;
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
 
  x1cout = (double)(pattern->size[0]/2);
  x2cout = (double)(pattern->size[1]/2);

  switch(pattern->type)
    {
    case PATTERN_QUADRUPOLE:
    case PATTERN_OCTOPOLE:
      nvec = pattern->size[2]/2;
      npix = pattern->size[0]*pattern->size[1];
      r2min = fabs(cd11*cd22-cd12*cd21)/10.0;
      r2max = BIG;
      cpix = pattern->modpix;
      spix = pattern->modpix+npix;
      angcoeff = (pattern->type==PATTERN_OCTOPOLE)? 4.0 : 2.0;
      for (p=0; p<nvec; p++, cpix+=npix, spix+=npix)
        {
        x1 = -x1cout;
        x2 = -x2cout;
        lr0 = log(3.0*(p+1)/nvec);
        for (ix2=pattern->size[1]; ix2--; x2+=1.0)
          {
          x1t = cd12*x2 + cd11*x1;
          x2t = cd22*x2 + cd21*x1;
          for (ix1=pattern->size[0]; ix1--;)
            {
            r2 = x1t*x1t+x2t*x2t;
            if (r2<r2max)
              {
              lr = 20.0*(0.5*log(r2 > r2min ? r2 : r2min)-lr0);
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

  return;
  }


