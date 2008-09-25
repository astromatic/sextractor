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
*	Last modify:	25/09/2008
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
#include	"fitswcs.h"
#include	"check.h"
#include	"pattern.h"
#include	"profit.h"

/*------------------------------- variables ---------------------------------*/

/****** pattern_init ***********************************************************
PROTO	patternstruct pattern_init(profitstruct *profit, pattern_type ptype,
				int ncomp)
PURPOSE	Allocate and initialize a new pattern structure.
INPUT	Pointer to a profit structure,
	Pattern type,
	Number of independent components.
OUTPUT	Pointer to the new pattern structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/09/2008
 ***/
patternstruct	*pattern_init(profitstruct *profit, pattypenum ptype, int ncomp)
  {
   patternstruct	*pattern;
   int			ninpix, noutpix;

  if (!ncomp)
    ncomp = PATTERN_NCOMP;
  QCALLOC(pattern, patternstruct, 1);
  pattern->type = ptype;
  pattern->ncomp = ncomp;
  pattern->size[0] = profit->modnaxisn[0];
  pattern->size[1] = profit->modnaxisn[1];
  switch(pattern->type)
    {
    case PATTERN_QUADRUPOLE:
    case PATTERN_OCTOPOLE:
      pattern->nmodes = 2;
      pattern->nfreq = 1;
      break;
    case PATTERN_POLARFOURIER:
      pattern->nfreq = PATTERN_FMAX+1;
      pattern->nmodes = 2*PATTERN_FMAX+1;
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown Pattern type","");
    }

  pattern->size[2] = ncomp*pattern->nmodes;
  ninpix = pattern->size[0]*pattern->size[1] * pattern->size[2];
  noutpix = profit->objnaxisn[0]*profit->objnaxisn[1] * pattern->size[2];
  pattern->aspect = *profit->paramlist[PARAM_DISK_ASPECT];
  pattern->posangle = fmod_m90_p90(*profit->paramlist[PARAM_DISK_POSANG]);
  pattern->scale = *profit->paramlist[PARAM_DISK_SCALE]/profit->pixstep;
  QMALLOC(pattern->modpix, double, ninpix);
  QMALLOC(pattern->lmodpix, PIXTYPE, noutpix);
  QMALLOC(pattern->coeff, double, pattern->size[2]);
  QMALLOC(pattern->mcoeff, double, ncomp*pattern->nfreq);
  QMALLOC(pattern->acoeff, double, ncomp*pattern->nfreq);

  return pattern;
  }  


/****** pattern_end ***********************************************************
PROTO	void pattern_end(patternstruct *pattern)
PURPOSE	End (deallocate) a pattern structure.
INPUT	Pattern structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/09/2008
 ***/
void	pattern_end(patternstruct *pattern)
  {
  free(pattern->modpix);
  free(pattern->lmodpix);
  free(pattern->coeff);
  free(pattern->mcoeff);
  free(pattern->acoeff);
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
VERSION	22/09/2008
 ***/
void	pattern_fit(patternstruct *pattern, profitstruct *profit)
  {
catstruct *cat;
char	name[MAXCHAR];
static int number;
   checkstruct	*check;
   double	*inpix, *doutpix1, *alpha,*beta,
		dval, dprod;
   PIXTYPE	*outpix,*outpix1,*outpix2;
   PIXTYPE	*weightpix;
   int		n,p,p2, nvec, ninpix, noutpix,nout;

  nvec = pattern->size[2];
  pattern_create(pattern, profit);
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
    alpha[p*(nvec+1)] += 1.0;
    beta[p] = dval;
    inpix += ninpix;
    outpix += noutpix;
    }

/* Solve the system */
  clapack_dpotrf(CblasRowMajor,CblasUpper,nvec,alpha,nvec);
  clapack_dpotrs(CblasRowMajor,CblasUpper,nvec,1,alpha,nvec,beta,nvec);

  pattern_compmodarg(pattern);

  free(alpha);

  if ((check = prefs.check[CHECK_PATTERNS]))
    {
    QCALLOC(outpix, PIXTYPE, noutpix);
    outpix2 = pattern->lmodpix;
    for (p=0; p<nvec; p++)
      {
      dval = pattern->coeff[p];
      outpix1 = outpix;
      for (n=noutpix; n--;)
        *(outpix1++) += dval**(outpix2++);
      }
    addcheck(check, outpix, profit->objnaxisn[0],profit->objnaxisn[1],
		profit->ix, profit->iy, 1.0);
    free(outpix);
    }
/*
nout = pattern->ncomp*pattern->nfreq;
QCALLOC(outpix, PIXTYPE, noutpix*nout);
outpix1 = outpix;
outpix2 = pattern->lmodpix;
for (p=0; p<nvec; p++)
{
dval = pattern->coeff[p];
for (n=noutpix; n--; )
*(outpix1++) += dval**(outpix2++);
if (pattern->type==PATTERN_POLARFOURIER)
  {
  if ((p%pattern->nmodes)%2)
    outpix1 -= noutpix;
  }
else if (!(p%2))
  outpix1 -= noutpix;
}
cat=new_cat(1);
init_cat(cat);
cat->tab->naxis=3;
QMALLOC(cat->tab->naxisn, int, 3);
cat->tab->naxisn[0]=profit->objnaxisn[0];
cat->tab->naxisn[1]=profit->objnaxisn[1];
cat->tab->naxisn[2]=nout;
cat->tab->bitpix=BP_FLOAT;
cat->tab->bytepix=4;
cat->tab->bodybuf=outpix;
cat->tab->tabsize=cat->tab->naxisn[0]*cat->tab->naxisn[1]*cat->tab->naxisn[2]*sizeof(PIXTYPE);
sprintf(name, "tata_%02d.fits", ++number);
save_cat(cat, name);
cat->tab->bodybuf=NULL;
free_cat(&cat, 1);
free(outpix);
*/
  return;
  }


/****** pattern_compmodarg ****************************************************
PROTO	void pattern_comparg(patternstruct *pattern)
PURPOSE	Compute modulus and argument for each pair of Fourier components.
INPUT	Pointer to pattern structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	19/09/2008
 ***/
void	pattern_compmodarg(patternstruct *pattern)
  {
   double	*coeff,*mcoeff,*acoeff,
		arg,argo,darg, ima,rea;
   int		f,r, nfreq;

  coeff = pattern->coeff;
  mcoeff = pattern->mcoeff;
  acoeff = pattern->acoeff;
  nfreq = pattern->nfreq;
  for (r=0; r<pattern->ncomp; r++)
    {
    for (f=0; f<pattern->nfreq; f++)
      {
      if (pattern->type == PATTERN_POLARFOURIER && !f)
        {
        *(mcoeff++) = fabs(*coeff);
        *(acoeff++) = *(coeff++)<0.0? 180.0 : 0.0;
        }
      else
        {
        rea = *(coeff++);
        ima = *(coeff++);
        *(mcoeff++) = sqrt(rea*rea + ima*ima);
        arg = atan2(ima, rea)/DEG;
        if (r>0)
          {
          darg = arg - argo;
/*-------- disambiguate increasing or decreasing phase angles */
          if (darg > 180.0)
            darg -= 360.0;
          else if (darg < -180.0)
            darg += 360.0;
          *acoeff = *(acoeff-nfreq) + darg;
          acoeff++;
          }
        else
          *(acoeff++) = arg;
        argo = arg;
        }
      }
    }

  return;
  }


/****** pattern_spiral ******************************************************
PROTO	float pattern_spiral(patternstruct *pattern)
PURPOSE	Compute a pattern spiral index.
INPUT	Pointer to pattern structure.
OUTPUT	Spiral index.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/09/2008
 ***/
float	pattern_spiral(patternstruct *pattern)
  {
   double	w,x,y, s,sx,sy,sxx,sxy;
   int		i,r, freq, rstart;

  if (pattern->ncomp<2)
    return 0.0;

  rstart = (int)(pattern->ncomp/pattern->rmax+0.4999) - 1;
  if (rstart<0)
    rstart = 0;
  else if (rstart>pattern->ncomp-2)
    rstart = pattern->ncomp-2;

  freq = (pattern->type == PATTERN_POLARFOURIER) ? 2 : 0;
  s = sx = sy = sxx = sxy = 0.0;
  for (r=rstart; r<pattern->ncomp; r++)
    {
    i = r*pattern->nfreq + freq;
    w = pattern->mcoeff[i];
    x = (double)(r - rstart);
    y = pattern->acoeff[i];
    s += w;
    sx += w*x;
    sy += w*y;
    sxx += w*x*x;
    sxy += w*x*y;
    }

  return (s*sxy - sx*sy)/(s*sxx - sx*sx);
  }


/****** pattern_create ******************************************************
PROTO	void pattern_create(patternstruct *pattern, profitstruct *profit)
PURPOSE create a pattern basis.
INPUT	Pointer to pattern structure,
	pointer to the profit structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/09/2008
 ***/
void	pattern_create(patternstruct *pattern, profitstruct *profit)
  {
   double		x1,x2, x1t,x2t, r2,r2min,r2max, lr, lr0, 
			mod,ang,ang0, cosang,sinang, angcoeff,
			ctheta,stheta, saspect,xscale,yscale,
			cd11,cd12,cd21,cd22, x1cout,x2cout, cmod,smod,
			cnorm,snorm,norm, dval, det, rad, dnrad, rscale2;
   double		*scbuf[PATTERN_FMAX],*scpix[PATTERN_FMAX],
			*scpixt,*cpix,*spix, *pix, *r2buf,*r2pix,*modpix;
   int			f,i,p, ix1,ix2, nrad, npix;

/* Compute Profile CD matrix */
  ctheta = cos(pattern->posangle*DEG);
  stheta = sin(pattern->posangle*DEG);
  saspect = fabs(pattern->aspect);
  xscale = (pattern->scale==0.0)? 0.0 : 1.0/fabs(pattern->scale);
  yscale = (pattern->scale*saspect == 0.0)?
			0.0 : 1.0/fabs(pattern->scale*saspect);
  cd11 = xscale*ctheta;
  cd12 = xscale*stheta;
  cd21 = -yscale*stheta;
  cd22 = yscale*ctheta;
 
  x1cout = (double)(pattern->size[0]/2);
  x2cout = (double)(pattern->size[1]/2);
/* Determinant of the change of coordinate system */
  det = xscale*yscale;
  r2min = det/10.0;
/* Stay within an ellipse contained in the pattern raster, both in x and y */
  r2max = x1cout*x1cout * det*det / (cd12*cd12+cd22*cd22);
  if (r2max > (dval = x2cout*x2cout * det*det / (cd21*cd21+cd11*cd11)))
    r2max = dval;
/* Set the limit of the pattern extent */
//  rad = 4.0*profit->obj->a*xscale;
/* The pattern limit does not exceed 90% of the mapped ellipse "radius" */
//  if (rad*rad > 0.9*0.9*r2max)
  nrad = pattern->ncomp;
  pattern->rmax = rad = 0.8*sqrt(r2max);/* Keep a margin using a fudge factor */
  if (!nrad)
      error(EXIT_FAILURE,
		"*Error*: insufficient number of vector elements",
		" for generating the pattern basis"); 
  dnrad = (double)nrad;
  npix = pattern->size[0]*pattern->size[1];
  switch(pattern->type)
    {
    case PATTERN_QUADRUPOLE:
    case PATTERN_OCTOPOLE:
      cpix = pattern->modpix;
      spix = pattern->modpix+npix;
      angcoeff = (pattern->type==PATTERN_OCTOPOLE)? 4.0 : 2.0;
      for (p=0; p<nrad; p++, cpix+=npix, spix+=npix)
        {
        rscale2 = (p+1)*dnrad;
        x1 = -x1cout;
        x2 = -x2cout;
        lr0 = log(rad*(p+1)/dnrad);
        cnorm = snorm = 0.0;
        for (ix2=pattern->size[1]; ix2--; x2+=1.0)
          {
          x1t = cd12*x2 + cd11*x1;
          x2t = cd22*x2 + cd21*x1;
          for (ix1=pattern->size[0]; ix1--;)
            {
            r2 = x1t*x1t+x2t*x2t;
            if (r2<r2max)
              {
              lr = 0.5*log(r2 > r2min ? r2 : r2min)-lr0;
              mod = exp(-0.5*rscale2*lr*lr);
              ang = angcoeff*atan2(x2t,x1t);
#ifdef HAVE_SINCOS
              sincos(ang, &sinang, &cosang);
#else
              sinang = sin(ang);
              cosang = cos(ang);
#endif
              *(cpix++) = cmod = mod*cosang;
              *(spix++) = smod = mod*sinang;
              cnorm += cmod*cmod;
              snorm += smod*smod;
              }
            else
              *(cpix++) = *(spix++) = 0.0;
            x1t += cd11;
            x2t += cd21;
            }
          }
        cpix -= npix;
        cnorm = cnorm > 0.0? 1.0/sqrt(cnorm) : 0.0;
        for (i=npix; i--;)
          *(cpix++) *= cnorm;
        spix -= npix;
        snorm = snorm > 0.0? 1.0/sqrt(snorm) : 0.0;
        for (i=npix; i--;)
          *(spix++) *= snorm;
        }
      break;
    case PATTERN_POLARFOURIER:
/*---- Pre-compute radii and quadrupoles to speed up computations later */
      QMALLOC(r2buf, double, npix);
      r2pix = r2buf;
      for (f=0; f<PATTERN_FMAX; f++)
        {
        QMALLOC(scbuf[f], double, 2*npix);
        scpix[f] = scbuf[f];
        }
      x1 = -x1cout;
      x2 = -x2cout;
      for (ix2=pattern->size[1]; ix2--; x2+=1.0)
        {
        x1t = cd12*x2 + cd11*x1;
        x2t = cd22*x2 + cd21*x1;
        for (ix1=pattern->size[0]; ix1--;)
          {
          *(r2pix++) = x1t*x1t+x2t*x2t;
          ang = ang0 = atan2(x2t,x1t);
          for (f=0; f<PATTERN_FMAX; f++)
            {
#ifdef HAVE_SINCOS
            sincos(ang, scpix[f]+npix, scpix[f]);
            scpix[f]++;
#else
            *(scpix[f]) = cos(ang);
            *(scpix[f]+++npix) = sin(ang);
#endif
            ang+=ang0;
            }
          x1t += cd11;
          x2t += cd21;
          }
        }
      modpix = NULL;		/* To avoid gcc -Wall warnings */
      pix = pattern->modpix;
      for (p=0; p<nrad; p++)
        {
        for (f=0; f<=PATTERN_FMAX; f++)
          {
          rscale2 = (p+1)*dnrad;
          norm = 0.0;
          lr0 = log(rad*(p+1)/dnrad);
          r2pix = r2buf;
          if (!f)
            {
            for (i=npix; i--;)
              {
              r2 = *(r2pix++);
              if (r2<r2max)
                {
                lr = 0.5*log(r2 > r2min ? r2 : r2min)-lr0;
                *(pix++) = dval = exp(-0.5*rscale2*lr*lr);
                norm += dval*dval;
                }
              else
                *(pix++) = 0.0;
              }
            pix -= npix;
            norm = norm > 1.0/BIG? 1.0/sqrt(norm) : 0.0;
            for (i=npix; i--;)
              *(pix++) *= norm;
            modpix = pix;
            }
          else
            {
            modpix -= npix;
            scpixt = scbuf[f-1];
            for (i=npix; i--;)
              {
              *(pix++) = dval = *(modpix++)**(scpixt++);
              norm += dval*dval;
              }
            pix -= npix;
            norm = norm > 0.0? 1.0/sqrt(norm) : 0.0;
            for (i=npix; i--;)
              *(pix++) *= norm;
            modpix -= npix;
            norm = 0.0;
            for (i=npix; i--;)
              {
              *(pix++) = dval = *(modpix++)**(scpixt++);
              norm += dval*dval;
              }
            pix -= npix;
            norm = norm > 0.0? 1.0/sqrt(norm) : 0.0;
            for (i=npix; i--;)
              *(pix++) *= norm;
            }
          }
        }
      free(r2buf);
      for (f=0; f<PATTERN_FMAX; f++)
        free(scbuf[f]);
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown Pattern type","");
    }

  return;
  }


