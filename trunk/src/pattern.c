/*
*				pattern.c
*
* Manage galaxy image patterns.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2007-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		09/07/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
#include	"fitswcs.h"
#include	"check.h"
#include	"pattern.h"
#include	"profit.h"

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
#endif

static double	psf_laguerre(double x, int p, int q);

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
VERSION	18/11/2009
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
      pattern->size[2] = ncomp*pattern->nmodes;
      break;
    case PATTERN_POLARFOURIER:
      pattern->nfreq = PATTERN_FMAX+1;
      pattern->nmodes = 2*PATTERN_FMAX+1;
      pattern->size[2] = ncomp*pattern->nmodes;
      break;
    case PATTERN_POLARSHAPELETS:
      pattern->nfreq = 0;
      pattern->nmodes = 0;
      pattern->size[2] = (ncomp+1)*(ncomp+2)/2;
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown Pattern type","");
    }

  ninpix = pattern->size[0]*pattern->size[1] * pattern->size[2];
  noutpix = profit->objnaxisn[0]*profit->objnaxisn[1] * pattern->size[2];
  QMALLOC(pattern->coeff, float, pattern->size[2]);
  QMALLOC(pattern->norm, float, pattern->size[2]);
  QMALLOC(pattern->modpix, float, ninpix);
  QMALLOC(pattern->lmodpix, PIXTYPE, noutpix);
  if (pattern->ncomp)
    {
    QMALLOC(pattern->r, float, pattern->ncomp);
    }
  if (pattern->nfreq)
    {
    QMALLOC(pattern->mcoeff, float, ncomp*pattern->nfreq);
    QMALLOC(pattern->acoeff, float, ncomp*pattern->nfreq);
    }

  return pattern;
  }  


/****** pattern_end ***********************************************************
PROTO	void pattern_end(patternstruct *pattern)
PURPOSE	End (deallocate) a pattern structure.
INPUT	Pattern structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/10/2008
 ***/
void	pattern_end(patternstruct *pattern)
  {
  free(pattern->norm);
  free(pattern->r);
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
VERSION	09/07/2012
 ***/
void	pattern_fit(patternstruct *pattern, profitstruct *profit)
  {
   checkstruct	*check;
   double	*alpha,*beta,*betat,
		dval, dprod;
   float	*inpix, *doutpix1, *coefft;
   PIXTYPE	*outpix,*outpix1,*outpix2;
   PIXTYPE	*weightpix;
   int		n,p,p2, nvec, ninpix, noutpix;

  nvec = pattern->size[2];
  pattern_create(pattern, profit);
  QMALLOC(alpha, double, nvec*nvec);
  QMALLOC(beta, double, nvec);
  inpix = pattern->modpix;
  ninpix = pattern->size[0]*pattern->size[1];
  outpix = pattern->lmodpix;
  noutpix = profit->objnaxisn[0]*profit->objnaxisn[1];
  for (p=0; p<nvec; p++)
    {
    profit_convolve(profit, inpix);
    profit_resample(profit, inpix, outpix, 1.0);
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
    alpha[p*(nvec+1)] += 0.2;
    beta[p] = dval;
    inpix += ninpix;
    outpix += noutpix;
    }

/* Solve the system */
#if defined(HAVE_LAPACKE)
  LAPACKE_dposv(LAPACK_COL_MAJOR, 'L', nvec, 1, alpha, nvec, beta, nvec);
#else
  clapack_dposv(CblasRowMajor, CblasUpper, nvec, 1, alpha, nvec, beta, nvec);
#endif

  betat = beta;
  coefft = pattern->coeff;
  for (p=nvec; p--;)
    *(coefft++) = (float)*(betat++);

  pattern_compmodarg(pattern, profit);

  free(alpha);
  free(beta);

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
{
catstruct *cat;
char	name[MAXCHAR];
static int number;
int	nout;

nout = nvec;
QCALLOC(outpix, PIXTYPE, ninpix*nout);
outpix1 = outpix;
doutpix1 = pattern->modpix;
for (p=0; p<nvec; p++)
{
dval = pattern->coeff[p];
for (n=ninpix; n--; )
*(outpix1++) += dval**(doutpix1++);
if (pattern->type==PATTERN_POLARFOURIER)
  {
  if ((p%pattern->nmodes)%2)
    outpix1 -= ninpix;
  }
else if (pattern->type==PATTERN_POLARSHAPELETS)
  {
  }
else if (!(p%2))
  outpix1 -= noutpix;
}
cat=new_cat(1);
init_cat(cat);
cat->tab->naxis=3;
QMALLOC(cat->tab->naxisn, int, 3);
cat->tab->naxisn[0]=profit->modnaxisn[0];
cat->tab->naxisn[1]=profit->modnaxisn[1];
cat->tab->naxisn[2]=nout;
cat->tab->bitpix=BP_FLOAT;
cat->tab->bytepix=4;
cat->tab->bodybuf=(char *)outpix;
cat->tab->tabsize=cat->tab->naxisn[0]*cat->tab->naxisn[1]*cat->tab->naxisn[2]*sizeof(PIXTYPE);
sprintf(name, "tata_%02d.fits", ++number);
save_cat(cat, name);
cat->tab->bodybuf=NULL;
free_cat(&cat, 1);
free(outpix);
}
*/
  return;
  }


/****** pattern_compmodarg ****************************************************
PROTO	void pattern_comparg(patternstruct *pattern, profitstruct *profit)
PURPOSE	Compute modulus and argument for each pair of Fourier components.
INPUT	Pointer to pattern structure,
	pointer to profit structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/11/2008
 ***/
void	pattern_compmodarg(patternstruct *pattern, profitstruct *profit)
  {
   float	*coeff,*mcoeff,*acoeff, *normt,
		arg,argo,darg, ima,rea;
   int		f,p, nfreq;

  if (pattern->type == PATTERN_POLARSHAPELETS)
    return;

  coeff = pattern->coeff;
  mcoeff = pattern->mcoeff;
  acoeff = pattern->acoeff;
  nfreq = pattern->nfreq;
  normt = pattern->norm;
  argo = 0.0;			/* To avoid gcc -Wall warnings */
  for (p=0; p<pattern->ncomp; p++)
    {
    for (f=0; f<nfreq; f++)
      {
      if (pattern->type == PATTERN_POLARFOURIER && !f)
        {
        *(mcoeff++) = fabs(*coeff) * *(normt++);
        *(acoeff++) = *(coeff++)<0.0? 180.0 : 0.0;
        }
      else
        {
        rea = *(coeff++) * *(normt++);
        ima = *(coeff++) * *(normt++);
        *(mcoeff++) = sqrt(rea*rea + ima*ima);
        arg = atan2(ima, rea)/DEG;
        if (p>0)
          {
          argo = *(acoeff-nfreq);
          darg = arg - fmod(argo+180.0, 360.0) + 180.0;;
/*-------- disambiguate increasing or decreasing phase angles */
          if (darg > 180.0)
            darg -= 360.0;
          else if (darg < -180.0)
            darg += 360.0;
          *acoeff = argo + darg;
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
VERSION	14/10/2008
 ***/
float	pattern_spiral(patternstruct *pattern)
  {
   double	w,x,y, s,sx,sy,sxx,sxy;
   int		f,i,p, pstart;

  if (pattern->ncomp<2)
    return 0.0;

  pstart = (int)(pattern->ncomp/pattern->rmax+0.4999) - 1;
  if (pstart<0)
    pstart = 0;
  else if (pstart>pattern->ncomp-2)
    pstart = pattern->ncomp-2;

  s = sx = sy = sxx = sxy = 0.0;
  for (p=pstart; p<pattern->ncomp; p++)
    {
    w = y = 0.0;
    for (f=0; f<pattern->nfreq; f++)
      {
      if (pattern->type == PATTERN_POLARFOURIER && (!f || f==1 || f==3))
        continue;
      i = p*pattern->nfreq + f;
      w += pattern->mcoeff[i];
      y += pattern->mcoeff[i]*pattern->acoeff[i]/f;
      }
    x = (double)(p - pstart);
    if (w>0.0)
      y /= w;
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
VERSION	01/12/2009
 ***/
void	pattern_create(patternstruct *pattern, profitstruct *profit)
  {
   double		*scbuf[PATTERN_FMAX],*scpix[PATTERN_FMAX],
			*scpixt,*r2buf,*r2pix,
			x1,x2, x1t,x2t, r,r2,r2min,r2max,
			mod,ang,ang0, cosang,sinang, angcoeff, posangle,flux,
			ctheta,stheta, saspect,xscale,yscale, scale, aspect,
			cd11,cd12,cd21,cd22, x1cout,x2cout, cmod,smod,
			cnorm,snorm,norm,norm0, dval, det, rad, dnrad,
			cpnorm,spnorm,pnorm, rl,rl2,rh,rh2,r0,r02, sbd,
			bt, wb, omwb, bflux, margin2, dposangle;
   int			f,i,p, ix1,ix2, nrad, npix;

   double		*fr2,*fr2t,*fexpr2,*fexpr2t,*ftheta,*fthetat,
			dm,fac, beta, invbeta2;
   float		*normt, *cpix,*spix, *pmodpix,*pix,*modpix,
			fnorm;
   int			m,n, nmax, kmax,hnmm;

/* Compute Profile CD matrix */
  aspect = fabs(*profit->paramlist[PARAM_DISK_ASPECT]);
  posangle = fmod_m90_p90(*profit->paramlist[PARAM_DISK_POSANG])*DEG;
  scale = fabs(*profit->paramlist[PARAM_DISK_SCALE]/profit->pixstep);
  flux = fabs(*profit->paramlist[PARAM_DISK_FLUX])*1.67835;
  bflux = fabs(*profit->paramlist[PARAM_SPHEROID_FLUX]);
  bt = bflux / (bflux+flux);
  if (bt > PATTERN_BTMAX)
    {
    wb = (bt - PATTERN_BTMAX) / (1.0 - PATTERN_BTMAX);
    if (wb > 1.0)
      wb = 1.0;
    omwb = 1.0 - wb;
    flux = wb*bflux + omwb*flux;
    scale = omwb*scale
	+ wb*fabs(*profit->paramlist[PARAM_SPHEROID_REFF]/profit->pixstep)*1.5;
    aspect = omwb*aspect
	+ wb*fabs(*profit->paramlist[PARAM_SPHEROID_ASPECT]);
    posangle /= DEG;
    dposangle = fmod_m90_p90(*profit->paramlist[PARAM_SPHEROID_POSANG])
		- posangle;
    if (dposangle > 90.0)
      dposangle -= 180.0;
    else if (dposangle < -90.0)
      dposangle += 180.0;
    posangle = fmod_m90_p90(posangle + wb*dposangle)*DEG;
    }

  ctheta = cos(posangle);
  stheta = sin(posangle);
  saspect = fabs(aspect);
  xscale = (scale==0.0)? 0.0 : 1.0/scale;
  yscale = (scale*saspect == 0.0)? 0.0 : 1.0/(scale*saspect);
  cd11 = xscale*ctheta;
  cd12 = xscale*stheta;
  cd21 = -yscale*stheta;
  cd22 = yscale*ctheta;
 
  x1cout = (double)(pattern->size[0]/2);
  x2cout = (double)(pattern->size[1]/2);

/* Determinant of the change of coordinate system */
  det = xscale*yscale;
  sbd = fabs(flux)*det/(2.0*PI);
  r2min = det/10.0;
/* Stay within an ellipse contained in the pattern raster, both in x and y */
  r2max = PATTERN_SCALE*PATTERN_SCALE;
  margin2 = (1.0-PATTERN_MARGIN)*(1.0-PATTERN_MARGIN);
  if (r2max > (dval = margin2
		* x1cout*x1cout * det*det / (cd12*cd12+cd22*cd22)))
    r2max = dval;
  if (r2max > (dval = margin2
		* x2cout*x2cout * det*det / (cd21*cd21+cd11*cd11)))
    r2max = dval;
/* Set the limit of the pattern extent */
//  rad = 4.0*profit->obj->a*xscale;
/* The pattern limit does not exceed 90% of the mapped ellipse "radius" */
//  if (rad*rad > 0.9*0.9*r2max)
  nrad = pattern->ncomp;
  pattern->rmax = rad = sqrt(r2max);
  if (!nrad)
      error(EXIT_FAILURE,
		"*Error*: insufficient number of vector elements",
		" for generating the pattern basis"); 
  dnrad = (double)nrad;
  npix = pattern->size[0]*pattern->size[1];
  normt = pattern->norm;
  switch(pattern->type)
    {
    case PATTERN_QUADRUPOLE:
    case PATTERN_OCTOPOLE:
      cpix = pattern->modpix;
      spix = pattern->modpix+npix;
      angcoeff = (pattern->type==PATTERN_OCTOPOLE)? 4.0 : 2.0;
      for (p=0; p<nrad; p++, cpix+=npix, spix+=npix)
        {
        x1 = -x1cout;
        x2 = -x2cout;       
        cnorm = snorm = cpnorm = spnorm = 0.0;
        rl = p*rad/dnrad;
        rl2 = rl*rl;
        pattern->r[p] = r0 = (p+1)*rad/dnrad;
        r02 = r0*r0;
        rh = (p+2)*rad/dnrad;
        rh2 = rh*rh;
        pmodpix = profit->modpix;
        for (ix2=pattern->size[1]; ix2--; x2+=1.0)
          {
          x1t = cd12*x2 + cd11*x1;
          x2t = cd22*x2 + cd21*x1;
          for (ix1=pattern->size[0]; ix1--; pmodpix++)
            {
            r2 = x1t*x1t+x2t*x2t;
            if (r2>rl2 && r2<rh2)
              {
              r = sqrt(r2);
              dval = (r<r0) ? (r-rl)/(r0-rl) : (rh-r)/(rh-r0);
              mod = (dval<0.5)? 2.0*dval*dval : 1.0-2.0*(1.0-dval)*(1.0-dval);
              ang = angcoeff*atan2(x2t,x1t);
#ifdef HAVE_SINCOS
              sincos(ang, &sinang, &cosang);
#else
              sinang = sin(ang);
              cosang = cos(ang);
#endif
              *(cpix++) = (float)(cmod = mod*cosang);
              *(spix++) = (float)(smod = mod*sinang);
              cnorm += cmod*cmod;
              snorm += smod*smod;
              cpnorm += cmod*cmod*(double)(*pmodpix**pmodpix);
              spnorm += smod*smod*(double)(*pmodpix**pmodpix);
              }
            else
              *(cpix++) = *(spix++) = 0.0;
            x1t += cd11;
            x2t += cd21;
            }
          }
        cpix -= npix;
        cnorm = (cnorm > 0.0? 1.0/sqrt(cnorm) : 1.0);
        *(normt++) = (float)cnorm*sqrt(cpnorm);
        for (i=npix; i--;)
          *(cpix++) *= cnorm;
        spix -= npix;
        snorm = (snorm > 0.0? 1.0/sqrt(snorm) : 1.0);
        *(normt++) = snorm*sqrt(spnorm);
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
        rl = p*rad/dnrad;
        rl2 = rl*rl;
        pattern->r[p] = r0 = (p+1)*rad/dnrad;
        r02 = r0*r0;
        rh = (p+2)*rad/dnrad;
        rh2 = rh*rh;
        for (f=0; f<=PATTERN_FMAX; f++)
          {
          norm = pnorm = 0.0;
          r2pix = r2buf;
          pmodpix = profit->modpix;
          if (!f)
            {
            for (i=npix; i--; pmodpix++)
              {
              r2 = *(r2pix++);
              if (r2>rl2 && r2<rh2)
                {
                r = sqrt(r2);
                dval = (r<r0) ? (r-rl)/(r0-rl) : (rh-r)/(rh-r0);
                *(pix++) = (float)(dval = (dval<0.5)?
			2.0*dval*dval : 1.0-2.0*(1.0-dval)*(1.0-dval));
                norm += dval*dval;
                }
              else
                *(pix++) = 0.0;
              }
            pix -= npix;
            pnorm = norm*sbd*sbd;
            norm0 = norm = (norm > 1.0/BIG? 1.0/sqrt(norm) : 1.0);
            *(normt++) = pnorm > 1.0/BIG? 1.0/sqrt(pnorm) : 0.0;
            fnorm = (float)norm;
            for (i=npix; i--;)
              *(pix++) *= fnorm;
            modpix = pix;
            }
          else
            {
            modpix -= npix;
            scpixt = scbuf[f-1];
            for (i=npix; i--; pmodpix++)
              {
              *(pix++) = (float)(dval = *(modpix++)**(scpixt++));
              norm += dval*dval;
              pnorm += dval*dval**pmodpix**pmodpix;
              }
            pix -= npix;
            pnorm = norm*sbd*sbd;
            norm = (norm > 0.0? 1.0/sqrt(norm) : 1.0);
            *(normt++) = (float)(pnorm > 1.0/BIG? norm0/sqrt(pnorm) : 0.0);
            fnorm = (float)norm;
            for (i=npix; i--;)
              *(pix++) *= fnorm;
            modpix -= npix;
            norm = pnorm = 0.0;
            pmodpix = profit->modpix;
            for (i=npix; i--; pmodpix++)
              {
              *(pix++) = (float)(dval = *(modpix++)**(scpixt++));
              norm += dval*dval;
              pnorm += dval*dval**pmodpix**pmodpix;
              }
            pix -= npix;
            pnorm = norm*sbd*sbd;
            norm = (norm > 0.0? 1.0/sqrt(norm) : 1.0);
            *(normt++) = pnorm > 1.0/BIG? norm0/sqrt(pnorm) : 0.0;
            fnorm = (float)norm;
            for (i=npix; i--;)
              *(pix++) *= fnorm;
            }
          }
        }
      free(r2buf);
      for (f=0; f<PATTERN_FMAX; f++)
        free(scbuf[f]);
      break;

    case PATTERN_POLARSHAPELETS:
      nmax = pattern->ncomp;
      kmax = (nmax+1)*(nmax+2)/2;
      beta = 0.667;

      invbeta2 = 1.0/(beta*beta);

/*---- Precompute some slow functions */
      QMALLOC(fr2, double, npix);
      QMALLOC(fexpr2, double, npix);
      QMALLOC(ftheta, double, npix);
      fr2t = fr2;
      fexpr2t = fexpr2;
      fthetat = ftheta;
      x1 = -x1cout;
      x2 = -x2cout;       
      for (ix2=pattern->size[1]; ix2--; x2+=1.0)
        {
        x1t = cd12*x2 + cd11*x1;
        x2t = cd22*x2 + cd21*x1;
        for (ix1=pattern->size[0]; ix1--;)
          {      
          *(fr2t++) = r2 = (x1t*x1t+x2t*x2t)*invbeta2;
          *(fexpr2t++) = exp(-r2/2.0);
          *(fthetat++) = atan2(x2t,x1t);
          x1t += cd11;
          x2t += cd21;
          }
        }

      pix = pattern->modpix;
      for (n=0; n<=nmax; n++)
        {
        for (m=n%2; m<=n; m+=2)
          {
          dm = (double)m;
/*-------- Compute ((n+m)/2)!/((n-m)/2)! */
          hnmm = (n-m)/2;
          fac = 1.0;
//          for (p=(n+m)/2; p>=hnmm; p--)
//            if (p)
//              fac *= (double)p;
//          fac = sqrt(1.0/(PI*fac))/beta;
          if ((hnmm%2))
            fac = -fac;
          fr2t = fr2;
          fexpr2t = fexpr2;
          fthetat = ftheta;
          norm = 0.0;
          for (i=npix; i--;fr2t++)
            {
            *(pix++) = (float)(dval = fac*pow(*fr2t, dm/2.0)
			*psf_laguerre(*fr2t, hnmm, m)
			**(fexpr2t++)*cos(dm**(fthetat++)));
            norm += dval*dval;
            }
          pix -= npix;
          pnorm = norm*sbd*sbd;
          norm = (norm > 0.0? 1.0/sqrt(norm) : 1.0);
          *(normt++) = pnorm > 1.0/BIG? norm/sqrt(pnorm) : 0.0;
          fnorm = (float)norm;
          for (i=npix; i--;)
            *(pix++) *= fnorm;
          if (m!=0)
            {
            fr2t = fr2;
            fexpr2t = fexpr2;
            fthetat = ftheta;
            norm = 0.0;
            for (i=npix; i--; fr2t++)
              {
              *(pix++) = (float)(dval = fac*pow(*fr2t, dm/2.0)
			*psf_laguerre(*fr2t, hnmm, m)
			**(fexpr2t++)*sin(dm**(fthetat++)));
              norm += dval*dval;
              }
            pix -= npix;
            pnorm = norm*sbd*sbd;
            norm = (norm > 0.0? 1.0/sqrt(norm) : 1.0);
            *(normt++) = pnorm > 1.0/BIG? norm/sqrt(pnorm) : 0.0;
            fnorm = (float)norm;
            for (i=npix; i--;)
              *(pix++) *= fnorm;
            }
          }
        }

      free(fr2);
      free(fexpr2);
      free(ftheta);
      break;
    default:
      error(EXIT_FAILURE, "*Internal Error*: Unknown Pattern type","");
    }

  return;
  }


/****** psf_laguerre **********************************************************
PROTO	double	psf_laguerre(double x, int p, int q)
PURPOSE	Return Laguerre polynomial value.
INPUT	x,
	p,
	q.
OUTPUT  Value of the Laguerre polynomial.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 12/11/2007
 ***/
static double	psf_laguerre(double x, int p, int q)
  {
   double	dn,dq, lpm1,lpm2, l;
   int		n;

  dq = q - 1.0;
  if (p==0)
    return 1.0;
  else if (p==1)
    return (2.0 - x + dq);
  else
    {
    l = 0.0;
    lpm2 = 1.0;
    lpm1 = 2.0 - x + dq;
    dn = 2.0;
    for (n=p-1; n--; dn+=1.0)
      {
      l = (2.0+(dq-x)/dn)*lpm1 - (1.0+dq/dn)*lpm2;
      lpm2 = lpm1;
      lpm1 = l;
      }
    }

  return l;
  }



