/*
*				photom.c
*
* Compute magnitudes and other photometric parameters.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
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
#include	<stdlib.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"photom.h"
#include	"plist.h"
#include	"subimage.h"

/****** photom_aper **********************************************************
PROTO	void photom_aper(fieldstruct *field, fieldstruct *wfield,
			obj2struct *obj2, int aper)
PURPOSE	Measure the flux within a circular aperture.
INPUT	Pointer to the image structure,
	pointer to the weight-map structure,
	pointer to the obj2 structure,
	aperture number.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/09/2014
 ***/
void  photom_aper(fieldstruct *field, fieldstruct *wfield, obj2struct *obj2,
		int aper)

  {
   subimagestruct	*subimage;
   float		r2, raper,raper2, rintlim,rintlim2,rextlim2,
			mx,my,dx,dx1,dy,dy2,
			offsetx,offsety,scalex,scaley,scale2, ngamma, locarea;
   double		tv, sigtv, area, pix, var, backnoise2, gain;
   int			x,y, x2,y2, xmin,xmax,ymin,ymax, sx,sy, w,h,
			pflag,corrflag, gainflag;
   long			pos;
   PIXTYPE		*image,*imaget, *weight,*weightt,
			wthresh = 0.0;

  if (wfield)
    wthresh = wfield->weight_thresh;
  weight = weightt = NULL;
  subimage = obj2->subimage;
  mx = subimage->dpos[0] - subimage->xmin[0];
  my = subimage->dpos[1] - subimage->xmin[1];
  w = subimage->size[0];
  h = subimage->size[1];
  ngamma = field->ngamma;
  pflag = (field->detector_type==DETECTOR_PHOTO);
  corrflag = (prefs.mask_type==MASK_CORRECT);
  gainflag = wfield && prefs.weightgain_flag;
  var = backnoise2 = field->backsig*field->backsig;
  gain = field->gain;
/* Integration radius */
  raper = prefs.apert[aper]/2.0;
  raper2 = raper*raper;
/* Internal radius of the oversampled annulus (<r-sqrt(2)/2) */
  rintlim = raper - 0.75;
  rintlim2 = (rintlim>0.0)? rintlim*rintlim: 0.0;
/* External radius of the oversampled annulus (>r+sqrt(2)/2) */
  rextlim2 = (raper + 0.75)*(raper + 0.75);
  tv = sigtv = area = 0.0;
  scaley = scalex = 1.0/APER_OVERSAMP;
  scale2 = scalex*scaley;
  offsetx = 0.5*(scalex-1.0);
  offsety = 0.5*(scaley-1.0);

  xmin = (int)(mx-raper+0.499999);
  xmax = (int)(mx+raper+1.499999);
  ymin = (int)(my-raper+0.499999);
  ymax = (int)(my+raper+1.499999);

  if (xmin < 0)
    {
    xmin = 0;
    obj2->flags |= OBJ_APERT_PB;
    }
  if (xmax > w)
    {
    xmax = w;
    obj2->flags |= OBJ_APERT_PB;
    }
  if (ymin < 0)
    {
    ymin = 0;
    obj2->flags |= OBJ_APERT_PB;
    }
  if (ymax > h)
    {
    ymax = h;
    obj2->flags |= OBJ_APERT_PB;
    }

  image = subimage->image;
  weight = subimage->imvar;
  for (y=ymin; y<ymax; y++)
    {
    imaget = image + (pos = y*w + xmin);
    if (wfield)
      weightt = weight + pos;
    for (x=xmin; x<xmax; x++, imaget++, weightt++)
      {
      dx = x - mx;
      dy = y - my;
      if ((r2=dx*dx+dy*dy) < rextlim2)
        {
        if (r2> rintlim2)
          {
          dx += offsetx;
          dy += offsety;
          locarea = 0.0;
          for (sy=APER_OVERSAMP; sy--; dy+=scaley)
            {
            dx1 = dx;
            dy2 = dy*dy;
            for (sx=APER_OVERSAMP; sx--; dx1+=scalex)
              if (dx1*dx1+dy2<raper2)
                locarea += scale2;
            }
          }
        else
          locarea = 1.0;
        area += locarea;
/*------ Here begin tests for pixel and/or weight overflows. Things are a */
/*------ bit intricated to have it running as fast as possible in the most */
/*------ common cases */
        if ((pix=*imaget)<=-BIG || (wfield && (var=*weightt)>=wthresh))
          {
          if (corrflag
		&& (x2=(int)(2*mx+0.49999-x))>=0 && x2<w
		&& (y2=(int)(2*my+0.49999-y))>=0 && y2<h
		&& (pix=*(image + (pos = y2*w + x2)))>-BIG)
            {
            if (wfield)
              {
              var = *(weight + pos);
              if (var>=wthresh)
                pix = var = 0.0;
              }
            }
          else
            {
            pix = 0.0;
            if (wfield)
              var = 0.0;
            }
          }
        if (pflag)
          {
          pix=exp(pix/ngamma);
          sigtv += var*locarea*pix*pix;
          }
        else
          sigtv += var*locarea;
        tv += locarea*pix;
        if (gainflag && pix>0.0 && gain>0.0)
          sigtv += pix/gain*var/backnoise2;
        }
      }
    }

  if (pflag)
    {
    tv = ngamma*(tv-area*exp(obj2->dbkg[0]/ngamma));
    sigtv /= ngamma*ngamma;
    }
  else
    {
    tv -= area*obj2->dbkg[0];
    if (!gainflag && gain > 0.0 && tv>0.0)
      sigtv += tv/gain;
    }

  if (aper<prefs.aper_size[0])
    obj2->flux_aper[aper] = tv;
  if (aper<prefs.aper_size[0])
    obj2->fluxerr_aper[aper] = sqrt(sigtv);
  if (aper<prefs.aper_size[0])
    obj2->mag_aper[aper] = tv>0.0? -2.5*log10(tv) + field->mag_zeropoint : 99.0;
  if (aper<prefs.aper_size[0])
    obj2->magerr_aper[aper] = tv>0.0? 1.086*sqrt(sigtv)/tv:99.0;

  return;
  }


/****** photom_petro *********************************************************
PROTO	void photom_petro(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2struct *obj2)
PURPOSE	Measure the flux within a Petrosian elliptical aperture
INPUT	Pointer to the image structure,
	pointer to the detection image structure,
	pointer to the weight-map structure,
	pointer to the detection weight-map structure,
	pointer to the obj2 structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/09/2014
 ***/
void  photom_petro(fieldstruct *field, fieldstruct *dfield,
	fieldstruct *wfield, fieldstruct *dwfield, obj2struct *obj2)
  {
   subimagestruct	*subimage;
   double		sigtv, tv, r1, v1,var,gain,backnoise2, muden,munum;
   float		bkg, ngamma, mx,my, dx,dy, cx2,cy2,cxy, r2,
			klim, klim2,kmin,kmin2,kmax,kmax2,kstep,kmea,kmea2,
			dxlim, dylim;
   int			area,areab, areaden, areanum,
			x,y, x2,y2, xmin,xmax,ymin,ymax, w,h,
			pflag, corrflag, gainflag, pos;
   PIXTYPE		*image,*imaget, *dimage,*dimaget, *weight,*weightt,
			*dweight,*dweightt,
			pix, wthresh=0.0, dwthresh=0.0;


/* Let's initialize some variables */
  if (!dfield)
    dfield = field;
  if (dwfield)
    dwthresh = dwfield->weight_thresh;
  weight = dweight = NULL;
  if (wfield)
    wthresh = wfield->weight_thresh;
  weightt = dweightt = NULL;
  subimage = obj2->subimage;
  w = subimage->size[0];
  h = subimage->size[1];
  ngamma = field->ngamma;
  bkg = (double)obj2->dbkg[0];
  mx = subimage->dpos[0] - (double)subimage->xmin[0];
  my = subimage->dpos[1] - (double)subimage->xmin[1];
  var = backnoise2 = field->backsig*field->backsig;
  gain = field->gain;
  pflag = (field->detector_type==DETECTOR_PHOTO);
  corrflag = (prefs.mask_type==MASK_CORRECT);
  gainflag = wfield && prefs.weightgain_flag;

/* First step: find the extent of the ellipse (the Petrosian factor) */
/* Clip boundaries in x and y */
/* We first check that the search ellipse is large enough... */
  if (PETRO_NSIG*sqrt(obj2->a*obj2->b)>prefs.autoaper[0]/2.0)
    {
    cx2 = obj2->cxx;
    cy2 = obj2->cyy;
    cxy = obj2->cxy;
    dxlim = cx2 - cxy*cxy/(4.0*cy2);
    dxlim = dxlim>0.0 ? PETRO_NSIG/sqrt(dxlim) : 0.0;
    dylim = cy2 - cxy*cxy/(4.0*cx2);
    dylim = dylim > 0.0 ? PETRO_NSIG/sqrt(dylim) : 0.0;
    klim2 = PETRO_NSIG*PETRO_NSIG;
    }
  else
/*-- ...if not, use the circular aperture provided by the user */
    {
    cx2 = cy2 = 1.0;
    cxy = 0.0;
    dxlim = dylim = prefs.autoaper[0]/2.0;
    klim2 =  dxlim*dxlim;
    }

  if ((xmin = RINT(mx-dxlim)) < 0)
    {
    xmin = 0;
    obj2->flags |= OBJ_APERT_PB;
    }
  if ((xmax = RINT(mx+dxlim)+1) > w)
    {
    xmax = w;
    obj2->flags |= OBJ_APERT_PB;
    }
  if ((ymin = RINT(my-dylim)) < 0)
    {
    ymin = 0;
    obj2->flags |= OBJ_APERT_PB;
    }
  if ((ymax = RINT(my+dylim)+1) > h)
    {
    ymax = h;
    obj2->flags |= OBJ_APERT_PB;
    }

  dimage = subimage->image;
  if (dwfield)
    dweight = subimage->imvar;
  klim = sqrt(klim2);
  kstep = klim/20.0;
  area = areab = areanum = areaden = 0;
  munum = muden = 0.0;
  kmea = 0.0;
  for (kmin=kstep; (kmax=kmin*1.2)<klim; kmin += kstep)
    {
    kmea = (kmin+kmax)/2.0;
    kmea2 = kmea*kmea;
    kmin2 = kmin*kmin;
    kmax2 = kmax*kmax;
    v1 = r1 = 0.0;
    area = areab = areanum = areaden = 0;
    munum = muden = 0.0;
    for (y=ymin; y<ymax; y++)
      {
      dimaget = dimage + (pos = y*w + xmin);
      if (dwfield)
        dweightt = dweight + pos;
      for (x=xmin; x<xmax; x++, dimaget++, dweightt++)
      {
      dx = x - mx;
      dy = y - my;
      if ((r2=cx2*dx*dx + cy2*dy*dy + cxy*dx*dy) <= kmax2)
        {
        if ((pix=*dimaget)>-BIG && (!dwfield || (dwfield&&*dweightt<dwthresh)))
          {
          area++;
          if (r2>=kmin2)
	    {
            munum += pix;
            areanum++;
            }
          if (r2<kmea2)
	    {
            muden += pix;
            areaden++;
            }
          }
        else
          areab++;
        }
      }
    }
  if (areanum && areaden)
    {
    munum /= (double)areanum;
    muden /= (double)areaden;
    if (munum<muden*0.2)
      break;
    }
  }

  area += areab;
  if (area)
    {
/*-- Go further only if some pixels are available !! */
    if (areanum && areaden && munum && muden)
      {
      obj2->petrofactor = prefs.petroparam[0]*kmea;
      if (obj2->petrofactor < prefs.petroparam[1])
        obj2->petrofactor = prefs.petroparam[1];
      }
    else
      obj2->petrofactor = prefs.petroparam[1];

/*-- Flag if the Petrosian photometry can be strongly affected by neighhours */
    if ((float)areab/area > CROWD_THRESHOLD)
      obj2->flags |= OBJ_CROWDED;

/*-- Second step: integrate within the ellipse */
/*-- Clip boundaries in x and y (bis) */
/*-- We first check that the derived ellipse is large enough... */
    if (obj2->petrofactor*sqrt(obj2->a*obj2->b)>prefs.autoaper[1]/2.0)
      {
      cx2 = obj2->cxx;
      cy2 = obj2->cyy;
      cxy = obj2->cxy;
      dxlim = cx2 - cxy*cxy/(4.0*cy2);
      dxlim = dxlim>0.0 ? obj2->petrofactor/sqrt(dxlim) : 0.0;
      dylim = cy2 - cxy*cxy/(4.0*cx2);
      dylim = dylim > 0.0 ? obj2->petrofactor/sqrt(dylim) : 0.0;
      klim2 = obj2->petrofactor*obj2->petrofactor;
      }
    else
/*---- ...if not, use the circular aperture provided by the user */
      {
      cx2 = cy2 = 1.0;
      cxy = 0.0;
      dxlim = dylim = prefs.autoaper[1]/2.0;
      klim2 =  dxlim*dxlim;
      obj2->petrofactor = 0.0;
      }

    if ((xmin = RINT(mx-dxlim)) < 0)
      {
      xmin = 0;
      obj2->flags |= OBJ_APERT_PB;
      }
    if ((xmax = RINT(mx+dxlim)+1) > w)
      {
      xmax = w;
      obj2->flags |= OBJ_APERT_PB;
      }
    if ((ymin = RINT(my-dylim)) < 0)
      {
      ymin = 0;
      obj2->flags |= OBJ_APERT_PB;
      }
    if ((ymax = RINT(my+dylim)+1) > w)
      {
      ymax = w;
      obj2->flags |= OBJ_APERT_PB;
      }

    area = areab = 0;
    tv = sigtv = 0.0;
    image = subimage->image;
    if (wfield)
      weight = subimage->imvar;
    for (y=ymin; y<ymax; y++)
      {
      imaget = image + (pos = y*w + xmin);
      if (wfield)
        weightt = weight + pos;
      for (x=xmin; x<xmax; x++, imaget++, weightt++)
        {
        dx = x - mx;
        dy = y - my;
        if ((cx2*dx*dx + cy2*dy*dy + cxy*dx*dy) <= klim2)
          {
          area++;
/*-------- Here begin tests for pixel and/or weight overflows. Things are a */
/*-------- bit intricated to have it running as fast as possible in the most */
/*-------- common cases */
          if ((pix=*imaget)<=-BIG || (wfield && (var=*weightt)>=wthresh))
            {
            areab++;
            if (corrflag
		&& (x2=(int)(2*mx+0.49999-x))>=0 && x2<w
		&& (y2=(int)(2*my+0.49999-y))>=0 && y2<h
		&& (pix=*(image + (pos = y2*w + x2)))>-BIG)
              {
              if (wfield)
                {
                var = *(weight + pos);
                if (var>=wthresh)
                  pix = var = 0.0;
                }
              }
            else
              {
              pix = 0.0;
              if (wfield)
                var = 0.0;
              }
            }
          if (pflag)
            {
            pix = exp(pix/ngamma);
            sigtv += var*pix*pix;
            }
          else
            sigtv += var;
          tv += pix;
          if (gainflag && pix>0.0 && gain>0.0)
            sigtv += pix/gain*var/backnoise2;
          }
        }
      }

/*-- Flag if the Petrosian photometry can be strongly affected by neighhours */
    if ((float)areab > CROWD_THRESHOLD*area)
      obj2->flags |= OBJ_CROWDED;

    if (pflag)
      {
      tv = ngamma*(tv-area*exp(bkg/ngamma));
      sigtv /= ngamma*ngamma;
      }
    else
      {
      tv -= area*bkg;
      if (!gainflag && gain > 0.0 && tv>0.0)
        sigtv += tv/gain;
      }
    }
  else
/*-- No available pixels: set the flux to zero */
    tv = sigtv = 0.0;


  obj2->flux_petro = tv;
  obj2->fluxerr_petro = sqrt(sigtv);

  if (FLAG(obj2.mag_petro))
    obj2->mag_petro = obj2->flux_petro>0.0?
			 -2.5*log10(obj2->flux_petro) + field->mag_zeropoint
			:99.0;
  if (FLAG(obj2.magerr_petro))
    obj2->magerr_petro = obj2->flux_petro>0.0?
			 1.086*obj2->fluxerr_petro/obj2->flux_petro
			:99.0;
  if (tv<=0.0)
    obj2->petrofactor = 0.0;

  return;
  }


/****** photom_auto *********************************************************
PROTO	void photom_auto(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2struct *obj2)
PURPOSE	Measure the flux within a Kron elliptical aperture
INPUT	Pointer to an array of image field pointers,
	pointer to an array of weight-map field pointers,
	number of images,
	pointer to the obj2 structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/05/2012
 ***/
void  photom_auto(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2struct *obj2)

  {
   fieldstruct		*field, *wfield;
   subimagestruct	*subimage;
   double		*jac,
			sigtv, tv, wv, r1, v1,var,gain,backnoise2;
   float		ngamma, mx,my, dx,dy, cx2,cy2,cxy, r2,klim2,
			dxlim, dylim, ftv,fwv;
   PIXTYPE		*image,*imaget, *weight,*weightt,
			pix, wthresh=0.0;
   int			f,s, area,areab, x,y, x2,y2, xmin,xmax,ymin,ymax,
			fymin,fymax, w,h, band, nband, nsubimage, pos,
			pflag, corrflag, gainflag, autoflag, cfluxflag;

/* field is the detection field */
  field = fields[0];
  wfield = (wfields && wfields[0])? wfields[0] : NULL;

/* Let's initialize some variables */
  weight = NULL;
  if (wfield)
    wthresh = wfield->weight_thresh;
  weightt = NULL;
  subimage = obj2->subimage;
  w = subimage->size[0];
  h = subimage->size[1];
  ngamma = field->ngamma;
  mx = subimage->dpos[0] - (double)subimage->xmin[0];
  my = subimage->dpos[1] - (double)subimage->xmin[1];
  var = backnoise2 = field->backsig*field->backsig;
  gain = field->gain;
  corrflag = (prefs.mask_type==MASK_CORRECT);
  gainflag = wfield && prefs.weightgain_flag;
  cfluxflag = FLAG(obj2.cflux_auto) | FLAG(obj2.cfluxerr_auto)
	| FLAG(obj2.cmag_auto) | FLAG(obj2.cmagerr_auto);

/* First step: find the extent of the ellipse (the kron factor r1) */
/* Clip boundaries in x and y */
/* We first check that the search ellipse is large enough... */
  if (KRON_NSIG*sqrt(obj2->a*obj2->b)>prefs.autoaper[0]/2.0)
    {
    cx2 = obj2->cxx;
    cy2 = obj2->cyy;
    cxy = obj2->cxy;
    dxlim = cx2 - cxy*cxy/(4.0*cy2);
    dxlim = dxlim>0.0 ? KRON_NSIG/sqrt(dxlim) : 0.0;
    dylim = cy2 - cxy*cxy/(4.0*cx2);
    dylim = dylim > 0.0 ? KRON_NSIG/sqrt(dylim) : 0.0;
    klim2 = KRON_NSIG*KRON_NSIG;
    }
  else
/*-- ...if not, use the circular aperture provided by the user */
    {
    cx2 = cy2 = 1.0;
    cxy = 0.0;
    dxlim = dylim = prefs.autoaper[0]/2.0;
    klim2 =  dxlim*dxlim;
    }

  if ((xmin = RINT(mx-dxlim)) < 0)
    {
    xmin = 0;
    obj2->flags |= OBJ_APERT_PB;
    }
  if ((xmax = RINT(mx+dxlim)+1) > w)
    {
    xmax = w;
    obj2->flags |= OBJ_APERT_PB;
    }
  if ((ymin = RINT(my-dylim)) < 0)
    {
    ymin = 0;
    obj2->flags |= OBJ_APERT_PB;
    }
  if ((ymax = RINT(my+dylim)+1) > h)
    {
    ymax = h;
    obj2->flags |= OBJ_APERT_PB;
    }

  v1 = r1 = 0.0;
  area = areab = 0;
  image = subimage->image;
  if (wfield)
    weight = subimage->imvar;
  for (y=ymin; y<ymax; y++)
    {
    imaget = image + (pos = y*w + xmin);
    if (wfield)
      weightt = weight + pos;
    for (x=xmin; x<xmax; x++, imaget++, weightt++)
      {
      dx = x - mx;
      dy = y - my;
      if ((r2=cx2*dx*dx + cy2*dy*dy + cxy*dx*dy) <= klim2)
        {
        if ((pix=*imaget)>-BIG && (!wfield || (wfield&&*weightt<wthresh)))
          {
          area++;
          r1 += sqrt(r2)*pix;
          v1 += pix;
          }
        else
          areab++;
        }
      }
    }

  nband = prefs.nphotinstru;
  for (band=0; band<nband; band++)
    obj2->cflux[band] = obj2->cfluxw[band] = 0.0;

  area += areab;
  if (area)
    {
/*-- Go further only if some pixels are available !! */
    if (r1>0.0 && v1>0.0)
      {
      obj2->auto_kronfactor = prefs.autoparam[0]*r1/v1;
      if (obj2->auto_kronfactor < prefs.autoparam[1])
        obj2->auto_kronfactor = prefs.autoparam[1];
      }
    else
      obj2->auto_kronfactor = prefs.autoparam[1];

/*-- Flag if the Kron photometry can be strongly affected by neighhours */
    if ((float)areab/area > CROWD_THRESHOLD)
      {
      obj2->flags |= OBJ_CROWDED;
      if (FLAG(obj2.auto_flags))
        for (f=0; f<nfield; f++)
          obj2->auto_flags[f] |= OBJ_CROWDED;
      }

/*-- Second step: integrate within the ellipse */
/*-- Clip boundaries in x and y (bis) */
/*-- We first check that the derived ellipse is large enough... */
    if (obj2->auto_kronfactor*sqrt(obj2->a*obj2->b)>prefs.autoaper[1]/2.0)
      {
      cx2 = obj2->cxx;
      cy2 = obj2->cyy;
      cxy = obj2->cxy;
      dxlim = cx2 - cxy*cxy/(4.0*cy2);
      dxlim = dxlim>0.0 ? obj2->auto_kronfactor/sqrt(dxlim) : 0.0;
      dylim = cy2 - cxy*cxy/(4.0*cx2);
      dylim = dylim > 0.0 ? obj2->auto_kronfactor/sqrt(dylim) : 0.0;
      klim2 = obj2->auto_kronfactor*obj2->auto_kronfactor;
      }
    else
/*---- ...if not, use the circular aperture provided by the user */
      {
      cx2 = cy2 = 1.0;
      cxy = 0.0;
      dxlim = dylim = prefs.autoaper[1]/2.0;
      klim2 =  dxlim*dxlim;
      obj2->auto_kronfactor = 0.0;
      }

    nsubimage = obj2->nsubimage;
    subimage = obj2->subimage;
    for (s=0; s<nsubimage; s++, subimage++)
      {
      field = subimage->field;
      wfield = subimage->wfield;
      f = field->imindex;
      pflag = (field->detector_type==DETECTOR_PHOTO);

      w = subimage->size[0];
      h = subimage->size[1];
      mx = subimage->dpos[0] - (double)subimage->xmin[0];
      my = subimage->dpos[1] - (double)subimage->xmin[1];

      jac = subimage->dinvjacob;
      cx2 = jac[0]*jac[0]*obj2->cxx
	+ jac[2]*jac[2]*obj2->cyy
	+ jac[0]*jac[2]*obj2->cxy;
      cy2 = jac[1]*jac[1]*obj2->cxx
	+ jac[3]*jac[3]*obj2->cyy
	+ jac[1]*jac[3]*obj2->cxy;
      cxy = 2.0*jac[0]*jac[1]*obj2->cxx
	+ 2.0*jac[2]*jac[3]*obj2->cyy
	+ (jac[0]*jac[3]+jac[1]*jac[2])*obj2->cxy;
      dxlim = cx2 - cxy*cxy/(4.0*cy2);
      dxlim = dxlim>0.0 ? obj2->auto_kronfactor/sqrt(dxlim) : 0.0;
      dylim = cy2 - cxy*cxy/(4.0*cx2);
      dylim = dylim > 0.0 ? obj2->auto_kronfactor/sqrt(dylim) : 0.0;

      autoflag = 0;
      if ((xmin = RINT(mx-dxlim)) < 0)
        {
        xmin = 0;
        autoflag = AUTOFLAG_APERT_PB;
        }
      if ((xmax = RINT(mx+dxlim)+1) > w)
        {
        xmax = w;
        autoflag = AUTOFLAG_APERT_PB;
        }
      if ((ymin = RINT(my-dylim)) < 0)
        {
        ymin = 0;
        autoflag = AUTOFLAG_APERT_PB;
        }
      if ((ymax = RINT(my+dylim)+1) > h)
        {
        ymax = h;
        autoflag = AUTOFLAG_APERT_PB;
        }

      area = areab = 0;
      tv = sigtv = 0.0;
      image = subimage->image;
      weight = subimage->imvar;
      for (y=ymin; y<ymax; y++)
        {
        imaget = image + (pos = xmin + (y%h)*w);
        if (wfield)
          weightt = weight + pos;
        for (x=xmin; x<xmax; x++, imaget++, weightt++)
          {
          dx = x - mx;
          dy = y - my;
          if ((cx2*dx*dx + cy2*dy*dy + cxy*dx*dy) <= klim2)
            {
            area++;
/*---------- Here begin tests for pixel and/or weight overflows. Things are a */
/*---------- bit intricated to have it running as fast as possible in the most */
/*---------- common cases */
            if ((pix=*imaget)<=-BIG || (wfield && (var=*weightt)>=wthresh))
              {
              areab++;
              if (corrflag
		&& (x2=(int)(2*mx+0.49999-x))>=0 && x2<w
		&& (y2=(int)(2*my+0.49999-y))>=0 && y2<h
		&& (pix=*(image + (pos = y2*w + x2)))>-BIG)
                {
                if (wfield)
                  {
                  var = *(weight + pos);
                  if (var>=wthresh)
                    pix = var = 0.0;
                  }
                }
              else
                {
                pix = 0.0;
                if (wfield)
                  var = 0.0;
                }
              }
            if (pflag)
              {
              pix = exp(pix/ngamma);
              sigtv += var*pix*pix;
              }
            else
              sigtv += var;
            tv += pix;
            if (gainflag && pix>0.0 && gain>0.0)
              sigtv += pix/gain*var/backnoise2;
            }
          }
        }

/*---- Flag if the Kron photometry can be strongly affected by neighhours */
      if ((float)areab > CROWD_THRESHOLD*area)
        autoflag |= AUTOFLAG_CROWDED;

      if (pflag)
        {
        tv = ngamma*(tv-area*exp(obj2->dbkg[f]/ngamma));
        sigtv /= ngamma*ngamma;
        }
      else
        {
        tv -= area*obj2->dbkg[f];
        if (!gainflag && gain > 0.0 && tv>0.0)
          sigtv += tv/gain;
        }

      if (FLAG(obj2.auto_flags))
        obj2->auto_flags[f] |= autoflag;
      if (FLAG(obj2.flux_auto))
        obj2->flux_auto[f] = tv;
      if (FLAG(obj2.mag_auto))
        obj2->mag_auto[f] = tv>0.0?
			-FDMAG * log(tv) + prefs.mag_zeropoint[f] : 99.0f;
      if (FLAG(obj2.fluxerr_auto))
        obj2->fluxerr_auto[f] = sqrt(sigtv);
      if (FLAG(obj2.magerr_auto))
        obj2->magerr_auto[f] = tv>0.0? FDMAG * sqrt(sigtv) / tv : 99.0f;

      if (FLAG(obj2.cflux_auto))
        {
        tv *= field->flux_factor;
        wv /= field->flux_factor*field->flux_factor;
        band = fields[f]->photomlabel;
        wv = sigtv>0.0? 1.0/sigtv : 0.0;
        obj2->cflux[band] += wv*tv;
        obj2->cfluxw[band] += wv;
        }
      }
    }
  else
/*---- No available pixels */
    {
    for (f=0; f<nfield; f++)
      {
      if (FLAG(obj2.auto_flags))
        obj2->auto_flags[f] |= AUTOFLAG_DETECT_PB;
      if (FLAG(obj2.flux_auto))
        obj2->flux_auto[f] = 0.0f;
      if (FLAG(obj2.mag_auto))
        obj2->mag_auto[f] = 99.0f;
      if (FLAG(obj2.fluxerr_auto))
        obj2->fluxerr_auto[f] = 0.0f;
      if (FLAG(obj2.magerr_auto))
        obj2->magerr_auto[f] = 99.0f;
      }
    if ((cfluxflag))
      for (band=0; band<nband; band++)
        {
        if (FLAG(obj2.cflux_auto))
          obj2->cflux_auto[band] = 0.0f;
        if (FLAG(obj2.cmag_auto))
          obj2->cmag_auto[band] = 99.0f;
        if (FLAG(obj2.cfluxerr_auto))
          obj2->cfluxerr_auto[band] = 0.0f;
        if (FLAG(obj2.cmagerr_auto))
          obj2->cmagerr_auto[band] = 99.0f;
        }
    }

/* Combined fluxes */

  if ((cfluxflag))
    for (band=0; band<nband; band++)
      {
      fwv = (float)obj2->cfluxw[band];
      ftv = fwv>0.0f? (float)(obj2->cflux[band]/obj2->cfluxw[band]) : 0.0f;
      if (FLAG(obj2.cflux_auto))
        obj2->cflux_auto[band] = ftv;
      if (FLAG(obj2.cmag_auto))
        obj2->cmag_auto[band] = ftv>0.0f ? -FDMAG * logf(ftv) : 99.0f;
      if (FLAG(obj2.cfluxerr_auto))
        obj2->cfluxerr_auto[band] = fwv>0.0? sqrtf(1.0/fwv) : 0.0f;
      if (FLAG(obj2.cmagerr_auto))
        obj2->cmagerr_auto[band] = obj2->cflux[band]>0.0?
			FDMAG * sqrtf(fwv)/(float)obj2->cflux[band] : 99.0f;
      }

  return;
  }


/****** photom_isocor ********************************************************
PROTO	void photom_isocor(fieldstruct *field, obj2struct *obj2)
PURPOSE	Correct isophotal flux as in Maddox et al. 1990.
INPUT	Pointer to the image structure,
	pointer to the obj2 structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	15/02/2012
 ***/
void  photom_isocor(fieldstruct *field, obj2struct *obj2)
  {
   double	ati;

  ati = (obj2->flux_iso>0.0)? (obj2->fdnpix*obj2->dthresh/obj2->flux_iso) : 0.0;
  if (ati>1.0)
    ati = 1.0;
  else if (ati<0.0)
    ati = 0.0;
  obj2->flux_isocor = obj2->flux_iso/(1.0-0.196099*ati-0.751208*ati*ati);
  if (FLAG(obj2.fluxerr_isocor))
    {
    if (obj2->flux_iso>0.0)
      {
       double	dati, sigtv;

      sigtv = obj2->fluxerr_iso/(obj2->flux_iso*obj2->flux_iso);
      dati = obj2->fdnpix?ati*sqrt(sigtv+1.0/obj2->fdnpix): 0.0;
      dati = 0.196099*dati + 0.751208*2*ati*dati;
      obj2->fluxerr_isocor = sqrt(sigtv+dati*dati)*obj2->flux_iso;
      }
    else
      obj2->fluxerr_isocor = sqrt(obj2->fluxerr_iso);
    }

  return;
  }


/****** photom_mags *********************************************************
PROTO	void photom_mags(fieldstruct *field, obj2struct *obj2)
PURPOSE	Compute magnitudes from flux measurements.
INPUT	Pointer to the image structure,
	pointer to the obj2 structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/05/2012
 ***/
void	photom_mags(fieldstruct *field, obj2struct *obj2)
  {
/* Mag. isophotal */
  if (FLAG(obj2.mag_iso))
    obj2->mag_iso = obj2->flux_iso>0.0?
			 -2.5*log10(obj2->flux_iso) + prefs.mag_zeropoint[0]
			:99.0;
  if (FLAG(obj2.magerr_iso))
    obj2->magerr_iso = obj2->flux_iso>0.0?
			 1.086*obj2->fluxerr_iso/obj2->flux_iso
			:99.0;

/* Mag. isophotal corrected */
  if (FLAG(obj2.mag_isocor))
    obj2->mag_isocor = obj2->flux_isocor>0.0?
			 -2.5*log10(obj2->flux_isocor) + prefs.mag_zeropoint[0]
			:99.0;
  if (FLAG(obj2.magerr_isocor))
    obj2->magerr_isocor = obj2->flux_isocor>0.0?
			 1.086*obj2->fluxerr_isocor/obj2->flux_isocor
			:99.0;

/* Choose the ``best'' flux according to the local crowding */

  if (FLAG(obj2.flux_best))
    {
    if (obj2->flags&OBJ_CROWDED)
      {
      obj2->flux_best = obj2->flux_isocor;
      obj2->fluxerr_best = obj2->fluxerr_isocor;
      }
    else
      {
      obj2->flux_best = obj2->flux_auto[0];
      obj2->fluxerr_best = obj2->fluxerr_auto[0];
      }
    }

/* Mag. Best */
  if (FLAG(obj2.mag_best))
    obj2->mag_best = obj2->flux_best>0.0?
			 -2.5*log10(obj2->flux_best) + prefs.mag_zeropoint[0]
			:99.0;
  if (FLAG(obj2.magerr_best))
    obj2->magerr_best = obj2->flux_best>0.0?
			 1.086*obj2->fluxerr_best/obj2->flux_best
			:99.0;

/* Mag. SOM-fit */
  if (FLAG(obj2.mag_somfit))
    obj2->mag_somfit = obj2->flux_somfit>0.0?
			 -2.5*log10(obj2->flux_somfit) + prefs.mag_zeropoint[0]
			:99.0;
  if (FLAG(obj2.magerr_somfit))
    obj2->magerr_somfit = obj2->flux_somfit>0.0?
			 1.086*obj2->fluxerr_somfit/obj2->flux_somfit
			:99.0;

/* Mag. models */
  if (FLAG(obj2.magcor_prof))
    obj2->magcor_prof = obj2->fluxcor_prof>0.0?
			 -2.5*log10(obj2->fluxcor_prof) + prefs.mag_zeropoint[0]
			:99.0;
  if (FLAG(obj2.magcorerr_prof))
    obj2->magcorerr_prof = obj2->fluxcor_prof>0.0?
			 1.086*obj2->fluxcorerr_prof/obj2->fluxcor_prof
			:99.0;

  if (FLAG(obj2.mumax_prof))
    obj2->mumax_prof = obj2->peak_prof > 0.0 ?
		-2.5*log10((obj2->peak_prof)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.mueff_prof))
    obj2->mueff_prof = obj2->fluxeff_prof > 0.0 ?
		-2.5*log10((obj2->fluxeff_prof)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.mumean_prof))
    obj2->mumean_prof = obj2->fluxmean_prof > 0.0 ?
		-2.5*log10((obj2->fluxmean_prof)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.prof_spheroid_mumax))
    obj2->prof_spheroid_mumax = obj2->prof_spheroid_peak > 0.0 ?
		-2.5*log10((obj2->prof_spheroid_peak)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.prof_spheroid_mueff))
    obj2->prof_spheroid_mueff = obj2->prof_spheroid_fluxeff > 0.0 ?
		-2.5*log10((obj2->prof_spheroid_fluxeff)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.prof_spheroid_mumean))
    obj2->prof_spheroid_mumean = obj2->prof_spheroid_fluxmean > 0.0 ?
		-2.5*log10((obj2->prof_spheroid_fluxmean)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.prof_disk_mumax))
    obj2->prof_disk_mumax = obj2->prof_disk_peak > 0.0 ?
		-2.5*log10((obj2->prof_disk_peak)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.prof_disk_mueff))
    obj2->prof_disk_mueff = obj2->prof_disk_fluxeff > 0.0 ?
		-2.5*log10((obj2->prof_disk_fluxeff)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.prof_disk_mumean))
    obj2->prof_disk_mumean = obj2->prof_disk_fluxmean > 0.0 ?
		-2.5*log10((obj2->prof_disk_fluxmean)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.prof_bar_mag))
    obj2->prof_bar_mag = obj2->prof_bar_flux>0.0?
			 -2.5*log10(obj2->prof_bar_flux)
			+ prefs.mag_zeropoint[0]
			:99.0;
  if (FLAG(obj2.prof_bar_magerr))
    obj2->prof_bar_magerr = obj2->prof_bar_flux>0.0?
			 1.086*obj2->prof_bar_fluxerr
				/obj2->prof_bar_flux
			:99.0;

  if (FLAG(obj2.prof_arms_mag))
    obj2->prof_arms_mag = obj2->prof_arms_flux>0.0?
			 -2.5*log10(obj2->prof_arms_flux)
			+ prefs.mag_zeropoint[0]
			:99.0;
  if (FLAG(obj2.prof_arms_magerr))
    obj2->prof_arms_magerr = obj2->prof_arms_flux>0.0?
			 1.086*obj2->prof_arms_fluxerr
				/obj2->prof_arms_flux
			:99.0;


/* Other surface brightnesses */
  if (FLAG(obj2.maxmu))
    obj2->maxmu = obj2->peak > 0.0 ?
		-2.5*log10((obj2->peak)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;

  if (FLAG(obj2.threshmu))
    obj2->threshmu = obj2->thresh > 0.0 ?
		-2.5*log10((obj2->thresh)
		 / (prefs.pixel_scale? field->pixscale*field->pixscale
				: obj2->pixscale2 * 3600.0*3600.0))
		+ prefs.mag_zeropoint[0]
		: 99.0;


  return;
  }


/****** photom_printinstruinfo **********************************************
PROTO	void photom_printinstruinfo(void)
PURPOSE	Print info about photometric instruments found
INPUT	-.
OUTPUT	-.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	21/03/2012
 ***/
void	photom_printinstruinfo(void)
  {
   int		i,l,len;

  QPRINTF(OUTPUT, "----- %d %s found for photometry:\n",
	prefs.nphotinstru, prefs.nphotinstru>1? "instruments":"instrument");
  for (i=0; i<prefs.nphotinstru; i++)
    {
    QPRINTF(OUTPUT, "\nInstrument P%-2d:\n", i+1);
    len = fitsfind(prefs.photinstrustr[i], "END     ");
    for (l=0; l<len; l++)
      QPRINTF(OUTPUT, "%.80s\n", prefs.photinstrustr[i]+l*80);
    }

  return;
  }


