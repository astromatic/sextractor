/*
*				winpos.c
*
* Compute iteratively windowed parameters.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1994,1997 ESO
*	          		(C) 1995,1996 Leiden Observatory 
*	          		(C) 1998-2021 IAP/CNRS/SorbonneU
*	          		(C)	2021-2023 CFHT/CNRS
*	          		(C) 2023-2025 CEA/AIM/UParisSaclay
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
*	Last modified:		19/03/2025
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
#include	"winpos.h"

static  obj2struct	*obj2 = &outobj2;

/****** compute_winpos ********************************************************
PROTO	void compute_winpos(picstruct *field, picstruct *wfield,
			objstruct *obj)
PURPOSE	Compute windowed source barycenter.
INPUT	Picture structure pointer,
	Weight-map structure pointer,
	object structure.
OUTPUT  -.
NOTES   obj->posx and obj->posy are taken as initial centroid guesses.
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	19/03/2025
 ***/
void	compute_winpos(picstruct *field, picstruct *wfield, 
			picstruct *dgeofield, objstruct *obj)

  {
   float		r2,invtwosig2, raper,raper2, rintlim,rintlim2,rextlim2,
			dx,dx1,dy,dy2, ddx, ddy, sig, invngamma,
                        offsetx,offsety,scalex,scaley,scale2, locarea;
   double               tv, norm, pix, var, backnoise2, invgain, locpix,
			dxpos,dypos, ddxpos, ddypos, err,err2, emx2,emy2,emxy,
			esum, temp,temp2, mx2, my2,mxy,pmx2, theta, mx,my,
			dmx, dmy, mx2ph, my2ph, cx2, cy2, cxy;
   int                  i,x,y, x2,y2, xmin,xmax,ymin,ymax, sx,sy, w,h,
                        fymin,fymax, pflag,corrflag, gainflag, errflag,
			momentflag;
   long                 pos;
   PIXTYPE              *strip,*stript, *wstrip, *wstript,
			*dgeoxstrip, *dgeoystrip, *dgeoxstript, *dgeoystript,
                        wthresh = 0.0;

  if (wfield)
    wthresh = wfield->weight_thresh;
  wstrip = wstript = NULL;
  w = field->width;
  h = field->stripheight;
  fymin = field->ymin;
  fymax = field->ymax;
  pflag = (prefs.detect_type==PHOTO)? 1:0;
  corrflag = (prefs.mask_type==MASK_CORRECT);
  gainflag = wfield && prefs.weightgain_flag;
  errflag = FLAG(obj2.winposerr_mx2) | FLAG(obj2.fluxerr_win);
  momentflag = FLAG(obj2.win_mx2) | FLAG(obj2.winposerr_mx2);
  var = backnoise2 = field->backsig*field->backsig;
  invgain = field->gain>0.0? 1.0/field->gain : 0.0;
  sig = obj2->hl_radius*2.0/2.35; /* From half-FWHM to sigma */
  invtwosig2 = 1.0/(2.0*sig*sig);

/* Integration radius */
  raper = WINPOS_NSIG*sig;

/* For photographic data */
  invngamma = pflag? 1./field->ngamma : 0.;

  raper2 = raper*raper;
/* Internal radius of the oversampled annulus (<r-sqrt(2)/2) */
  rintlim = raper - 0.75;
  rintlim2 = (rintlim>0.0)? rintlim*rintlim: 0.0;
/* External radius of the oversampled annulus (>r+sqrt(2)/2) */
  rextlim2 = (raper + 0.75)*(raper + 0.75);
  scaley = scalex = 1.0/WINPOS_OVERSAMP;
  scale2 = scalex*scaley;
  offsetx = 0.5*(scalex-1.0);
  offsety = 0.5*(scaley-1.0);
/* Use isophotal centroid as a first guess */
  mx = obj2->posx - 1.0;
  my = obj2->posy - 1.0;
  dmx = obj2->pos_dgeox;
  dmy = obj2->pos_dgeoy;

  for (i=0; i<WINPOS_NITERMAX; i++)
    {
    xmin = (int)(mx-raper+0.499999);
    xmax = (int)(mx+raper+1.499999);
    ymin = (int)(my-raper+0.499999);
    ymax = (int)(my+raper+1.499999);
    mx2ph = (mx+dmx)*2.0 + 0.49999;	// Pixel position, not the physical one
    my2ph = (my+dmy)*2.0 + 0.49999;

    if (xmin < 0)
      {
      xmin = 0;
      obj->flag |= OBJ_APERT_PB;
      }
    if (xmax > w)
      {
      xmax = w;
      obj->flag |= OBJ_APERT_PB;
      }
    if (ymin < fymin)
      {
      ymin = fymin;
      obj->flag |= OBJ_APERT_PB;
      }
    if (ymax > fymax)
      {
      ymax = fymax;
      obj->flag |= OBJ_APERT_PB;
      }

    tv = esum = emxy = emx2 = emy2 = mx2 = my2 = mxy = 0.0;
    dxpos = dypos = ddxpos = ddypos = 0.0;
    strip = field->strip;
    wstrip = wstript = NULL;		/* To avoid gcc -Wall warnings */
    if (wfield)
      wstrip = wfield->strip;
    if (dgeofield) {			// Differential geometry map is present
      dgeoxstrip = dgeofield->dgeostrip[0];
      dgeoystrip = dgeofield->dgeostrip[1];
    }
    for (y=ymin; y<ymax; y++)
      {
      stript = strip + (pos = (y%h)*w + xmin);
      if (wfield)
        wstript = wstrip + pos;
      if (dgeofield) {
        dgeoxstript = dgeoxstrip + pos;
        dgeoystript = dgeoystrip + pos;
      }
      for (x=xmin; x<xmax; x++)
        {
        dx = x - mx;
        dy = y - my;
        if (dgeofield) {
          dx -= (ddx = *(dgeoxstript++));
          dy -= (ddy = *(dgeoystript++));
        }
        if ((r2=dx*dx+dy*dy)<rextlim2)
          {
          if (WINPOS_OVERSAMP > 1 && r2 > rintlim2)
            {
            dx += offsetx;
            dy += offsety;
            locarea = 0.0;
            for (sy=WINPOS_OVERSAMP; sy--; dy+=scaley)
              {
              dx1 = dx;
              dy2 = dy*dy;
              for (sx=WINPOS_OVERSAMP; sx--; dx1+=scalex)
                if (dx1*dx1+dy2<raper2)
                  locarea += scale2;
              }
            }
          else
            locarea = 1.0;
          locarea *= expf(-r2*invtwosig2);
/*-------- Here begin tests for pixel and/or weight overflows. Things are a */
/*-------- bit intricated to have it running as fast as possible in the most */
/*-------- common cases */
          if ((pix=*stript)<=-BIG || (wfield && (var=*wstript)>=wthresh))
            {
            if (corrflag
		&& (x2=(int)(mx2ph-x))>=0 && x2<w
		&& (y2=(int)(my2ph-y))>=fymin && y2<fymax
		&& (pix=*(strip + (pos = (y2%h)*w + x2)))>-BIG)
              {
              if (wfield)
                {
                var = *(wstrip + pos);
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
            pix = expf(pix*invngamma);
          locpix = locarea*pix;
          tv += locpix;
          dx = x - mx;
          dy = y - my;
          if (dgeofield) {
            dx -= ddx;
            dy -= ddy;
            ddxpos += locpix*ddx;
            ddypos += locpix*ddy;
          }
          dxpos += locpix*dx;
          dypos += locpix*dy;
          if (errflag)
            {
            err = var;
            if (pflag)
              err *= locpix*pix*invngamma*invngamma;
            else if (invgain>0.0 && pix>0.0)
              {
              if (gainflag)
                err += pix*invgain*var/backnoise2;
              else
                err += pix*invgain;
              }
            err2 = locarea*locarea*err;
            esum += err2;
            emx2 += err2*(dx*dx+0.0833);	/* Finite pixel size */
            emy2 += err2*(dy*dy+0.0833);	/* Finite pixel size */
            emxy += err2*dx*dy;
            }
          mx2 += locpix*dx*dx;
          my2 += locpix*dy*dy;
          mxy += locpix*dx*dy;
          }
        stript++;
        if (wfield)
          wstript++;
        }
      }

    if (tv>0.0)
      {
      mx += (dxpos /= tv) * WINPOS_FAC;
      my += (dypos /= tv) * WINPOS_FAC;
      }
    else
      break;

/*-- Stop here if position does not change */
    if (dxpos*dxpos+dypos*dypos < WINPOS_STEPMIN*WINPOS_STEPMIN)
      break;
    }

  mx2 = mx2/tv - dxpos*dxpos;
  my2 = my2/tv - dypos*dypos;
  mxy = mxy/tv - dxpos*dypos;
  obj2->winpos_x = mx + 1.0;	/* The dreaded 1.0 FITS offset */
  obj2->winpos_y = my + 1.0; 	/* The dreaded 1.0 FITS offset */
  if (dgeofield) {
    obj2->winpos_dgeox = (ddxpos / tv) * WINPOS_FAC;
    obj2->winpos_dgeoy = (ddypos / tv) * WINPOS_FAC;
  }

  obj2->winpos_niter = i+1;

/* WINdowed flux */
  if (FLAG(obj2.flux_win))
    {
    obj2->flux_win = tv;
    obj2->fluxerr_win = sqrt(esum);
    obj2->snr_win = esum>(1.0/BIG)? obj2->flux_win / obj2->fluxerr_win: BIG;
    }
  temp2=mx2*my2-mxy*mxy;
  cx2 = obj->cxx;
  cy2 = obj->cyy;
  cxy = obj->cxy;
  dx = obj2->winpos_x - obj2->posx;
  dy = obj2->winpos_y - obj2->posy;
  obj2->win_flag = 8 * (dx*dx > 1.0 && dy*dy > 1.0
			&& (cx2*dx*dx + cy2*dy*dy + cxy*dx*dy) > 1.0)
			// Check that the difference vector with isophotal
			// centroid does not lie outside of the 1-sigma
			// footprint ellipse
		 | 4 * (tv <= 0.0)	// Check negative flux
		 | 2 * (mx2 < 0.0 || my2 < 0.0)	// Negative 2nd order moments
		 | (temp2<0.0);
  if (obj2->win_flag)
    {
/*--- Negative values: revert to isophotal estimates */
    obj2->winpos_x = obj2->posx;
    obj2->winpos_y = obj2->posy;
    if (FLAG(obj2.winposerr_mx2))
      {
      obj2->winposerr_mx2 = obj->poserr_mx2;
      obj2->winposerr_my2 = obj->poserr_my2;
      obj2->winposerr_mxy = obj->poserr_mxy;
      if (FLAG(obj2.winposerr_a))
        {
        obj2->winposerr_a = obj2->poserr_a;
        obj2->winposerr_b = obj2->poserr_b;
        obj2->winposerr_theta = obj2->poserr_theta;
        }
      if (FLAG(obj2.winposerr_cxx))
        {
        obj2->winposerr_cxx = obj2->poserr_cxx;
        obj2->winposerr_cyy = obj2->poserr_cyy;
        obj2->winposerr_cxy = obj2->poserr_cxy;
        }
      }
    if (momentflag)
      {
      obj2->win_mx2 = obj->mx2;
      obj2->win_my2 = obj->my2;
      obj2->win_mxy = obj->mxy;
      if (FLAG(obj2.win_cxx))
        {
        obj2->win_cxx = obj->cxx;
        obj2->win_cyy = obj->cyy;
        obj2->win_cxy = obj->cxy;
        }
      if (FLAG(obj2.win_a))
        {
        obj2->win_a = obj->a;
        obj2->win_b = obj->b;
        obj2->win_polar = obj2->polar;
        obj2->win_theta = obj->theta;
        }
      }
    }
  else
    {
    if (FLAG(obj2.winposerr_mx2))
      {
      norm = WINPOS_FAC*WINPOS_FAC/(tv*tv);
      emx2 *= norm;
      emy2 *= norm;
      emxy *= norm;
/*-- Handle fully correlated profiles (which cause a singularity...) */
      esum *= 0.08333*norm;
      if (obj->singuflag && (emx2*emy2-emxy*emxy) < esum*esum)
        {
        emx2 += esum;
        emy2 += esum;
        }

      obj2->winposerr_mx2 = emx2;
      obj2->winposerr_my2 = emy2;
      obj2->winposerr_mxy = emxy;
/*---- Error ellipse parameters */
      if (FLAG(obj2.winposerr_a))
        {
         double	pmx2,pmy2,temp,theta;

        if (fabs(temp=emx2-emy2) > 0.0)
          theta = atan2(2.0 * emxy,temp) / 2.0;
        else
          theta = PI/4.0;

        temp = sqrt(0.25*temp*temp+ emxy*emxy);
        pmy2 = pmx2 = 0.5*(emx2+emy2);
        pmx2+=temp;
        pmy2-=temp;

        obj2->winposerr_a = (float)sqrt(pmx2);
        obj2->winposerr_b = (float)sqrt(pmy2);
        obj2->winposerr_theta = theta*180.0/PI;
        }

      if (FLAG(obj2.winposerr_cxx))
        {
         double	temp;

        obj2->winposerr_cxx = (float)(emy2/(temp=emx2*emy2-emxy*emxy));
        obj2->winposerr_cyy = (float)(emx2/temp);
        obj2->winposerr_cxy = (float)(-2*emxy/temp);
        }
      }

    if (momentflag)
      {
/*-- Handle fully correlated profiles (which cause a singularity...) */
      if ((temp2=mx2*my2-mxy*mxy)<0.00694)
        {
        mx2 += 0.0833333;
        my2 += 0.0833333;
        temp2 = mx2*my2-mxy*mxy;
        }
      obj2->win_mx2 = mx2;
      obj2->win_my2 = my2;
      obj2->win_mxy = mxy;

      if (FLAG(obj2.win_cxx))
        {
        obj2->win_cxx = (float)(my2/temp2);
        obj2->win_cyy = (float)(mx2/temp2);
        obj2->win_cxy = (float)(-2*mxy/temp2);
        }

      if (FLAG(obj2.win_a))
        {
        if ((fabs(temp=mx2-my2)) > 0.0)
          theta = atan2(2.0 * mxy,temp) / 2.0;
        else
          theta = PI/4.0;

        temp = sqrt(0.25*temp*temp+mxy*mxy);
        pmx2 = 0.5*(mx2+my2);
        obj2->win_a = (float)sqrt(pmx2 + temp);
        obj2->win_b = (float)sqrt(pmx2 - temp);
        if (FLAG(obj2.win_polar))
          obj2->win_polar = temp / pmx2;
        obj2->win_theta = theta*180.0/PI;
        }
      }
    }

  return;
  }


