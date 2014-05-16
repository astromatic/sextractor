/*
*				winpos.c
*
* Compute iteratively windowed parameters.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2005-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		02/04/2012
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
#include	"subimage.h"
#include	"winpos.h"

/****** win_pos **************************************************************
PROTO	void win_pos(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2struct *obj2)
PURPOSE	Compute windowed source barycenter.
INPUT	Pointer to an array of image field pointers,
	pointer to an array of weight-map field pointers,
	number of images,
	obj2 structure.
OUTPUT  -.
NOTES   obj2->mx and obj2->my are taken as initial centroid guesses.
AUTHOR  E. Bertin (IAP)
VERSION 02/04/2012
 ***/
void	win_pos(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2struct *obj2)
  {
   fieldstruct		*field, *wfield;
   subimagestruct	*subimage;
   float		*sposx,*sposy,
			r2,invtwosig2, raper,raper2, rintlim2, rextlim2,
			dx,dx1,dy,dy2, sig, invngamma, pdbkg,
                        offsetx,offsety,scalex,scaley,scale2, locarea;
   double               tv, norm, pix, var, backnoise2, invgain, locpix,
			dxpos,dypos, err,err2, emx2,emy2,emxy,
			esum, temp,temp2, mx2, my2,mxy,pmx2, theta, mx,my,
			mx2ph, my2ph, ftv, fdxpos,fdypos, femx2,femy2,femxy,
			fesum, fmx2,fmy2,fmxy, wsum, fw;
   int                  *sokflag,
			f,i,s,x,y, x2,y2, xmin,xmax,ymin,ymax, sx,sy, w,h,
                        pflag,corrflag, gainflag, errflag, momentflag, winflag,
			band,nband, nsubimage;
   long                 pos;
   PIXTYPE              *image, *imaget, *weight,*weightt,
                        wthresh = 0.0;

  corrflag = (prefs.mask_type==MASK_CORRECT);
  errflag = FLAG(obj2.winposerr_mx2) | FLAG(obj2.fluxerr_win);
  momentflag = FLAG(obj2.win_mx2) | FLAG(obj2.winposerr_mx2);
  sig = obj2->hl_radius*2.0/2.35482;	/* From half FWHM to sigma */
  invtwosig2 = 1.0/(2.0*sig*sig);

/* Use isophotal centroid as a first guess */
  nsubimage = obj2->nsubimage;
  QMALLOC(sposx, float, nsubimage);
  QMALLOC(sposy, float, nsubimage);
  QMALLOC(sokflag, int, nsubimage);
  subimage = obj2->subimage;
  for (s=0; s<nsubimage; s++, subimage++)
    {
/*-- Store positions in sub-images */
    sposx[s] = subimage->dpos[0] - 1.0 - subimage->immin[0];
    sposy[s] = subimage->dpos[1] - 1.0 - subimage->immin[1];
    sokflag[s] = 1;
    }

/* Initialize WIN measurements to "N/A" values */
  for (f=0; f<nfield; f++)
    {
    if (FLAG(obj2.winpos_x))
     obj2->winpos_x[f] = 0.0;
    if (FLAG(obj2.winpos_y))
     obj2->winpos_y[f] = 0.0;
    if (FLAG(obj2.winpos_niter))
      obj2->winpos_niter[f] = 0;
    if (FLAG(obj2.flux_win))
      obj2->flux_win[f] = 0.0;
    if (FLAG(obj2.fluxerr_win))
      obj2->fluxerr_win[f] = 0.0;
    if (FLAG(obj2.snr_win))
      obj2->snr_win[f] = 0.0;
    if (FLAG(obj2.win_flags))
      obj2->win_flags[f] = WINFLAG_NOTCOVERED;
    if (errflag)
      {
      if (FLAG(obj2.winposerr_mx2))
        obj2->winposerr_mx2[f] = 0.0;
      if (FLAG(obj2.winposerr_my2))
        obj2->winposerr_my2[f] = 0.0;
      if (FLAG(obj2.winposerr_mxy))
        obj2->winposerr_mxy[f] = 0.0;
      if (FLAG(obj2.winposerr_a))
        obj2->winposerr_a[f] = 0.0;
      if (FLAG(obj2.winposerr_b))
        obj2->winposerr_b[f] = 0.0;
      if (FLAG(obj2.winposerr_theta))
        obj2->winposerr_theta[f] = 0.0;
      if (FLAG(obj2.winposerr_cxx))
        obj2->winposerr_cxx[f] = 0.0;
      if (FLAG(obj2.winposerr_cyy))
        obj2->winposerr_cyy[f] = 0.0;
      if (FLAG(obj2.winposerr_cxy))
        obj2->winposerr_cxy[f] = 0.0;
      }
    if (momentflag)
      {
      if (FLAG(obj2.win_mx2))
        obj2->win_mx2[f] = 0.0;
      if (FLAG(obj2.win_my2))
        obj2->win_my2[f] = 0.0;
      if (FLAG(obj2.win_mxy))
        obj2->win_mxy[f] = 0.0;
      if (FLAG(obj2.win_cxx))
        obj2->win_cxx[f] = 0.0;
      if (FLAG(obj2.win_cyy))
        obj2->win_cyy[f] = 0.0;
      if (FLAG(obj2.win_cxy))
        obj2->win_cxy[f] = 0.0;
      if (FLAG(obj2.win_a))
        obj2->win_a[f] = 0.0;
      if (FLAG(obj2.win_b))
        obj2->win_b[f] = 0.0;
      if (FLAG(obj2.win_theta))
        obj2->win_theta[f] = 0.0;
      }
    }

  nband = prefs.nphotinstru;
  nsubimage = obj2->nsubimage;
  for (band=0; band<nband; band++)
    {
    for (i=0; i<WINPOS_NITERMAX; i++)
      {
      tv = wsum = esum = emxy = emx2 = emy2 = mx2 = my2 = mxy = 0.0;
      dxpos = dypos = 0.0;
      subimage = obj2->subimage;
      for (s=0; s<nsubimage; s++, subimage++)
        {
        field = subimage->field;
        if (!sokflag[s] || field->photomlabel != band)
          continue;
        ftv = fesum = femxy = femx2 = femy2 = fmx2 = fmy2 = fmxy = 0.0;
        fdxpos = fdypos = 0.0;
        wfield = subimage->wfield;
        if (wfield)
          wthresh = wfield->weight_thresh;
/*------ For photographic data */
        pflag = (field->detector_type==DETECTOR_PHOTO);
        if (pflag)
          {
          invngamma = 1.0/field->ngamma;
          pdbkg = expf(obj2->dbkg[0]*invngamma);
          }
        else
          {
          invngamma = 0.0;
          pdbkg = 0.0;
          }
        gainflag = wfield && field->weightgain_flag;
        var = backnoise2 = field->backsig*field->backsig;
        invgain = field->gain>0.0? 1.0/field->gain : 0.0;

        mx = sposx[s];
        my = sposy[s];
        raper = WINPOS_NSIG*sig*subimage->dscale;
        raper2 = raper*raper;
/*------ Internal radius of the oversampled annulus (<r-sqrt(2)/2) */
        rintlim2 = raper - 0.75;
        rintlim2 = (rintlim2>0.0)? rintlim2*rintlim2: 0.0;

/*------ External radius of the oversampled annulus (>r+sqrt(2)/2) */
        rextlim2 = (raper + 0.75)*(raper + 0.75);
        scaley = scalex = 1.0/WINPOS_OVERSAMP;
        scale2 = scalex*scaley;
        offsetx = 0.5*(scalex-1.0);
        offsety = 0.5*(scaley-1.0);

        xmin = (int)(mx-raper+0.499999);
        xmax = (int)(mx+raper+1.499999);
        ymin = (int)(my-raper+0.499999);
        ymax = (int)(my+raper+1.499999);
        mx2ph = mx*2.0 + 0.49999;
        my2ph = my*2.0 + 0.49999;

        w = subimage->imsize[0];
        h = subimage->imsize[1];
        if (xmin < 0)
          {
          xmin = 0;
          winflag = WINFLAG_APERT_PB;
          }
        if (xmax > w)
          {
          xmax = w;
          winflag = WINFLAG_APERT_PB;
          }
        if (ymin < 0)
          {
          ymin = 0;
          winflag = WINFLAG_APERT_PB;
          }
        if (ymax > h)
          {
          ymax = h;
          winflag = WINFLAG_APERT_PB;
          }

        image = subimage->image;
        weight = weightt = NULL;	/* To avoid gcc -Wall warnings */
        if (wfield)
          weight = subimage->weight;
        for (y=ymin; y<ymax; y++)
          {
          imaget = image + (pos = y*w + xmin);
          if (wfield)
            weightt = weight + pos;
          for (x=xmin; x<xmax; x++, imaget++, weightt++)
            {
            dx = x - mx;
            dy = y - my;
            if ((r2=dx*dx+dy*dy)<rextlim2)
              {
              if (WINPOS_OVERSAMP>1 && r2>rintlim2)
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
/*------------ Here begin tests for pixel/weight overflows. Things are a */
/*------------ bit intricated to have it running as fast as possible in the */
/*------------ most common cases */
              if ((pix=*imaget)<=-BIG || (wfield && (var=*weightt)>=wthresh))
                {
                if (corrflag
			&& (x2=(int)(mx2ph-x))>=0 && x2<w
			&& (y2=(int)(my2ph-y))>=0 && y2<h
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
                pix = expf(pix*invngamma);
              dx = x - mx;
              dy = y - my;
              locpix = locarea*pix;
              ftv += locpix;
              fdxpos += locpix*dx;
              fdypos += locpix*dy;
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
              fesum += err2;
              femx2 += err2*(dx*dx+0.0833);	/* Finite pixel size */
              femy2 += err2*(dy*dy+0.0833);	/* Finite pixel size */
              femxy += err2*dx*dy;
              if (momentflag)
                {
                fmx2 += locpix*dx*dx;
                fmy2 += locpix*dy*dy;
                fmxy += locpix*dx*dy;
                }
              }
            }
          }
        if (fesum <= 1.0/BIG)
          continue;
        fw = 1.0/fesum;
        wsum += fw;
        esum += 1.0;
        tv += fw*ftv;
        dxpos += fw*fdxpos;
        dypos += fw*fdypos;
        emx2 += fw*femx2;
        emy2 += fw*femy2;
        emxy += fw*femxy;

        if (momentflag)
          {
          mx2 += fw*fmx2;
          my2 += fw*fmy2;
          mxy += fw*fmxy;
          }
        }

      if (wsum>0.0)
        {
        wsum = 1.0/wsum;
        esum *= wsum;
        tv *= wsum;
        dxpos *= wsum;
        dypos *= wsum;
        emx2 *= wsum;
        emy2 *= wsum;
        emxy *= wsum;
        if (momentflag)
          {
          mx2 *= wsum;
          my2 *= wsum;
          mxy *= wsum;
          }
        }

/*---- Exit here if the summed flux is negative */
      if (tv<=0.0)
        break;

      dxpos /= tv;
      dypos /= tv;

/*---- Exit here if the incremental displacement is too small */
      if (dxpos*dxpos+dypos*dypos < WINPOS_STEPMIN*WINPOS_STEPMIN)
        break;

/*---- Introduce a small dampening factor before iterating */
      dxpos *= WINPOS_FAC;
      dypos *= WINPOS_FAC;
      for (s=0; s<obj2->nsubimage; s++)
        {
        sposx[s] += dxpos;
        sposy[s] += dypos;
        }
      }

    mx2 = mx2/tv - dxpos*dxpos;
    my2 = my2/tv - dypos*dypos;
    mxy = mxy/tv - dxpos*dypos;

    temp2=mx2*my2-mxy*mxy;

    winflag =  (tv <= 0.0)*WINFLAG_NEGFLUX
		+ (mx2 < 0.0 || my2 < 0.0)*WINFLAG_NEGMOMENT
		+ (temp2<0.0)*WINFLAG_SINGULAR;

/*-- Write final parameters */
    subimage = obj2->subimage;
    for (s=0; s<nsubimage; s++, subimage++)
      {
      field = subimage->field;
      if (!sokflag[s] || field->photomlabel != band)
        continue;

      f = field->imindex;
      if (FLAG(obj2.winpos_x))
        obj2->winpos_x[f] = sposx[s] + subimage->immin[0]
				+ 1.0;	/* Mind the 1 pixel FITS offset */
      if (FLAG(obj2.winpos_y))
        obj2->winpos_y[f] = sposy[s] + subimage->immin[1]
				+ 1.0;	/* Mind the 1 pixel FITS offset */
      if (FLAG(obj2.winpos_niter))
        obj2->winpos_niter[f] = i+1;

/*---- WINdowed flux */
      if (FLAG(obj2.flux_win))
        obj2->flux_win[f] = tv;
      if (FLAG(obj2.fluxerr_win))
        obj2->fluxerr_win[f] = sqrt(esum);
      if (FLAG(obj2.snr_win))
        obj2->snr_win[f] = esum>(1.0/BIG)? tv / sqrt(esum) : BIG;
      
      if (FLAG(obj2.win_flags))
        obj2->win_flags[f] = winflag;

      if (winflag)
        {      
/*------- Negative values: revert to isophotal estimates */
        if (FLAG(obj2.winposerr_mx2))
          {
          obj2->winposerr_mx2[f] = obj2->poserr_mx2;
          obj2->winposerr_my2[f] = obj2->poserr_my2;
          obj2->winposerr_mxy[f] = obj2->poserr_mxy;
          if (FLAG(obj2.winposerr_a))
            {
            obj2->winposerr_a[f] = obj2->poserr_a;
            obj2->winposerr_b[f] = obj2->poserr_b;
            obj2->winposerr_theta[f] = obj2->poserr_theta;
            }
          if (FLAG(obj2.winposerr_cxx))
            {
            obj2->winposerr_cxx[f] = obj2->poserr_cxx;
            obj2->winposerr_cyy[f] = obj2->poserr_cyy;
            obj2->winposerr_cxy[f] = obj2->poserr_cxy;
            }
          }
        if (momentflag)
          {
          obj2->win_mx2[f] = obj2->mx2;
          obj2->win_my2[f] = obj2->my2;
          obj2->win_mxy[f] = obj2->mxy;
          if (FLAG(obj2.win_cxx))
            {
            obj2->win_cxx[f] = obj2->cxx;
            obj2->win_cyy[f] = obj2->cyy;
            obj2->win_cxy[f] = obj2->cxy;
            }
          if (FLAG(obj2.win_a))
            {
            obj2->win_a[f] = obj2->a;
            obj2->win_b[f] = obj2->b;
            obj2->win_polar[f] = obj2->polar;
            obj2->win_theta[f] = obj2->theta;
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
/*-------- Handle fully correlated profiles (which cause a singularity...) */
          esum *= 0.08333*norm;
          if (obj2->singuflag && (emx2*emy2-emxy*emxy) < esum*esum)
            {
            emx2 += esum;
            emy2 += esum;
            }

          obj2->winposerr_mx2[f] = emx2;
          obj2->winposerr_my2[f] = emy2;
          obj2->winposerr_mxy[f] = emxy;
/*-------- Error ellipse parameters */
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

            obj2->winposerr_a[f] = (float)sqrt(pmx2);
            obj2->winposerr_b[f] = (float)sqrt(pmy2);
            obj2->winposerr_theta[f] = theta*180.0/PI;
            }

          if (FLAG(obj2.winposerr_cxx))
            {
             double	temp;

            obj2->winposerr_cxx[f] = (float)(emy2/(temp=emx2*emy2-emxy*emxy));
            obj2->winposerr_cyy[f] = (float)(emx2/temp);
            obj2->winposerr_cxy[f] = (float)(-2*emxy/temp);
            }
          }

        if (momentflag)
          {
/*------ Handle fully correlated profiles (which cause a singularity...) */
          if ((temp2=mx2*my2-mxy*mxy)<0.00694)
            {
            mx2 += 0.0833333;
            my2 += 0.0833333;
            temp2 = mx2*my2-mxy*mxy;
            }
          obj2->win_mx2[f] = mx2;
          obj2->win_my2[f] = my2;
          obj2->win_mxy[f] = mxy;

          if (FLAG(obj2.win_cxx))
            {
            obj2->win_cxx[f] = (float)(my2/temp2);
            obj2->win_cyy[f] = (float)(mx2/temp2);
            obj2->win_cxy[f] = (float)(-2*mxy/temp2);
            }

          if (FLAG(obj2.win_a))
            {
            if ((fabs(temp=mx2-my2)) > 0.0)
              theta = atan2(2.0 * mxy,temp) / 2.0;
            else
              theta = PI/4.0;

            temp = sqrt(0.25*temp*temp+mxy*mxy);
            pmx2 = 0.5*(mx2+my2);
            obj2->win_a[f] = (float)sqrt(pmx2 + temp);
            obj2->win_b[f] = (float)sqrt(pmx2 - temp);
            if (FLAG(obj2.win_polar))
              obj2->win_polar[f] = temp / pmx2;
            obj2->win_theta[f] = theta*180.0/PI;
            }
          }
        }
      }
    }

  free(sposx);
  free(sposy);
  free(sokflag);

  return;
  }



