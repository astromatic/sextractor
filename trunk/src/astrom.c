/*
*				astrom.c
*
* High level astrometric computations.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		11/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdlib.h>
#include	<string.h>

#include 	"wcs/wcs.h"
#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"astrom.h"
#include	"fitswcs.h"
#include	"wcs/tnx.h"

static obj2struct	*obj2 = &outobj2;

/****************************** initastrom **********************************/
/*
Initialize astrometrical structures.
*/
void	initastrom(picstruct *field)

  {
   wcsstruct	*wcs;

  wcs = field->wcs;

/* Test if the WCS is in use */
  if (wcs->lng != wcs->lat)
    {
    if (FLAG(obj2.dtheta2000) || FLAG(obj2.dtheta1950))
      {
      if (fabs(wcs->equinox-2000.0)>0.003)
        precess(wcs->equinox, 0.0, 90.0, 2000.0, &wcs->ap2000, &wcs->dp2000);
      else
        {
        wcs->ap2000 = 0.0;
        wcs->dp2000 = 90.0;
        }

      if (FLAG(obj2.theta1950) || FLAG(obj2.poserr_theta1950))
        j2b(wcs->equinox, wcs->ap2000, wcs->dp2000, &wcs->ap1950, &wcs->dp1950);
      }
    }

/* Override astrometric definitions only if user supplies a pixel-scale */
  if (prefs.pixel_scale == 0.0)
    field->pixscale = wcs->pixscale*3600.0;	/* in arcsec */
  else
    field->pixscale = prefs.pixel_scale;

  return;
  }


/***************************** astrom_pos **********************************/
/*
Compute real FOCAL and WORLD coordinates according to FITS info.
*/
void	astrom_pos(picstruct *field, objstruct *obj)

  {
   wcsstruct	*wcs;
   double	rawpos[NAXIS], wcspos[NAXIS],
		da,dd;
   int		lng,lat;

  wcs = field->wcs;
  lng = wcs->lng;
  lat = wcs->lat;

/* If working with WCS, compute FOCAL coordinates and local matrix */
  if (FLAG(obj2.mxf))
    {
    rawpos[0] = obj2->posx;
    rawpos[1] = obj2->posy;
    raw_to_red(wcs, rawpos, wcspos);
    obj2->mxf = wcspos[0];
    obj2->myf = wcspos[1];
    }

/* If working with WCS, compute WORLD coordinates and local matrix */
  if (FLAG(obj2.mxw))
    {
    rawpos[0] = obj2->posx;
    rawpos[1] = obj2->posy;
    raw_to_wcs(wcs, rawpos, wcspos);
    obj2->mxw = wcspos[0];
    obj2->myw = wcspos[1];
    if (lng != lat)
      {
      obj2->alphas = lng<lat? obj2->mxw : obj2->myw;
      obj2->deltas = lng<lat? obj2->myw : obj2->mxw;
      if (FLAG(obj2.alpha2000))
        {
        if (fabs(wcs->equinox-2000.0)>0.003)
          precess(wcs->equinox, wcspos[lng<lat?0:1], wcspos[lng<lat?1:0],
		2000.0, &obj2->alpha2000, &obj2->delta2000);
        else
          {
          obj2->alpha2000 = lng<lat? obj2->mxw : obj2->myw;
          obj2->delta2000 = lng<lat? obj2->myw : obj2->mxw;
          }
        if (FLAG(obj2.dtheta2000))
          {
          da = wcs->ap2000 - obj2->alpha2000;
          dd = (sin(wcs->dp2000*DEG)
		-sin(obj2->delta2000*DEG)*sin(obj2->deltas*DEG))
		/(cos(obj2->delta2000*DEG)*cos(obj2->deltas*DEG));
          dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
          obj2->dtheta2000 = (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
          }
        if (FLAG(obj2.alpha1950))
          {
          j2b(wcs->equinox, obj2->alpha2000, obj2->delta2000,
		&obj2->alpha1950, &obj2->delta1950);
          if (FLAG(obj2.dtheta1950))
            {
            da = wcs->ap1950 - obj2->alpha1950;
            dd = (sin(wcs->dp1950*DEG)
		-sin(obj2->delta1950*DEG)*sin(obj2->deltas*DEG))
		/(cos(obj2->delta1950*DEG)*cos(obj2->deltas*DEG));
            dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
            obj2->dtheta1950 = (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
           }
          }
        }
      }
    }

/* Custom coordinate system for the MAMA machine */
  if (FLAG(obj2.mamaposx))
    {
    rawpos[0] = obj2->posx - 0.5;
    rawpos[1] = obj2->posy - 0.5;
    raw_to_wcs(wcs, rawpos, wcspos);
    obj2->mamaposx = wcspos[1]*(MAMA_CORFLEX+1.0);
    obj2->mamaposy = wcspos[0]*(MAMA_CORFLEX+1.0);
    }

  return;
  }


/***************************** astrom_peakpos *******************************/
/*
Compute real FOCAL and WORLD peak coordinates according to FITS info.
*/
void	astrom_peakpos(picstruct *field, objstruct *obj)

  {
   wcsstruct	*wcs;
   double	rawpos[NAXIS], wcspos[NAXIS];
   int		lng,lat;

  wcs = field->wcs;
  lng = wcs->lng;
  lat = wcs->lat;

  if (FLAG(obj2.peakxf))
    {
    rawpos[0] = obj->peakx;
    rawpos[1] = obj->peaky;
    raw_to_red(wcs, rawpos, wcspos);
    obj2->peakxf = wcspos[0];
    obj2->peakyf = wcspos[1];
    }
  if (FLAG(obj2.peakxw))
    {
    rawpos[0] = obj->peakx;
    rawpos[1] = obj->peaky;
    raw_to_wcs(wcs, rawpos, wcspos);
    obj2->peakxw = wcspos[0];
    obj2->peakyw = wcspos[1];
    if (lng != lat)
      {
      obj2->peakalphas = lng<lat? obj2->peakxw : obj2->peakyw;
      obj2->peakdeltas = lng<lat? obj2->peakyw : obj2->peakxw;
      if (FLAG(obj2.peakalpha2000))
        {
        if (fabs(wcs->equinox-2000.0)>0.003)
          precess(wcs->equinox, wcspos[lng<lat?0:1], wcspos[lng<lat?1:0],
		2000.0, &obj2->peakalpha2000, &obj2->peakdelta2000);
        else
          {
          obj2->peakalpha2000 = lng<lat? obj2->peakxw : obj2->peakyw;
          obj2->peakdelta2000 = lng<lat? obj2->peakyw : obj2->peakxw;
          }
        if (FLAG(obj2.peakalpha1950))
          j2b(wcs->equinox, obj2->peakalpha2000, obj2->peakdelta2000,
		&obj2->peakalpha1950, &obj2->peakdelta1950);
        }
      }
    }

  return;
  }


/****************************** astrom_winpos *******************************/
/*
Compute real FOCAL and WORLD windowed coordinates according to FITS info.
*/
void	astrom_winpos(picstruct *field, objstruct *obj)

  {
   wcsstruct	*wcs;
   double	rawpos[NAXIS], wcspos[NAXIS];
   int		lng,lat;

  wcs = field->wcs;
  lng = wcs->lng;
  lat = wcs->lat;

  if (FLAG(obj2.winpos_xf))
    {
    rawpos[0] = obj2->winpos_x;
    rawpos[1] = obj2->winpos_y;
    raw_to_red(wcs, rawpos, wcspos);
    obj2->winpos_xf = wcspos[0];
    obj2->winpos_yf = wcspos[1];
    }
  if (FLAG(obj2.winpos_xw))
    {
    rawpos[0] = obj2->winpos_x;
    rawpos[1] = obj2->winpos_y;
    raw_to_wcs(wcs, rawpos, wcspos);
    obj2->winpos_xw = wcspos[0];
    obj2->winpos_yw = wcspos[1];
    if (lng != lat)
      {
      obj2->winpos_alphas = lng<lat? obj2->winpos_xw : obj2->winpos_yw;
      obj2->winpos_deltas = lng<lat? obj2->winpos_yw : obj2->winpos_xw;
      if (FLAG(obj2.winpos_alpha2000))
        {
        if (fabs(wcs->equinox-2000.0)>0.003)
          precess(wcs->equinox, wcspos[0], wcspos[1],
		2000.0, &obj2->winpos_alpha2000, &obj2->winpos_delta2000);
        else
          {
          obj2->winpos_alpha2000 = lng<lat? obj2->winpos_xw : obj2->winpos_yw;
          obj2->winpos_delta2000 = lng<lat? obj2->winpos_yw : obj2->winpos_xw;
          }
        if (FLAG(obj2.winpos_alpha1950))
          j2b(wcs->equinox, obj2->winpos_alpha2000, obj2->winpos_delta2000,
		&obj2->winpos_alpha1950, &obj2->winpos_delta1950);
        }
      }
    }

  return;
  }


/****************************** astrom_psfpos *******************************/
/*
Compute real FOCAL and WORLD PSF coordinates according to FITS info.
*/
void	astrom_psfpos(picstruct *field, objstruct *obj)

  {
   wcsstruct	*wcs;
   double	rawpos[NAXIS], wcspos[NAXIS];
   int		lng,lat;

  wcs = field->wcs;
  lng = wcs->lng;
  lat = wcs->lat;

  if (FLAG(obj2.xf_psf))
    {
    rawpos[0] = obj2->x_psf;
    rawpos[1] = obj2->y_psf;
    raw_to_red(wcs, rawpos, wcspos);
    obj2->xf_psf = wcspos[0];
    obj2->yf_psf = wcspos[1];
    }
  if (FLAG(obj2.xw_psf))
    {
    rawpos[0] = obj2->x_psf;
    rawpos[1] = obj2->y_psf;
    raw_to_wcs(wcs, rawpos, wcspos);
    obj2->xw_psf = wcspos[0];
    obj2->yw_psf = wcspos[1];
    if (lng != lat)
      {
      obj2->alphas_psf = lng<lat? obj2->xw_psf : obj2->yw_psf;
      obj2->deltas_psf = lng<lat? obj2->yw_psf : obj2->xw_psf;
      if (FLAG(obj2.alpha2000_psf))
        {
        if (fabs(wcs->equinox-2000.0)>0.003)
          precess(wcs->equinox, wcspos[0], wcspos[1],
		2000.0, &obj2->alpha2000_psf, &obj2->delta2000_psf);
        else
          {
          obj2->alpha2000_psf = lng<lat? obj2->xw_psf : obj2->yw_psf;
          obj2->delta2000_psf = lng<lat? obj2->yw_psf : obj2->xw_psf;
          }
        if (FLAG(obj2.alpha1950_psf))
          j2b(wcs->equinox, obj2->alpha2000_psf, obj2->delta2000_psf,
		&obj2->alpha1950_psf, &obj2->delta1950_psf);
        }
      }
    }

  return;
  }


/****************************** astrom_profpos *******************************/
/*
Compute real FOCAL and WORLD profit coordinates according to FITS info.
*/
void	astrom_profpos(picstruct *field, objstruct *obj)

  {
   wcsstruct	*wcs;
   double	rawpos[NAXIS], wcspos[NAXIS];
   int		lng,lat;

  wcs = field->wcs;
  lng = wcs->lng;
  lat = wcs->lat;

  if (FLAG(obj2.xf_prof))
    {
    rawpos[0] = obj2->x_prof;
    rawpos[1] = obj2->y_prof;
    raw_to_red(wcs, rawpos, wcspos);
    obj2->xf_prof = wcspos[0];
    obj2->yf_prof = wcspos[1];
    }
  if (FLAG(obj2.xw_prof))
    {
    rawpos[0] = obj2->x_prof;
    rawpos[1] = obj2->y_prof;
    raw_to_wcs(wcs, rawpos, wcspos);
    obj2->xw_prof = wcspos[0];
    obj2->yw_prof = wcspos[1];
    if (lng != lat)
      {
      obj2->alphas_prof = lng<lat? obj2->xw_prof : obj2->yw_prof;
      obj2->deltas_prof = lng<lat? obj2->yw_prof : obj2->xw_prof;
      if (FLAG(obj2.alpha2000_prof))
        {
        if (fabs(wcs->equinox-2000.0)>0.003)
          precess(wcs->equinox, wcspos[0], wcspos[1],
		2000.0, &obj2->alpha2000_prof, &obj2->delta2000_prof);
        else
          {
          obj2->alpha2000_prof = lng<lat? obj2->xw_prof : obj2->yw_prof;
          obj2->delta2000_prof = lng<lat? obj2->yw_prof : obj2->xw_prof;
          }
        if (FLAG(obj2.alpha1950_prof))
          j2b(wcs->equinox, obj2->alpha2000_prof, obj2->delta2000_prof,
		&obj2->alpha1950_prof, &obj2->delta1950_prof);
        }
      }
    }

  return;
  }


/****************************** astrom_shapeparam ****************************/
/*
Compute shape parameters in WORLD and SKY coordinates.
*/
void	astrom_shapeparam(picstruct *field, objstruct *obj)
  {
   wcsstruct	*wcs;
   double	dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  wcs = field->wcs;
  naxis = wcs->naxis;
  lng = wcs->lng;
  lat = wcs->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = obj2->jacob[lng+naxis*lng];
  lm1 = obj2->jacob[lat+naxis*lng];
  lm2 = obj2->jacob[lng+naxis*lat];
  lm3 = obj2->jacob[lat+naxis*lat];


/* All WORLD params based on 2nd order moments have to pass through here */
  dx2 = obj->mx2;
  dy2 = obj->my2;
  dxy = obj->mxy;
  obj2->mx2w = xm2 = lm0*lm0*dx2 + lm1*lm1*dy2 + lm0*lm1*dxy;
  obj2->my2w = ym2 = lm2*lm2*dx2 + lm3*lm3*dy2 + lm2*lm3*dxy;
  obj2->mxyw = xym = lm0*lm2*dx2 + lm1*lm3*dy2 + (lm0*lm3+lm1*lm2)*dxy;
  temp=xm2-ym2;
  if (FLAG(obj2.thetaw))
    {
    obj2->thetaw = fmod_m90_p90((temp == 0.0)?
				(45.0) : (0.5*atan2(2.0 * xym,temp)/DEG));

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
      if (FLAG(obj2.thetas))
        obj2->thetas = fmod_m90_p90(lng<lat?
				((obj2->thetaw>0.0?90:-90.0) - obj2->thetaw)
				: obj2->thetaw);
      if (FLAG(obj2.theta2000))
        obj2->theta2000 = fmod_m90_p90(obj2->thetas + obj2->dtheta2000);
      if (FLAG(obj2.theta1950))
        obj2->theta1950 = fmod_m90_p90(obj2->thetas + obj2->dtheta1950);
      }
    }

  if (FLAG(obj2.aw))
    {
    temp = sqrt(0.25*temp*temp+xym*xym);
    pm2 = 0.5*(xm2+ym2);
    obj2->aw = (float)sqrt(pm2+temp);
    obj2->bw = (float)sqrt(pm2-temp);
    obj2->polarw = temp / pm2;
    }

  if (FLAG(obj2.cxxw))
    {
/*-- Handle large, fully correlated profiles (can cause a singularity...) */
    if ((temp=xm2*ym2-xym*xym)<1e-6)
      {
      temp = 1e-6;
      xym *= 0.99999;
      }
    obj2->cxxw = (float)(ym2/temp);
    obj2->cyyw = (float)(xm2/temp);
    obj2->cxyw = (float)(-2*xym/temp);
    }

  return;
  }


/**************************** astrom_winshapeparam ***************************/
/*
Compute shape parameters in WORLD and SKY coordinates.
*/
void	astrom_winshapeparam(picstruct *field, objstruct *obj)
  {
   wcsstruct	*wcs;
   double	dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  wcs = field->wcs;
  naxis = wcs->naxis;
  lng = wcs->lng;
  lat = wcs->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = obj2->jacob[lng+naxis*lng];
  lm1 = obj2->jacob[lat+naxis*lng];
  lm2 = obj2->jacob[lng+naxis*lat];
  lm3 = obj2->jacob[lat+naxis*lat];

/* All WORLD params based on 2nd order moments have to pass through here */
  dx2 = obj2->win_mx2;
  dy2 = obj2->win_my2;
  dxy = obj2->win_mxy;
  obj2->win_mx2w = xm2 = lm0*lm0*dx2 + lm1*lm1*dy2 + lm0*lm1*dxy;
  obj2->win_my2w = ym2 = lm2*lm2*dx2 + lm3*lm3*dy2 + lm2*lm3*dxy;
  obj2->win_mxyw = xym = lm0*lm2*dx2 + lm1*lm3*dy2 + (lm0*lm3+lm1*lm2)*dxy;
  temp=xm2-ym2;
  if (FLAG(obj2.win_thetaw))
    {
    obj2->win_thetaw = fmod_m90_p90((temp == 0.0)?
			(45.0) : (0.5*atan2(2.0*xym,temp)/DEG));

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
      if (FLAG(obj2.win_thetas))
        obj2->win_thetas = fmod_m90_p90(lng<lat?
			((obj2->win_thetaw>0.0?90:-90.0) - obj2->win_thetaw)
			: obj2->win_thetaw);
      if (FLAG(obj2.win_theta2000))
        obj2->win_theta2000 = fmod_m90_p90(obj2->win_thetas + obj2->dtheta2000);
      if (FLAG(obj2.win_theta1950))
        obj2->win_theta1950 = fmod_m90_p90(obj2->win_thetas + obj2->dtheta1950);
      }
    }

  if (FLAG(obj2.win_aw))
    {
    temp = sqrt(0.25*temp*temp+xym*xym);
    pm2 = 0.5*(xm2+ym2);
    obj2->win_aw = (float)sqrt(pm2+temp);
    obj2->win_bw = (float)sqrt(pm2-temp);
    obj2->win_polarw = temp / pm2;
    }

  if (FLAG(obj2.win_cxxw))
    {
/*-- Handle large, fully correlated profiles (can cause a singularity...) */
    if ((temp=xm2*ym2-xym*xym)<1e-6)
      {
      temp = 1e-6;
      xym *= 0.99999;
      }
    obj2->win_cxxw = (float)(ym2/temp);
    obj2->win_cyyw = (float)(xm2/temp);
    obj2->win_cxyw = (float)(-2*xym/temp);
    }

  return;
  }


/******************************* astrom_errparam *****************************/
/*
Compute error ellipse parameters in WORLD and SKY coordinates.
*/
void	astrom_errparam(picstruct *field, objstruct *obj)
  {
   wcsstruct	*wcs;
   double	dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  wcs = field->wcs;
  naxis = wcs->naxis;
  lng = wcs->lng;
  lat = wcs->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = obj2->jacob[lng+naxis*lng];
  lm1 = obj2->jacob[lat+naxis*lng];
  lm2 = obj2->jacob[lng+naxis*lat];
  lm3 = obj2->jacob[lat+naxis*lat];

/* All WORLD params based on 2nd order moments have to pass through here */
  dx2 = obj->poserr_mx2;
  dy2 = obj->poserr_my2;
  dxy = obj->poserr_mxy;
  obj2->poserr_mx2w = xm2 = lm0*lm0*dx2 + lm1*lm1*dy2 + lm0*lm1*dxy;
  obj2->poserr_my2w = ym2 = lm2*lm2*dx2 + lm3*lm3*dy2 + lm2*lm3*dxy;
  obj2->poserr_mxyw = xym = lm0*lm2*dx2 + lm1*lm3*dy2 + (lm0*lm3+lm1*lm2)*dxy;
  temp=xm2-ym2;
  if (FLAG(obj2.poserr_thetaw))
    {
    obj2->poserr_thetaw = fmod_m90_p90((temp==0.0)?
				(45.0):(0.5*atan2(2.0*xym,temp)/DEG));

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
      if (FLAG(obj2.poserr_thetas))
        obj2->poserr_thetas = fmod_m90_p90(lng<lat?
		((obj2->poserr_thetaw>0.0?90:-90.0) - obj2->poserr_thetaw)
		: obj2->poserr_thetaw);
      if (FLAG(obj2.poserr_theta2000))
        obj2->poserr_theta2000 = fmod_m90_p90(obj2->poserr_thetas
				+ obj2->dtheta2000);
      if (FLAG(obj2.poserr_theta1950))
        obj2->poserr_theta1950 = fmod_m90_p90(obj2->poserr_thetas
				+ obj2->dtheta1950);
      }
    }

  if (FLAG(obj2.poserr_aw))
    {
    temp = sqrt(0.25*temp*temp+xym*xym);
    pm2 = 0.5*(xm2+ym2);
    obj2->poserr_aw = (float)sqrt(pm2+temp);
    obj2->poserr_bw = (float)sqrt(pm2-temp);
    }

  if (FLAG(obj2.poserr_cxxw))
    {
/*-- Handle large, fully correlated profiles (can cause a singularity...) */
    if ((temp=xm2*ym2-xym*xym)<1e-6)
      {
      temp = 1e-6;
      xym *= 0.99999;
      }
    obj2->poserr_cxxw = (float)(ym2/temp);
    obj2->poserr_cyyw = (float)(xm2/temp);
    obj2->poserr_cxyw = (float)(-2*xym/temp);
    }

  return;
  }


/***************************** astrom_winerrparam ***************************/
/*
Compute error ellipse parameters in WORLD and SKY coordinates.
*/
void	astrom_winerrparam(picstruct *field, objstruct *obj)
  {
   wcsstruct	*wcs;
   double	dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  wcs = field->wcs;
  naxis = wcs->naxis;
  lng = wcs->lng;
  lat = wcs->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = obj2->jacob[lng+naxis*lng];
  lm1 = obj2->jacob[lat+naxis*lng];
  lm2 = obj2->jacob[lng+naxis*lat];
  lm3 = obj2->jacob[lat+naxis*lat];

/* All WORLD params based on 2nd order moments have to pass through here */
  dx2 = obj2->winposerr_mx2;
  dy2 = obj2->winposerr_my2;
  dxy = obj2->winposerr_mxy;
  obj2->winposerr_mx2w = xm2 = lm0*lm0*dx2 + lm1*lm1*dy2 + lm0*lm1*dxy;
  obj2->winposerr_my2w = ym2 = lm2*lm2*dx2 + lm3*lm3*dy2 + lm2*lm3*dxy;
  obj2->winposerr_mxyw = xym = lm0*lm2*dx2 + lm1*lm3*dy2 + (lm0*lm3+lm1*lm2)*dxy;
  temp=xm2-ym2;
  if (FLAG(obj2.winposerr_thetaw))
    {
    obj2->winposerr_thetaw = (fmod_m90_p90(temp==0.0)?
				(45.0):(0.5*atan2(2.0*xym,temp)/DEG));

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
      if (FLAG(obj2.winposerr_thetas))
        obj2->winposerr_thetas = fmod_m90_p90(lng<lat?
		((obj2->winposerr_thetaw>0.0?90:-90.0) - obj2->winposerr_thetaw)
		: obj2->winposerr_thetaw);
      if (FLAG(obj2.winposerr_theta2000))
        obj2->winposerr_theta2000 = fmod_m90_p90(obj2->winposerr_thetas
				+ obj2->dtheta2000);
      if (FLAG(obj2.winposerr_theta1950))
        obj2->winposerr_theta1950 = fmod_m90_p90(obj2->winposerr_thetas
				+ obj2->dtheta1950);
      }
    }

  if (FLAG(obj2.winposerr_aw))
    {
    temp = sqrt(0.25*temp*temp+xym*xym);
    pm2 = 0.5*(xm2+ym2);
    obj2->winposerr_aw = (float)sqrt(pm2+temp);
    obj2->winposerr_bw = (float)sqrt(pm2-temp);
    }

  if (FLAG(obj2.winposerr_cxxw))
    {
/*-- Handle large, fully correlated profiles (can cause a singularity...) */
    if ((temp=xm2*ym2-xym*xym)<1e-6)
      {
      temp = 1e-6;
      xym *= 0.99999;
      }
    obj2->winposerr_cxxw = (float)(ym2/temp);
    obj2->winposerr_cyyw = (float)(xm2/temp);
    obj2->winposerr_cxyw = (float)(-2*xym/temp);
    }

  return;
  }


/***************************** astrom_psferrparam ***************************/
/*
Compute error ellipse parameters in WORLD and SKY coordinates.
*/
void	astrom_psferrparam(picstruct *field, objstruct *obj)
  {
   wcsstruct	*wcs;
   double	dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  wcs = field->wcs;
  naxis = wcs->naxis;
  lng = wcs->lng;
  lat = wcs->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = obj2->jacob[lng+naxis*lng];
  lm1 = obj2->jacob[lat+naxis*lng];
  lm2 = obj2->jacob[lng+naxis*lat];
  lm3 = obj2->jacob[lat+naxis*lat];

/* All WORLD params based on 2nd order moments have to pass through here */
  dx2 = obj2->poserrmx2_psf;
  dy2 = obj2->poserrmy2_psf;
  dxy = obj2->poserrmxy_psf;
  obj2->poserrmx2w_psf = xm2 = lm0*lm0*dx2+lm1*lm1*dy2+lm0*lm1*dxy;
  obj2->poserrmy2w_psf = ym2 = lm2*lm2*dx2+lm3*lm3*dy2+lm2*lm3*dxy;
  obj2->poserrmxyw_psf = xym = lm0*lm2*dx2+lm1*lm3*dy2+(lm0*lm3+lm1*lm2)*dxy;
  temp=xm2-ym2;
  if (FLAG(obj2.poserrthetaw_psf))
    {
    obj2->poserrthetaw_psf = (fmod_m90_p90(temp==0.0)?
				(45.0):(0.5*atan2(2.0*xym,temp)/DEG));

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
      if (FLAG(obj2.poserrthetas_psf))
        obj2->poserrthetas_psf = fmod_m90_p90(lng<lat?
		((obj2->poserrthetaw_psf>0.0?90:-90.0)-obj2->poserrthetaw_psf)
		: obj2->poserrthetaw_psf);
      if (FLAG(obj2.poserrtheta2000_psf))
        obj2->poserrtheta2000_psf = fmod_m90_p90(obj2->poserrthetas_psf
				+ obj2->dtheta2000);
      if (FLAG(obj2.poserrtheta1950_psf))
        obj2->poserrtheta1950_psf = fmod_m90_p90(obj2->poserrthetas_psf
				+ obj2->dtheta1950);
      }
    }

  if (FLAG(obj2.poserraw_psf))
    {
    temp = sqrt(0.25*temp*temp+xym*xym);
    pm2 = 0.5*(xm2+ym2);
    obj2->poserraw_psf = (float)sqrt(pm2+temp);
    obj2->poserrbw_psf = (float)sqrt(pm2-temp);
    }

  if (FLAG(obj2.poserrcxxw_psf))
    {
/*-- Handle large, fully correlated profiles (can cause a singularity...) */
    if ((temp=xm2*ym2-xym*xym)<1e-6)
      {
      temp = 1e-6;
      xym *= 0.99999;
      }
    obj2->poserrcxxw_psf = (float)(ym2/temp);
    obj2->poserrcyyw_psf = (float)(xm2/temp);
    obj2->poserrcxyw_psf = (float)(-2*xym/temp);
    }

  return;
  }


/***************************** astrom_proferrparam ***************************/
/*
Compute error ellipse parameters in WORLD and SKY coordinates.
*/
void	astrom_proferrparam(picstruct *field, objstruct *obj)
  {
   wcsstruct	*wcs;
   double	dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  wcs = field->wcs;
  naxis = wcs->naxis;
  lng = wcs->lng;
  lat = wcs->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = obj2->jacob[lng+naxis*lng];
  lm1 = obj2->jacob[lat+naxis*lng];
  lm2 = obj2->jacob[lng+naxis*lat];
  lm3 = obj2->jacob[lat+naxis*lat];

/* All WORLD params based on 2nd order moments have to pass through here */
  dx2 = obj2->poserrmx2_prof;
  dy2 = obj2->poserrmy2_prof;
  dxy = obj2->poserrmxy_prof;
  obj2->poserrmx2w_prof = xm2 = lm0*lm0*dx2+lm1*lm1*dy2+lm0*lm1*dxy;
  obj2->poserrmy2w_prof = ym2 = lm2*lm2*dx2+lm3*lm3*dy2+lm2*lm3*dxy;
  obj2->poserrmxyw_prof = xym = lm0*lm2*dx2+lm1*lm3*dy2+(lm0*lm3+lm1*lm2)*dxy;
  temp=xm2-ym2;
  if (FLAG(obj2.poserrthetaw_prof))
    {
    obj2->poserrthetaw_prof = (fmod_m90_p90(temp==0.0)?
				(45.0):(0.5*atan2(2.0*xym,temp)/DEG));

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
      if (FLAG(obj2.poserrthetas_prof))
        obj2->poserrthetas_prof = fmod_m90_p90(lng<lat?
		((obj2->poserrthetaw_prof>0.0?90:-90.0)-obj2->poserrthetaw_prof)
		: obj2->poserrthetaw_prof);
      if (FLAG(obj2.poserrtheta2000_prof))
        obj2->poserrtheta2000_prof = fmod_m90_p90(obj2->poserrthetas_prof
				+ obj2->dtheta2000);
      if (FLAG(obj2.poserrtheta1950_prof))
        obj2->poserrtheta1950_prof = fmod_m90_p90(obj2->poserrthetas_prof
				+ obj2->dtheta1950);
      }
    }

  if (FLAG(obj2.poserraw_prof))
    {
    temp = sqrt(0.25*temp*temp+xym*xym);
    pm2 = 0.5*(xm2+ym2);
    obj2->poserraw_prof = (float)sqrt(pm2+temp);
    obj2->poserrbw_prof = (float)sqrt(pm2-temp);
    }

  if (FLAG(obj2.poserrcxxw_prof))
    {
/*-- Handle large, fully correlated profiles (can cause a singularity...) */
    if ((temp=xm2*ym2-xym*xym)<1e-6)
      {
      temp = 1e-6;
      xym *= 0.99999;
      }
    obj2->poserrcxxw_prof = (float)(ym2/temp);
    obj2->poserrcyyw_prof = (float)(xm2/temp);
    obj2->poserrcxyw_prof = (float)(-2*xym/temp);
    }

  return;
  }


/*************************** astrom_profshapeparam ***************************/
/*
Compute profile-fitting shape parameters in WORLD and SKY coordinates.
*/
void	astrom_profshapeparam(picstruct *field, objstruct *obj)
  {
   wcsstruct	*wcs;
   double	mat[9], tempmat[9], mx2wcov[9], dpdmx2[6], cov[4],
		dx2,dy2,dxy, xm2,ym2,xym, pm2, lm0,lm1,lm2,lm3, ct,st,
		temp, invstemp, den, invden, dval;
   int		lng,lat, naxis;

  wcs = field->wcs;
  naxis = wcs->naxis;
  lng = wcs->lng;
  lat = wcs->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = obj2->jacob[lng+naxis*lng];
  lm1 = obj2->jacob[lat+naxis*lng];
  lm2 = obj2->jacob[lng+naxis*lat];
  lm3 = obj2->jacob[lat+naxis*lat];

/* Spheroid World coordinates */
  obj2->prof_spheroid_thetaerrw=obj2->prof_spheroid_thetaerr; /* quick & dirty*/
  if (FLAG(obj2.prof_spheroid_reffw))
    {
    ct = cos(obj2->prof_spheroid_theta*DEG);
    st = sin(obj2->prof_spheroid_theta*DEG);
    dx2 = obj2->prof_spheroid_reff*obj2->prof_spheroid_reff * (ct*ct
	+ st*st * obj2->prof_spheroid_aspect*obj2->prof_spheroid_aspect);
    dy2 = obj2->prof_spheroid_reff*obj2->prof_spheroid_reff * (st*st
	+ ct*ct * obj2->prof_spheroid_aspect*obj2->prof_spheroid_aspect);
    dxy = ct*st * obj2->prof_spheroid_reff*obj2->prof_spheroid_reff
	*(1.0 - obj2->prof_spheroid_aspect*obj2->prof_spheroid_aspect);
    xm2 = lm0*lm0*dx2 + lm1*lm1*dy2 + lm0*lm1*dxy;
    ym2 = lm2*lm2*dx2 + lm3*lm3*dy2 + lm2*lm3*dxy;
    xym = lm0*lm2*dx2 + lm1*lm3*dy2 + (lm0*lm3+lm1*lm2)*dxy;
    temp=xm2-ym2;
    if (FLAG(obj2.prof_spheroid_thetaw))
      {
      obj2->prof_spheroid_thetaw = fmod_m90_p90((temp == 0.0)?
		(45.0) : (0.5*atan2(2.0 * xym,temp)/DEG));

      if (wcs->lng != wcs->lat)
        {
        if (FLAG(obj2.prof_spheroid_thetas))
          obj2->prof_spheroid_thetas = fmod_m90_p90(lng<lat?
			((obj2->prof_spheroid_thetaw>0.0?90:-90.0)
				- obj2->prof_spheroid_thetaw)
			: obj2->prof_spheroid_thetaw);
        if (FLAG(obj2.prof_spheroid_theta2000))
          obj2->prof_spheroid_theta2000=fmod_m90_p90(obj2->prof_spheroid_thetas
					+ obj2->dtheta2000);
        if (FLAG(obj2.prof_spheroid_theta1950))
          obj2->prof_spheroid_theta1950=fmod_m90_p90(obj2->prof_spheroid_thetas
					+ obj2->dtheta1950);
        }
      }
    temp = sqrt(0.25*temp*temp+xym*xym);
    pm2 = 0.5*(xm2+ym2);
    obj2->prof_spheroid_reffw = sqrt(pm2+temp);
    obj2->prof_spheroid_refferrw = obj2->prof_spheroid_reff > 0.0?
	obj2->prof_spheroid_refferr/obj2->prof_spheroid_reff
				*obj2->prof_spheroid_reffw
	: 0.0;
    obj2->prof_spheroid_aspectw = (obj2->prof_spheroid_reffw>0.0 && pm2>temp)?
				  sqrt(pm2-temp) / obj2->prof_spheroid_reffw
				: obj2->prof_spheroid_aspect;
    obj2->prof_spheroid_aspecterrw = obj2->prof_spheroid_aspect > 0.0?
	obj2->prof_spheroid_aspecterr/obj2->prof_spheroid_aspect
				*obj2->prof_spheroid_aspectw
	: 0.0;
    }

/* Disk World coordinates */
  obj2->prof_disk_thetaerrw = obj2->prof_disk_thetaerr; /* quick & dirty*/
  if (FLAG(obj2.prof_disk_scalew))
    {
    ct = cos(obj2->prof_disk_theta*DEG);
    st = sin(obj2->prof_disk_theta*DEG);
    dx2 = obj2->prof_disk_scale*obj2->prof_disk_scale * (ct*ct
	+ st*st * obj2->prof_disk_aspect*obj2->prof_disk_aspect);
    dy2 = obj2->prof_disk_scale*obj2->prof_disk_scale * (st*st
	+ ct*ct * obj2->prof_disk_aspect*obj2->prof_disk_aspect);
    dxy = ct*st * obj2->prof_disk_scale*obj2->prof_disk_scale
	*(1.0 - obj2->prof_disk_aspect*obj2->prof_disk_aspect);
    xm2 = lm0*lm0*dx2 + lm1*lm1*dy2 + lm0*lm1*dxy;
    ym2 = lm2*lm2*dx2 + lm3*lm3*dy2 + lm2*lm3*dxy;
    xym = lm0*lm2*dx2 + lm1*lm3*dy2 + (lm0*lm3+lm1*lm2)*dxy;
    temp=xm2-ym2;
    if (FLAG(obj2.prof_disk_thetaw))
      {
      obj2->prof_disk_thetaw = fmod_m90_p90((temp == 0.0)?
		(45.0) : (0.5*atan2(2.0 * xym,temp)/DEG));

/*---- Compute position angles in J2000 or B1950 reference frame */
      if (wcs->lng != wcs->lat)
        {
        if (FLAG(obj2.prof_disk_thetas))
          obj2->prof_disk_thetas = fmod_m90_p90(lng<lat?
			((obj2->prof_disk_thetaw>0.0?90:-90.0)
				- obj2->prof_disk_thetaw)
			: obj2->prof_disk_thetaw);
        if (FLAG(obj2.prof_disk_theta2000))
          obj2->prof_disk_theta2000 = fmod_m90_p90(obj2->prof_disk_thetas
					+ obj2->dtheta2000);
        if (FLAG(obj2.prof_disk_theta1950))
          obj2->prof_disk_theta1950 = fmod_m90_p90(obj2->prof_disk_thetas
					+ obj2->dtheta1950);
        }
      }
    temp = sqrt(0.25*temp*temp+xym*xym);
    pm2 = 0.5*(xm2+ym2);
    obj2->prof_disk_scalew = sqrt(pm2+temp);
    obj2->prof_disk_scaleerrw = obj2->prof_disk_scale > 0.0?
	obj2->prof_disk_scaleerr/obj2->prof_disk_scale*obj2->prof_disk_scalew
	: 0.0;
    obj2->prof_disk_aspectw = (obj2->prof_disk_scalew>0.0 && pm2>temp)?
				  sqrt(pm2-temp) / obj2->prof_disk_scalew
				: obj2->prof_disk_aspect;
    obj2->prof_disk_aspecterrw = obj2->prof_disk_aspect > 0.0?
	obj2->prof_disk_aspecterr/obj2->prof_disk_aspect*obj2->prof_disk_aspectw
	: 0.0;
/*-- Arms World coordinates */
    if (FLAG(obj2.prof_arms_scalew))
      {
      obj2->prof_arms_scalew = (obj2->prof_disk_scale > 0.0) ?
	  obj2->prof_arms_scale*obj2->prof_disk_scalew/obj2->prof_disk_scale
	: 0.0;
      obj2->prof_arms_scaleerrw = (obj2->prof_arms_scale > 0.0) ?
	obj2->prof_arms_scaleerr/obj2->prof_arms_scale*obj2->prof_arms_scalew
	: 0.0;
      obj2->prof_arms_startw = (obj2->prof_disk_scale > 0.0) ?
	  obj2->prof_arms_start*obj2->prof_disk_scalew/obj2->prof_disk_scale
	: 0.0;
      obj2->prof_arms_starterrw = (obj2->prof_arms_start > 0.0) ?
	obj2->prof_arms_starterr/obj2->prof_arms_start*obj2->prof_arms_startw
	: 0.0;
      }
    }

/* Bar World coordinates */
  obj2->prof_bar_thetaerrw = obj2->prof_bar_thetaerr;
  if (FLAG(obj2.prof_bar_lengthw))
    {
    ct = cos(obj2->prof_bar_theta*DEG);
    st = sin(obj2->prof_bar_theta*DEG);
    dx2 = obj2->prof_bar_length*obj2->prof_bar_length * (ct*ct
	+ st*st * obj2->prof_bar_aspect*obj2->prof_bar_aspect);
    dy2 = obj2->prof_bar_length*obj2->prof_bar_length * (st*st
	+ ct*ct * obj2->prof_bar_aspect*obj2->prof_bar_aspect);
    dxy = ct*st * obj2->prof_bar_length*obj2->prof_bar_length
	*(1.0 - obj2->prof_bar_aspect*obj2->prof_bar_aspect);
    xm2 = lm0*lm0*dx2 + lm1*lm1*dy2 + lm0*lm1*dxy;
    ym2 = lm2*lm2*dx2 + lm3*lm3*dy2 + lm2*lm3*dxy;
    xym = lm0*lm2*dx2 + lm1*lm3*dy2 + (lm0*lm3+lm1*lm2)*dxy;
    temp=xm2-ym2;
    if (FLAG(obj2.prof_bar_thetaw))
      {
      obj2->prof_bar_thetaw = fmod_m90_p90((temp == 0.0)?
		(45.0) : (0.5*atan2(2.0 * xym,temp)/DEG));

/*---- Compute position angles in J2000 or B1950 reference frame */
      if (wcs->lng != wcs->lat)
        {
        if (FLAG(obj2.prof_bar_thetas))
          obj2->prof_bar_thetas = fmod_m90_p90(lng<lat?
			((obj2->prof_bar_thetaw>0.0?90:-90.0)
				- obj2->prof_bar_thetaw)
			: obj2->prof_bar_thetaw);
        if (FLAG(obj2.prof_bar_theta2000))
          obj2->prof_bar_theta2000 = fmod_m90_p90(obj2->prof_bar_thetas
					+ obj2->dtheta2000);
        if (FLAG(obj2.prof_bar_theta1950))
          obj2->prof_bar_theta1950 = fmod_m90_p90(obj2->prof_bar_thetas
					+ obj2->dtheta1950);
        }
      }
    temp = sqrt(0.25*temp*temp+xym*xym);
    pm2 = 0.5*(xm2+ym2);
    obj2->prof_bar_lengthw = sqrt(pm2+temp);
    obj2->prof_bar_lengtherrw = obj2->prof_bar_length > 0.0?
	obj2->prof_bar_lengtherr/obj2->prof_bar_length*obj2->prof_bar_lengthw
	: 0.0;
    obj2->prof_bar_aspectw = obj2->prof_bar_lengthw>0.0?
				  sqrt(pm2-temp) / obj2->prof_bar_lengthw
				: obj2->prof_bar_aspect;
    obj2->prof_bar_aspecterrw = obj2->prof_bar_aspect > 0.0?
	obj2->prof_bar_aspecterr/obj2->prof_bar_aspect*obj2->prof_bar_aspectw
	: 0.0;
    }

/* Global 2nd-order moments */
  if (FLAG(obj2.prof_mx2w))
    {
    dx2 = obj2->prof_mx2;
    dy2 = obj2->prof_my2;
    dxy = obj2->prof_mxy;
    mat[0] = lm0*lm0;
    mat[3] = lm1*lm1;
    mat[6] = lm0*lm1;
    mat[1] = lm2*lm2;
    mat[4] = lm3*lm3;
    mat[7] = lm2*lm3;
    mat[2] = lm0*lm2;
    mat[5] = lm1*lm3;
    mat[8] = lm0*lm3+lm1*lm2;
    obj2->prof_mx2w = xm2 = mat[0]*dx2 + mat[3]*dy2 + mat[6]*dxy;
    obj2->prof_my2w = ym2 = mat[1]*dx2 + mat[4]*dy2 + mat[7]*dxy;
    obj2->prof_mxyw = xym = mat[2]*dx2 + mat[5]*dy2 + mat[8]*dxy;
    temp=xm2-ym2;
    if (FLAG(obj2.prof_thetaw))
      {
      obj2->prof_thetaw = fmod_m90_p90((temp == 0.0)?
			(45.0) : (0.5*atan2(2.0*xym,temp)/DEG));

/*---- Compute position angles in J2000 or B1950 reference frame */
      if (wcs->lng != wcs->lat)
        {
        if (FLAG(obj2.prof_thetas))
          obj2->prof_thetas = fmod_m90_p90(lng<lat?
			((obj2->prof_thetaw>0.0?90:-90.0) - obj2->prof_thetaw)
			: obj2->prof_thetaw);
        if (FLAG(obj2.prof_theta2000))
          obj2->prof_theta2000 = fmod_m90_p90(obj2->prof_thetas
		+ obj2->dtheta2000);
        if (FLAG(obj2.prof_theta1950))
          obj2->prof_theta1950 = fmod_m90_p90(obj2->prof_thetas
		+ obj2->dtheta1950);
        }
      }

    if (FLAG(obj2.prof_aw))
      {
      temp = sqrt(0.25*temp*temp+xym*xym);
      pm2 = 0.5*(xm2+ym2);
      obj2->prof_aw = (float)sqrt(pm2+temp);
      obj2->prof_bw = (float)sqrt(pm2-temp);
      }

    if (FLAG(obj2.prof_cxxw))
      {
/*---- Handle large, fully correlated profiles (can cause a singularity...) */
      if ((temp=xm2*ym2-xym*xym)<1e-6)
        {
        temp = 1e-6;
        xym *= 0.99999;
        }
      obj2->prof_cxxw = (float)(ym2/temp);
      obj2->prof_cyyw = (float)(xm2/temp);
      obj2->prof_cxyw = (float)(-2*xym/temp);
      }
  
/*-- Use the Jacobians to compute the moment covariance matrix */
    if (FLAG(obj2.prof_pol1errw) || FLAG(obj2.prof_e1errw))
      propagate_covar(obj2->prof_mx2cov, mat, mx2wcov, 3, 3, tempmat);

    if (FLAG(obj2.prof_pol1w))
      {
      if (xm2+ym2 > 1.0/BIG)
        {
        obj2->prof_pol1w = (xm2 - ym2) / (xm2+ym2);
        obj2->prof_pol2w = 2.0*xym / (xm2 + ym2);
        if (FLAG(obj2.prof_pol1errw))
          {
/*-------- Compute the Jacobian of polarisation */
          invden = 1.0/(xm2+ym2);
          dpdmx2[0] =  2.0*ym2*invden*invden;
          dpdmx2[1] = -2.0*xm2*invden*invden;
          dpdmx2[2] =  0.0;
          dpdmx2[3] = -2.0*xym*invden*invden;
          dpdmx2[4] = -2.0*xym*invden*invden;
          dpdmx2[5] =  2.0*invden;
          propagate_covar(mx2wcov, dpdmx2, cov, 3, 2, tempmat);
          obj2->prof_pol1errw = (float)sqrt(cov[0]<0.0? 0.0: cov[0]);
          obj2->prof_pol2errw = (float)sqrt(cov[3]<0.0? 0.0: cov[3]);
          obj2->prof_pol12corrw = (dval=cov[0]*cov[3]) > 0.0?
					(float)(cov[1]/sqrt(dval)) : 0.0;
          }
        }
      else
        obj2->prof_pol1w = obj2->prof_pol2w = obj2->prof_pol1errw
		= obj2->prof_pol2errw = obj2->prof_pol12corrw = 0.0;
      }

    if (FLAG(obj2.prof_e1w))
      {
      if (xm2+ym2 > 1.0/BIG)
        {
        temp = xm2*ym2 - xym*xym;
        den = (temp>=0.0) ? xm2+ym2+2.0*sqrt(temp) : xm2+ym2;
        invden = 1.0/den;
        obj2->prof_e1w = (float)(invden*(xm2 - ym2));
        obj2->prof_e2w = (float)(2.0 * invden * xym);
        if (FLAG(obj2.prof_e1errw))
        {
/*------ Compute the Jacobian of ellipticity */
        invstemp = (temp>=0.0) ? 1.0/sqrt(temp) : 0.0;
        dpdmx2[0] = ( den - (1.0+ym2*invstemp)*(xm2-ym2))*invden*invden;
        dpdmx2[1] = (-den - (1.0+xm2*invstemp)*(xm2-ym2))*invden*invden;
        dpdmx2[2] = 2.0*xym*invstemp*(xm2-ym2)*invden*invden;
        dpdmx2[3] = -2.0*xym*(1.0+ym2*invstemp)*invden*invden;
        dpdmx2[4] = -2.0*xym*(1.0+xm2*invstemp)*invden*invden;
        dpdmx2[5] =  (2.0*den+4.0*xym*xym*invstemp)*invden*invden;

/*------ Use the Jacobian to compute the ellipticity covariance matrix */
        propagate_covar(mx2wcov, dpdmx2, cov, 3, 2, tempmat);
        obj2->prof_e1errw = (float)sqrt(cov[0]<0.0? 0.0: cov[0]);
        obj2->prof_e2errw = (float)sqrt(cov[3]<0.0? 0.0: cov[3]);
        obj2->prof_e12corrw = (dval=cov[0]*cov[3]) > 0.0?
					(float)(cov[1]/sqrt(dval)) : 0.0;
        }
        }
      else
        obj2->prof_e1w = obj2->prof_e2w = obj2->prof_e1errw
		= obj2->prof_e2errw = obj2->prof_e12corrw = 0.0;
      }
    }

  return;
  }

