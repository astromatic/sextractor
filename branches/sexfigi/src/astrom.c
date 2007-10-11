 /*
 				astrom.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Astrometrical computations.
*
*	Last modify:	11/10/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
    if (FLAG(obj2.theta2000) || FLAG(obj2.theta1950)
	|| FLAG(obj2.poserr_theta2000) || FLAG(obj2.poserr_theta1950)
	|| FLAG(obj2.win_theta2000) || FLAG(obj2.win_theta1950)
	|| FLAG(obj2.winposerr_theta2000) || FLAG(obj2.winposerr_theta1950))
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


/**************************** computeastrom *********************************/
/*
Compute real WORLD coordinates and dimensions according to FITS info.
*/
void	computeastrom(picstruct *field, objstruct *obj)

  {
  wcsstruct	*wcs;
  double	rawpos[NAXIS], wcspos[NAXIS],
		pixscale2;

  wcs = field->wcs;

/* If working with WCS, compute WORLD coordinates and local matrix */
  if (FLAG(obj2.mxw))
    {
    rawpos[0] = obj2->posx;
    rawpos[1] = obj2->posy;
    raw_to_wcs(wcs, rawpos, wcspos);
    obj2->mxw = wcspos[0];
    obj2->myw = wcspos[1];
    if (wcs->lng != wcs->lat)
      {
      obj2->alphas = obj2->mxw;
      obj2->deltas = obj2->myw;
      if (FLAG(obj2.alpha2000))
        {
        if (fabs(wcs->equinox-2000.0)>0.003)
          precess(wcs->equinox, wcspos[0], wcspos[1],
		2000.0, &obj2->alpha2000, &obj2->delta2000);
        else
          {
          obj2->alpha2000 = obj2->mxw;
          obj2->delta2000 = obj2->myw;
          }
        if (FLAG(obj2.alpha1950))
          j2b(wcs->equinox, obj2->alpha2000, obj2->delta2000,
		&obj2->alpha1950, &obj2->delta1950);
        }
      }
    }

/* Idem for peak-flux positions */
  if (FLAG(obj2.peakxw))
    {
    rawpos[0] = obj->peakx;
    rawpos[1] = obj->peaky;
    raw_to_wcs(wcs, rawpos, wcspos);
    obj2->peakxw = wcspos[0];
    obj2->peakyw = wcspos[1];
    if (wcs->lng != wcs->lat)
      {
      obj2->peakalphas = obj2->peakxw;
      obj2->peakdeltas = obj2->peakyw;
      if (FLAG(obj2.peakalpha2000))
        {
        if (fabs(wcs->equinox-2000.0)>0.003)
          precess(wcs->equinox, wcspos[0], wcspos[1],
		2000.0, &obj2->peakalpha2000, &obj2->peakdelta2000);
        else
          {
          obj2->peakalpha2000 = obj2->peakxw;
          obj2->peakdelta2000 = obj2->peakyw;
          }
        if (FLAG(obj2.peakalpha1950))
          j2b(wcs->equinox, obj2->peakalpha2000, obj2->peakdelta2000,
		&obj2->peakalpha1950, &obj2->peakdelta1950);
        }
      }
    }

/* Idem for Windowed positions */
  if (FLAG(obj2.winpos_xw))
    {
    rawpos[0] = obj2->winpos_x;
    rawpos[1] = obj2->winpos_y;
    raw_to_wcs(wcs, rawpos, wcspos);
    obj2->winpos_xw = wcspos[0];
    obj2->winpos_yw = wcspos[1];
    if (wcs->lng != wcs->lat)
      {
      obj2->winpos_alphas = obj2->winpos_xw;
      obj2->winpos_deltas = obj2->winpos_yw;
      if (FLAG(obj2.winpos_alpha2000))
        {
        if (fabs(wcs->equinox-2000.0)>0.003)
          precess(wcs->equinox, wcspos[0], wcspos[1],
		2000.0, &obj2->winpos_alpha2000, &obj2->winpos_delta2000);
        else
          {
          obj2->winpos_alpha2000 = obj2->winpos_xw;
          obj2->winpos_delta2000 = obj2->winpos_yw;
          }
        if (FLAG(obj2.winpos_alpha1950))
          j2b(wcs->equinox, obj2->winpos_alpha2000, obj2->winpos_delta2000,
		&obj2->winpos_alpha1950, &obj2->winpos_delta1950);
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

  if (FLAG(obj2.mx2w)
	|| FLAG(obj2.win_mx2w)
	|| FLAG(obj2.poserr_mx2w)
	|| FLAG(obj2.winposerr_mx2w)
	|| ((!prefs.pixel_scale) && (FLAG(obj2.npixw)
		|| FLAG(obj2.fdnpixw)
		|| FLAG(obj2.fwhmw))))
    {
    rawpos[0] = obj2->posx;
    rawpos[1] = obj2->posy;
    pixscale2 = wcs_jacobian(wcs, rawpos, obj2->jacob);
    }

/* Express shape parameters in WORLD frame */
  if (FLAG(obj2.mx2w))
    astrom_shapeparam(field, obj);
  if (FLAG(obj2.win_mx2w))
    astrom_winshapeparam(field, obj);

/* Express position error parameters in WORLD frame */
  if (FLAG(obj2.poserr_mx2w))
    astrom_errparam(field, obj);
  if (FLAG(obj2.winposerr_mx2w))
    astrom_winerrparam(field, obj);

  if (FLAG(obj2.npixw))
    obj2->npixw = obj->npix * (prefs.pixel_scale?
	field->pixscale/3600.0*field->pixscale/3600.0 : pixscale2);
  if (FLAG(obj2.fdnpixw))
    obj2->fdnpixw = obj->fdnpix * (prefs.pixel_scale?
	field->pixscale/3600.0*field->pixscale/3600.0 :pixscale2);

  if (FLAG(obj2.fwhmw))
    obj2->fwhmw = obj->fwhm * (prefs.pixel_scale?
	field->pixscale/3600.0 : sqrt(pixscale2));

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
    obj2->thetaw = (temp == 0.0)? (45.0) : (0.5*atan2(2.0 * xym,temp)/DEG);

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
       double	da,dd;

      if (FLAG(obj2.thetas))
        obj2->thetas = obj2->thetaw + (obj2->thetaw>0.0?-90:90.0);
      if (FLAG(obj2.theta2000))
        {
        da = wcs->ap2000 - obj2->alpha2000;
        dd = (sin(wcs->dp2000*DEG)
		-sin(obj2->delta2000*DEG)*sin(obj2->deltas*DEG))
		/(cos(obj2->delta2000*DEG)*cos(obj2->deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->theta2000 = obj2->thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }

      if (FLAG(obj2.theta1950))
        {
        da = wcs->ap1950 - obj2->alpha1950;
        dd = (sin(wcs->dp1950*DEG)
		-sin(obj2->delta1950*DEG)*sin(obj2->deltas*DEG))
		/(cos(obj2->delta1950*DEG)*cos(obj2->deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->theta1950 = obj2->thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }
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
    obj2->win_thetaw = (temp == 0.0)? (45.0) : (0.5*atan2(2.0*xym,temp)/DEG);

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
       double	da,dd;

      if (FLAG(obj2.win_thetas))
        obj2->win_thetas = obj2->win_thetaw + (obj2->win_thetaw>0.0?-90:90.0);
      if (FLAG(obj2.win_theta2000))
        {
        da = wcs->ap2000 - obj2->winpos_alpha2000;
        dd = (sin(wcs->dp2000*DEG)
		-sin(obj2->winpos_delta2000*DEG)*sin(obj2->winpos_deltas*DEG))
		/(cos(obj2->winpos_delta2000*DEG)*cos(obj2->winpos_deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->win_theta2000 = obj2->win_thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }

      if (FLAG(obj2.win_theta1950))
        {
        da = wcs->ap1950 - obj2->winpos_alpha1950;
        dd = (sin(wcs->dp1950*DEG)
		-sin(obj2->winpos_delta1950*DEG)*sin(obj2->winpos_deltas*DEG))
		/(cos(obj2->winpos_delta1950*DEG)*cos(obj2->winpos_deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->win_theta1950 = obj2->win_thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }
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
    obj2->poserr_thetaw = (temp==0.0)? (45.0):(0.5*atan2(2.0*xym,temp)/DEG);

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
       double	da,dd;

      if (FLAG(obj2.poserr_thetas))
        obj2->poserr_thetas = obj2->poserr_thetaw
				+ (obj2->poserr_thetaw>0.0? -90:90.0);
      if (FLAG(obj2.poserr_theta2000))
        {
        da = wcs->ap2000 - obj2->alpha2000;
        dd = (sin(wcs->dp2000*DEG)
		-sin(obj2->delta2000*DEG)*sin(obj2->deltas*DEG))
		/(cos(obj2->delta2000*DEG)*cos(obj2->deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->poserr_theta2000 = obj2->poserr_thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }

      if (FLAG(obj2.poserr_theta1950))
        {
        da = wcs->ap1950 - obj2->alpha1950;
        dd = (sin(wcs->dp1950*DEG)
		-sin(obj2->delta1950*DEG)*sin(obj2->deltas*DEG))
		/(cos(obj2->delta1950*DEG)*cos(obj2->deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->poserr_theta1950 = obj2->poserr_thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }
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
    obj2->winposerr_thetaw = (temp==0.0)? (45.0):(0.5*atan2(2.0*xym,temp)/DEG);

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (wcs->lng != wcs->lat)
      {
       double	da,dd;

      if (FLAG(obj2.winposerr_thetas))
        obj2->winposerr_thetas = obj2->winposerr_thetaw
				+ (obj2->winposerr_thetaw>0.0? -90:90.0);
      if (FLAG(obj2.winposerr_theta2000))
        {
        da = wcs->ap2000 - obj2->winpos_alpha2000;
        dd = (sin(wcs->dp2000*DEG)
		-sin(obj2->winpos_delta2000*DEG)*sin(obj2->winpos_deltas*DEG))
		/(cos(obj2->winpos_delta2000*DEG)*cos(obj2->winpos_deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->winposerr_theta2000 = obj2->winposerr_thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }

      if (FLAG(obj2.winposerr_theta1950))
        {
        da = wcs->ap1950 - obj2->winpos_alpha1950;
        dd = (sin(wcs->dp1950*DEG)
		-sin(obj2->winpos_delta1950*DEG)*sin(obj2->winpos_deltas*DEG))
		/(cos(obj2->winpos_delta1950*DEG)*cos(obj2->winpos_deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->winposerr_theta1950 = obj2->winposerr_thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }
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

