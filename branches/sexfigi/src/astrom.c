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
*	Last modify:	13/07/2006
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
#include	"wcs/tnx.h"

static obj2struct	*obj2 = &outobj2;

/****************************** initastrom **********************************/
/*
Initialize astrometrical structures.
*/
void	initastrom(picstruct *field)

  {
   astromstruct	*as;
   double	*lm;
   int		l,n, lng,lat, naxis;

  as = field->astrom;

  naxis = as->naxis;

/* Test if the WCS is in use */
  if (as->wcs_flag)
/*-- ...Yes: launch the WCS stuff! */
    {
    QCALLOC(as->lin, struct linprm, 1);
    QMALLOC(as->cel, struct celprm, 1);
    QMALLOC(as->prj, struct prjprm, 1);
    QMALLOC(as->lin->cdelt, double, naxis);
    QMALLOC(as->lin->crpix, double, naxis);
    QMALLOC(as->lin->pc, double, naxis*naxis);
/* Set WCS flags to 0: structures will be reinitialized by the WCS library */
    as->wcs->flag = as->lin->flag = as->cel->flag = as->prj->flag = 0;
    as->lin->naxis = as->naxis;

/* wcsprm structure */
    lng = as->lng = as->wcs->lng;
    lat = as->lat = as->wcs->lat;

/*-- linprm structure */
    for (l=0; l<naxis; l++)
      {
      as->lin->crpix[l] = as->crpix[l];
      as->lin->cdelt[l] = as->cdelt[l];
      as->cel->ref[l] = as->crval[l];
      }
    for (l=0; l<naxis*naxis; l++)
      as->lin->pc[l] = as->pc[l];

/*-- celprm structure */
    if (lng>=0)
      {
      as->cel->ref[0] = as->crval[lng];
      as->cel->ref[1] = as->crval[lat];
      }
    else
      {
      as->cel->ref[0] = as->crval[0];
      as->cel->ref[1] = as->crval[1];
      }
    as->cel->ref[2] = as->longpole;
    as->cel->ref[3] = as->latpole;

/* prjprm structure */
    as->prj->r0 = as->r0;
    as->prj->tnx_lngcor = as->tnx_lngcor;
    as->prj->tnx_latcor = as->tnx_latcor;
    if (lng>=0)
      {
      n = 0;
      for (l=100; l--;)
        {
        as->prj->p[l] = as->projp[l+lng*100];
        as->prj->p[l+100] = as->projp[l+lat*100];
        if (!n && (as->prj->p[l] || as->prj->p[l+100]))
          n = l+1;
        }
    }

/*-- Compute an "average linear matrix" (at field center) */
    compute_wcs(field, (field->width+1)/2.0, (field->height+1)/2.0);

/*---- Compute Pole coordinates in J2000 and/or B1950 for THETAs */
    if (FLAG(obj2.theta2000) || FLAG(obj2.theta1950)
	|| FLAG(obj2.poserr_theta2000) || FLAG(obj2.poserr_theta1950)
	|| FLAG(obj2.win_theta2000) || FLAG(obj2.win_theta1950)
	|| FLAG(obj2.winposerr_theta2000) || FLAG(obj2.winposerr_theta1950))
      {
      if (fabs(as->equinox-2000.0)>0.003)
        precess(as->equinox, 0.0, 90.0, 2000.0, &as->ap2000, &as->dp2000);
      else
        {
        as->ap2000 = 0.0;
        as->dp2000 = 90.0;
        }

      if (FLAG(obj2.theta1950) || FLAG(obj2.poserr_theta1950))
        j2b(as->equinox, as->ap2000, as->dp2000, &as->ap1950, &as->dp1950);
      }
    }
  else
/*-- ...No: compute only the determinant */
    {
/*-- Simplify the original FITS PC matrix */
    lm = as->linmat;
    for (l=0; l<naxis*naxis; l++)
      lm[l] = as->pc[l]*as->cdelt[l/naxis];
/*-- Check valid only in 2D */
    if ((as->lindet = lm[0]*lm[3] - lm[1]*lm[2]) == 0.0)
      warning ("Null determinant in the global distortion matrix:\n",
	"         Some WORLD-parameters will be incorrect");
    }

/* Override astrometric definitions only if user supplies a pixel-scale */
  if (prefs.pixel_scale == 0.0)
    {
    as->pixscale = sqrt(fabs(as->lindet));
    field->pixscale = 3600.0*as->pixscale;	/* in arcsec2 */
    }
  else
    as->pixscale = (field->pixscale=prefs.pixel_scale)/3600.0;

  return;
  }


/**************************** computeastrom *********************************/
/*
Compute real WORLD coordinates and dimensions according to FITS info.
*/
void	computeastrom(picstruct *field, objstruct *obj)

  {
  astromstruct	*as;
  double	*lm, *wcspos;

  as = field->astrom;
  lm = as->linmat;

/* If working with WCS, compute WORLD coordinates and local matrix */
  if (FLAG(obj2.mxw))
    {
    if (as->wcs_flag)
      {
      wcspos = compute_wcs(field, obj2->posx, obj2->posy);
      obj2->alphas = obj2->mxw = wcspos[0];
      obj2->deltas = obj2->myw = wcspos[1];
      if (FLAG(obj2.alpha2000))
        {
        if (fabs(as->equinox-2000.0)>0.003)
          precess(as->equinox, wcspos[0], wcspos[1],
		2000.0, &obj2->alpha2000, &obj2->delta2000);
        else
          {
          obj2->alpha2000 = obj2->mxw;
          obj2->delta2000 = obj2->myw;
          }
        if (FLAG(obj2.alpha1950))
          j2b(as->equinox, obj2->alpha2000, obj2->delta2000,
		&obj2->alpha1950, &obj2->delta1950);
        }
      }
    else
      {
       double	dx,dy;

      dx = obj2->posx - as->crpix[0];
      dy = obj2->posy - as->crpix[1];
      obj2->mxw = as->crval[0]+ lm[0]*dx + lm[1]*dy;	/* CDELT included! */
      obj2->myw = as->crval[1]+ lm[2]*dx + lm[3]*dy;	/* CDELT included! */
      }
    }

/* Idem for peak-flux positions */
  if (FLAG(obj2.peakxw))
    {
    if (as->wcs_flag)
      {
      wcspos = compute_wcs(field, (double)obj->peakx, (double)obj->peaky);
      obj2->peakalphas = obj2->peakxw = wcspos[0];
      obj2->peakdeltas = obj2->peakyw = wcspos[1];
      if (FLAG(obj2.peakalpha2000))
        {
        if (fabs(as->equinox-2000.0)>0.003)
          precess(as->equinox, wcspos[0], wcspos[1],
		2000.0, &obj2->peakalpha2000, &obj2->peakdelta2000);
        else
          {
          obj2->peakalpha2000 = obj2->peakxw;
          obj2->peakdelta2000 = obj2->peakyw;
          }
        if (FLAG(obj2.peakalpha1950))
          j2b(as->equinox, obj2->peakalpha2000, obj2->peakdelta2000,
		&obj2->peakalpha1950, &obj2->peakdelta1950);
        }
      }
    else 
      {
       double	dx,dy;

      dx = obj->peakx - as->crpix[0];
      dy = obj->peaky - as->crpix[1];
      obj2->peakxw = as->crval[0]+ lm[0]*dx + lm[1]*dy;	/* CDELT included! */
      obj2->peakyw = as->crval[1]+ lm[2]*dx + lm[3]*dy;	/* CDELT included! */
      }
    }

/* Idem for Windowed positions */
  if (FLAG(obj2.winpos_xw))
    {
    if (as->wcs_flag)
      {
      wcspos = compute_wcs(field, obj2->winpos_x, obj2->winpos_y);
      obj2->winpos_alphas = obj2->winpos_xw = wcspos[0];
      obj2->winpos_deltas = obj2->winpos_yw = wcspos[1];
      if (FLAG(obj2.winpos_alpha2000))
        {
        if (fabs(as->equinox-2000.0)>0.003)
          precess(as->equinox, wcspos[0], wcspos[1],
		2000.0, &obj2->winpos_alpha2000, &obj2->winpos_delta2000);
        else
          {
          obj2->winpos_alpha2000 = obj2->winpos_xw;
          obj2->winpos_delta2000 = obj2->winpos_yw;
          }
        if (FLAG(obj2.winpos_alpha1950))
          j2b(as->equinox, obj2->winpos_alpha2000, obj2->winpos_delta2000,
		&obj2->winpos_alpha1950, &obj2->winpos_delta1950);
        }
      }
    else
      {
       double	dx,dy;

      dx = obj2->winpos_x - as->crpix[0];
      dy = obj2->winpos_y - as->crpix[1];
      obj2->winpos_xw = as->crval[0]+ lm[0]*dx + lm[1]*dy;/* CDELT included! */
      obj2->winpos_yw = as->crval[1]+ lm[2]*dx + lm[3]*dy;/* CDELT included! */
      }
    }

/* Custom coordinate system for the MAMA machine */
  if (FLAG(obj2.mamaposx))
    {
     double	dx,dy;

    dx = obj2->posx - 0.5;
    dy = obj2->posy - 0.5;
    obj2->mamaposx = (as->crval[1]+lm[2]*dx+lm[3]*dy)
			*(MAMA_CORFLEX+1.0);		/* CDELT included! */
    obj2->mamaposy = (as->crval[0]+lm[0]*dx+lm[1]*dy);	/* CDELT included! */
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
    obj2->npixw = obj->npix*as->pixscale*as->pixscale;
  if (FLAG(obj2.fdnpixw))
    obj2->fdnpixw = obj->fdnpix*as->pixscale*as->pixscale;

  if (FLAG(obj2.fwhmw))
    obj2->fwhmw = obj->fwhm*as->pixscale;

  return;
  }


/****************************** compute_wcs *********************************/
/*
Compute real WORLD coordinates and local distortion matrix according to the
WCS info.
*/
double	*compute_wcs(picstruct *field, double mx, double my)

  {
   astromstruct	*as;
   static double	pixpos[NAXIS], wcspos[NAXIS],wcspos0[2], imgcrd[NAXIS],
			phi,theta;
   double	*lm, al,da,de,cde;
   int		rcode, lng,lat, naxis;

  as = field->astrom;
  lm = as->linmat;

  naxis = as->naxis;
  lng = as->lng;
  lat = as->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }

  pixpos[lng] = mx;
  pixpos[lat] = my;

  if ((rcode=wcsrev((const char(*)[9])as->ctype, as->wcs, pixpos, as->lin,
	imgcrd, as->prj, &phi, &theta, as->crval, as->cel, wcspos)))
    error(EXIT_FAILURE, "*Error* in WCSlib: ", (char *)wcsrev_errmsg[rcode]);

/* Compute the local distortion matrix */
  al = wcspos0[lng<lat?0:1] = wcspos[lng];
  de = wcspos0[lng<lat?1:0] = wcspos[lat];

/* Get world coordinates for vector 1,0 */
  pixpos[lng] = mx + 1.0;
  pixpos[lat] = my;
  if ((rcode=wcsrev((const char(*)[9])as->ctype, as->wcs, pixpos, as->lin,
	imgcrd, as->prj, &phi, &theta, as->crval, as->cel, wcspos)))
    error(EXIT_FAILURE, "*Error* in WCSlib: ", (char *)wcsrev_errmsg[rcode]);

  da = wcspos[lng]-al;
  if (da>180.0)
    da -= 360.0;
  else if (da<-180.0)
    da += 360.0;

  lm[lng] = da*(cde=cos(de*DEG));
  lm[lat] = wcspos[lat] - de;

/* Get world coordinates for vector 0,1 */
/* Second one */
  pixpos[lng] = mx;
  pixpos[lat] = my + 1.0;
  if ((rcode=wcsrev((const char(*)[9])as->ctype, as->wcs, pixpos, as->lin,
	imgcrd, as->prj, &phi, &theta, as->crval, as->cel, wcspos)))
    error(EXIT_FAILURE, "*Error* in WCSlib: ", (char *)wcsrev_errmsg[rcode]);

  da = wcspos[lng]-al;
  if (da>180.0)
    da -= 360.0;
  else if (da<-180.0)
    da += 360.0;

  lm[2] = da*cde;
  lm[3] = wcspos[lat] - de;

  as->lindet = lm[lng+lng*naxis]*lm[lat+lat*naxis]
	- lm[lng+lat*naxis]*lm[lat+lng*naxis];
  if (as->lindet == 0.0)
    warning ("Null determinant in the local distortion matrix:\n",
	"         Some WORLD-parameters will be incorrect");

  if (prefs.pixel_scale == 0.0)
    as->pixscale = sqrt(fabs(as->lindet));

  return wcspos0;
  }


/****************************** astrom_shapeparam ****************************/
/*
Compute shape parameters in WORLD and SKY coordinates.
*/
void	astrom_shapeparam(picstruct *field, objstruct *obj)
  {
   astromstruct	*as;
   double	*lm,
		dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  as = field->astrom;
  lm = as->linmat;

  naxis = as->naxis;
  lng = as->lng;
  lat = as->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = lm[lng+naxis*lng];
  lm1 = lm[lat+naxis*lng];
  lm2 = lm[lng+naxis*lat];
  lm3 = lm[lat+naxis*lat];


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
    if (as->wcs_flag && FLAG(obj2.thetas))
      obj2->thetas = obj2->thetaw + (obj2->thetaw>0.0?-90:90.0);

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (as->wcs_flag)
      {
       double	da,dd;

      if (FLAG(obj2.theta2000))
        {
        da = as->ap2000 - obj2->alpha2000;
        dd = (sin(as->dp2000*DEG)
		-sin(obj2->delta2000*DEG)*sin(obj2->deltas*DEG))
		/(cos(obj2->delta2000*DEG)*cos(obj2->deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->theta2000 = obj2->thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }

      if (FLAG(obj2.theta1950))
        {
        da = as->ap1950 - obj2->alpha1950;
        dd = (sin(as->dp1950*DEG)
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
   astromstruct	*as;
   double	*lm,
		dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  as = field->astrom;
  lm = as->linmat;

  naxis = as->naxis;
  lng = as->lng;
  lat = as->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = lm[lng+naxis*lng];
  lm1 = lm[lat+naxis*lng];
  lm2 = lm[lng+naxis*lat];
  lm3 = lm[lat+naxis*lat];

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
    if (as->wcs_flag && FLAG(obj2.win_thetas))
      obj2->win_thetas = obj2->win_thetaw +
	(obj2->win_thetaw>0.0?-90:90.0);

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (as->wcs_flag)
      {
       double	da,dd;

      if (FLAG(obj2.win_theta2000))
        {
        da = as->ap2000 - obj2->winpos_alpha2000;
        dd = (sin(as->dp2000*DEG)
		-sin(obj2->winpos_delta2000*DEG)*sin(obj2->winpos_deltas*DEG))
		/(cos(obj2->winpos_delta2000*DEG)*cos(obj2->winpos_deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->win_theta2000 = obj2->win_thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }

      if (FLAG(obj2.win_theta1950))
        {
        da = as->ap1950 - obj2->winpos_alpha1950;
        dd = (sin(as->dp1950*DEG)
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
   astromstruct	*as;
   double	*lm,
		dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  as = field->astrom;
  lm = as->linmat;

  naxis = as->naxis;
  lng = as->lng;
  lat = as->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = lm[lng+naxis*lng];
  lm1 = lm[lat+naxis*lng];
  lm2 = lm[lng+naxis*lat];
  lm3 = lm[lat+naxis*lat];

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
    if (as->wcs_flag && FLAG(obj2.poserr_thetas))
      obj2->poserr_thetas = obj2->poserr_thetaw
				+ (obj2->poserr_thetaw>0.0? -90:90.0);

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (as->wcs_flag)
      {
       double	da,dd;

      if (FLAG(obj2.poserr_theta2000))
        {
        da = as->ap2000 - obj2->alpha2000;
        dd = (sin(as->dp2000*DEG)
		-sin(obj2->delta2000*DEG)*sin(obj2->deltas*DEG))
		/(cos(obj2->delta2000*DEG)*cos(obj2->deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->poserr_theta2000 = obj2->poserr_thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }

      if (FLAG(obj2.poserr_theta1950))
        {
        da = as->ap1950 - obj2->alpha1950;
        dd = (sin(as->dp1950*DEG)
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
   astromstruct	*as;
   double	*lm,
		dx2,dy2,dxy, xm2,ym2,xym, temp,pm2, lm0,lm1,lm2,lm3;
   int		lng,lat, naxis;

  as = field->astrom;
  lm = as->linmat;

  naxis = as->naxis;
  lng = as->lng;
  lat = as->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }
  lm0 = lm[lng+naxis*lng];
  lm1 = lm[lat+naxis*lng];
  lm2 = lm[lng+naxis*lat];
  lm3 = lm[lat+naxis*lat];

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
    if (as->wcs_flag && FLAG(obj2.winposerr_thetas))
      obj2->winposerr_thetas = obj2->winposerr_thetaw
				+ (obj2->winposerr_thetaw>0.0? -90:90.0);

/*-- Compute position angles in J2000 or B1950 reference frame */
    if (as->wcs_flag)
      {
       double	da,dd;

      if (FLAG(obj2.winposerr_theta2000))
        {
        da = as->ap2000 - obj2->winpos_alpha2000;
        dd = (sin(as->dp2000*DEG)
		-sin(obj2->winpos_delta2000*DEG)*sin(obj2->winpos_deltas*DEG))
		/(cos(obj2->winpos_delta2000*DEG)*cos(obj2->winpos_deltas*DEG));
        dd = dd<1.0? (dd>-1.0?acos(dd)/DEG:180.0) : 0.0;
        obj2->winposerr_theta2000 = obj2->winposerr_thetas
		+ (((da>0.0 && da<180.0) || da<-180.0)?-dd:dd);
        }

      if (FLAG(obj2.winposerr_theta1950))
        {
        da = as->ap1950 - obj2->winpos_alpha1950;
        dd = (sin(as->dp1950*DEG)
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


/******************************* copyastrom *********************************/
/*
Copy astrometrical structures.
*/
void	copyastrom(picstruct *infield, picstruct *outfield)

  {
   astromstruct	*inas, *outas;
   int		naxis;

  if (infield->astrom)
    {
    QMEMCPY(infield->astrom, outfield->astrom, astromstruct, 1);
    inas = infield->astrom;
    outas = outfield->astrom;
    naxis = inas->naxis;
    if (inas->wcs_flag)
      {
      QMEMCPY(inas->wcs, outas->wcs, struct wcsprm, 1);
      QMEMCPY(inas->lin, outas->lin, struct linprm, 1);
      QMEMCPY(inas->cel, outas->cel, struct celprm, 1);
      QMEMCPY(inas->prj, outas->prj, struct prjprm, 1);
      QMEMCPY(inas->lin->cdelt, outas->lin->cdelt, double, naxis);
      QMEMCPY(inas->lin->crpix, outas->lin->crpix, double, naxis);
      QMEMCPY(inas->lin->pc, outas->lin->pc, double, naxis*naxis);
      outas->tnx_lngcor = copy_tnxaxis(inas->tnx_lngcor);
      outas->tnx_latcor = copy_tnxaxis(inas->tnx_latcor);
      }
    }

  return;
  }


/******************************* endastrom ***********************************/
/*
Free astrometrical structures.
*/
void	endastrom(picstruct *field)

  {
   astromstruct	*as;

  as = field->astrom;
  if (as->wcs_flag)
    {
    free(as->lin->cdelt);
    free(as->lin->crpix);
    free(as->lin->pc);
    free(as->wcs);
    free(as->lin);
    free(as->cel);
    free(as->prj);
    free_tnxaxis(as->tnx_lngcor);
    free_tnxaxis(as->tnx_latcor);
    }

  free(as);

  return;
  }


/********************************* precess ***********************************/
/*
precess equatorial coordinates according to the equinox (from Ephemerides du
Bureau des Longitudes 1992). Epoch for coordinates should be J2000
(FK5 system).
*/
void	precess(double yearin, double alphain, double deltain,
		double yearout, double *alphaout, double *deltaout)

  {
   double	dzeta,theta,z, t1,t1t1, t2,t2t2,t2t2t2,
		cddsadz, cddcadz, cdd, sdd, adz, cdin,sdin,ct,st,caindz;

  alphain *= DEG;
  deltain *= DEG;

  t1 = (yearin - 2000.0)/1000.0;
  t2 = (yearout - yearin)/1000.0;
  t1t1 = t1*t1;
  t2t2t2 = (t2t2 = t2*t2)*t2;
  theta = (97171.735e-06 - 413.691e-06*t1 - 1.052e-06 * t1t1) * t2
	+ (-206.846e-06 - 1.052e-06*t1) * t2t2 - 202.812e-06 * t2t2t2;
  dzeta = (111808.609e-06 + 677.071e-06*t1 - 0.674e-06 * t1t1) * t2
	+ (146.356e-06 - 1.673e-06*t1) * t2t2 + 87.257e-06 * t2t2t2;
  z = (111808.609e-06 +677.071e-06*t1 - 0.674e-06 * t1t1) * t2
	+ (530.716e-06 + 0.320e-06*t1) * t2t2 + 88.251e-06 * t2t2t2;
  cddsadz = (cdin=cos(deltain)) * sin(alphain+dzeta);
  cddcadz = -(sdin=sin(deltain))*(st=sin(theta))
	+cdin*(ct=cos(theta))*(caindz=cos(alphain+dzeta));
  sdd = sdin*ct + cdin*st*caindz;
  cdd = cos(*deltaout = asin(sdd));
  adz = asin(cddsadz/cdd);
  if (cddcadz<0.0)
    adz = PI - adz;
  if (adz<0.0)
    adz += 2.0*PI;
  adz += z;
  *alphaout = adz/DEG;
  *deltaout /= DEG;

  return;
  }


/*********************************** j2b *************************************/
/*
conver equatorial coordinates from equinox and epoch J2000 to equinox and
epoch B1950 for extragalactic sources (from Aoki et al. 1983, after
inversion of their matrix and some custom arrangements).
*/
void    j2b(double yearobs, double alphain, double deltain,
	double *alphaout, double *deltaout)
  {
   int			i,j;
   static double	a[3] = {-1.62557e-6, -0.31919e-6, -0.13843e-6},
                        ap[3] = {1.245e-3, -1.580e-3, -0.659e-3},
                        m[6][6] = {
  { 0.9999256794678425,    0.01118148281196562,   0.004859003848996022,
   -2.423898417033081e-06,-2.710547600126671e-08,-1.177738063266745e-08},
  {-0.01118148272969232,   0.9999374849247641,   -2.717708936468247e-05,
    2.710547578707874e-08,-2.423927042585208e-06, 6.588254898401055e-11},
  {-0.00485900399622881,  -2.715579322970546e-05, 0.999988194643078,
    1.177738102358923e-08, 6.582788892816657e-11,-2.424049920613325e-06},
  {-0.0005508458576414713, 0.2384844384742432,   -0.4356144527773499,
    0.9999043171308133,    0.01118145410120206,   0.004858518651645554},
  {-0.2385354433560954,   -0.002664266996872802,  0.01225282765749546,
   -0.01118145417187502,   0.9999161290795875,   -2.717034576263522e-05},
  { 0.4357269351676567,   -0.008536768476441086,  0.002113420799663768,
   -0.004858518477064975, -2.715994547222661e-05, 0.9999668385070383}},
			a1[3], r[3], ro[3], r1[3], r2[3], v1[3], v[3];
   static double	cai, sai, cdi, sdi, dotp, rmod, alpha, delta, t1;

/* Convert Julian years from J2000.0 to tropic centuries from B1950.0 */
  t1 = ((yearobs - 2000.0) + (MJD2000 - MJD1950)/365.25)*JU2TROP/100.0;
  alphain *= DEG;
  deltain *= DEG;
  cai = cos(alphain);
  sai = sin(alphain);
  cdi = cos(deltain);
  sdi = sin(deltain);
  r[0] = cdi*cai;
  r[1] = cdi*sai;
  r[2] = sdi;
  for (i=0; i<3; i++)
    v[i] = r2[i] = v1[i] = 0.0;
  for (j=0; j<6; j++)
    for (i=0; i<6; i++)
      if (j<3)
        r2[j] += m[j][i]*(i<3?r[i]:v[i-3]);
      else
        v1[j-3] += m[j][i]*(i<3?r[i]:v[i-3]);

  for (i=0; i<3; i++)
    r1[i] = r2[i]+v1[i]*ARCSEC*t1;

  dotp = 0.0;
  for (i=0; i<3; i++)
    {
    a1[i] = a[i]+ap[i]*ARCSEC*t1;
    dotp += a1[i]*(r1[i]+a1[i]);
    }
  dotp = 2.0/(sqrt(1+4.0*dotp)+1.0);
  rmod = 0.0;
  for (i=0; i<3; i++)
    {
    ro[i] = dotp*(r1[i]+a1[i]);
    rmod += ro[i]*ro[i];
    }
  rmod = sqrt(rmod);
  delta = asin(ro[2]/rmod);
  alpha = acos(ro[0]/cos(delta)/rmod);
  if (ro[1]<0)
    alpha = 2.0*PI - alpha;
  *alphaout = alpha/DEG;
  *deltaout = delta/DEG;

  return;
  }

