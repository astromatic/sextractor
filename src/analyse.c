/*
*				analyse.c
*
* Do measurements on detections.
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
*	Last modified:		26/06/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"analyse.h"
#include	"back.h"
#include	"catout.h"
#include	"clean.h"
#include	"check.h"
#include	"assoc.h"
#include	"astrom.h"
#include	"plist.h"
#include	"flag.h"
#include	"graph.h"
#include	"growth.h"
#include	"image.h"
#include	"misc.h"
#include	"neurro.h"
#include	"objlist.h"
#include	"photom.h"
#include	"psf.h"
#include	"profit.h"
#include	"retina.h"
#include	"som.h"
#include	"subimage.h"
#include	"weight.h"
#include	"winpos.h"

#ifdef USE_THREADS
#include	"threads.h"

extern pthread_t	*pthread_thread;
extern pthread_attr_t	pthread_attr;
extern pthread_mutex_t	pthread_objlistmutex, pthread_freeobjmutex;
extern pthread_cond_t	pthread_objlistaddcond, pthread_objsavecond;
extern fieldstruct	**pthread_fields,**pthread_wfields;
extern objliststruct	*pthread_objlist;
extern int		pthread_nobjlist, pthread_objlistaddindex,
			pthread_objlistprocindex, pthread_objlistsaveindex,
			pthread_nfield, pthread_nthreads, pthread_endflag;

#endif

/****** analyse_iso *********************************************************
PROTO	void analyse_iso(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist, int n)
PURPOSE	Do (isophotal) measurements on pixel lists.
INPUT	Pointer to an array of image field pointers,
	pointer to an array of weight-map field pointers,
	number of images,
	pointer to the objlist,
	object index in the objlist.
OUTPUT	-.
NOTES	Requires access to the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	03/10/2012
 ***/
void  analyse_iso(fieldstruct **fields, fieldstruct **wfields, int nfield,
			objliststruct *objlist, int n)
  {
   fieldstruct		*field, *wfield;
   checkstruct		*check;
   objstruct		*obj;
   pliststruct		*pixel, *pixt;
   PIXTYPE		threshs[NISO],
			pix, cdpix, tpix, peak,cdpeak, thresh,dthresh,minthresh;
   double		tv,sigtv, ngamma,
			esum, emx2,emy2,emxy, err,gain,backnoise2,dbacknoise2,
			xm,ym, x,y,var,var2, threshfac;
   float		*heap,*heapt,*heapj,*heapk, swap;
   int			i,j,k,h, photoflag,area,errflag, cleanflag,
			pospeakflag, minarea, gainflag;


  obj = objlist->obj+n;
  pixel = objlist->plist;
  field = fields[0];
/* Prepare computation of positional error */
  esum = emx2 = emy2 = emxy = 0.0;
  if ((errflag=FLAG(obj2.poserr_mx2)))
    {
    dbacknoise2 = field->backsig*field->backsig;
    xm = obj->mx;
    ym = obj->my;
    }
  else
    xm = ym = dbacknoise2 = 0.0;	/* to avoid gcc -Wall warnings */

  pospeakflag = FLAG(obj2.peakx);
  gain = field->gain;
  ngamma = field->ngamma;
  photoflag = (field->detector_type==DETECTOR_PHOTO);
  gainflag = PLISTEXIST(var) && prefs.weightgain_flag;

  h = minarea = prefs.ext_minarea;

/* Prepare selection of the heap selection for CLEANing */
  if ((cleanflag = prefs.clean_flag))
    {
    if (obj->fdnpix < minarea)
      {
      obj->mthresh = 0.0;
      cleanflag = 0;
      heapt = heap = NULL;		/* to avoid gcc -Wall warnings */
      }
    else
      {
      QMALLOC(heap, float, minarea);
      heapt = heap;
      }
    }
  else
    {
    obj->mthresh = 0.0;
    heapt = heap = NULL;		/* to avoid gcc -Wall warnings */
    }


/* Measure essential isophotal parameters in the measurement image... */
  tv = sigtv = 0.0;
  var = backnoise2 = field->backsig*field->backsig;
  peak = -BIG;
  cdpeak = -BIG;
  thresh = field->thresh;
  minthresh = (PLISTEXIST(var))? BIG : thresh;
  threshfac = field->backsig > 0.0 ? field->thresh / field->backsig : 1.0;
dthresh = field->dthresh;
  area = 0;
  for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
    {
    pix = PLIST(pixt,value);
    if (pix>peak)
      peak = pix;

    cdpix=PLISTPIX(pixt,cvalue);
    if (pospeakflag && cdpix>cdpeak)
      {
      cdpeak=cdpix;
      obj->dpeakx =  PLIST(pixt,x) + 1;
      obj->dpeaky =  PLIST(pixt,y) + 1;
      }
    if (PLISTEXIST(var))
      {
      var = PLISTPIX(pixt, var);
      thresh = threshfac*sqrt(var);
      if (thresh < minthresh)
        minthresh = thresh;
      }

    if (photoflag)
      {
      pix = exp(pix/ngamma);
      var2 = pix*pix*var;
      }
    else
      var2 = var;

    if (gainflag && pix>0.0 && gain>0.0)
      var2 += pix/gain*var/backnoise2;

    sigtv += var2;
    if (pix>thresh)
      area++;
    tv += pix;
    if (errflag)
      {
      err = dbacknoise2;
      if (gain>0.0 && cdpix>0.0)
        err += cdpix/gain;
      x = PLIST(pixt,x) - xm;
      y = PLIST(pixt,y) - ym;
      esum += err;
      emx2 += err*x*x;
      emy2 += err*y*y;
      emxy += err*x*y;
      }

/*-- Find the minareath pixel in decreasing intensity for CLEANing */
    if (cleanflag)
      {
      tpix = PLISTPIX(pixt, cvalue) - (PLISTEXIST(dthresh)?
		PLISTPIX(pixt, dthresh):dthresh);
      if (h>0)
        *(heapt++) = (float)tpix;
      else if (h)
        {
        if ((float)tpix>*heap)
          {
          *heap = (float)tpix;
          for (j=0; (k=(j+1)<<1)<=minarea; j=k)
            {
            heapk = heap+k;
            heapj = heap+j;
            if (k != minarea && *(heapk-1) > *heapk)
              {
              heapk++;
              k++;
              }
            if (*heapj <= *(--heapk))
              break;
            swap = *heapk;
            *heapk = *heapj;
            *heapj = swap;
            }
          }
        }
      else
        fqmedian(heap, minarea);
      h--;
      }
    }

/* Flagging from the flag-map */
  if (PLISTEXIST(flag))
    flag_get(obj, pixel);

/* Flag and count pixels with a low weight */
  if (PLISTEXIST(wflag))
    weight_count(obj, pixel);

  if (cleanflag)
    {
    obj->mthresh = *heap;
    free(heap);
    }

  if (errflag)
    {
     double	flux2;

    flux2 = obj->fdflux*obj->fdflux;
/*-- Estimation of the error ellipse moments: we re-use previous variables */
    emx2 /= flux2;	/* variance of xm */
    emy2 /= flux2;	/* variance of ym */
    emxy /= flux2;	/* covariance */

/*-- Handle fully correlated profiles (which cause a singularity...) */
    esum *= 0.08333/flux2;
    if (obj->singuflag && (emx2*emy2-emxy*emxy) < esum*esum)
      {
      emx2 += esum;
      emy2 += esum;
      }

    obj->poserr_mx2 = emx2;
    obj->poserr_my2 = emy2;
    obj->poserr_mxy = emxy;
    }
 
  if (photoflag)
    {
    tv = ngamma*(tv-obj->fdnpix*exp(obj->dbkg/ngamma));
    sigtv /= ngamma*ngamma;
    }
  else
    {
    tv -= obj->fdnpix*obj->dbkg;
    if (!gainflag && gain > 0.0 && tv>0.0)
      sigtv += tv/gain;
    }

  obj->npix = area;
  obj->dflux = tv;
  obj->dfluxerr = sigtv;
  obj->dpeak = peak;
  obj->thresh = minthresh - obj->dbkg;
  obj->dpeak -= obj->dbkg;

/* Initialize isophotal thresholds so as to sample optimally the full profile*/

  if (FLAG(obj2.iso[0]))
    {
     int	*iso;
     PIXTYPE	*thresht;

    memset(obj->iso, 0, NISO*sizeof(int));
    if (field->detector_type == DETECTOR_PHOTO)
      for (i=0; i<NISO; i++)
        threshs[i] = obj->thresh + (obj->dpeak-obj->thresh)*i/NISO;
    else
      {
      if (obj->dpeak>0.0 && obj->thresh>0.0)
        for (i=0; i<NISO; i++)
          threshs[i] = obj->thresh*pow(obj->dpeak/obj->thresh, (double)i/NISO);
      else
        for (i=0; i<NISO; i++)
          threshs[i] = 0.0;
      }
    for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
      for (i=NISO,iso=obj->iso,thresht=threshs;
		i-- && PLIST(pixt,value)>*(thresht++);)
        (*(iso++))++;
    }

/* Put objects in "segmentation check-image" */

  if ((check = prefs.check[CHECK_SEGMENTATION]))
    for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
      ((ULONG *)check->pix)[check->width*PLIST(pixt,y)+PLIST(pixt,x)]
		= -(ULONG)obj->number;

  if ((check = prefs.check[CHECK_OBJECTS]))
    for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
      ((PIXTYPE *)check->pix)[check->width*PLIST(pixt,y)+PLIST(pixt,x)]
		= PLIST(pixt,value);

/* Compute the FWHM of the object */
  if (FLAG(obj2.fwhm))
    {
     PIXTYPE	thresh0;

    thresh0 = obj->dpeak/5.0;
    if (thresh0<obj->thresh)
      thresh0 = obj->thresh;
    if (thresh0>0.0)
      {
       double	mx,my, s,sx,sy,sxx,sxy, dx,dy,d2, lpix,pix, b, inverr2, sat,
		dbkg, d, bmax;

      mx = obj->mx;
      my = obj->my;
      dbkg = obj->dbkg;
      sat = (double)(field->satur_level - obj->bkg);
      s = sx = sy = sxx = sxy = 0.0;
      for (pixt=pixel+obj->firstpix;pixt>=pixel;pixt=pixel+PLIST(pixt,nextpix))
        {
        pix = PLIST(pixt,value)-dbkg;
        if (pix>thresh0 && pix<sat)
          {
          dx = PLIST(pixt,x) - mx;
          dy = PLIST(pixt,y) - my;
          lpix = log(pix);
          inverr2 = pix*pix;
          s += inverr2;
          d2 = dx*dx+dy*dy;
          sx += d2*inverr2;
          sxx += d2*d2*inverr2;
          sy += lpix*inverr2;
          sxy += lpix*d2*inverr2;
          }        
        }
      d = s*sxx-sx*sx;
      if (fabs(d)>0.0)
        {
        b = -(s*sxy-sx*sy)/d;
        if (b<(bmax = 1/(13*obj->a*obj->b)))	/* to have FWHM # 6 sigma */
          b = bmax;
        obj->fwhm = (float)(1.6651/sqrt(b));
/*----- correction for undersampling effects (established from simulations) */
        if (obj->fwhm>0.0)
          obj->fwhm -= 1/(4*obj->fwhm);
        }
      else
        obj->fwhm = 0.0;
      }
    else
      obj->fwhm = 0.0;
    }

  return;
  }


/****** analyse_final *******************************************************
PROTO	void analyse_final(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist)
PURPOSE Do the final analysis based on a list of detections and a detection
	index.
INPUT	Pointer to an array of image field pointers,
	pointer to an array of weight-map field pointers,
	number of images,
	objlist pointer.
OUTPUT  -.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 15/05/2014
 ***/
void analyse_final(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist)
  {
   fieldstruct		*field, *wfield;
   objstruct		*obj;
   int			o;

// field is the detection field */
  field = fields[0];			// field is the detection field
  wfield = wfields? wfields[0]:NULL;	// wfield is the detection weight map

  obj = objlist->obj;
  for (o=objlist->nobj; o--; obj++)
    if (!(obj->obj2 = analyse_obj2obj2(fields, wfields, nfield, obj)))
      error(EXIT_FAILURE, "*Error*: ", "obj2 stack full");

#ifdef USE_THREADS
  if (prefs.nthreads>1)
/*-- Push in objlist */
    pthread_objlist_add(*objlist);
  else
    {
/*-- 1 single thread: don't bother with independent measurement threads */
/*-- Analyse the objects of objlist */
    obj = objlist->obj;
    for (o=objlist->nobj; o--; obj++)
      obj->writable_flag
	= (analyse_full(fields, wfields, nfield, obj) == RETURN_OK);
/*-- Write to catalogue */
    analyse_end(fields, wfields, nfield, objlist);
    }
#else
    obj = objlist->obj;
    for (o=objlist->nobj; o--; obj++)
      obj->writable_flag
	= (analyse_full(fields, wfields, nfield, obj) == RETURN_OK);
/*-- Write to catalogue */
    analyse_end(fields, wfields, nfield, objlist);
#endif

  return;
  }



/****** analyse_obj2obj2 ******************************************************
PROTO	obj2struct *analyse_obj2obj2(fieldstruct **fields,
			fieldstruct **wfields, int nfield, objstruct *obj)
PURPOSE Move object data from obj to obj2 structure.
INPUT	Pointer to an array of image field pointers,
	pointer to an array of weight-map field pointers,
	number of images,
	obj pointer.
OUTPUT  New obj2 pointer.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 26/06/2014
TODO	Free memory allocated for bkg, dbkg, etc.
 ***/
obj2struct	*analyse_obj2obj2(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objstruct *obj)
  {
   checkstruct		*check;
   fieldstruct		*field;
   obj2struct		*obj2;
   float		sigbkg;
   static int		number;
   int			f, idx,idy;

  QMALLOC(obj2, obj2struct, 1);
  obj->obj2 = obj2;

/* Allocate memory for multidimensional catalog parameter arrays */
  catout_allocobjparams(obj);

/* Allocate memory for other arrays (not catalogue measurements) */
/* Sky background */
//  catout_alloobjcother(obj, &flagobj2.diffbkg, prefs.nimage*sizeof(float));
  catout_allocobjother(obj, &flagobj2.bkg, prefs.nimage*sizeof(float));
  catout_allocobjother(obj, &flagobj2.dbkg, prefs.nimage*sizeof(float));
  catout_allocobjother(obj, &flagobj2.sigbkg, prefs.nimage*sizeof(float));
/* Flux combination */
  catout_allocobjother(obj, &flagobj2.cflux, prefs.nphotinstru*sizeof(double));
  catout_allocobjother(obj, &flagobj2.cfluxw, prefs.nphotinstru*sizeof(double));
/* Position combination */
  catout_allocobjother(obj, &flagobj2.cposx, prefs.nphotinstru*sizeof(double));
  catout_allocobjother(obj, &flagobj2.cposy, prefs.nphotinstru*sizeof(double));
  catout_allocobjother(obj, &flagobj2.cposw, prefs.nphotinstru*sizeof(double));

/*-- Local backgrounds */
  for (f=0; f<nfield; f++)
    {
    field = fields[f];
    if (FLAG(obj2.bkg))
      obj2->bkg[f] = (float)back_interpolate(field, obj->mx, obj->my);
    obj2->dbkg[f] = 0.0;
    if (prefs.pback_type == LOCAL)
      {
      obj2->dbkg[f] = back_local(field, obj, &sigbkg);
      if (FLAG(obj2.bkg))
        obj2->bkg[f] += obj2->dbkg[f];
      obj2->sigbkg[f] = sigbkg<0.0? field->backsig : sigbkg;          
      }
    else
      obj2->sigbkg[f] = field->backsig;
    if (FLAG(obj2.wflag))
      obj2->wflag[f] = obj->wflag;
    }

/* Copy main data */
  obj2->number = ++number;
  obj2->fdnpix = obj->fdnpix;
  obj2->dnpix = obj->dnpix;
  obj2->npix = obj->npix;
  obj2->nzdwpix = obj->nzdwpix;
  obj2->nzwpix = obj->nzwpix;
  obj2->dflux = obj->dflux;
  obj2->dfluxerr = obj->dfluxerr;
  obj2->fdpeak = obj->fdpeak;
  obj2->peak = obj->dpeak;
  obj2->peakx = obj->dpeakx;
  obj2->peaky = obj->dpeaky;
  obj2->mx = obj->mx;
  obj2->my = obj->my;
/* Integer coordinates */
  obj2->ix=(int)(obj2->mx+0.49999);		/* Integer coordinates */
  obj2->iy=(int)(obj2->my+0.49999);
  obj2->posx[0] = obj2->mx+1.0;			/* That's standard FITS */
  obj2->posy[0] = obj2->my+1.0;
  obj2->poserr_mx2 = obj->poserr_mx2;
  obj2->poserr_my2 = obj->poserr_my2;
  obj2->poserr_mxy = obj->poserr_mxy;
  obj2->xmin = obj->xmin;
  obj2->xmax = obj->xmax;
  obj2->ymin = obj->ymin;
  obj2->ymax = obj->ymax;
  obj2->flags = obj->flag;
  obj2->singuflag = obj->singuflag;
  if (FLAG(obj2.imaflags))
    {
    memcpy(obj2->imaflags, obj->imaflags, prefs.nfimage*sizeof(FLAGTYPE));
    memcpy(obj2->imanflags, obj->imanflags, prefs.nfimage*sizeof(int));
    }
  obj2->mx2 = obj->mx2;
  obj2->my2 = obj->my2;
  obj2->mxy = obj->mxy;
  obj2->a = obj->a;
  obj2->b = obj->b;
  obj2->theta = obj->theta;
  obj2->abcor = obj->abcor;
  obj2->cxx = obj->cxx;
  obj2->cyy = obj->cyy;
  obj2->cxy = obj->cxy;
  obj2->thresh = obj->thresh;
  obj2->dthresh = obj->dthresh;
  obj2->mthresh = obj->mthresh;
  memcpy(obj2->iso, obj->iso, NISO*sizeof(int));
  obj2->fwhm = obj->fwhm;

/* Copy image data around current object */
  subimage_getall(fields, wfields, nfield, obj2);

/* if BLANKing is on, paste back the object pixels in the first subimage*/
  if (prefs.blank_flag && obj->isoimage)
    subimage_fill(obj2->subimage, obj->isoimage);

  if ((check=prefs.check[CHECK_SEGMENTATION]))
    {
/*-- Re-number segmentation map */
     ULONG	*pix;
     ULONG	oldsnumber=-obj->number, newsnumber=obj2->number;
     int	dx,dx0,dy,dpix;

    pix = (ULONG *)check->pix + check->width*obj->ymin + obj->xmin;
    dx0 = obj->xmax-obj->xmin+1;
    dpix = check->width-dx0;
    for (dy=obj->ymax-obj->ymin+1; dy--; pix += dpix)
      for (dx=dx0; dx--; pix++)
        if (*pix==oldsnumber)
          *pix = newsnumber;
    }

  return obj2;
  }


/****** analyse_full *********************************************************
PROTO	int analyse_full(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objstruct *obj)
PURPOSE Final analysis of object data.
INPUT   Pointer to an array of image field pointers,
	pointer to an array of weight-map field pointers,
	number of images,
	objstruct pointer.
OUTPUT  RETURN_OK if the object has been processed, RETURN_ERROR otherwise.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 19/05/2014
 ***/
int	analyse_full(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objstruct *obj)
  {
   obj2struct		*obj2;
   fieldstruct		*field, *dfield, *dwfield, *wfield;
   checkstruct		*check;
   double		rawpos[NAXIS],
			analtime1;
   int			i,j, ix,iy, idx,idy, selecflag, newnumber,nsub;

  if (FLAG(obj2.analtime))
    analtime1 = counter_seconds();
  else
    analtime1 = 0.0;		/* To avoid gcc -Wall warnings */

  obj2 = obj->obj2;

/* field is the detection field */
  field = fields[0];

dfield = dwfield = wfield = NULL;
   if (prefs.psf_flag)
     obj2->psf_flag = 0;	/* Reset PSF building flag */
   if (prefs.dpsf_flag)
     obj2->dpsf_flag = 0;	/* Reset PSF building flag */

/*------------------------------ Association ------------------------------*/
  if (prefs.assoc_flag)
    {
    obj2->assoc_number = do_assoc(field, obj2->posx[0],obj2->posy[0],
			obj2->assoc);
    if ((prefs.assocselec_type==ASSOCSELEC_MATCHED && !(obj2->assoc_number))
	||
	(prefs.assocselec_type==ASSOCSELEC_NOMATCHED && (obj2->assoc_number)))
      {
/*---- Treatment of discarded detections */
/*---- update segmentation map  and exit */
      if ((check=prefs.check[CHECK_SEGMENTATION]))
        {
         ULONG	*pix;
         ULONG	oldsnumber = obj2->number;
         int	dx,dx0,dy,dpix;

        pix = (ULONG *)check->pix + check->width*obj2->ymin + obj2->xmin;
        dx0 = obj2->xmax-obj2->xmin+1;
        dpix = check->width-dx0;
        for (dy=obj2->ymax-obj2->ymin+1; dy--; pix += dpix)
          for (dx=dx0; dx--; pix++)
            if (*pix==oldsnumber)
              *pix = 0;
        }
      return RETURN_ERROR;
      }
    }

/*------------------------- Error ellipse parameters ------------------------*/
  if (FLAG(obj2.poserr_a))
    {
     double	pmx2,pmy2,temp,theta;

    if (fabs(temp=obj2->poserr_mx2-obj2->poserr_my2) > 0.0)
      theta = atan2(2.0 * obj2->poserr_mxy,temp) / 2.0;
    else
      theta = PI/4.0;

    temp = sqrt(0.25*temp*temp+obj2->poserr_mxy*obj2->poserr_mxy);
    pmy2 = pmx2 = 0.5*(obj2->poserr_mx2+obj2->poserr_my2);
    pmx2+=temp;
    pmy2-=temp;

    obj2->poserr_a = (float)sqrt(pmx2);
    obj2->poserr_b = (float)sqrt(pmy2);
    obj2->poserr_theta = theta*180.0/PI;
    }

  if (FLAG(obj2.poserr_cxx))
    {
     double	xm2,ym2, xym, temp;

    xm2 = obj2->poserr_mx2;
    ym2 = obj2->poserr_my2;
    xym = obj2->poserr_mxy;
    obj2->poserr_cxx = (float)(ym2/(temp=xm2*ym2-xym*xym));
    obj2->poserr_cyy = (float)(xm2/temp);
    obj2->poserr_cxy = (float)(-2*xym/temp);
    }

/* Aspect ratio */
  if (FLAG(obj2.elong))
    obj2->elong = obj2->a/obj2->b;

  if (FLAG(obj2.ellip))
    obj2->ellip = 1-obj2->b/obj2->a;

  if (FLAG(obj2.polar))
    obj2->polar = (obj2->a*obj2->a-obj2->b*obj2->b)
		/ (obj2->a*obj2->a+obj2->b*obj2->b);

/* Express positions in FOCAL or WORLD coordinates */
  if (FLAG(obj2.posxf) || FLAG(obj2.posxw))
    astrom_pos(fields, nfield, obj2);

  obj2->pixscale2 = 0.0;	/* To avoid gcc -Wall warnings */
  if (FLAG(obj2.mx2w)
	|| FLAG(obj2.win_mx2w)
	|| FLAG(obj2.poserr_mx2w)
	|| FLAG(obj2.winposerr_mx2w)
	|| FLAG(obj2.poserrmx2w_psf)
	|| FLAG(obj2.poserrmx2w_prof)
	|| FLAG(obj2.prof_flagw)
	|| ((!prefs.pixel_scale) && FLAG(obj2.area_flagw))
	|| ((!prefs.pixel_scale) && FLAG(obj2.fwhmw_psf)))
    {
    rawpos[0] = obj2->posx[0];
    rawpos[1] = obj2->posy[0];
    obj2->pixscale2 = wcs_jacobian(field->wcs, rawpos, obj2->jacob);
    }

/* Express shape parameters in the FOCAL or WORLD frame */
  if (FLAG(obj2.mx2w))
    astrom_shapeparam(field, obj2);
/* Express position error parameters in the FOCAL or WORLD frame */
  if (FLAG(obj2.poserr_mx2w))
    astrom_errparam(field, obj2);

  if (FLAG(obj2.npixw))
    obj2->npixw = obj2->npix * (prefs.pixel_scale?
	field->pixscale/3600.0*field->pixscale/3600.0 : obj2->pixscale2);
  if (FLAG(obj2.fdnpixw))
    obj2->fdnpixw = obj2->fdnpix * (prefs.pixel_scale?
	field->pixscale/3600.0*field->pixscale/3600.0 : obj2->pixscale2);

  if (FLAG(obj2.fwhmw))
    obj2->fwhmw = obj2->fwhm * (prefs.pixel_scale?
	field->pixscale/3600.0 : sqrt(obj2->pixscale2));

/*------------------------------- Photometry -------------------------------*/

/* Convert the father of photom. error estimates from variance to RMS */
  obj2->flux_iso = obj2->dflux;
  obj2->fluxerr_iso = sqrt(obj2->dfluxerr);

  if (FLAG(obj2.flux_isocor))
    photom_isocor(field, obj2);

  if (FLAG(obj2.flux_aper))
    for (i=0; i<prefs.naper; i++)
      photom_aper(field, wfield, obj2, i);

  if ((prefs.auto_flag) && !prefs.prof_flag)
    photom_auto(fields, wfields, nfield, obj2);

  if (FLAG(obj2.flux_petro))
    photom_petro(field, dfield, wfield, dwfield, obj2);

/*------------------------------ Astrometry -------------------------------*/
/* Express positions in FOCAL or WORLD coordinates */
  if (FLAG(obj2.peakxf) || FLAG(obj2.peakxw))
    astrom_peakpos(field, obj2);
  if (obj2->peak+obj2->bkg[0] >= field->satur_level)
      obj2->flags |= OBJ_SATUR;

/* Estimate of shape */
  growth_aver(fields, wfields, nfield, obj2);
/* Get a good estimate of position and flux */
  if ((prefs.win_flag))
    win_pos(fields, wfields, nfield, obj2);
/* Express positions in FOCAL or WORLD coordinates */
  if (FLAG(obj2.winpos_xf) || FLAG(obj2.winpos_xw))
    astrom_winpos(fields, nfield, obj2);
/* Express shape parameters in the FOCAL or WORLD frame */
  if (FLAG(obj2.win_mx2w))
    astrom_winshapeparam(fields, nfield, obj2);
/* Express position error parameters in the FOCAL or WORLD frame */
  if (FLAG(obj2.winposerr_mx2w))
    astrom_winerrparam(fields, nfield, obj2);

/* Check-image CHECK_APERTURES option */
  if ((check = prefs.check[CHECK_APERTURES]))
    {
    if (FLAG(obj2.flux_aper))
      for (i=0; i<prefs.naper; i++)
        sexcircle(check->pix, check->width, check->height,
		obj2->mx, obj2->my, prefs.apert[i]/2.0, check->overlay);

    if (FLAG(obj2.flux_auto))
      sexellipse(check->pix, check->width, check->height,
		obj2->mx, obj2->my, obj2->a*obj2->auto_kronfactor,
		obj2->b*obj2->auto_kronfactor, obj2->theta,
		check->overlay, obj2->flags&OBJ_CROWDED);

    if (FLAG(obj2.flux_petro))
      sexellipse(check->pix, check->width, check->height,
		obj2->mx, obj2->my, obj2->a*obj2->petrofactor,
		obj2->b*obj2->petrofactor, obj2->theta,
		check->overlay, obj2->flags&OBJ_CROWDED);
    }

/*----------------------------- Model fitting -----------------------------*/
  if (prefs.prof_flag)
    {
/*-- Perform measurements on the fitted models */
    profit_measure(obj2->profit, obj2);
/*-- Express positions in FOCAL or WORLD coordinates */
    if (FLAG(obj2.xf_prof) || FLAG(obj2.xw_prof))
      astrom_profpos(field, obj2);
/*-- Express shape parameters in the FOCAL or WORLD frame */
    if (FLAG(obj2.prof_flagw))
      astrom_profshapeparam(field, obj2);
/*-- Express position error parameters in the FOCAL or WORLD frame */
    if (FLAG(obj2.poserrmx2w_prof))
      astrom_proferrparam(field, obj2);
/*-- PSF- and model-guided star/galaxy separation */
    if (FLAG(obj2.prof_class_star) || FLAG(obj2.prof_concentration))
      profit_spread(obj2->profit, field, wfield, obj);
    }

/*------------------------------- PSF fitting ------------------------------
  nsub = 1;
  if (prefs.psffit_flag)
      {
      if (prefs.dpsffit_flag)
        double_psf_fit(ppsf, field, wfield, obj2, thepsf, dfield, dwfield);
      else
        psf_fit(thepsf, field, wfield, obj2);
      obj2->npsf = thepsfit->npsf;
      nsub = thepsfit->npsf;
      if (nsub<1)
        nsub = 1;
      }

    for (j=0; j<nsub; j++)
      {
      if (prefs.psffit_flag)
        {
        obj2->x_psf = thepsfit->x[j];
        obj2->y_psf = thepsfit->y[j];
        if (FLAG(obj2.xf_psf) || FLAG(obj2.xw_psf))
          astrom_psfpos(field, obj2);
*------ Express position error parameters in the FOCAL or WORLD frame *
        if (FLAG(obj2.poserrmx2w_psf))
          astrom_psferrparam(field, obj2);
        if (FLAG(obj2.flux_psf))
          obj2->flux_psf = thepsfit->flux[j]>0.0? thepsfit->flux[j]:0.0; *?*
        if (FLAG(obj2.mag_psf))
          obj2->mag_psf = thepsfit->flux[j]>0.0?
		prefs.mag_zeropoint -2.5*log10(thepsfit->flux[j]) : 99.0;
        if (FLAG(obj2.fluxerr_psf))
          obj2->fluxerr_psf= thepsfit->fluxerr[j];
        if (FLAG(obj2.magerr_psf))
          obj2->magerr_psf =
		(thepsfit->flux[j]>0.0 && thepsfit->fluxerr[j]>0.0) ? *?*
			1.086*thepsfit->fluxerr[j]/thepsfit->flux[j] : 99.0;
        }
      }
*/

/* Express everything in magnitude units */
  photom_mags(field, obj2);

/*--------------------------------- "Vignets" ------------------------------*/
  if (FLAG(obj2.vignet))
    copyimage(field,obj2->vignet,prefs.vignet_size[0],prefs.vignet_size[1],
	obj2->ix, obj2->iy, -BIG);

  if (FLAG(obj2.vigshift))
    copyimage_center(field, obj2->vigshift, prefs.vigshift_size[0],
		prefs.vigshift_size[1], obj2->mx, obj2->my);

/* Edit min and max coordinates to follow the FITS conventions */
  obj2->xmin += 1;
  obj2->ymin += 1;
  obj2->xmax += 1;
  obj2->ymax += 1;

/* Count source */
#ifdef USE_THREADS
  if (prefs.nthreads>1)
    {
    QPTHREAD_MUTEX_LOCK(&pthread_countobjmutex);
    thecat.ntotal++;
     if (((prefs.prof_flag) && !(thecat.ntotal%(10*pthread_nthreads)))
		|| !(thecat.ntotal%(100*pthread_nthreads)))
      NPRINTF(OUTPUT, "\33[1M> Line:%5d  "
	"Objects: %8d detected / %8d sextracted\n\33[1A",
	thecat.nline, thecat.ndetect, thecat.ntotal);
    QPTHREAD_MUTEX_UNLOCK(&pthread_countobjmutex);
    }
  else
#endif
    {
    thecat.ntotal++;
    if (((prefs.prof_flag) && !(thecat.ntotal%10)) || !(thecat.ntotal%100))
      NPRINTF(OUTPUT, "\33[1M> Line:%5d  "
	"Objects: %8d detected / %8d sextracted\n\33[1A",
	thecat.nline, thecat.ndetect, thecat.ntotal);
    }

/* Processing time */
  obj2->analtime = (float)(counter_seconds() - analtime1);

  return RETURN_OK;
  }


/****** analyse_end *********************************************************
PROTO	void analyse_end(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist)
PURPOSE Write to catalogue measurements made on a list of detections and
	release objects.
INPUT   Pointer to an array of image field pointers,
	pointer to an array of weight-map field pointers,
	number of images,
	objlist pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 09/06/2014
 ***/
void	analyse_end(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist)
  {
   objstruct	*obj;
   obj2struct	*obj2;
   int		o;

  obj=objlist->obj;
  for (o=objlist->nobj; o--; obj++)
    {
    obj2 = obj->obj2;
    if ((obj->writable_flag))
      {
/*---- Catalogue output */
      FPRINTF(OUTPUT, "%8d %6.1f %6.1f %5.1f %5.1f %12g "
			"%c%c%c%c%c%c%c%c\n",
	obj2->number, obj2->mx+1.0, obj2->my+1.0,
	obj2->a, obj2->b,
	obj2->dflux,
	obj2->flags&OBJ_CROWDED?'C':'_',
	obj2->flags&OBJ_MERGED?'M':'_',
	obj2->flags&OBJ_SATUR?'S':'_',
	obj2->flags&OBJ_TRUNC?'T':'_',
	obj2->flags&OBJ_APERT_PB?'A':'_',
	obj2->flags&OBJ_ISO_PB?'I':'_',
	obj2->flags&OBJ_DOVERFLOW?'D':'_',
	obj2->flags&OBJ_OVERFLOW?'O':'_');
      catout_writeobj(obj);
      }
    }

  objlist->done_flag = 0;

#ifdef USE_THREADS
  if (prefs.nthreads>1)
    QPTHREAD_MUTEX_LOCK(&pthread_freeobjmutex);
  QPTHREAD_COND_BROADCAST(&pthread_objsavecond);
  if (prefs.nthreads>1)
    QPTHREAD_MUTEX_UNLOCK(&pthread_freeobjmutex);
#endif

  return;
  }


