/*
*				analyse.c
*
* Do measurements on detections.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		05/10/2011
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
#include	"growth.h"
#include	"image.h"
#include	"photom.h"
#include	"psf.h"
#include	"profit.h"
#include	"retina.h"
#include	"som.h"
#include	"weight.h"
#include	"winpos.h"

/****** analyse_iso *********************************************************
PROTO	void analyse_iso(picstruct *field, picstruct *dfield,
			objliststruct *objlist, int n)
PURPOSE	Do (isophotal) on pixel lists in the MEASUREMENT image.
INPUT	Pointer to the measurement image,
	pointer to the detection image (if different from measurement),
	pointer to the objlist,
	object index in the objlist.
OUTPUT	-.
NOTES	Requires access to the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	18/07/2011
 ***/
void  analyse_iso(picstruct *field, picstruct *dfield,
			objliststruct *objlist, int n)
  {
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


  if (!dfield)
    dfield = field;

  obj = objlist->obj+n;
  pixel = objlist->plist;

/* Prepare computation of positional error */
  esum = emx2 = emy2 = emxy = 0.0;
  if ((errflag=FLAG(obj.poserr_mx2)))
    {
    dbacknoise2 = dfield->backsig*dfield->backsig;
    xm = obj->mx;
    ym = obj->my;
    }
  else
    xm = ym = dbacknoise2 = 0.0;	/* to avoid gcc -Wall warnings */

  pospeakflag = FLAG(obj.peakx);
  gain = field->gain;
  ngamma = field->ngamma;
  photoflag = (prefs.detect_type==PHOTO);
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
  dthresh = dfield->dthresh;
  area = 0;
  for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
    {
    pix = PLIST(pixt,value);
    if (pix>peak)
      peak = pix;

    cdpix=PLISTPIX(pixt,cdvalue);
    if (pospeakflag && cdpix>cdpeak)
      {
      cdpeak=cdpix;
      obj->peakx =  PLIST(pixt,x) + 1;
      obj->peaky =  PLIST(pixt,y) + 1;
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
      tpix = PLISTPIX(pixt, cdvalue) - (PLISTEXIST(dthresh)?
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
    getflags(obj, pixel);

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
  obj->flux = tv;
  obj->fluxerr = sigtv;
  obj->peak = peak;
  obj->thresh = minthresh - obj->dbkg;
  obj->peak -= obj->dbkg;

/* Initialize isophotal thresholds so as to sample optimally the full profile*/

  if (FLAG(obj.iso[0]))
    {
     int	*iso;
     PIXTYPE	*thresht;

    memset(obj->iso, 0, NISO*sizeof(int));
    if (prefs.detect_type == PHOTO)
      for (i=0; i<NISO; i++)
        threshs[i] = obj->thresh + (obj->peak-obj->thresh)*i/NISO;
    else
      {
      if (obj->peak>0.0 && obj->thresh>0.0)
        for (i=0; i<NISO; i++)
          threshs[i] = obj->thresh*pow(obj->peak/obj->thresh, (double)i/NISO);
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
		= (ULONG)obj->number;

  if ((check = prefs.check[CHECK_OBJECTS]))
    for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
      ((PIXTYPE *)check->pix)[check->width*PLIST(pixt,y)+PLIST(pixt,x)]
		= PLIST(pixt,value);

/* Compute the FWHM of the object */
  if (FLAG(obj.fwhm))
    {
     PIXTYPE	thresh0;

    thresh0 = obj->peak/5.0;
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
PROTO	void analyse_final(picstruct *field, picstruct *dfield,
		picstruct *wfield, picstruct *dwfield,
		objliststruct *objlist, int iobj)
PURPOSE Do the final analysis based on a list of detections and a detection
	index.
INPUT   Measurement field pointer,
        detection field pointer,
        measurement weight-map field pointer,
        detection weight-map field pointer,
	objlist pointer,
	obj index.
OUTPUT  -.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 07/10/2011
 ***/
void analyse_final(picstruct *field, picstruct *dfield,
		picstruct *wfield, picstruct *dwfield,
		objliststruct *objlist, int iobj)
  {
   picstruct		*cfield;
   obj2liststruct	*obj2list;
   objstruct		*obj;
   obj2struct		*obj2, *prevobj2, *firstobj2;
   int			i, ymax, nextiobj, noverlap, iobj2;

  obj2list = thecat.obj2list;

/* cfield is the detection field in any case */
  cfield = dfield? dfield:field;
/* Find overlapping detections and link them */
  noverlap = analyse_overlapness(cleanobjlist, iobj);
/* Convert every linked detection to a linked obj2 */
  prevobj2 = NULL;
  for (i=noverlap; i--;)
    {
    obj = &objlist->obj[iobj];
/*-- Warn if there is a possibility for any aperture to be truncated */
    if ((ymax=obj->ycmax) > cfield->ymax)
      {
      sprintf(gstr, "Object at position %.0f,%.0f ", obj->mx+1, obj->my+1);
      QWARNING(gstr, "may have some apertures truncated:\n"
		"          You might want to increase MEMORY_BUFSIZE");
      }
    else if (ymax>cfield->yblank && prefs.blank_flag)
      {
      sprintf(gstr, "Object at position %.0f,%.0f ", obj->mx+1, obj->my+1);
      QWARNING(gstr, "may have some unBLANKed neighbours:\n"
		"          You might want to increase MEMORY_PIXSTACK");
      }
    obj2 = analyse_obj2obj2(field,dfield, wfield,dwfield, obj, obj2list);
    if (!obj2)
      error(EXIT_FAILURE, "*Error*: ", "obj2 stack full");
    obj2->nextobj2 = NULL;
    if (prevobj2)
      {
      prevobj2->nextobj2 = obj2;
      obj2->prevobj2 = prevobj2;
      }
    else
      {
      obj2->prevobj2 = NULL;
      firstobj2 = obj2;
      }
    prevobj2 = obj2;
    nextiobj = obj->next;
/*-- Take care of next obj that might be swapped by subcleanobj! */
    if (nextiobj == objlist->nobj-1)
      nextiobj = iobj;
    subcleanobj(iobj);
    iobj = nextiobj;
    }

/* Analyse the group of obj2s and write out catalogue */
  analyse_group(field, dfield, wfield, dwfield, firstobj2);

/* Free the group of obj2s */
  for (obj2=firstobj2; obj2; obj2=obj2->nextobj2)
    {
    QFREE(obj2->image);
    if (dfield)
      QFREE(obj2->dimage);
    if (wfield)
      QFREE(obj2->weight);
    if (dwfield)
      QFREE(obj2->dweight);
    }

  for (obj2=firstobj2; obj2->nextobj2; obj2=obj2->nextobj2);
  obj2->nextobj2 = obj2list->freeobj2;
  obj2list->freeobj2->prevobj2 = obj2->nextobj2;
  obj2list->freeobj2 = firstobj2;

  return;
  }


/****** analyse_overlapness ******************************************************
PROTO	obj2struct analyse_overlapness(objliststruct *objlist, int iobj)
PURPOSE Link together overlapping detections.
INPUT   objliststruct pointer,
	obj index.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 07/10/2011
 ***/
int analyse_overlapness(objliststruct *objlist, int iobj)
  {
   objstruct	*obj,*cobj,*fobj;
   int		i, blend, nblend,nobj;

  nblend = 1;
  nobj = objlist->nobj;
  obj = objlist->obj;
  fobj = &obj[iobj];
  fobj->prev = -1;
  blend = fobj->blend;
  cobj = fobj;
  for (i=0; i<nobj; i++, obj++)
    if (obj->blend == blend && obj!=fobj)
      {
      cobj->next = i;
      obj->prev = iobj;
      cobj = obj;
      iobj = i;
      nblend++;
      }

  cobj->next = -1;

  return nblend;
  }


/****** analyse_obj2obj2 ******************************************************
PROTO	void analyse_obj2obj2(picstruct *field, picstruct *dfield,
		picstruct *wfield, picstruct *dwfield,
		objstruct *obj, obj2liststruct *obj2list)
PURPOSE Move object data from obj to obj2 structure.
INPUT   Measurement field pointer,
        detection field pointer,
        measurement weight-map field pointer,
        detection weight-map field pointer,
	obj pointer,
	obj2list pointer,
OUTPUT  New obj2 pointer.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 06/10/2011
 ***/
obj2struct	*analyse_obj2obj2(picstruct *field, picstruct *dfield,
		picstruct *wfield, picstruct *dwfield,
		objstruct *obj, obj2liststruct *obj2list)
  {
   obj2struct	*obj2;
   int		idx,idy;

  obj2 = obj2list->freeobj2;
  if (obj2->nextobj2)
    {
    obj2list->freeobj2 = obj2->nextobj2;
    obj2->nextobj2 = obj2->prevobj2 = NULL;
    }
  else
    return NULL;

/* Copy main data */
  obj2->number = obj->number;
  obj2->fdnpix = obj->fdnpix;
  obj2->dnpix = obj->dnpix;
  obj2->npix = obj->npix;
  obj2->nzdwpix = obj->nzdwpix;
  obj2->nzwpix = obj->nzwpix;
  obj2->fdflux = obj->fdflux;
  obj2->dflux = obj->dflux;
  obj2->flux = obj->flux;
  obj2->fluxerr = obj->fluxerr;
  obj2->fdpeak = obj->fdpeak;
  obj2->dpeak = obj->dpeak;
  obj2->peak = obj->peak;
  obj2->peakx = obj->peakx;
  obj2->peaky = obj->peaky;
  obj2->mx = obj->mx;
  obj2->my = obj->my;
/* Integer coordinates */
  obj2->ix=(int)(obj2->mx+0.49999);		    /* Integer coordinates */
  obj2->iy=(int)(obj2->my+0.49999);
  obj2->sposx = (float)(obj2->posx = obj2->mx+1.0); /* That's standard FITS */
  obj2->sposy = (float)(obj2->posy = obj2->my+1.0);
  obj2->poserr_mx2 = obj->poserr_mx2;
  obj2->poserr_my2 = obj->poserr_my2;
  obj2->poserr_mxy = obj->poserr_mxy;
  obj2->xmin = obj->xmin;
  obj2->xmax = obj->xmax;
  obj2->ymin = obj->ymin;
  obj2->ymax = obj->ymax;
  obj2->flag = obj->flag;
  obj2->wflag = obj->wflag;
  memcpy(obj2->imaflag, obj->imaflag, prefs.imaflag_size*sizeof(FLAGTYPE));
  obj2->singuflag = obj->singuflag;
  memcpy(obj2->imanflag, obj->imanflag, prefs.imanflag_size*sizeof(int));
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
  obj2->bkg = obj->bkg;
  obj2->dbkg = obj->dbkg;
  obj2->sigbkg = obj->sigbkg;
  obj2->thresh = obj->thresh;
  obj2->dthresh = obj->dthresh;
  obj2->mthresh = obj->mthresh;
  memcpy(obj2->iso, obj->iso, NISO*sizeof(int));
  obj2->fwhm = obj->fwhm;

/* Copy image data around current object */
  obj2->imsize[0] = 2.0*(obj2->xmax-obj2->xmin)+1+2*field->stripmargin;
  obj2->imsize[1] = 2.0*(obj2->ymax-obj2->ymin)+1+2*field->stripmargin;
  obj2->immin[0] = obj2->ix - obj2->imsize[0]/2;
  obj2->immin[1] = obj2->iy - obj2->imsize[1]/2;
  obj2->immax[0] = obj2->immin[0] + obj2->imsize[0];
  obj2->immax[1] = obj2->immin[1] + obj2->imsize[1];
  QMALLOC(obj2->image, PIXTYPE, obj2->imsize[0]*obj2->imsize[1]);
  copyimage(field, obj2->image, obj2->imsize[0],obj2->imsize[1],
	obj2->ix,obj2->iy);
  if (dfield)
    {
    QMALLOC(obj2->dimage, PIXTYPE, obj2->imsize[0]*obj2->imsize[1]);
    copyimage(dfield, obj2->dimage, obj2->imsize[0],obj2->imsize[1],
	obj2->ix,obj2->iy);
    }
  if (wfield)
    {
    QMALLOC(obj2->weight, PIXTYPE, obj2->imsize[0]*obj2->imsize[1]);
    copyimage(wfield, obj2->weight, obj2->imsize[0],obj2->imsize[1],
	obj2->ix,obj2->iy);
    }
  if (dwfield)
    {
    QMALLOC(obj2->dweight, PIXTYPE, obj2->imsize[0]*obj2->imsize[1]);
    copyimage(dwfield, obj2->dweight, obj2->imsize[0],obj2->imsize[1],
	obj2->ix,obj2->iy);
    }
/* if BLANKing is on, paste back the object pixels in the image*/
  if (prefs.blank_flag)
    {
/*-- Compute coordinates of blank start in object image */
    idx = obj->subx - obj2->immin[0];
    idy = obj->suby - obj2->immin[1];
    if (obj->blank)
      {
      deblankimage(obj->blank, obj->subw, obj->subh,
		obj2->image, obj2->imsize[0],obj2->imsize[1], idx,idy);
      free(obj->blank);
      }
    if (obj->dblank)
      {
      deblankimage(obj->dblank, obj->subw, obj->subh,
		obj2->dimage, obj2->imsize[0],obj2->imsize[1], idx,idy);
      free(obj->dblank);
      }
    }

  return obj2;
  }


/****** analyse_full *********************************************************
PROTO	void analyse_full(picstruct *field, picstruct *dfield,
		picstruct *wfield, picstruct *dwfield, obj2struct *obj2)
PURPOSE Final analysis of object data.
INPUT   Measurement field pointer,
        Detection field pointer,
        Measurement weight-map field pointer,
        Detection weight-map field pointer,
	obj2struct pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 06/10/2011
 ***/
void	analyse_full(picstruct *field, picstruct *dfield,
		picstruct *wfield, picstruct *dwfield, obj2struct *obj2)
  {
   checkstruct		*check;
   double		rawpos[NAXIS],
			analtime1;
   int			i,j, ix,iy, idx,idy, selecflag, newnumber,nsub;

  if (FLAG(obj2.analtime))
    analtime1 = counter_seconds();
  else
    analtime1 = 0.0;		/* To avoid gcc -Wall warnings */

   if (prefs.psf_flag)
     obj2->psf_flag = 0;	/* Reset PSF building flag */
   if (prefs.dpsf_flag)
     obj2->dpsf_flag = 0;	/* Reset PSF building flag */

/*------------------------------ Association ------------------------------*/
  if (prefs.assoc_flag)
    {
    obj2->assoc_number = do_assoc(field, obj2->posx, obj2->posy, obj2->assoc);

    if (prefs.assocselec_type!=ASSOCSELEC_ALL
		&& ((prefs.assocselec_type==ASSOCSELEC_MATCHED)?
		obj2->assoc_number:!obj2->assoc_number))
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
      return;
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
  if (FLAG(obj2.mxf) || FLAG(obj2.mxw))
    astrom_pos(field, obj2);

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
    rawpos[0] = obj2->posx;
    rawpos[1] = obj2->posy;
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
  obj2->flux_iso = obj2->flux;
  obj2->fluxerr_iso = sqrt(obj2->fluxerr);

  if (FLAG(obj2.flux_isocor))
    photom_isocor(field, obj2);

  if (FLAG(obj2.flux_aper))
    for (i=0; i<prefs.naper; i++)
      photom_aper(field, wfield, obj2, i);

  if (FLAG(obj2.flux_auto))
    photom_auto(field, dfield, wfield, dwfield, obj2);

  if (FLAG(obj2.flux_petro))
    photom_petro(field, dfield, wfield, dwfield, obj2);

/*------------------------------ Astrometry -------------------------------*/
/* Express positions in FOCAL or WORLD coordinates */
  if (FLAG(obj2.peakxf) || FLAG(obj2.peakxw))
    astrom_peakpos(field, obj2);
  if (obj2->peak+obj2->bkg >= field->satur_level)
    obj2->flag |= OBJ_SATUR;
/* Estimate of shape */
  growth_aver(field, wfield, obj2);
/* Get a good estimate of position and flux */
  compute_winpos(field, wfield, obj2);
/* Express positions in FOCAL or WORLD coordinates */
  if (FLAG(obj2.winpos_xf) || FLAG(obj2.winpos_xw))
    astrom_winpos(field, obj2);
/* Express shape parameters in the FOCAL or WORLD frame */
  if (FLAG(obj2.win_mx2w))
    astrom_winshapeparam(field, obj2);
/* Express position error parameters in the FOCAL or WORLD frame */
  if (FLAG(obj2.winposerr_mx2w))
    astrom_winerrparam(field, obj2);

/* Check-image CHECK_APERTURES option */
  if ((check = prefs.check[CHECK_APERTURES]))
    {
    if (FLAG(obj2.flux_aper))
      for (i=0; i<prefs.naper; i++)
        sexcircle(check->pix, check->width, check->height,
		obj2->mx, obj2->my, prefs.apert[i]/2.0, check->overlay);

    if (FLAG(obj2.flux_auto))
      sexellips(check->pix, check->width, check->height,
		obj2->mx, obj2->my, obj2->a*obj2->kronfactor,
		obj2->b*obj2->kronfactor, obj2->theta,
		check->overlay, obj2->flag&OBJ_CROWDED);

    if (FLAG(obj2.flux_petro))
      sexellips(check->pix, check->width, check->height,
		obj2->mx, obj2->my, obj2->a*obj2->petrofactor,
		obj2->b*obj2->petrofactor, obj2->theta,
		check->overlay, obj2->flag&OBJ_CROWDED);
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
      profit_spread(obj2->profit, field, wfield, obj2);
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
    copyimage(field,obj2->vignet,prefs.vignetsize[0],prefs.vignetsize[1],
	obj2->ix, obj2->iy);

  if (FLAG(obj2.vigshift))
    copyimage_center(field, obj2->vigshift, prefs.vigshiftsize[0],
		prefs.vigshiftsize[1], obj2->mx, obj2->my);


/*------------------------------- Source index -----------------------------*/
  newnumber = ++thecat.ntotal;
/*-- update segmentation map */
  if ((check=prefs.check[CHECK_SEGMENTATION]))
    {
     ULONG	*pix;
     ULONG	newsnumber = newnumber,
		oldsnumber = obj2->number;
     int	dx,dx0,dy,dpix;

    pix = (ULONG *)check->pix + check->width*obj2->ymin + obj2->xmin;
    dx0 = obj2->xmax-obj2->xmin+1;
    dpix = check->width-dx0;
    for (dy=obj2->ymax-obj2->ymin+1; dy--; pix += dpix)
      for (dx=dx0; dx--; pix++)
        if (*pix==oldsnumber)
          *pix = newsnumber;
    }
  obj2->number = newnumber;

/* Edit min and max coordinates to follow the FITS conventions */
  obj2->xmin += 1;
  obj2->ymin += 1;
  obj2->xmax += 1;
  obj2->ymax += 1;

/* Processing time */
  obj2->analtime = (float)(counter_seconds() - analtime1);

  return;
  }


/****** analyse_group ********************************************************
PROTO	void analyse_group(picstruct *field, picstruct *dfield,
		picstruct *wfield, picstruct *dwfield, obj2struct *fobj2)
PURPOSE Perform measurements on a group of detections.
INPUT   Measurement field pointer,
        Detection field pointer,
        Measurement weight-map field pointer,
        Detection weight-map field pointer,
	obj2struct pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 07/10/2011
 ***/
void	analyse_group(picstruct *field, picstruct *dfield,
		picstruct *wfield, picstruct *dwfield, obj2struct *fobj2)
  {
   obj2struct		*obj2, *modobj2;
   int			i;

  if (prefs.prof_flag)
    {
/*-- Setup model fitting for this group */
    for (obj2=fobj2; obj2; obj2=obj2->nextobj2)
      {
      photom_auto(field, dfield, wfield, dwfield, obj2);
      growth_aver(field, wfield, obj2);
      obj2->profit = profit_init(field, wfield,
		obj2, thepsf, prefs.prof_modelflags);
      }
    for (obj2=fobj2; obj2; obj2=obj2->nextobj2)
      {
      photom_auto(field, dfield, wfield, dwfield, obj2);
      growth_aver(field, wfield, obj2);
      }
    if (fobj2->nextobj2)
      {
/*---- Iterative multiple fit if several sources overlap */
      for (i=0; i<ANALYSE_NMULTITER; i++)
        {
        for (obj2=fobj2; obj2; obj2=obj2->nextobj2)
          profit_fit(obj2->profit, obj2);
        for (obj2=fobj2; obj2; obj2=obj2->nextobj2)
          {
          if (i)
            profit_copyobjpix(obj2->profit, field, wfield, obj2);
          for (modobj2=fobj2; modobj2; modobj2=modobj2->nextobj2)
            if (modobj2 != obj2)
              profit_submodpix(obj2->profit, modobj2->profit, 0.9);
          }
        }
      for (obj2=fobj2; obj2; obj2=obj2->nextobj2)
        for (modobj2=fobj2; modobj2; modobj2=modobj2->nextobj2)
          if (modobj2 != obj2)
            addtobig(modobj2->profit->lmodpix,
		modobj2->profit->objnaxisn[0],modobj2->profit->objnaxisn[1],
		obj2->image, obj2->imsize[0], obj2->imsize[1],
		modobj2->profit->ix-obj2->ix, modobj2->profit->iy-obj2->iy,
		-1.0);

      }
    else
/*---- One single source */
      profit_fit(fobj2->profit, fobj2);
    }

/* Full source analysis and catalogue output */
  for (obj2=fobj2; obj2; obj2=obj2->nextobj2)
    {
    analyse_full(field, dfield, wfield, dwfield, obj2);
/*-- Catalogue output */
    FPRINTF(OUTPUT, "%8d %6.1f %6.1f %5.1f %5.1f %12g "
			"%c%c%c%c%c%c%c%c\n",
	obj2->number, obj2->mx+1.0, obj2->my+1.0,
	obj2->a, obj2->b,
	obj2->flux,
	obj2->flag&OBJ_CROWDED?'C':'_',
	obj2->flag&OBJ_MERGED?'M':'_',
	obj2->flag&OBJ_SATUR?'S':'_',
	obj2->flag&OBJ_TRUNC?'T':'_',
	obj2->flag&OBJ_APERT_PB?'A':'_',
	obj2->flag&OBJ_ISO_PB?'I':'_',
	obj2->flag&OBJ_DOVERFLOW?'D':'_',
	obj2->flag&OBJ_OVERFLOW?'O':'_');
    catout_writeobj(obj2);
    }

/* Deallocate memory used for model-fitting */
  if (prefs.prof_flag)
    for (obj2=fobj2; obj2; obj2=obj2->nextobj2)
      profit_end(obj2->profit);

  return;
  }


