/**
* @file		objlist.c
* @brief	Manage object lists (e.g., for advanced deblending)
* @date		11/07/2014
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2014 IAP/CNRS/UPMC
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
#include	"catout.h"
#include	"clean.h"
#include	"check.h"
#include	"plist.h"
#include	"image.h"
#include	"lutz.h"
#include	"misc.h"
#include	"objlist.h"
#include	"profit.h"
#include	"scan.h"
#include	"subimage.h"
#include	"weight.h"

#ifdef USE_THREADS
#include	"threads.h"

pthread_t	*pthread_thread;
pthread_attr_t	pthread_attr;
pthread_mutex_t	pthread_objlistmutex, pthread_freeobjmutex;
pthread_cond_t	pthread_objlistaddcond, pthread_objsavecond;
fieldstruct	**pthread_fields,**pthread_wfields;
objliststruct	*pthread_objlist;
int		pthread_nobjlist, pthread_objlistaddindex,
		pthread_objlistprocindex, pthread_objlistsaveindex,
		pthread_nfield, pthread_nthreads, pthread_endflag;

#endif
extern checkstruct	*thecheck;

/****** objlist_new *******************************************************//**
Create a new, empty objlist.
@param[out] 		Pointer to a new empty objlist.

@author 		E. Bertin (IAP)
@date			16/05/2014
 ***/
objliststruct	*objlist_new(void) {

   objliststruct *objlist;

  QCALLOC(objlist, objliststruct, 1);

  return objlist;
}


/****** objlist_end *******************************************************//**
Free resources allocated for an objlist
@param[in] objlist	Pointer to the objlist.

@author 		E. Bertin (IAP)
@date			26/06/2014
 ***/
void	objlist_end(objliststruct *objlist) {

  objstruct	*obj;
  int		o;

  obj = objlist->obj;
  for (o=objlist->nobj; o--; obj++)
    obj_end(obj);

  free(objlist->obj);
  free(objlist->objindex);
  free(objlist->plist);

  if (objlist->subimage) {
    subimage_end(objlist->subimage);
    free(objlist->subimage);
  }
  free(objlist);

  return;
}


/****** obj_end **********************************************************//**
Free resources allocated for an object
@param[in] obj	Pointer to the object.

@author 		E. Bertin (IAP)
@date			09/07/2014
 ***/
void	obj_end(objstruct *obj) {

  if (obj->isoimage) {
    subimage_end(obj->isoimage);
    QFREE(obj->isoimage);
  }
  if (obj->fullimage) {
    subimage_end(obj->fullimage);
    QFREE(obj->fullimage);
  }

  if (!obj->obj2)
    return;

  catout_freeobjparams(obj);
//  catout_freeother(obj, &flagobj2.diffbkg);
  catout_freeother(obj, &flagobj2.sigbkg);
  catout_freeother(obj, &flagobj2.cflux);
  catout_freeother(obj, &flagobj2.cfluxw);
  catout_freeother(obj, &flagobj2.cposx);
  catout_freeother(obj, &flagobj2.cposy);
  catout_freeother(obj, &flagobj2.cposw);
  QFREE(obj->obj2);

  return;
}


/****** objlist_addobj ****************************************************//**
Add an object to an objlist.
@param[in] objlist	Pointer to the objlist
@param[in] wfields	Pointer to the object to be added
@param[out] 		RETURN_FATAL_ERROR if a memory (re-)allocation issue
			happened,
			RETURN_OK otherwise.
@warning		The ordering of the objlist and the obj pointer may be
			modified.
@author 		E. Bertin (IAP)
@date			26/06/2014
 ***/
int	objlist_addobj(objliststruct *objlist, objstruct *obj,
		pliststruct *plist) {

   objstruct	*objt;
   pliststruct	*plistt;
   int		i, j, nobj, npix;

// Check if memory (re-)allocation is necessary for objects
  if (objlist->nobj >= objlist->nobjmax) {
    nobj = objlist->nobjmax + OBJLIST_NOBJMAXINC;
    if (objlist->nobj >= nobj)
      nobj = objlist->nobj;

    if (!(objt = (objlist->nobjmax) ?
		  (objstruct *)realloc(objlist->obj, nobj * sizeof(objstruct))
		: (objstruct *)malloc(nobj * sizeof(objstruct))))
      return RETURN_FATAL_ERROR;

    objlist->obj = objt;
    objlist->nobjmax = nobj;
  }

// Copy the content of the input object
  objt = objlist->obj + objlist->nobj;
  *objt = *obj;

// Copy new pixels if requested
  if (plist && obj->fdnpix) {
    npix = objlist->npix + obj->fdnpix;
//-- (Re-)allocate memory for pixels

    if (!(plistt = (objlist->npix) ?
		  (pliststruct *)realloc(objlist->plist, npix * plistsize)
		: (pliststruct *)malloc(npix * plistsize)))
      return RETURN_FATAL_ERROR;

    objlist->plist = plistt;	// Update plist pointer
    objt->firstpix = j = objlist->npix * plistsize;
    plistt += j;		// Move to end of current plist
    for (i = obj->firstpix; i != -1; i = PLIST(plist + i, nextpix)) {
      memcpy(plistt, plist + i, (size_t)plistsize);
      PLIST(plistt, nextpix) = (j += plistsize);
      plistt += plistsize;
    }

    PLIST(plistt -= plistsize, nextpix) = -1;
    objt->lastpix = j - plistsize;
    objlist->npix = npix;
  }

// Nothing wrong can happen anymore: update number of objects
  objlist->nobj++;

  return RETURN_OK;
}


/****** objlist_subobj ****************************************************//**
Remove an object from an objlist.
@param[in] objlist	Pointer to the objlist
@param[in] objindex	Index of the object to be removed
@param[out] 		RETURN_ERROR if the objlist is empty,
			RETURN_FATAL_ERROR if a memory (re-)allocation issue
			happened,
			RETURN_OK otherwise.
@warning		The ordering of the list and the obj pointer may be
			modified.
@author 		E. Bertin (IAP)
@date			09/06/2014
 ***/
int	objlist_subobj(objliststruct *objlist, int objindex) {

   int nobjmax;

  if (!(objlist->nobj))
    return RETURN_ERROR;

  if (--objlist->nobj) {
    if (objindex != objlist->nobj)
      objlist->obj[objindex] = objlist->obj[objlist->nobj];
    nobjmax = objlist->nobjmax;
    while (nobjmax - OBJLIST_NOBJMAXINC > objlist->nobj)
      nobjmax -= OBJLIST_NOBJMAXINC;
    if (nobjmax != objlist->nobjmax
	&& !(objlist->obj = (objstruct *)realloc(objlist->obj,
	    (objlist->nobjmax = nobjmax) * sizeof(objstruct))))
      return RETURN_FATAL_ERROR;
  } else {
    QFREE(objlist->obj);
    objlist->nobjmax = 0;
  }

  return RETURN_OK;
}


/****** objlist_movobj ****************************************************//**
Move an object from an objlist to another.
@param[in] objlistin	Pointer to the source objlist
@param[in] objindex	Index of the object to be removed
@param[in] objlistout	Pointer to the destination objlist
@param[out] 		RETURN_ERROR if the objlist is empty,
			RETURN_FATAL_ERROR if a memory (re-)allocation issue
			happened,
			RETURN_OK otherwise.
@warning		The ordering of the source objlist and the source obj
			pointer may be modified. The destination obj pointer may
			may be modified.
@author 		E. Bertin (IAP)
@date			24/06/2014
 ***/
int	objlist_movobj(objliststruct *objlistin, int objindex,
		objliststruct *objlistout) {

  if (objlist_addobj(objlistout, objlistin->obj + objindex, NULL) == RETURN_OK)
    return objlist_subobj(objlistin, objindex);
  else
    RETURN_ERROR;
}


/****** objlist_deblend ***************************************************//**
Perform model-fitting and deblending on a list of objects.
@param[in] fields	Pointer to an array of image field pointers
@param[in] wfields	Pointer to an array of weight-map field pointers
@param[in] nfield	Number of images
@param[in] objlist	Pointer to objlist
@param[in] objindex	Index of blended object in the objlist
@author 		E. Bertin (IAP)
@date			11/09/2014
 ***/
objliststruct	*objlist_deblend(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist, int objindex) {
   fieldstruct		*field, *wfield;
   subprofitstruct	*subprofit, *modsubprofit;
   subimagestruct	*subimage;
   objliststruct	*newobjlist, *overobjlist;
   objstruct		*obj, *modobj, *fobj;
   PIXTYPE		*imagecopy, *imvarcopy, *imwork, *pix, *pixwork,
			val, tol;
   int			xmin,ymin, xmax,ymax, blend, nobj, npix,
			i, j,k, o, o2, p;

  field = fields[0];			// field is the detection field
  wfield = wfields? wfields[0]:NULL;	// wfield is the detection weight map

// Find overlapping detections and link them
  overobjlist = objlist_overlap(objlist, objlist->obj + objindex);

  xmax = ymax = -(xmin = ymin = 0x7FFFFFFF);	// largest signed 32-bit int
  obj = overobjlist->obj;
  blend = obj->blend;				// Keep note of the blend index
  for (o=overobjlist->nobj; o--; obj++) {	// find boundaries of the list
    if (obj->xmin < xmin)			// of overlapping objects
      xmin = obj->xmin;
    if (obj->xmax > xmax)
      xmax = obj->xmax;
    if (obj->ymin < ymin)
      ymin = obj->ymin;
    if (obj->ymax > ymax)
      ymax = obj->ymax;
  }

// TODO: add margin
// Extract large subimage covering the whole objlist
  overobjlist->subimage = subimage = subimage_fromfield(field, wfield,
				xmin, xmax + 1, ymin, ymax + 1);
  npix = subimage->size[0]*subimage->size[1];
  if (prefs.filter_flag) {
//-- Allocate memory for storing the convolved model-subtracted image
    QCALLOC(subimage->fimage, PIXTYPE, npix);
    if (subimage->imvar)
//---- Allocate memory for storing the convolved variance map
      QCALLOC(subimage->fimvar, PIXTYPE, npix);
  }

  obj = overobjlist->obj;
  for (o=overobjlist->nobj; o--; obj++) {
//-- if BLANKing is on, paste back the object pixels in the sub-image
    if (prefs.blank_flag)
      subimage_fill(subimage, obj->isoimage);
  }


// Backup the large image and the variance map
  QMEMCPY(subimage->image, imagecopy, PIXTYPE, npix);
  if (subimage->imvar) {
    QMEMCPY(subimage->imvar, imvarcopy, PIXTYPE, npix);
  } else
    imvarcopy = NULL;

// Create a work buffer of the large image
  QMALLOC(imwork, PIXTYPE, npix);

  obj = overobjlist->obj;
  for (i=overobjlist->nobj; i--; obj++) {
//-- Create individual object subimages
    obj->fullimage = subimage_fromfield(field, wfield,
		obj->xmin, obj->xmax + 1, obj->ymin, obj->ymax + 1);
//-- if BLANKing is on, paste back the object pixels in the sub-images
    if (prefs.blank_flag) {
      subimage_fill(obj->fullimage, obj->isoimage);
    }
    obj->profit = profit_init(obj, obj->fullimage, 1,
			MODEL_MOFFAT, PROFIT_NOCONV);
  }

  for (j=0; j<GROUP_NDEBLENDITER; j++) {
    nobj = overobjlist->nobj;
//-- Iterative multiple fit if several sources overlap
    for (i=0; i<GROUP_NMULTITER; i++) {
      obj = overobjlist->obj;
      for (o=nobj; o--; obj++)
        profit_fit(obj->profit, obj);

      obj = overobjlist->obj;
      for (o=nobj; o--; obj++) {
        subprofit = obj->profit->subprofit;
//------ Subtract the contribution from existing overlapping neighbors (models)
        if (i)
          subprofit_copyobjpix(subprofit, obj->fullimage);
        if (subprofit->nobjpix) {
          modobj = overobjlist->obj;
          for (o2=nobj; o2--; modobj++)
          if (modobj != obj) {
            modsubprofit = modobj->profit->subprofit;
            if (modsubprofit->lmodpix)
              subprofit_addmodpix(modsubprofit, subprofit->objpix,
			subprofit->ix, subprofit->iy,
			subprofit->objnaxisn[0], subprofit->objnaxisn[1],
			modsubprofit->subsamp, nobj>1 ? -0.95: -1.0);
          }
        }
      }
      if (nobj <= 1)
        break;
    }

//-- Reset work buffer
    memset(imwork, 0, npix*sizeof(PIXTYPE));

    obj = overobjlist->obj;
    for (o=nobj; o--; obj++)
      if ((subprofit = obj->profit->subprofit) && subprofit->lmodpix)
//------ Stack up models from all the sources
        subprofit_addmodpix(subprofit, imwork,
			subimage->ipos[0], subimage->ipos[1],
			subimage->size[0], subimage->size[1],
			subprofit->subsamp, 1.0);

//-- Subtract the result from the subimage
    pix = subimage->image;
    pixwork = imwork;
    for (p=npix; p--;)
      *(pix++) -= *(pixwork++);

//-- Update the variance map
    if (subimage->imvar) {
      pix = subimage->imvar;
      pixwork = imwork;
      tol = prefs.deblend_fittol;
      for (p=npix; p--;) {
        val = tol * *(pixwork++);
        *(pix++) += val*val;
      }
    }

    if (prefs.filter_flag) {
      convolve_image(field, subimage->image, subimage->fimage,
		subimage->size[0], subimage->size[1]);
      if (subimage->imvar)
        convolve_image(field, subimage->imvar, subimage->fimvar,
		subimage->size[0], subimage->size[1]);
    }

//-- Extract residual sources (if they exist)
    newobjlist = lutz_subextract(subimage, field->dthresh*2, xmin, xmax + 1,
			ymin, ymax + 1);

//-- Recover (sub) image and variance map
    memcpy(subimage->image, imagecopy, npix*sizeof(PIXTYPE));
    if (subimage->imvar)
      memcpy(subimage->imvar, imvarcopy, npix*sizeof(PIXTYPE));

//-- Loop over input objlist; we keep indexing independent of loop counter to
//-- accommodate changes in the ordering of the list
    for (o=0, k=newobjlist->nobj; k--; o++) {
      obj = newobjlist->obj + o;
      scan_preanalyse(obj, newobjlist->plist, ANALYSE_FULL|ANALYSE_ROBUST);
      if (prefs.ext_maxarea && obj->fdnpix > prefs.ext_maxarea)
        continue; 
      obj->number = ++thecat.ndetect;
      obj->blend = blend;
      obj->deblend_npass = j + 1;
//---- Isophotal measurements
      analyse_iso(fields, wfields, nfield, newobjlist, o);
      if (prefs.blank_flag) {
        if (!(obj->isoimage = subimage_fromplist(field, wfield, obj,
		newobjlist->plist))) {
//-------- Not enough memory for the BLANKing subimage: flag the object now
          obj->flag |= OBJ_OVERFLOW;
          sprintf(gstr, "%.0f,%.0f", obj->mx+1, obj->my+1);
          warning("Memory overflow during masking for detection at ", gstr);
        }
      }
//---- Create individual object subimages
      obj->fullimage = subimage_fromfield(field, wfield,
		obj->xmin, obj->xmax + 1, obj->ymin, obj->ymax + 1);
//---- if BLANKing is on, paste back the object pixels in the sub-images
      if (prefs.blank_flag)
        subimage_fill(obj->fullimage, obj->isoimage);
      obj->profit = profit_init(obj, obj->fullimage, 1,
			MODEL_MOFFAT, PROFIT_NOCONV);
      objlist_movobj(newobjlist, o--, overobjlist);
    }
    objlist_end(newobjlist);
  }

  obj = overobjlist->obj;
  for (o=overobjlist->nobj; o--; obj++)
/*
{
subprofit = obj->profit->subprofit;
profit_residuals(obj->profit, 0.0, obj->profit->paraminit, NULL);
if (subprofit->lmodpix) {
check_add(thecheck, subprofit->lmodpix,
subprofit->objnaxisn[0],subprofit->objnaxisn[1],
subprofit->ix,subprofit->iy, 1.0);
}
*/
    profit_end(obj->profit);
/*
}
*/
  free(imagecopy);
  free(imvarcopy);
  free(imwork);

  return overobjlist;
}


/****** objlist_overlap ***************************************************//**
Move objects overlapping a given object from the input list to a new list.
from the input list.
@param[in] objlist	Pointer to objlist
@param[in] fobj		Pointer to object

@author 		E. Bertin (IAP)
@date			10/07/2014
@todo	The selection algorithm is currently very basic and inefficient.
 ***/
objliststruct *objlist_overlap(objliststruct *objlist, objstruct *fobj) {
   objliststruct	*overobjlist;
   int			i,o, blend;

  overobjlist = objlist_new();
  blend = fobj->blend;
// Loop over input objlist; we keep indexing independent of loop counter to
// accommodate changes in the ordering of the list
  for (o=0, i=objlist->nobj; i--; o++)
    if ( objlist->obj[o].blend == blend)
      objlist_movobj(objlist, o--, overobjlist);

  return overobjlist;
}


#ifdef USE_THREADS

/****** pthread_objlist_init **********************************************//**
Setup threads, mutexes and semaphores for multithreaded objlist processing
@param[in] fields	Pointer to an array of image field pointers
@param[in] wfields	Pointer to an array of weight-map field pointers
@param[in] nfield	Number of images
@param[in] nthreads	Number of threads

@author 		E. Bertin (IAP)
@date			09/06/2014
 ***/
void	pthread_objlist_init(fieldstruct **fields, fieldstruct **wfields,
			int nfield, int nthreads)
  {
   int	n, p;

  pthread_fields = fields;
  pthread_wfields = wfields;
  pthread_nfield = nfield;
  pthread_nthreads = nthreads;
  pthread_nobjlist = prefs.obj2_stacksize;
  QMALLOC(pthread_thread, pthread_t, nthreads);
  QCALLOC(pthread_objlist, objliststruct, pthread_nobjlist);
  QPTHREAD_COND_INIT(&pthread_objlistaddcond, NULL);
  QPTHREAD_COND_INIT(&pthread_objsavecond, NULL);
  QPTHREAD_MUTEX_INIT(&pthread_objlistmutex, NULL);
  QPTHREAD_MUTEX_INIT(&pthread_freeobjmutex, NULL);
  QPTHREAD_MUTEX_INIT(&pthread_countobjmutex, NULL);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_objlistaddindex = pthread_objlistprocindex
	= pthread_objlistsaveindex= 0;
  pthread_endflag = 0;

/* Start the measurement/write_to_catalog threads */
  for (p=0; p<nthreads; p++)
    QPTHREAD_CREATE(&pthread_thread[p], &pthread_attr,
		&pthread_objlist_analyse, (void *)p);

  return;
  }


/****** pthread_objlist_end ***********************************************//**
Terminate threads, mutexes and semaphores for multithreaded objlist processing

@author 	E. Bertin (IAP)
@date		09/06/2014
 ***/
void	pthread_objlist_end(void)
  {
   int	p;

  QPTHREAD_MUTEX_LOCK(&pthread_objlistmutex);
/* Call all threads to exit */
  pthread_endflag = 1;
  QPTHREAD_COND_BROADCAST(&pthread_objlistaddcond);
  QPTHREAD_COND_BROADCAST(&pthread_objlistaddcond);
  QPTHREAD_MUTEX_UNLOCK(&pthread_objlistmutex);

  for (p=0; p<pthread_nthreads; p++)
    QPTHREAD_JOIN(pthread_thread[p], NULL);

  QPTHREAD_MUTEX_DESTROY(&pthread_objlistmutex);
  QPTHREAD_MUTEX_DESTROY(&pthread_freeobjmutex);
  QPTHREAD_MUTEX_DESTROY(&pthread_countobjmutex);
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  QPTHREAD_COND_DESTROY(&pthread_objlistaddcond);
  QPTHREAD_COND_DESTROY(&pthread_objsavecond);

  free(pthread_thread);
  free(pthread_objlist);

  return;
  }


/****** pthread_objlist_add **********************************************//**
Add an objlist to the list of objlists that need to be processed.
@param[in] objlist	objlist to be added (note: not a pointer!)

@author 	E. Bertin (IAP)
@date		09/06/2014
 ***/
void	pthread_objlist_add(objliststruct objlist)
  {

  QPTHREAD_MUTEX_LOCK(&pthread_objlistmutex);
  while (pthread_objlistaddindex>=pthread_objlistsaveindex+pthread_nobjlist)
/*-- Wait for stack to flush if limit on the number of stored objs is reached*/
    QPTHREAD_COND_WAIT(&pthread_objsavecond, &pthread_objlistmutex);
  pthread_objlist[pthread_objlistaddindex++%pthread_nobjlist] = objlist;
  QPTHREAD_COND_BROADCAST(&pthread_objlistaddcond);
  QPTHREAD_MUTEX_UNLOCK(&pthread_objlistmutex);

  return;
  }


/****** pthread_objlist_analyse ******************************************//**
Thread that takes care of measuring and saving objlists.
@paran[in] arg	unused (here only for compliancy with POSIX threads)

@author 	E. Bertin (IAP)
@date		09/06/2014
 ***/
void	*pthread_objlist_analyse(void *arg)
  {
   objliststruct	*objlist;
   objstruct		*obj;
   int			o;

  while (1)
    {
    QPTHREAD_MUTEX_LOCK(&pthread_objlistmutex);
/*-- Flush objects for which measurements have been completed */
    while (pthread_objlistsaveindex<pthread_objlistprocindex
		&& (pthread_objlist[pthread_objlistsaveindex].done_flag))
      analyse_end(pthread_fields, pthread_wfields, pthread_nfield,
		&pthread_objlist[pthread_objlistsaveindex++]);
    while (pthread_objlistprocindex>=pthread_objlistaddindex)
/*---- Wait for more objects to be pushed in stack */
      {
      if ((pthread_endflag))
        {
        QPTHREAD_MUTEX_UNLOCK(&pthread_objlistmutex);
        pthread_exit(NULL);
        }
      QPTHREAD_COND_WAIT(&pthread_objlistaddcond, &pthread_objlistmutex);
      }
    objlist = &pthread_objlist[pthread_objlistprocindex++%pthread_nobjlist];
    QPTHREAD_MUTEX_UNLOCK(&pthread_objlistmutex);
    obj = objlist->obj;
    for (o=objlist->nobj; o--; obj++)
      analyse_full(pthread_fields, pthread_wfields, pthread_nfield, obj);
/*-- Flag groups as done */
    QPTHREAD_MUTEX_LOCK(&pthread_objlistmutex);
    objlist->done_flag = 1;
    QPTHREAD_MUTEX_UNLOCK(&pthread_objlistmutex);
    }

  return (void *)NULL;
  }

#endif



