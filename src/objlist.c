/**
* @file		objlist.c
* @brief	Manage groups of detection (advanced deblending)
* @date		04/06/2014
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

pthread_t	*pthread_thread;
pthread_attr_t	pthread_attr;
pthread_mutex_t	pthread_groupmutex, pthread_freeobjmutex;
pthread_cond_t	pthread_groupaddcond, pthread_objsavecond;
fieldstruct	**pthread_fields,**pthread_wfields;
objgroupstruct	*pthread_group;
int		pthread_ngroup, pthread_groupaddindex,pthread_groupprocindex,
		pthread_groupsaveindex, pthread_nfield, pthread_nthreads,
		pthread_endflag;

#endif

/****** objlist_new *******************************************************//**
Create a new, empty objlist.
@param[out] 		Pointer to a new empty objlist.

@author 		E. Bertin (IAP)
@date			16/05/2014
 ***/
objliststruct	*objlist_new(void) {

  QCALLOC(objlist, objliststruct, 1);

  return objlist;
}


/****** objlist_end *******************************************************//**
Free resources allocated to an objlist
@param[in] objlist	Pointer to the objlist.

@author 		E. Bertin (IAP)
@date			04/06/2014
 ***/
void	objlist_end(objliststruct *objlist) {

  free(objlist->obj);
  free(objlist->plist);
  free(objlist->subimage);
  free(objlist);
}


/****** objlist_addobj ****************************************************//**
Add an object to an objlist.
@param[in] objlist	Pointer to the objlist
@param[in] wfields	Pointer to the object to be added
@param[out] 		RETURN_FATAL_ERROR if a memory (re-)allocation issue
			happened,
			RETURN_OK otherwise.

@author 		E. Bertin (IAP)
@date			16/05/2014
 ***/
int	objlist_addobj(objliststruct *objlist, *objstruct obj) {

// Check if memory (re-)allocation is necessary.
  if (objlist->nobj >= objlist->nobjmax) {
    objlist->nobjmax += OBJLIST_NOBJMAXINC;
    if (objlist->nobjmax && !(objlist->obj = (objstruct *)realloc(objlist->obj,
		objlist->nobjmax * sizeof(objstruct))))
      return RETURN_FATAL_ERROR;
    else if (!(objlist->obj = (objstruct *)malloc(objlist->nobjmax
	 * sizeof(objstruct))))
      return RETURN_FATAL_ERROR;
  }

// Copy the content of the input object
  objlist->obj[objlist->nobj] = *obj;
  objlist->nobj++;

  return RETURN_OK
}


/****** objlist_subobj ****************************************************//**
Remove an object from an objlist.
@param[in] objlist	Pointer to the objlist
@param[in] objindex	Index of the object to be removed
@param[out] 		RETURN_ERROR if the objlist is empty,
			RETURN_FATAL_ERROR if a memory (re-)allocation issue
			happened,
			RETURN_OK otherwise.

@author 		E. Bertin (IAP)
@date			16/05/2014
 ***/
int	objlist_subobj(objliststruct *objlist, int objindex) {
  if (--objlist->nobj) {
    if (objlist->nobj != objindex)
      objlist->obj[objindex] = objlist->obj[objlist->nobj];
    nobjmax = objlist->nobjmax;
    while (nobjmax - OBJLIST_NOBJMAXINC > objlist->nobj)
      nobjmax -= OBJLIST_NOBJMAXINC;
    if (nobjmax != objlist->nobjmax
	&& !(objlist->obj = (objstruct *)realloc(objlist->obj,
	    (objlist->nobjmax = nobjmax) * sizeof(objstruct))))
      return RETURN_ERROR;
  } else {
    QFREE(objlist->obj);
    objlist->nobjmax = 0;
  }

  return RETURN_OK
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

@author 		E. Bertin (IAP)
@date			05/06/2014
 ***/
int	objlist_movobj(objliststruct *objlistin, int objindex,
		objliststruct *objlistout) {

  if (objlist_addobj(objlistout, objlistin->obj[objindex]) == RETURN_OK)
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

@author 		E. Bertin (IAP)
@date			05/06/2014
 ***/
void	objlist_deblend(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist, int objindex) {
   fieldstruct		*field, *wfield;
   subprofitstruct	*subprofit,*modsubprofit;
   subimagestruct	*subimage;
   objstruct		*obj, *modobj, *fobj;
   int			i,o,s;

  field = fields[0];			// field is the detection field
  wfield = wfields? wfields[0]:NULL;	// wfield is the detection weight map

// Find overlapping detections and link them
  overobjlist = objlist_overlap(objlist, objlist->obj[objindex]);

  xmax = ymax = -(xmin = ymin = 0x7FFFFFFF);	// largest signed 32-bit int
  obj = overobjlist->obj;
  blend = obj->blend;				// Keep note of the blend index
  for (o=overobjlist->nobj; o--; obj++) {	// find boundaries of the group
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
// Extract subimage covering the whole group
  overobjlist->subimage = subimage_fromfield(field, wfield,
				xmin, xmax, ymin, ymax);
  obj = overobjlist->obj;
  for (o=overobjlist->nobj; o--; obj++) {
//-- if BLANKing is on, paste back the object pixels in the sub-image
    if (prefs.blank_flag && obj->isoimage) {
      subimage_fill(overobjlist->subimage, obj->isoimage);
  }

  obj = overobjlist->obj;
  for (i=overobjlist->nobj; i--; obj++) {
//-- Create individual object subimages
    obj->fullimage = subimage_fromfield(field, wfield,
		obj->xmin, obj->xmax, obj->ymin, obj->ymax);
//-- if BLANKing is on, paste back the object pixels in the sub-images
    if (prefs.blank_flag && obj->isoimage) {
      subimage_fill(obj->fullimage, obj->isoimage);
    }
    obj->profit = profit_init(obj, obj->fullimage, 1 ,
		MODEL_MOFFAT, PROFIT_NOCONV);
  }

  subimage = overobjlist->subimage;
  for (j=0; j<GROUP_NDEBLENDITER; j++) {
    nobj = overobjlist->nobj;
//-- Iterative multiple fit if several sources overlap
    for (i=0; i<GROUP_NMULTITER; i++) {
      obj = overobjlist->obj;
      for (o=nobj; o--; obj++)
        profit_fit(obj->profit);
      obj = overobjlist->obj;
      for (o=nobj; o--; obj++) {
//------ Subtract the contribution from all overlapping neighbors (models)
        if (i)
          subprofit_copyobjpix(obj->profit->subprofit, obj->fullimage);
        modobj = overobjlist->obj;
        for (o2=nobj; o2--; modobj++)
          if (modobj != obj) {
            subprofit = obj->profit->subprofit;
            modsubprofit = modobj->profit->subprofit;
            subprofit_submodpix(modsubprofit, subprofit->objpix,
			subprofit->ix, subprofit->iy,
			subprofit->objnaxisn[0], subprofit->objnaxisn[1],
			subprofit->subsamp, nobj>1 ? 0.95: 1.0);
          }
      }
      if (nobj <= 1)
        break;
    }

    obj = overobjlist->obj;
    for (o=nobj; o--; obj++)
      subprofit_submodpix(modsubprofit, subimage->image,
			subimage->xmin[0], subimage->xmin[1],
			subimage->size[0], subimage->size[1],
			subprofit->subsamp, nobj>1 ? 0.95: 1.0);
    newobjlist = lutz_subextract(subimage, overobjlist->dthresh, xmin, xmax,
			ymin, ymax);
    obj = newobjlist->obj;
    for (o=0; o<newobjlist->nobj; o++, obj++) {
      scan_preanalyse(newobjlist, o, ANALYSE_FULL|ANALYSE_ROBUST);
      if (prefs.ext_maxarea && obj->fdnpix > prefs.ext_maxarea)
        continue; 
      obj->number = ++thecat.ndetect;
      obj->blend = blend;
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
		obj->xmin, obj->xmax, obj->ymin, obj->ymax);
//---- if BLANKing is on, paste back the object pixels in the sub-images
      if (prefs.blank_flag)
        subimage_fill(obj->fullimage, obj->isoimage);
      objlist_movobj(newobjlist, o, overobjlist);
    }
    objlist_end(newobjlist);
  }


  return overobjlist;

/* Full source analysis and decide if detection should be written to catalogue*/
  for (obj=fobj; obj; obj=obj->nextobj)
    obj->writable_flag =
		(analyse_full(fields, wfields, nfield, obj) == RETURN_OK);
/* Deallocate memory used for model-fitting */
  if (prefs.prof_flag)
    for (obj=fobj; obj; obj=obj->nextobj)
      profit_end(obj->profit);

/* Free the group of objs */
  for (obj=fobj; obj; obj=obj->nextobj)
    subimage_endall(obj);

  subimage_end(group->subimage);

#ifdef USE_THREADS
  if (prefs.nthreads>1) {
/*-- Flag groups as done */
    QPTHREAD_MUTEX_LOCK(&pthread_groupmutex);
    group->done_flag = 1;
    QPTHREAD_MUTEX_UNLOCK(&pthread_groupmutex);
  }
#endif

  return;
}


/****** objlist_overlap ***************************************************//**
Create a list of objects overlapping a given object.
@param[in] objlist	Pointer to objlist
@param[in] fobj		Pointer to object

@author 		E. Bertin (IAP)
@date			22/05/2014
@todo	The selection algorithm is currently very basic and inefficient.
 ***/
objliststruct *objlist_overlap(objliststruct *objlist, objstruct *fobj)
  {
   objlistruct	*overobjlist;
   objstruct	*obj;
   int		i, blend, nobj;

  QMALLOC(overobjlist, 1, objliststruct);
  QMALLOC(overobjlist->obj, ANALYSE_NOVERLAP, objstruct);
  overobjlist = overobjlist_new();
  nobj = objlist->nobj;
  obj = objlist->obj;
  blend = fobj->blend;
  for (i=0; i<nobj; i++, obj++)
    if (obj->blend == blend && obj!=fobj)
      objlist_add(overobjlist, obj);

  return nblend;
  }



#ifdef USE_THREADS

/****** pthread_group_init ***********************************************//**
Setup threads, mutexes and semaphores for multithreaded group processing
@param[in] fields	Pointer to an array of image field pointers
@param[in] wfields	Pointer to an array of weight-map field pointers
@param[in] nfield	Number of images
@param[in] nthreads	Number of threads

@author 		E. Bertin (IAP)
@date			03/01/2014
 ***/
void	pthread_group_init(fieldstruct **fields, fieldstruct **wfields,
			int nfield, int nthreads)
  {
   int	n,p;

  pthread_fields = fields;
  pthread_wfields = wfields;
  pthread_nfield = nfield;
  pthread_nthreads = nthreads;
  pthread_ngroup = prefs.obj_stacksize;
  QMALLOC(pthread_thread, pthread_t, nthreads);
  QCALLOC(pthread_group, objgroupstruct, pthread_ngroup);
  QPTHREAD_COND_INIT(&pthread_groupaddcond, NULL);
  QPTHREAD_COND_INIT(&pthread_objsavecond, NULL);
  QPTHREAD_MUTEX_INIT(&pthread_groupmutex, NULL);
  QPTHREAD_MUTEX_INIT(&pthread_freeobjmutex, NULL);
  QPTHREAD_MUTEX_INIT(&pthread_countobjmutex, NULL);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_groupaddindex = pthread_groupprocindex = pthread_groupsaveindex= 0;
  pthread_endflag = 0;

/* Start the measurement/write_to_catalog threads */
  for (p=0; p<nthreads; p++)
    QPTHREAD_CREATE(&pthread_thread[p], &pthread_attr,
		&pthread_group_analyse, (void *)p);

  return;
  }


/****** pthread_group_end ************************************************//**
Terminate threads, mutexes and semaphores for multithreaded group processing

@author 	E. Bertin (IAP)
@date		03/01/2014
 ***/
void	pthread_group_end(void)
  {
   int	p;

  QPTHREAD_MUTEX_LOCK(&pthread_groupmutex);
/* Call all threads to exit */
  pthread_endflag = 1;
  QPTHREAD_COND_BROADCAST(&pthread_groupaddcond);
  QPTHREAD_COND_BROADCAST(&pthread_groupaddcond);
  QPTHREAD_MUTEX_UNLOCK(&pthread_groupmutex);

  for (p=0; p<pthread_nthreads; p++)
    QPTHREAD_JOIN(pthread_thread[p], NULL);

  QPTHREAD_MUTEX_DESTROY(&pthread_groupmutex);
  QPTHREAD_MUTEX_DESTROY(&pthread_freeobjmutex);
  QPTHREAD_MUTEX_DESTROY(&pthread_countobjmutex);
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  QPTHREAD_COND_DESTROY(&pthread_groupaddcond);
  QPTHREAD_COND_DESTROY(&pthread_objsavecond);

  free(pthread_thread);
  free(pthread_group);

  return;
  }


/****** pthread_group_add ************************************************//**
Add a group to the list of groups that need to be processed.
@param[in] group	group to be added (note: not a pointer!)

@author 	E. Bertin (IAP)
@date		03/01/2014
 ***/
void	pthread_group_add(objgroupstruct group)
  {

  QPTHREAD_MUTEX_LOCK(&pthread_groupmutex);
  while (pthread_groupaddindex>=pthread_groupsaveindex+pthread_ngroup)
/*-- Wait for stack to flush if limit on the number of stored objs is reached*/
    QPTHREAD_COND_WAIT(&pthread_objsavecond, &pthread_groupmutex);
  pthread_group[pthread_groupaddindex++%pthread_ngroup] = group;
  QPTHREAD_COND_BROADCAST(&pthread_groupaddcond);
  QPTHREAD_MUTEX_UNLOCK(&pthread_groupmutex);

  return;
  }


/****** pthread_group_analyse ********************************************//**
Thread that takes care of measuring and saving obj groups.
@paran[in] arg	unused (here only for compliancy with POSIX threads)

@author 	E. Bertin (IAP)
@date		03/01/2014
 ***/
void	*pthread_analyse_objgroup(void *arg)
  {
   objgroupstruct	*group;

  while (1)
    {
    QPTHREAD_MUTEX_LOCK(&pthread_groupmutex);
/*-- Flush objects for which measurements have been completed */
    while (pthread_groupsaveindex<pthread_groupprocindex
		&& (pthread_group[pthread_groupsaveindex].done_flag))
      analyse_end(pthread_fields, pthread_wfields, pthread_nfield,
		&pthread_group[pthread_groupsaveindex++]);
    while (pthread_groupprocindex>=pthread_groupaddindex)
/*---- Wait for more objects to be pushed in stack */
      {
      if ((pthread_endflag))
        {
        QPTHREAD_MUTEX_UNLOCK(&pthread_groupmutex);
        pthread_exit(NULL);
        }
      QPTHREAD_COND_WAIT(&pthread_groupaddcond, &pthread_groupmutex);
      }
    group = &pthread_group[pthread_groupprocindex++%pthread_ngroup];
    QPTHREAD_MUTEX_UNLOCK(&pthread_groupmutex);
    analyse_group(pthread_fields, pthread_wfields, pthread_nfield, group);
    }

  return (void *)NULL;
  }

#endif



