/**
* @file		objlist.h
* @brief	Include file for objlist.c.
* @date		05/06/2014
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2011-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
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

#ifndef _SUBIMAGE_H_
#include "subimage.h"
#endif

#ifndef _PLIST_H_
#include "plist.h"
#endif

#ifndef _OBJLIST_H_
#define _OBJLIST_H_

/*----------------------------- Internal constants --------------------------*/

#define	OBJLIST_NOBJMAXINC	16	/// memory allocation increment */
#define	GROUP_NDEBLENDITER	10	/// number of deblending iterations
#define	GROUP_NMULTITER		10	/// number of multi-model iterations

/*--------------------------------- typedefs --------------------------------*/
typedef struct
  {
  int		nobj;			/* number of objects in list */
  int		nobjmax;		/* number of allocated slots in list */
  objstruct	*obj;			/* pointer to the object array */
  int		npix;			/* number of pixels in pixel-list */
  struct subimage	*subimage;	/* Array of sub-images */
  pliststruct	*plist;			/* pointer to the pixel-list */
  PIXTYPE	dthresh;		/* detection threshold */
  PIXTYPE	thresh;			/* analysis threshold */
  }	objliststruct;

/*------------------------------ Prototypes ---------------------------------*/

objliststruct	*objlist_new(void);

obj2struct	*analyse_obj2obj2(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objstruct *obj, obj2liststruct *obj2list);
int		objlist_addobj(objliststruct *objlist, objstruct *obj),
		objlist_movobj(objliststruct *objlistin, int objindex,
			objliststruct *objlistout),
		objlist_subobj(objliststruct *objlist, int objindex),
		analyse_full(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2struct *obj2),
		analyse_overlapness(objliststruct *objlist, int iobj);

void		objlist_end(objliststruct *objlist),
		analyse_end(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2groupstruct *group2),
		analyse_final(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist, int iobj),
		analyse_group(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2groupstruct *group2),
		analyse_iso(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist, int n);

#ifdef	USE_THREADS
void		*pthread_analyse_obj2group(void *arg),
		pthread_add_obj2group(obj2groupstruct group2),
		pthread_end_obj2group(void),
		pthread_init_obj2group(fieldstruct **fields,
			fieldstruct **wfields, int nfield, int nthreads);

pthread_mutex_t	pthread_countobj2mutex;
#endif

#endif

