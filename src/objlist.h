/**
* @file		objlist.h
* @brief	Include file for objlist.c.
* @date		25/06/2014
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

#ifndef _OBJLIST_H_
#define _OBJLIST_H_

/*----------------------------- Internal constants --------------------------*/

#define	OBJLIST_NOBJMAXINC	16	/// memory allocation increment */
#define	GROUP_NDEBLENDITER	3	/// number of deblending iterations
#define	GROUP_NMULTITER		7	/// number of multi-model iterations

/*--------------------------------- typedefs --------------------------------*/
typedef struct
  {
  objstruct	*obj;			/// pointer to the object array
  int		nobj;			/// number of objects in list
  int		nobjmax;		/// number of allocated slots in list
  int		*objindex;		/// pointer to the object index array
  pliststruct	*plist;			/// pointer to the pixel-list
  int		npix;			/// number of pixels in pixel-list
  struct subimage	*subimage;	/// Array of sub-images
  int		done_flag;		/// Set if objs are ready to save
  }	objliststruct;

/*------------------------------ Prototypes ---------------------------------*/

objliststruct	*objlist_deblend(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objliststruct *objlist, int objindex),
		*objlist_new(void),
		*objlist_overlap(objliststruct *objlist, objstruct *obj);

int		objlist_addobj(objliststruct *objlist, objstruct *obj,
			pliststruct *plist),
		objlist_movobj(objliststruct *objlistin, int objindex,
			objliststruct *objlistout),
		objlist_subobj(objliststruct *objlist, int objindex);

void		obj_end(objstruct *obj),
		objlist_end(objliststruct *objlist);

#ifdef	USE_THREADS
void		*pthread_objlist_analyse(void *arg),
		pthread_objlist_add(objliststruct objlist),
		pthread_objlist_end(void),
		pthread_objlist_init(fieldstruct **fields,
			fieldstruct **wfields, int nfield, int nthreads);

pthread_mutex_t	pthread_countobjmutex;
#endif

#endif

