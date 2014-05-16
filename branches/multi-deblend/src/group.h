/**
* @file		group.h
* @brief	Include file for group.c.
* @date		28/03/2014
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

#ifndef _GROUP_H_
#define _GROUP_H_

/*----------------------------- Internal constants --------------------------*/

#define	GROUP_NDEBLENDITER	10	/** number of deblending iterations */
#define	GROUP_NMULTITER		10	/** number of multi-model iterations */

/*--------------------------------- typedefs --------------------------------*/
/*------------------------------ Prototypes ---------------------------------*/

obj2struct	*analyse_obj2obj2(fieldstruct **fields, fieldstruct **wfields,
			int nfield, objstruct *obj, obj2liststruct *obj2list);
int		analyse_full(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2struct *obj2),
		analyse_overlapness(objliststruct *objlist, int iobj);

void		analyse_end(fieldstruct **fields, fieldstruct **wfields,
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

