/*
*				clean.h
*
* Include file for clean.c.
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

#ifndef _OBJLIST_H_
#include "objlist.h"
#endif

/*------------------------------ definitions --------------------------------*/

#define		CLEAN_ZONE		10.0	/* zone (in sigma) to */
						/* consider for processing */

/*------------------------------- functions ---------------------------------*/

objliststruct	*clean_init(void);


extern int	clean_process(objliststruct *cleanobjlist, fieldstruct *field,
			objstruct *objin);

extern void	clean_add(objliststruct *cleanobjlist, objstruct *objin),
		clean_end(objliststruct *cleanobjlist),
		clean_merge(objstruct *objin, objstruct *objout),
		clean_sub(objliststruct *cleanobjlist, int objindex);


