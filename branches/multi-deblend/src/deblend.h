/*
*				deblend.h
*
* Include file for deblend.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2012-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		25/06/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _OBJLIST_H_
#include "objlist.h"
#endif

#ifndef _SUBIMAGE_H_
#include "subimage.h"
#endif


/*------------------------------ definitions --------------------------------*/

#ifndef	RAND_MAX
#define	RAND_MAX		2147483647
#endif
#define	DEBLEND_NSONMAX		1024		/* max. number per level */
#define	DEBLEND_NBRANCH		16		/* starting number per branch */

/*------------------------------- functions ---------------------------------*/

objliststruct	*deblend_parcelout(objstruct *objin, subimagestruct *subimage,
			pliststruct *plist);

int		deblend_belong(int corenb, objliststruct *coreobjlist,
			int shellnb, objliststruct *shellobjlist),
		deblend_gatherup(objliststruct *objlistin,
			objliststruct *objlistout);
	
void		deblend_alloc(void),
		deblend_free(void);

