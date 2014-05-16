/*
*				lutz.h
*
* Include file for lutz.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		11/01/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _LUTZ_H_
#define _LUTZ_H_

/*------------------------------ definitions --------------------------------*/

#define	NOBJ			256		/* starting number of obj. */
#define	UNKNOWN			-1		/* flag for LUTZ */

/*--------------------------------- typedefs --------------------------------*/

typedef	enum		{COMPLETE, INCOMPLETE, NONOBJECT, OBJECT}
				status;	/* Extraction status */

/*--------------------------------- variables -------------------------------*/
PIXTYPE		*dumscan;

/*------------------------------- structures --------------------------------*/
/* Temporary object parameters during extraction */
typedef struct structinfo
  {
  LONG		pixnb;			/* Number of pixels included */
  LONG		firstpix;		/* Pointer to first pixel of pixlist */
  LONG		lastpix;		/* Pointer to last pixel of pixlist */
  short		flag;			/* Extraction flag */
  }       infostruct;


/*------------------------------- functions ---------------------------------*/
void		lutz_alloc(int width, int height),
		lutz_free(void),
		lutz_output(infostruct *info, objliststruct *objlist),
		lutz_update(infostruct *infoptr1, infostruct *infoptr2,
			pliststruct *pixel);

int		lutz_subextract(objliststruct *objlistroot, int nroot,
			objstruct *objparent, objliststruct *objlist);

#endif
