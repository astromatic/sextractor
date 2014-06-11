/*
*				lutz.h
*
* Include file for lutz.c.
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
*	Last modified:		11/06/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _OBJLIST_H_
#include "objlist.h"
#endif

#ifndef _LUTZ_H_
#define _LUTZ_H_

/*------------------------------ definitions --------------------------------*/

#define	NSUBOBJ_START		256		/// Nb of subobjects at start
#define	UNKNOWN			-1		/// Lutz algorithm flag code

/*--------------------------------- typedefs --------------------------------*/

typedef	enum		{COMPLETE, INCOMPLETE, NONOBJECT, OBJECT}
				status;	/* Extraction status */

/*--------------------------------- variables -------------------------------*/
PIXTYPE		*dumscan;

/*------------------------------- structures --------------------------------*/
 /// Temporary detection parameters during extraction
typedef struct structinfo
  {
  LONG		pixnb;			/// Nb of pixels included in detection
  LONG		firstpix;		/// Pointer to first pixel of pixel list
  LONG		lastpix;		/// Pointer to last pixel of pixel list
  short		flag;			/// Extraction flag
  }       infostruct;

/*------------------------------- functions ---------------------------------*/
objliststruct	*lutz_subextract(subimagestruct *subimage, PIXTYPE thresh,
			int xmin, int xmax, int ymin, int ymax);

int		lutz_output(infostruct *info, objliststruct *objlist);

void		lutz_update(infostruct *infoptr1, infostruct *infoptr2,
			pliststruct *pixel);

#endif
