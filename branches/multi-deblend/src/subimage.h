/*
*				subimage.h
*
* Include file for subimage.c.
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
*	Last modified:		09/06/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _SUBIMAGE_H_
#define _SUBIMAGE_H_

#ifndef _FIELD_H_
#include        "field.h"
#endif

/*------------------------------- structures --------------------------------*/
typedef struct subimage
  {
  struct field	*field;		/// pointer to the field hosting local data
  struct field	*wfield;	/// pointer to the field hosting local weights
  PIXTYPE	*image;			/// Pointer to subimage data
  PIXTYPE	*fimage;		/// Pointer to filtered subimage data
  PIXTYPE	*weight;		/// Pointer to subimage weights
  double	dpos[2];		/// Object coordinates (if applicable)
  int		ipos[2];		/// Coordinates of subimage centre
  int		size[2];		/// Subimage dimensions [pixels]
  int		xmin[2];		/// Subimage min x/y coordinates
  int		xmax[2];		/// Subimage max x/y coordinates - 1
  double	djacob[4];		/// Local astrometric Jacobian matrix
  double	dinvjacob[4];		/// Inverse Jacobian matrix
  double	dscale;			/// Local relative pixel scale
  PIXTYPE	bkg;			/// Subimage background level
  }	subimagestruct;

/*------------------------------- functions ---------------------------------*/

subimagestruct	*subimage_fromfield(fieldstruct *field, fieldstruct *wfield,
			int xmin, int xmax, int ymin, int ymax),
		*subimage_fromplist(fieldstruct *field, fieldstruct *wfield,
			objstruct *obj, pliststruct *plist),
		*subimage_getall(fieldstruct **fields, fieldstruct **wfields,
				int nfield, obj2struct *obj2),
		*subimage_new(fieldstruct *field, fieldstruct *wfield,
			int xmin, int xmax, int ymin, int ymax);

void		subimage_end(subimagestruct *subimage),
		subimage_endall(obj2struct *obj2),
		subimage_fill(subimagestruct *subimage,
			subimagestruct *submask),
		subimage_init(subimagestruct *subimage,
				fieldstruct *field, fieldstruct *wfield,
				obj2struct *obj2, subimagestruct *dsubimage);

#endif
