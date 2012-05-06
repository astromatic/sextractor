/*
*				subimage.h
*
* Include file for subimage.c.
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
*	Last modified:		03/04/2012
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
  struct field	*field;	/* pointer to the image field that hosts local data */
  struct field	*wfield;/* pointer to the field that hosts local weights */
  PIXTYPE	*image;				/* Copy of local image data */
  PIXTYPE	*weight;			/* Copy of local weight data */
  double	dpos[2];		/* Local object coordinates of centre */
  int		ipos[2];		/* Local coordinates of image centre */
  int		imsize[2];			/* Local image data size */
  int		immin[2];			/* Local image data min x/y */
  int		immax[2];			/* Local image data max x/y */
  double	djacob[4];			/* Local image Jacobian matrix*/
  double	dinvjacob[4];			/* Inverse Jacobian matrix */
  double	dscale;				/* Local relative pixel scale */
  PIXTYPE	bkg;				/* Local background level */
  }	subimagestruct;

/*------------------------------- functions ---------------------------------*/

subimagestruct	*subimage_getall(fieldstruct **fields, fieldstruct **wfields,
				int nfield, obj2struct *obj2);

void		subimage_endall(obj2struct *obj2),
		subimage_init(subimagestruct *subimage,
				fieldstruct *field, fieldstruct *wfield,
				obj2struct *obj2, subimagestruct *dsubimage);

#endif
