/*
*				subimage.c
*
* Manage subimage structures.
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
*	Last modified:		02/04/2012
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
#include	"field.h"
#include	"fitswcs.h"
#include	"image.h"
#include	"subimage.h"

/****** subimage_getall ******************************************************
PROTO	subimagestruct *subimage_getall(fieldstruct **fields,
			fieldstruct **wfields, int nfield, obj2struct *obj2)
PURPOSE	Create array of object sub-images (copy of local image data).
INPUT	Pointer to an array of image fields,
	pointer to an array of weight-map fields,
	number of fields,
	pointer to obj2 "object",
OUTPUT	Pointer to a new malloc'ed subimage structure.
NOTES	Global preferences are used. Input fields must have been through
	frame_wcs() with detection field as a reference.
AUTHOR	E. Bertin (IAP)
VERSION	02/04/2012
 ***/
subimagestruct	*subimage_getall(fieldstruct **fields, fieldstruct **wfields,
			int nfield, obj2struct *obj2)

  {
   subimagestruct	*subimage;
   fieldstruct		*field,*wfield;
   int			f,s, nsubimage;

  if (prefs.multigrids_flag)
    {
    nsubimage = 0;
    for (f=0; f<nfield; f++)
      {
      field = fields[f];
      if (!(field->flags&MULTIGRID_FIELD)
		|| (field->wcs->outmin[0] <= obj2->xmax
			&& field->wcs->outmax[0] >= obj2->xmin
			&& field->wcs->outmin[1] <= obj2->ymax
			&& field->wcs->outmax[1] >= obj2->ymin))
          nsubimage++;
      }
    obj2->nsubimage = nsubimage;
    }
  else
    obj2->nsubimage = nsubimage = nfield;

  QMALLOC(obj2->subimage, subimagestruct, obj2->nsubimage);

  subimage = obj2->subimage;
/* Initialize first (detection) sub-image */
  subimage_init(subimage, fields[0], wfields? wfields[0]:NULL, obj2, NULL);
/* Initialize remaining sub-images */
  for (f=1; f<nfield; f++)
    {
    field = fields[f];
    wfield = wfields? wfields[f] : NULL;
    if (field->flags&MULTIGRID_FIELD
	&& (field->wcs->outmin[0] > obj2->xmax
		|| field->wcs->outmax[0] < obj2->xmin
		|| field->wcs->outmin[1] > obj2->ymax
		|| field->wcs->outmax[1] < obj2->ymin))
      continue;
    subimage_init(++subimage, field, wfield, obj2, obj2->subimage);
    }

  subimage = obj2->subimage;
  for (s=nsubimage; s--; subimage++)
    {
    QMALLOC(subimage->image, PIXTYPE, subimage->imsize[0]*subimage->imsize[1]);
    copyimage(subimage->field, subimage->image,
	subimage->imsize[0],subimage->imsize[1],
	subimage->ipos[0],subimage->ipos[1]);
    if (subimage->wfield)
      {
      QMALLOC(subimage->weight, PIXTYPE,
		subimage->imsize[0]*subimage->imsize[1]);
      copyimage(subimage->wfield, subimage->weight,
	subimage->imsize[0],subimage->imsize[1],
	subimage->ipos[0],subimage->ipos[1]);
      }
    else
      subimage->weight = NULL;
    }

  return obj2->subimage;
  }


/****** subimage_init ********************************************************
PROTO	void subimage_init(subimagestruct *subimage,
		fieldstruct *field, fieldstruct *wfield,
		obj2struct *obj2, subimagestruct *dsubimage)
PURPOSE	Initialize a new sub-image (copy of local image data).
INPUT	Pointer to subimage,
	pointer to current image field,
	pointer to current weight-map field (or NULL if no weighting),
	pointer to obj2 "object",
	pointer to reference sub-image (or NULL to make one).
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/04/2012
 ***/
void	subimage_init(subimagestruct *subimage,
			fieldstruct *field, fieldstruct *wfield,
			obj2struct *obj2, subimagestruct *dsubimage)

  {
   fieldstruct	*dfield;
   double	dval,dmax, det;

  if (dsubimage)
    {
    dfield = dsubimage->field;
    if (field->flags&MULTIGRID_FIELD)
      {
      det = wcs_rawtoraw(field->wcs, dfield->wcs,
			dsubimage->dpos, subimage->dpos, subimage->djacob);
      subimage->dscale = sqrt(det);
      if (fabs(det) > 0.0 )
        det = 1.0/det;
      subimage->dinvjacob[0] = det*subimage->djacob[3];
      subimage->dinvjacob[1] = -det*subimage->djacob[1];
      subimage->dinvjacob[2] = -det*subimage->djacob[2];
      subimage->dinvjacob[3] = det*subimage->djacob[0];
      subimage->ipos[0] = (int)(subimage->dpos[0]-0.50001);/* Integer coords */
      subimage->ipos[1] = (int)(subimage->dpos[1]-0.50001);/* Integer coords */
      dmax = fabs(subimage->djacob[0]*dsubimage->imsize[0]);
      if ((dval=fabs(subimage->djacob[1]*dsubimage->imsize[1]))>dmax)
        dmax = dval;
      subimage->imsize[0] = (int)(dmax+0.49999);
      dmax = fabs(subimage->djacob[2]*dsubimage->imsize[0]);
      if ((dval=fabs(subimage->djacob[3]*dsubimage->imsize[1]))>dmax)
        dmax = dval;
      subimage->imsize[1] = (int)(dmax+0.49999);
      subimage->immin[0] = subimage->ipos[0] - subimage->imsize[0]/2;
      subimage->immin[1] = subimage->ipos[1] - subimage->imsize[1]/2;
      subimage->immax[0] = subimage->immin[0] + subimage->imsize[0];
      subimage->immax[1] = subimage->immin[1] + subimage->imsize[1];
      }
    else
      *subimage = *dsubimage;
    }
  else
    {
    subimage->dscale = 1.0;
    subimage->dpos[0] = obj2->mx + 1.0;
    subimage->dpos[1] = obj2->my + 1.0;
    subimage->ipos[0] = (int)(subimage->dpos[0]-0.50001); /* Integer coords */
    subimage->ipos[1] = (int)(subimage->dpos[1]-0.50001); /* Integer coords */
    subimage->djacob[0] = subimage->djacob[2]
		= subimage->dinvjacob[0] = subimage->dinvjacob[2] = 1.0;
    subimage->djacob[1] = subimage->djacob[3]
		= subimage->dinvjacob[1] = subimage->dinvjacob[3] = 0.0;
    subimage->imsize[0] = 2.0*(obj2->xmax-obj2->xmin)+1+2*field->stripmargin;
    subimage->imsize[1] = 2.0*(obj2->ymax-obj2->ymin)+1+2*field->stripmargin;
    subimage->immin[0] = subimage->ipos[0] - subimage->imsize[0]/2;
    subimage->immin[1] = subimage->ipos[1] - subimage->imsize[1]/2;
    subimage->immax[0] = subimage->immin[0] + subimage->imsize[0];
    subimage->immax[1] = subimage->immin[1] + subimage->imsize[1];
    }

  subimage->field = field;
  subimage->wfield = wfield;

  return;
  }


/****** subimage_endall ******************************************************
PROTO	void subimage_endall(obj2struct *obj2)
PURPOSE	Free resources related to an array of sub-images.
INPUT	Pointer to an obj2 "object" structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	29/03/2012
 ***/
void	subimage_endall(obj2struct *obj2)

  {
   subimagestruct	*subimage;
   int			s;

  subimage = obj2->subimage;
  for (s=obj2->nsubimage; s--; subimage++)
    {
    free(subimage->image);
    free(subimage->weight);
    }

  QFREE(obj2->subimage);
  obj2->nsubimage = 0;

  return;
  }



