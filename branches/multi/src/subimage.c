/**
* @file		subimage.c
* @brief	Manage sub-images
* @date		21/10/2014
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2014 IAP/CNRS/UPMC
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
#include	"plist.h"
#include	"subimage.h"

/****** subimage_fromplist ************************************************//**
Create a sub-image from an object and a pixel list, without any margin.
@param[in] field	Pointer to the image field.
@param[in] wfield	Pointer to the weight field.
@param[in] obj		Pointer to the parent object.
@param[in] plist	Pixel list.
@param[out] 		Pointer to a new subimage.

@author 		E. Bertin (IAP)
@date			21/10/2014
 ***/
subimagestruct	*subimage_fromplist(fieldstruct *field, fieldstruct *wfield,
			objstruct *obj, pliststruct *plist) {

   subimagestruct	*subimage;
   pliststruct	*plistt;
   PIXTYPE	*pix, *fpix, *wpix, *pixt;
   long		pos, npix;
   int		i, n, xmin,ymin, w;

/* Initialize sub-image */
  QCALLOC(subimage, subimagestruct, 1);
  subimage->dscale = 1.0;
  subimage->xmin[0] = xmin = obj->xmin;
  subimage->xmin[1] = ymin = obj->ymin;
  subimage->xmax[0] = obj->xmax + 1;
  subimage->xmax[1] = obj->ymax + 1;
  subimage->size[0] = w = subimage->xmax[0] - subimage->xmin[0];
  subimage->size[1] = subimage->xmax[1] - subimage->xmin[1];
  subimage->ipos[0] = subimage->xmin[0] + subimage->size[0]/2;
  subimage->ipos[1] = subimage->xmin[1] + subimage->size[1]/2;
  subimage->field = field;
  subimage->wfield = wfield;
  subimage->bkg = 0.0;

  npix = subimage->size[0]*subimage->size[1];

  if (!(subimage->image = (PIXTYPE *)malloc(npix*sizeof(PIXTYPE)))) {
    free(subimage);
    return NULL;
  }

  pix = subimage->image;
  for (i=npix; i--;)
    *(pix++) = -BIG;

  if (prefs.filter_flag) {
    if (!(subimage->fimage = (PIXTYPE *)malloc(npix*sizeof(PIXTYPE)))) {
      free(subimage->image);
      free(subimage);
      return NULL;
    }

    pix = subimage->fimage;
    for (i=npix; i--;)
      *(pix++) = -BIG;
  }

  if (subimage->wfield) {
    if (!(subimage->imvar = (PIXTYPE *)malloc(npix*sizeof(PIXTYPE)))) {
      free(subimage->image);
      free(subimage->fimage);
      free(subimage);
      return NULL;
    }
    wpix = subimage->imvar;
    for (i=npix; i--;)
      *(wpix++) = BIG;
  }

  pix = subimage->image;
  fpix = subimage->fimage;
  wpix = subimage->imvar;
  for (i=obj->firstpix; i!=-1; i=PLIST(plistt,nextpix)) {
    plistt = plist+i;
    pos = (PLIST(plistt,x)-xmin) + (PLIST(plistt,y)-ymin)*w;
    *(pix+pos) = PLISTPIX(plistt, value);
    if (fpix)
      *(fpix+pos) = PLISTPIX(plistt, cvalue);
    if (wpix)
      *(wpix+pos) = PLISTPIX(plistt,var);
  }

  return subimage;
}


/****** subimage_fromfield ************************************************//**
Create a sub-image from an image field, without any margin.
@param[in] field	Pointer to the image field.
@param[in] wfield	Pointer to the weight field.
@param[in] xmin		Lower frame limit in x (included).
@param[in] xmax		Upper frame limit in x (excluded).
@param[in] ymin		Lower frame limit in y (included).
@param[in] ymax		Upper frame limit in y (excluded).
@param[out] 		Pointer to a new subimage.

@author 		E. Bertin (IAP)
@date			11/09/2014
 ***/
subimagestruct	*subimage_fromfield(fieldstruct *field, fieldstruct *wfield,
			int xmin, int xmax, int ymin, int ymax) {

   subimagestruct	*subimage;
   PIXTYPE	*pix,*wpix;
   long		npix, pos;
   int		i, n, w;

/* Initialize sub-image */
  QCALLOC(subimage, subimagestruct, 1);
  subimage->dscale = 1.0;
  subimage->xmin[0] = xmin;
  subimage->xmin[1] = ymin;
  subimage->xmax[0] = xmax;
  subimage->xmax[1] = ymax ;
  subimage->size[0] = w = xmax - xmin;
  subimage->size[1] = ymax - ymin;
  subimage->ipos[0] = subimage->xmin[0] + subimage->size[0]/2;
  subimage->ipos[1] = subimage->xmin[1] + subimage->size[1]/2;
  subimage->field = field;
  subimage->wfield = wfield;
  subimage->bkg = 0.0;

  npix = subimage->size[0]*subimage->size[1];

  if (!(subimage->image = (PIXTYPE *)calloc(npix, sizeof(PIXTYPE)))) {
    free(subimage);
    return NULL;
  }

  copyimage(field, subimage->image, subimage->size[0],subimage->size[1],
	subimage->ipos[0],subimage->ipos[1], -BIG);

  if (subimage->wfield) {
    if (!(subimage->imvar = (PIXTYPE *)malloc(npix*sizeof(PIXTYPE)))) {
      free(subimage->image);
      free(subimage);
      return NULL;
    }
    copyimage(wfield, subimage->imvar, subimage->size[0],subimage->size[1],
	subimage->ipos[0],subimage->ipos[1], BIG);
  }

  return subimage;
}


/****** subimage_fill ****************************************************//**
Replace pixels of the first subimage depending on the valid pixels from the
second subimage.
@param[in] subimage	Pointer to the first (destination) subimage.
@param[in] submask	Pointer to the second subimage.

@author 		E. Bertin (IAP)
@date			17/10/2014
 ***/
void	subimage_fill(subimagestruct *subimage, subimagestruct *submask,
			subimage_fillenum	fill_type) {

   int		x,y, xmin,xmax,ymin,ymax, w,dwima,dwmask;
   PIXTYPE	*pixima,*pixmask,
		val;

/* Don't go further if out of frame!! */
  if (subimage->xmax[0] <= submask->xmin[0]
	|| subimage->xmin[0] >= submask->xmax[0]
	|| subimage->xmax[1] <= submask->xmin[1]
	|| subimage->xmin[1] >= submask->xmax[1])
    return;

/* Set the image boundaries */

  pixmask = submask->image;
  xmin = submask->xmin[0] - subimage->xmin[0];
  xmax = submask->xmax[0] - subimage->xmin[0];

  if (xmin < 0) {
    pixmask -= xmin;
    xmin = 0;
  }
  if (xmax > subimage->size[0])
    xmax = subimage->size[0];

  ymin = submask->xmin[1] - subimage->xmin[1];
  ymax = submask->xmax[1] - subimage->xmin[1];
  if (ymin < 0) {
    pixmask -= ymin * submask->size[0];
    ymin = 0;
  }
  if (ymax > subimage->size[1])
    ymax = subimage->size[1];

/* Copy the right pixels to the destination */
  w = xmax - xmin;
  dwmask = submask->size[0] - w;
  dwima = subimage->size[0] - w;
  pixima = subimage->image + ymin*subimage->size[0] + xmin;
  for (y=ymax-ymin; y--; pixmask += dwmask, pixima += dwima)
    for (x=w; x--; pixima++)
      if ((val = *(pixmask++)) > -BIG)
        *pixima = (fill_type==SUBIMAGE_FILL_INPUT) ? val : -BIG;

  return;
}


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
VERSION	11/09/2014
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

  QCALLOC(obj2->subimage, subimagestruct, obj2->nsubimage);

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
  for (s=0; s<nsubimage; s++, subimage++)
    {
    QMALLOC(subimage->image, PIXTYPE, subimage->size[0]*subimage->size[1]);
    copyimage(subimage->field, subimage->image,
	subimage->size[0],subimage->size[1],
	subimage->ipos[0],subimage->ipos[1], -BIG);
    if (subimage->wfield)
      {
      QMALLOC(subimage->imvar, PIXTYPE,
		subimage->size[0]*subimage->size[1]);
      copyimage(subimage->wfield, subimage->imvar,
	subimage->size[0],subimage->size[1],
	subimage->ipos[0],subimage->ipos[1], BIG);
      }

    subimage->bkg = obj2->bkg[s];
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
VERSION	07/05/2012
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
      dmax = fabs(subimage->djacob[0]*dsubimage->size[0]);
      if ((dval=fabs(subimage->djacob[1]*dsubimage->size[1]))>dmax)
        dmax = dval;
      subimage->size[0] = (int)(dmax+0.49999);
      dmax = fabs(subimage->djacob[2]*dsubimage->size[0]);
      if ((dval=fabs(subimage->djacob[3]*dsubimage->size[1]))>dmax)
        dmax = dval;
      subimage->size[1] = (int)(dmax+0.49999);
      subimage->xmin[0] = subimage->ipos[0] - subimage->size[0]/2;
      subimage->xmin[1] = subimage->ipos[1] - subimage->size[1]/2;
      subimage->xmax[0] = subimage->xmin[0] + subimage->size[0];
      subimage->xmax[1] = subimage->xmin[1] + subimage->size[1];
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
		= subimage->dinvjacob[0] = subimage->dinvjacob[3] = 1.0;
    subimage->djacob[1] = subimage->djacob[2]
		= subimage->dinvjacob[1] = subimage->dinvjacob[2] = 0.0;
    subimage->size[0] = 3.0*(obj2->xmax-obj2->xmin)+1+2*field->stripmargin;
    subimage->size[1] = 3.0*(obj2->ymax-obj2->ymin)+1+2*field->stripmargin;
    subimage->xmin[0] = subimage->ipos[0] - subimage->size[0]/2;
    subimage->xmin[1] = subimage->ipos[1] - subimage->size[1]/2;
    subimage->xmax[0] = subimage->xmin[0] + subimage->size[0];
    subimage->xmax[1] = subimage->xmin[1] + subimage->size[1];
    }

  subimage->field = field;
  subimage->wfield = wfield;

  return;
  }


/****** subimage_end ********************************************************
PROTO	void subimage_endall(subimagestruct *subimage)
PURPOSE	Free resources related to asub-image.
INPUT	Pointer to a subimage structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/09/2014
 ***/
void	subimage_end(subimagestruct *subimage)

  {
  free(subimage->image);
  free(subimage->imvar);
  free(subimage->fimage);
  free(subimage->fimvar);

  return;
  }


/****** subimage_endall ******************************************************
PROTO	void subimage_endall(obj2struct *obj2)
PURPOSE	Free resources related to an array of sub-images.
INPUT	Pointer to an obj2 "object" structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2014
 ***/
void	subimage_endall(obj2struct *obj2)

  {
   subimagestruct	*subimage;
   int			s;

  subimage = obj2->subimage;
  for (s=obj2->nsubimage; s--; subimage++)
    subimage_end(subimage);

  QFREE(obj2->subimage);
  obj2->nsubimage = 0;

  return;
  }



