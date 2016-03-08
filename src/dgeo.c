/**
* @file         dgeo.c
* @brief        Manage differential geometry maps (to correct pixel positions)
* @date         12/02/2015
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       This file part of:      SExtractor
*
*       Copyright:              (C) 1993-2015 IAP/CNRS/UPMC
*
*       License:                GNU General Public License
*
*       SExtractor is free software: you can redistribute it and/or modify
*       it under the terms of the GNU General Public License as published by
*       the Free Software Foundation, either version 3 of the License, or
*       (at your option) any later version.
*       SExtractor is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*       GNU General Public License for more details.
*       You should have received a copy of the GNU General Public License
*       along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
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
#include	"dgeo.h"
#include	"field.h"

/****** dgeo_copy ************************************************************
PROTO	int dgeo_copy(picstruct *dgeofield, PIXTYPE *destx, PIXTYPE *desty,
		int w,int h, int ix,int iy)
PURPOSE	Copy a small part of the differential geometry map components to output
	rasters.
INPUT	Pointer to the dgeofield structure,
	Pointer to the output x map component,
	Pointer to the output y map component,
	width in pixels,
	height in pixels,
	offset in x,
	offset in y.
OUTPUT	RETURN_OK if successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/12/2015
 ***/
int	dgeo_copy(picstruct *dgeofield, PIXTYPE *destx, PIXTYPE *desty,
		int w,int h, int ix,int iy)
  {
   int		i,y, xmin,xmax,ymin,ymax,w2;

/* First set the output to 0 */
  if (destx)
    memset(destx, 0, w*h*sizeof(PIXTYPE));
  if (desty)
    memset(desty, 0, w*h*sizeof(PIXTYPE));

/* Don't go further if out of frame!! */
  if (!dgeofield || ix<0 || ix>=dgeofield->width
	|| iy<dgeofield->ymin || iy>=dgeofield->ymax)
    return RETURN_ERROR;

/* Set the image boundaries */
  w2 = w;
  ymin = iy-h/2;
  ymax = ymin + h;
  if (ymin<dgeofield->ymin)
    {
    if (destx)
      destx += (dgeofield->ymin-ymin)*w;
    if (desty)
      desty += (dgeofield->ymin-ymin)*w;
    ymin = dgeofield->ymin;
    }
  if (ymax>dgeofield->ymax)
    ymax = dgeofield->ymax;

  xmin = ix-w/2;
  xmax = xmin + w;
  if (xmax>dgeofield->width)
    {
    w2 -= xmax-dgeofield->width;
    xmax = dgeofield->width;
    }
  if (xmin<0)
    {
    if (destx)
      destx += -xmin;
    if (desty)
      desty += -xmin;
    w2 -= -xmin;
    xmin = 0;
    }

/* Copy the right pixels to the destination */
  if (destx)
    for (y=ymin; y<ymax; y++, destx += w)
      memcpy(destx, &DGEOPIXX(dgeofield, xmin, y), w2*sizeof(PIXTYPE));
  if (desty)
    for (y=ymin; y<ymax; y++, desty += w)
      memcpy(desty, &DGEOPIXY(dgeofield, xmin, y), w2*sizeof(PIXTYPE));

  return RETURN_OK;
  }

 

