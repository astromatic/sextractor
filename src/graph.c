/*
*				graph.c
*
* Add simple graphics to image rasters.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		23/11/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdlib.h>

#include	"define.h"
#include	"globals.h"

float	sexx1, sexy1;

/****** sexmove **************************************************************
PROTO	void sexmove(double x, double y)
PURPOSE	Move graphic pointer in a PIXTYPE image raster.
INPUT	Input pixel coordinate on FITS AXIS1,
	input pixel coordinate on FITS AXIS2.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/11/2011
 ***/
void	sexmove(float x, float y)
  {
  sexx1 = x;
  sexy1 = y;

  return;
  }


/****** sexdraw *************************************************************
PROTO	void sexdraw(PIXTYPE *raster, int w, int h, float x, float y,
		PIXTYPE val)
PURPOSE	Draw a line in a PIXTYPE image raster.
INPUT	Pointer to a PIXTYPE image raster,
	image width [pixels],
	image height [pixels],
	destination pixel coordinate on FITS AXIS1,
	destination pixel coordinate on FITS AXIS2,
	line grey level.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/11/2011
 ***/
void	sexdraw(PIXTYPE *raster, int w, int h, double x, double y, PIXTYPE val)
  {
   float	dx,dy, slope;
   int		ix1,iy1, ix2,iy2, ix,iy;

  dx = x - sexx1;
  dy = y - sexy1;
  if (fabs(dx) > fabs(dy))
    {
    slope = dy/dx;
    ix1 = RINT(sexx1);
    ix2 = RINT(x);
    if (ix2>ix1)
      {
      for (ix=ix1+1; ix<=ix2; ix++)
        if (ix>=0 && ix<w)
          {
          iy = RINT(sexy1+(ix-sexx1)*slope);
          if (iy>=0 && iy<h)
            raster[ix+w*iy] += val;
          }
      }
    else
      {
      for (ix=ix1-1; ix>=ix2; ix--)
        if (ix>=0 && ix<w)
          {
          iy = RINT(sexy1+(ix-sexx1)*slope);
          if (iy>=0 && iy<h)
            raster[ix+w*iy] += val;
          }
      }
    }
  else
    {
    slope = dx/(dy == 0.0? 1.0:dy);
    iy1 = RINT(sexy1);
    iy2 = RINT(y);
    if (iy2>iy1)
      {
      for (iy=iy1+1; iy<=iy2; iy++)
        if (iy>=0 && iy<h)
          {
          ix = RINT(sexx1+(iy-sexy1)*slope);
          if (ix>=0 && ix<w)
            raster[ix+w*iy] += val;
          }
      }
    else
      for (iy=iy1-1; iy>=iy2; iy--)
        {
        if (iy>=0 && iy<h)
          {
          ix = RINT(sexx1+(iy-sexy1)*slope);
          if (ix>=0 && ix<w)
            raster[ix+w*iy] += val;
          }
        }
    }

  sexx1 = x;
  sexy1 = y;

  return;
  }


/****** sexcircle ***********************************************************
PROTO	void sexcircle(PIXTYPE *raster, int w, int h, float x, float y,
		float radius, PIXTYPE val)
PURPOSE	Draw a circle in a PIXTYPE image raster.
INPUT	Pointer to a PIXTYPE image raster,
	image width [pixels],
	image height [pixels],
	destination pixel coordinate on FITS AXIS1,
	destination pixel coordinate on FITS AXIS2,
	circle radius [pixels],
	circle grey level.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/11/2011
 ***/
void	sexcircle(PIXTYPE *raster, int w,int h, float x, float y, float radius,
		PIXTYPE val)
  {
   float	alpha, cta,sta;
   int		i;

  sexmove(x+radius, y);

  for (i=1; i<37; i++)
    {
    alpha = i*PI/18.0;
#ifdef HAVE_SINCOSF
    sincosf(alpha, &sta, &cta);
#else
    sta = sinf(alpha);
    cta = cosf(alpha);
#endif
    sexdraw(raster,w,h, x + radius*cta, y + radius*sta, val);
    }

  return;
  }


/****** sexellipse ***********************************************************
PROTO	void sexellipse(PIXTYPE *raster, int w, int h, float x, float y,
		float a, float b, float theta, PIXTYPE val)
PURPOSE	Draw an ellipse in a PIXTYPE image raster.
INPUT	Pointer to a PIXTYPE image raster,
	image width [pixels],
	image height [pixels],
	destination pixel coordinate on FITS AXIS1,
	destination pixel coordinate on FITS AXIS2,
	ellipse semi-major axis length [pixels],
	ellipse semi-minor axis length [pixels],
	ellipse position angle from FITS AXIS1 towards FITS AXIS2 [degrees],
	ellipse grey level.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/11/2011
 ***/
void	sexellipse(PIXTYPE *raster, int w, int h, float x, float y, float a,
		float b, float theta, PIXTYPE val, int dotflag)

  {
   float	alpha, ct,st, cta,sta;
   int		i;

  ct = cosf(PI*theta/180);
  st = sinf(PI*theta/180);

  sexmove(x+a*ct, y+a*st);

  for (i=1; i<37; i++)
    {
    alpha = i*PI/18.0;
#ifdef HAVE_SINCOSF
    sincosf(alpha, &sta, &cta);
#else
    sta = sinf(alpha);
    cta = cosf(alpha);
#endif
    sta *= b;
    cta *= a;
    if (dotflag && !(i&1))
      sexmove(x + cta*ct - sta*st, y + cta*st + sta*ct);
    else
      sexdraw(raster,w,h, x + cta*ct - sta*st, y + cta*st + sta*ct, val);
    }

  return;
  }


