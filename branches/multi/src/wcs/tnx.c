/*
*				tnx.c
*
* Manage the TNX astrometric format (from IRAF)
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic WCS library
*
*	Copyright:		(C) 2000-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"tnx.h"

/******* read_tnxaxis *********************************************************
PROTO	tnxaxisstruct *read_tnxaxis(char *tnxstr)
PURPOSE	Read a TNX axis mapping structure.
INPUT	String containing the TNX info.
OUTPUT	TNXAXIS structure if OK, or NULL in case of error.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/07/2006
 ***/

tnxaxisstruct	*read_tnxaxis(char *tnxstr)

  {
   tnxaxisstruct	*tnxaxis;
   char		*pstr, *ptr;
   double	min, max;
   int		i, order;

  if ((pstr=strpbrk(tnxstr, "1234567890-+.")))
    {
    if (!(tnxaxis=malloc(sizeof(tnxaxisstruct))))
      return NULL;
    tnxaxis->type = (int)(atof(strtok_r(pstr, " ", &ptr))+0.5);
    tnxaxis->xorder = (pstr=strtok_r(NULL, " ", &ptr))?
			(int)(atof(pstr)+0.5) : 0;
    tnxaxis->yorder = (pstr=strtok_r(NULL, " ", &ptr))?
			(int)(atof(pstr)+0.5) : 0;
    tnxaxis->xterms = (pstr=strtok_r(NULL, " ", &ptr))?
			(int)(atof(pstr)+0.5) : 0;
    min = (pstr=strtok_r(NULL, " ", &ptr))? atof(pstr) : 0.0;
    max = (pstr=strtok_r(NULL, " ", &ptr))? atof(pstr) : 0.0;
    if (max <= min)
      return NULL;
    tnxaxis->xrange = 2.0 / (max - min);
    tnxaxis->xmaxmin =  - (max + min) / 2.0;
    min = (pstr=strtok_r(NULL, " ", &ptr))? atof(pstr) : 0.0;
    max = (pstr=strtok_r(NULL, " ", &ptr))? atof(pstr) : 0.0;
    if (max <= min)
      return NULL;
    tnxaxis->yrange = 2.0 / (max - min);
    tnxaxis->ymaxmin =  - (max + min) / 2.0;
    switch (tnxaxis->xterms)
      {
      case TNX_XNONE:
        tnxaxis->ncoeff = tnxaxis->xorder + tnxaxis->yorder - 1;
        break;
      case TNX_XHALF:
        order = tnxaxis->xorder<tnxaxis->yorder?
			tnxaxis->xorder : tnxaxis->yorder;
        tnxaxis->ncoeff = tnxaxis->xorder*tnxaxis->yorder - order*(order-1)/2;
        break;
      case TNX_XFULL:
        tnxaxis->ncoeff = tnxaxis->xorder * tnxaxis->yorder;
        break;
      default:
        return NULL;
      }
/*-- Now read the mapping coefficients */
    if (!(tnxaxis->coeff=malloc(tnxaxis->ncoeff*sizeof(double))))
      return NULL;
    for (i=0; i<tnxaxis->ncoeff && (pstr=strtok_r(NULL, " ", &ptr)); i++)
      tnxaxis->coeff[i] = atof(pstr);
    if (i!=tnxaxis->ncoeff)
      return NULL;
    if (!(tnxaxis->xbasis=malloc(tnxaxis->xorder*sizeof(double))))
      return NULL;
    if (!(tnxaxis->ybasis=malloc(tnxaxis->yorder*sizeof(double))))
      return NULL;
    return tnxaxis;
    }
  else
    return NULL;
  }


/******* copy_tnxaxis *********************************************************
PROTO	tnxaxisstruct *copy_tnxaxis(tnxaxisstruct *axis)
PURPOSE	Copy a TNX axis mapping structure.
INPUT	TNXAXIS structure pointer.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/11/2003
 ***/

tnxaxisstruct	*copy_tnxaxis(tnxaxisstruct *axis)

  {
   tnxaxisstruct	*tnxaxis;
   int			i;

  if (axis)
    {
    if (!axis->ncoeff)
      return NULL;
    if (!(tnxaxis=malloc(sizeof(tnxaxisstruct))))
      return NULL;
    *tnxaxis = *axis;
    if (!(tnxaxis->coeff=malloc(tnxaxis->ncoeff*sizeof(double))))
      return NULL;
    for (i=0; i<tnxaxis->ncoeff; i++)
      tnxaxis->coeff[i] = axis->coeff[i];
    if (!(tnxaxis->xbasis=malloc(tnxaxis->xorder*sizeof(double))))
      return NULL;
    for (i=0; i<tnxaxis->xorder; i++)
      tnxaxis->xbasis[i] = axis->xbasis[i];
    if (!(tnxaxis->ybasis=malloc(tnxaxis->yorder*sizeof(double))))
      return NULL;
    for (i=0; i<tnxaxis->yorder; i++)
      tnxaxis->ybasis[i] = axis->ybasis[i];
    return tnxaxis;
    }

  return NULL;
  }


/******* free_tnxaxis *********************************************************
PROTO	void free_tnxaxis(tnxaxisstruct *axis)
PURPOSE	Free a TNX axis mapping structure.
INPUT	TNXAXIS structure pointer.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/04/2000
 ***/

void	free_tnxaxis(tnxaxisstruct *axis)

  {
  if (axis)
    {
    
    free(axis->coeff);
    free(axis->xbasis);
    free(axis->ybasis);
    free(axis);
    }

  return;
  }


/******* raw_to_tnxaxis *******************************************************
PROTO	double raw_to_tnxaxis(tnxaxisstruct *axis, double x, double y)
PURPOSE	Compute the correction value on a TNX axis at current position.
INPUT	TNXAXIS structure pointer,
	x coordinate,
	y coordinate.
OUTPUT	Value on the TNXaxis.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/04/2000
 ***/

double	raw_to_tnxaxis(tnxaxisstruct *axis, double x, double y)

  {
   double	*xbasis, *ybasis,*coeff,
		norm, accum, val;
   int		i, j, xorder,xorder0,yorder,maxorder,xterms;

  xbasis = axis->xbasis;
  ybasis = axis->ybasis;
  xorder = axis->xorder;
  yorder = axis->yorder;
  xterms = axis->xterms;

  switch (axis->type)
    {
    case TNX_CHEBYSHEV:
      xbasis[0] = 1.0;
      if (xorder > 1)
        {
        xbasis[1] = norm = (x + axis->xmaxmin)*axis->xrange;
        if (xorder > 2)
          for (i = 2; i < xorder; i++)
	    xbasis[i] = 2.0*norm*xbasis[i-1] - xbasis[i-2];
        }
      ybasis[0] = 1.0;
      if (yorder > 1)
        {
        ybasis[1] = norm = (y + axis->ymaxmin)*axis->yrange;
        if (yorder > 2)
          for (i = 2; i < yorder; i++)
	    ybasis[i] = 2.0*norm*xbasis[i-1] - ybasis[i-2];
        }
      break;

    case TNX_LEGENDRE:
      xbasis[0] = 1.0;
      if (xorder > 1)
        {
        xbasis[1] = norm = (x + axis->xmaxmin)*axis->xrange;
        if (xorder > 2)
          for (i = 2; (j=i) < xorder; i++)
            xbasis[i] = ((2.0*j - 3.0) * norm * xbasis[i-1] -
                       (j - 2.0) * xbasis[i-2]) / (j - 1.0);
        }
      ybasis[0] = 1.0;
      if (yorder > 1)
        {
        ybasis[1] = norm = (y + axis->ymaxmin)*axis->yrange;
        if (yorder > 2)
          for (i = 2; (j=i) < xorder; i++)
            ybasis[i] = ((2.0*j - 3.0) * norm * ybasis[i-1] -
                       (j - 2.0) * ybasis[i-2]) / (j - 1.0);
        }
      break;

    case TNX_POLYNOMIAL:
      xbasis[0] = 1.0;
      if (xorder > 1)
        {
        xbasis[1] = x;
        if (xorder > 2)
          for (i = 2; i < xorder; i++)
            xbasis[i] = x * xbasis[i-1];
        }
      ybasis[0] = 1.0;
      if (yorder > 1)
        {
        ybasis[1] = y;
        if (yorder > 2)
          for (i = 2; i < yorder; i++)
            ybasis[i] = y * ybasis[i-1];
        }
      break;

    default:
      return 0.0;
    }

/* Loop over y basis functions */
  maxorder = xorder > yorder ? xorder : yorder;
  xorder0 = xorder;
  coeff = axis->coeff;
  val = 0.0;
  for (i = 0; i<yorder; i++)
    {
/*-- Loop over the x basis functions */
    accum = 0.0;
    xbasis = axis->xbasis;
    for (j = xorder; j--;)
      accum += *(coeff++) * *(xbasis++);
    val += accum**(ybasis++);

/*-- Elements of the coefficient vector where neither k = 1 or i = 1
           are not calculated if sf->xterms = no. */
    if (xterms == TNX_XNONE)
      xorder = 1;
    else if (xterms == TNX_XHALF && (i + 1 + xorder0) > maxorder)
      xorder--;
    }

  return val;
  }


