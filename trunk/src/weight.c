/*
*				weight.c
*
* Handle external weight maps.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1997-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		11/10/2010
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
#include	"field.h"
#include	"plist.h"
#include	"weight.h"

/******************************* newweight **********************************/
/*
Load a weight map and initialize relevant parameters.
*/
picstruct	*newweight(char *filename, picstruct *reffield,
			weightenum wtype, int nok)

  {
   picstruct	*wfield;

  switch(wtype)
    {
    case WEIGHT_FROMINTERP:
      wfield = inheritfield(reffield, INTERP_FIELD);
      break;

    case WEIGHT_FROMBACK:
      wfield = inheritfield(reffield, BACKRMS_FIELD);
      break;

    case WEIGHT_FROMRMSMAP:
      wfield = newfield(filename, RMS_FIELD, nok);
      if ((wfield->width!=reffield->width)||(wfield->height!=reffield->height))
        error(EXIT_FAILURE,
	"*Error*: measured frame and weight map have different sizes","");
      wfield->sigfac = 1.0;
      break;

    case WEIGHT_FROMVARMAP:
      wfield = newfield(filename, VAR_FIELD, nok);
      if ((wfield->width!=reffield->width)||(wfield->height!=reffield->height))
        error(EXIT_FAILURE,
	"*Error*: measured frame and weight map have different sizes","");
      break;

    case WEIGHT_FROMWEIGHTMAP:
      wfield = newfield(filename, WEIGHT_FIELD, nok);
      if ((wfield->width!=reffield->width)||(wfield->height!=reffield->height))
        error(EXIT_FAILURE,
	"*Error*: measured frame and weight map have different sizes","");
      break;
    default:
      wfield = NULL;			/* To avoid gcc -Wall warnings */
      error(EXIT_FAILURE,
	"*Internal Error*: Unknown weight-map type in ", "makeit()");
      break;
    }

/* Default normalization factor (will be changed if necessary later) */
  wfield->sigfac = 1.0;

  return wfield;
  }


/******************************* weight_to_var *******************************/
/*
Transform an array of possibily unnormalized weights into a calibrated
variance map.
*/
void	weight_to_var(picstruct *wfield, PIXTYPE *data, int npix)

  {
   float	sigfac2;
   int		i;

  switch(wfield->flags&(BACKRMS_FIELD|RMS_FIELD|VAR_FIELD|WEIGHT_FIELD))
    {
    case BACKRMS_FIELD:
    case RMS_FIELD:
      for (i=npix; i--; data++)
        if (*data<BIG)
          *data *= *data;
      break;
    case VAR_FIELD:
      sigfac2 = wfield->sigfac*wfield->sigfac;
      for (i=npix; i--;)
        *(data++) *= sigfac2;
      break;
    case WEIGHT_FIELD:
     sigfac2 = wfield->sigfac*wfield->sigfac;
     for (i=npix; i--; data++)
      if (*data > 0.0)
        *data = sigfac2/(*data);
      else
        *data = BIG;
      break;
    default:
      error(EXIT_FAILURE,
	"*Internal Error*: Unknown weight-map type in ", "weight_to_var()");
      break;
    }

  return;
  }


/******************************** weight_count *******************************
PROTO	void weight_count(objstruct *obj, pliststruct *pixel)
PURPOSE	Count pixels with zero weights.
INPUT   Objstruct pointer,
	pixel list pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 01/10/2009
 ***/
void	weight_count(objstruct *obj, pliststruct *pixel)

  {
   pliststruct	*pixt;
   int		i, nw,ndw, wflag;

  nw = ndw = wflag = 0;

  for (i=obj->firstpix; i!=-1; i=PLIST(pixt,nextpix))
    {
    pixt = pixel+i;
    if (PLISTFLAG(pixt, wflag) & OBJ_LOWWEIGHT)
      nw++;
    if (PLISTFLAG(pixt, wflag) & OBJ_LOWDWEIGHT)
      ndw++;
    }

  obj->nzwpix = nw;
  obj->nzdwpix = ndw;
  obj->wflag = nw? OBJ_LOWWEIGHT : 0;
  if (ndw)
    obj->wflag |= OBJ_LOWDWEIGHT;

  return;
  }



