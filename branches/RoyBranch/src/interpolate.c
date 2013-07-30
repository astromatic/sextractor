/*
*				interpolate.c
*
* Interpolate bad/missing image pixels.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1998-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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

#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"field.h"
#include	"interpolate.h"


/****************************** init_interpolate ****************************/
/*
Init resources required for data interpolation.
*/
void    init_interpolate(picstruct *field, int xtimeout, int ytimeout)

  {
  QMALLOC(field->interp_backup, PIXTYPE, field->width);
/* ytimeout < 0 means we don't need a timeout buffer (it won't be used) */
  if (ytimeout>=0)
    {
    QMALLOC(field->interp_ytimeoutbuf, int, field->width);
    memset(field->interp_ytimeoutbuf, 0, (size_t)(field->width*sizeof(int)));
    }

  field->interp_xtimeout = xtimeout;
  field->interp_ytimeout = ytimeout;

  field->interp_flag = 1;

  return;
  }


/******************************** interpolate *******************************/
/*
Interpolate (crudely) input data.
*/
void    interpolate(picstruct *field, picstruct *wfield,
		PIXTYPE *data, PIXTYPE *wdata)

  {
   PIXTYPE	*backup,*wbackup,
		thresh;
   int		*ytimeout,
		xtimeout,xtimeout0,ytimeout0, allflag, i;

  thresh = wfield->weight_thresh;
  backup = field->interp_backup;
  wbackup = wfield->interp_backup;
  ytimeout = wfield->interp_ytimeoutbuf;
  xtimeout0 = wfield->interp_xtimeout;
  ytimeout0 = wfield->interp_ytimeout;
  xtimeout = 0;	/* Start as if the previous pixel was already interpolated */
  allflag = field->interp_flag;
  for (i=field->width; i--; ytimeout++)
    {
/*-- Check if interpolation is needed */
    if (*wdata>=thresh)	/* It's a variance map: the bigger the worse */
      {
/*---- Check if the previous pixel was already interpolated */
      if (!xtimeout)
        {
        if (*ytimeout)
          {
          (*ytimeout)--;
          *wdata = *wbackup;
          if (allflag)
            *data = *backup;
          }
        }
      else
        {
        xtimeout--;
        *wdata = *(wdata-1);
        if (allflag)
          *data = *(data-1);
        }
      }
    else
      {
      xtimeout = xtimeout0;
     *ytimeout = ytimeout0;
      }
    *(wbackup++) = *(wdata++);
    if (allflag)
      *(backup++) = *(data++);
    }

  return;
  }


/******************************* end_interpolate ****************************/
/*
Free memory allocated for data interpolation.
*/
void    end_interpolate(picstruct *field)

  {
  free(field->interp_backup);
  free(field->interp_ytimeoutbuf);

  return;
  }

