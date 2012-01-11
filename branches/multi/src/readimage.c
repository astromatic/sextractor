/*
*				readimage.c.
*
* Read image data.
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
*	Last modified:		11/01/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"wcs/wcs.h"
#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"check.h"
#include	"field.h"
#include	"fits/fitscat.h"
#include	"fitswcs.h"
#include	"interpolate.h"
#include	"back.h"
#include	"astrom.h"
#include	"weight.h"
#include        "wcs/tnx.h"

/******************************* loadstrip ***********************************/
/*
Load a new strip of pixel data into the buffer.
*/
void	*loadstrip(fieldstruct *field, fieldstruct *wfield)

  {
   tabstruct	*tab;
   checkstruct	*check;
   int		y, w, flags, interpflag;
   PIXTYPE	*data, *wdata, *rmsdata;

  tab = field->tab;
  w = field->width;
  flags = field->flags;
  interpflag = (wfield && wfield->interp_flag);
  wdata = NULL;			/* To avoid gcc -Wall warnings */

  if (!field->y)
    {
/*- First strip */
     int	nbpix;

    nbpix = w*field->stripheight;

    if (flags ^ FLAG_FIELD)
      {
/*---- Allocate space for the frame-buffer */
      if (!(field->strip=(PIXTYPE *)malloc(field->stripheight*field->width
        *sizeof(PIXTYPE))))
        error(EXIT_FAILURE,"Not enough memory for the image buffer of ",
		field->rfilename);

      data = field->strip;
/*---- We assume weight data have been read just before */
      if (interpflag)
        wdata = wfield->strip;
      if (flags & BACKRMS_FIELD)
        for (y=0, rmsdata=data; y<field->stripheight; y++, rmsdata += w)
          back_rmsline(field, y, rmsdata);
      else if (flags & INTERP_FIELD)
        copydata(field, 0, nbpix);
      else
        read_body(tab, data, nbpix);
      if (flags & (WEIGHT_FIELD|RMS_FIELD|BACKRMS_FIELD|VAR_FIELD))
        weight_to_var(field, data, nbpix);
      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_IDENTICAL]))
        check_write(check, data, nbpix);
      for (y=0; y<field->stripheight; y++, data += w)
        {
/*------ This is the only place where one can pick-up safely the current bkg */
        if (flags & (MEASURE_FIELD|DETECT_FIELD))
          back_subline(field, y, data);
/*------ Go to interpolation process */
        if (interpflag)
          {
          interpolate(field,wfield, data, wdata);
          wdata += w;
          }
/*------ Check-image stuff */
        if (prefs.check_flag)
          {
          if (flags & MEASURE_FIELD)
            {
            if ((check = prefs.check[CHECK_BACKGROUND]))
              check_write(check, field->backline, w);
            if ((check = prefs.check[CHECK_SUBTRACTED]))
              check_write(check, data, w);
            if ((check = prefs.check[CHECK_APERTURES]))
              check_write(check, data, w);
            if ((check = prefs.check[CHECK_SUBPSFPROTOS]))
              check_write(check, data, w);
            if ((check = prefs.check[CHECK_SUBPROFILES]))
              check_write(check, data, w);
            if ((check = prefs.check[CHECK_SUBSPHEROIDS]))
              check_write(check, data, w);
            if ((check = prefs.check[CHECK_SUBDISKS]))
              check_write(check, data, w);
            }
          if ((flags&DETECT_FIELD) && (check=prefs.check[CHECK_BACKRMS]))
            {
            back_rmsline(field, y, (PIXTYPE *)check->pix);
            check_write(check, check->pix, w);
            }
          }
        }
      }
    else
      {
      if (!(field->fstrip=(FLAGTYPE *)malloc(field->stripheight*field->width
		*sizeof(FLAGTYPE))))
      error(EXIT_FAILURE,"Not enough memory for the flag buffer of ",
	field->rfilename);
      read_ibody(field->tab, field->fstrip, nbpix);
      }

    field->ymax = field->stripheight;
    if (field->ymax < field->height)
      field->stripysclim = field->stripheight - field->stripmargin;
    }
  else
    {
/*- other strips */
    if (flags ^ FLAG_FIELD)
      {
      data = field->strip + field->stripylim*w;
/*---- We assume weight data have been read just before */
      if (interpflag)
        wdata = wfield->strip + field->stripylim*w;

/*---- copy to Check-image the "oldest" line before it is replaced */
      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_SUBOBJECTS]))
        check_write(check, data, w);

      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_MASK]))
        check_write(check, data, w);

      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_SUBMASK]))
        check_write(check, data, w);

      if (flags & BACKRMS_FIELD)
        back_rmsline(field, field->ymax, data);
      else if (flags & INTERP_FIELD)
        copydata(field, field->stripylim*w, w);
      else
        read_body(tab, data, w);
      if (flags & (WEIGHT_FIELD|RMS_FIELD|BACKRMS_FIELD|VAR_FIELD))
        weight_to_var(field, data, w);

      if ((flags & MEASURE_FIELD) && (check=prefs.check[CHECK_IDENTICAL]))
        check_write(check, data, w);
/*---- Interpolate and subtract the background at current line */
      if (flags & (MEASURE_FIELD|DETECT_FIELD))
        back_subline(field, field->ymax, data);
      if (interpflag)
        interpolate(field,wfield, data, wdata);
/*---- Check-image stuff */
      if (prefs.check_flag)
        {
        if (flags & MEASURE_FIELD)
          {
          if ((check = prefs.check[CHECK_BACKGROUND]))
            check_write(check, field->backline, w);
          if ((check = prefs.check[CHECK_SUBTRACTED]))
            check_write(check, data, w);
          if ((check = prefs.check[CHECK_APERTURES]))
            check_write(check, data, w);
          if ((check = prefs.check[CHECK_SUBPSFPROTOS]))
            check_write(check, data, w);
          if ((check = prefs.check[CHECK_SUBPROFILES]))
            check_write(check, data, w);
          if ((check = prefs.check[CHECK_SUBSPHEROIDS]))
            check_write(check, data, w);
          if ((check = prefs.check[CHECK_SUBDISKS]))
            check_write(check, data, w);
          }
        if ((flags&DETECT_FIELD) && (check=prefs.check[CHECK_BACKRMS]))
          {
          back_rmsline(field, field->ymax, (PIXTYPE *)check->pix);
          check_write(check, check->pix, w);
          }
        }
      }
    else
      read_ibody(tab, field->fstrip + field->stripylim*w, w);

    field->stripylim = (++field->ymin)%field->stripheight;
    if ((++field->ymax)<field->height)
      field->stripysclim = (++field->stripysclim)%field->stripheight;
    }

  return (flags ^ FLAG_FIELD)?
		  (void *)(field->strip + field->stripy*w)
		: (void *)(field->fstrip + field->stripy*w);
  }


/******************************** copydata **********************************/
/*
Copy image data from one field to the other.
*/
void	copydata(fieldstruct *field, int offset, int size)
  {
  memcpy(field->strip+offset, field->reffield->strip+offset,
		size*sizeof(PIXTYPE));
  return;
  }


/******************************* readimagehead *******************************/
/*
extract some data from the FITS-file header
*/
void	readimagehead(fieldstruct *field)
  {
#define FITSREADS(buf, k, str, def) \
                {if (fitsread(buf,k,str, H_STRING,T_STRING) != RETURN_OK) \
                   strcpy(str, (def)); \
                }

   tabstruct	*tab;

  tab = field->tab;

  if(tab->naxis < 2)
    error(EXIT_FAILURE, field->filename, " does NOT contain 2D-data!");

/*---------------------------- Basic keywords ------------------------------*/
  if (tab->bitpix != BP_BYTE
	&& tab->bitpix != BP_SHORT
	&& tab->bitpix != BP_LONG
	&& tab->bitpix != BP_FLOAT
	&& tab->bitpix != BP_DOUBLE)
    error(EXIT_FAILURE, "Sorry, I don't know that kind of data.", "");

  field->width = tab->naxisn[0];
  field->height = tab->naxisn[1];
  field->npix = (KINGSIZE_T)field->width*field->height;
  field->bitpix = tab->bitpix;
  if (tab->bitsgn && prefs.fitsunsigned_flag)
    tab->bitsgn = 0;

  FITSREADS(tab->headbuf, "OBJECT  ", field->ident, "Unnamed");

/*----------------------------- Astrometry ---------------------------------*/
/* Presently, astrometry is done only on the measurement and detect images */
  if (field->flags&(MEASURE_FIELD|DETECT_FIELD))
    field->wcs = read_wcs(tab);

  QFSEEK(field->cat->file, tab->bodypos, SEEK_SET, field->filename);

  return;

#undef FITSREADS

  }


