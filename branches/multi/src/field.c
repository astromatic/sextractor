/*
*				field.c
*
* Handle field (image) structures.
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
*	Last modified:		10/02/2012
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
#include	"assoc.h"
#include	"astrom.h"
#include	"back.h"
#include	"field.h"
#include	"filter.h"
#include	"fitswcs.h"
#include	"header.h"
#include	"interpolate.h"

/****** field_init ***********************************************************
PROTO	fieldstruct *field_init(char *filename, int ext, int flags)
PURPOSE	Create and initialize a new field (image structure).
INPUT	Image filename,
	position among valid extensions in FITS file,
	image flags (e.g. DETECT_FIELD, MEASURE_FIELD...).
OUTPUT	Pointer to a new malloc'ed field structure.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	10/02/2012
 ***/
fieldstruct	*field_init(char *filename, int ext, int flags)

  {
   fieldstruct	*field;
   catstruct	*cat;
   tabstruct	*tab;
   h_type	htype;
   t_type	ttype;
   char		str[MAXCHAR], label[72],
		*photombuf, *pstr;
   int		j,s, line, nok, ntab, margin;

/* Move to ext'th valid FITS image extension */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: cannot open ", filename);

/* First allocate memory for the new field (and nullify pointers) */
  QCALLOC(field, fieldstruct, 1);
  field->flags = flags;
  field->cat = cat;
  tab = cat->tab;
  ext++;	/* At least one pass through the loop */
  nok = ext;
  for (ntab=cat->ntab; ext && ntab--; tab=tab->nexttab)
    {
    if ((tab->naxis < 2)
	|| !strncmp(tab->xtension, "BINTABLE", 8)
	|| !strncmp(tab->xtension, "ASCTABLE", 8))
      continue;
    field->tab = tab;
    ext--;
    }
  if (ntab<0)
    error(EXIT_FAILURE, "Not enough valid FITS image extensions in ",filename);

  strcpy (field->filename, filename);
/* A short, "relative" version of the filename */
  if (!(field->rfilename = strrchr(field->filename, '/')))
    field->rfilename = field->filename;
  else
    field->rfilename++;

/* Create a file name with a "header" extension */
  strcpy(field->hfilename, filename);
  if (!(pstr = strrchr(field->hfilename, '.')))
    pstr = field->hfilename+strlen(field->hfilename);
  sprintf(pstr, "%s", prefs.head_suffix);

  sprintf(gstr, "Looking for %s", field->rfilename);
  NFPRINTF(OUTPUT, gstr);
/* Check the image exists and read important info (image size, etc...) */
  field->headflag = !read_aschead(field->hfilename, nok - 1, field->tab);
  readimagehead(field);

  if (cat->ntab>1)
    sprintf(gstr, " [%d/%d]", nok, cat->tab->naxis<2? cat->ntab-1 : cat->ntab);
  QPRINTF(OUTPUT, "----- %s %s%s\n",
	flags&FLAG_FIELD?   "Flagging  from:" :
       (flags&(RMS_FIELD|VAR_FIELD|WEIGHT_FIELD)?
			     "Weighting from:" :
       (flags&MEASURE_FIELD? "Measuring from:" :
			     "Detecting from:")),
	field->rfilename,
        cat->ntab>1? gstr : "");
  QPRINTF(OUTPUT, "      \"%.20s\" / %s / %dx%d / %d bits %s\n",
	field->ident,
	field->headflag? "EXT. HEADER" : "no ext. header",
	field->width, field->height, field->tab->bytepix*8,
	field->tab->bitpix>0?
	(field->tab->compress_type!=COMPRESS_NONE?"(compressed)":"(integers)")
	:"(floats)");

/* Check the astrometric system and do the setup of the astrometric stuff */
  if (prefs.world_flag && (flags & (MEASURE_FIELD|DETECT_FIELD)))
    initastrom(field);
  else
    field->pixscale=prefs.pixel_scale;

/* Gain and Saturation */
  if (flags & (DETECT_FIELD|MEASURE_FIELD))
    {
    if (fitsread(field->tab->headbuf, prefs.gain_key, &field->gain,
	H_FLOAT, T_DOUBLE) != RETURN_OK)
      field->gain = prefs.gain;
    if (fitsread(field->tab->headbuf, prefs.satur_key, &field->satur_level,
	H_FLOAT, T_DOUBLE) !=RETURN_OK)
      field->satur_level = prefs.satur_level;
    }

  if ((flags & MEASURE_FIELD))
    {
/*-- Put a photometric label to the present field */
    field->photomlabel = 0;
/*-- Create dummy a FITS header to store all keyword values */
    QMALLOC(photombuf, char, FBSIZE);
    memset(photombuf, ' ', FBSIZE);
    strncpy(photombuf, "END     ", 8);
    for (s=0; s<prefs.nphotinstru_key; s++)
      {
      fitsadd(photombuf, prefs.photinstru_key[s], "");
/*---- Look in the main header */
      if ((line=fitsfind(cat->tab->headbuf,prefs.photinstru_key[s]))
		!= RETURN_ERROR)
        {
        fitspick(cat->tab->headbuf+line*80, str,(void *)label,&htype,&ttype, str);
        fitswrite(photombuf, prefs.photinstru_key[s], label, htype, ttype);
        }
/*---- Override with what is in the current header */
      if ((line=fitsfind(field->tab->headbuf,prefs.photinstru_key[s]))
		!= RETURN_ERROR)
        {
        fitspick(field->tab->headbuf+line*80,str,(void *)label,&htype,&ttype,
		str);
        fitswrite(photombuf, prefs.photinstru_key[s], label, htype, ttype);
        }
      }
/*-- Compare the dummy photometric FITS header to the ones previously stored */
    for (j=0; j<prefs.nphotinstru; j++)
      if (!strncmp((const char *)prefs.photinstrustr[j], photombuf,
		80*prefs.nphotinstru_key))
        {
        field->photomlabel = j;
        break;
        }
    if (j>=prefs.nphotinstru)
      {
      if (prefs.nphotinstru >= prefs.nphotinstrumax)
        {
        prefs.nphotinstrumax += prefs.nimage;
        QREALLOC(prefs.photinstrustr, char *, prefs.nphotinstrumax);
        }
      QMEMCPY(photombuf, prefs.photinstrustr[prefs.nphotinstru], char, FBSIZE);
      field->photomlabel = prefs.nphotinstru++;
      }
    free(photombuf);
    }

/* Background */
  if (flags & (DETECT_FIELD|MEASURE_FIELD|WEIGHT_FIELD|VAR_FIELD|RMS_FIELD))
    {
    field->ngamma = prefs.mag_gamma/log(10.0);

    field->backw = prefs.backsize[0]<field->width ? prefs.backsize[0]
						  : field->width;
    field->backh = prefs.backsize[1]<field->height ? prefs.backsize[1]
						   : field->height;
    field->nbackp = field->backw * field->backh;
    if ((field->nbackx = (field->width-1)/field->backw + 1) < 1)
      field->nbackx = 1;
    if ((field->nbacky = (field->height-1)/field->backh + 1) < 1)
      field->nbacky = 1;
    field->nback = field->nbackx * field->nbacky;
    field->nbackfx = field->nbackx>1 ? prefs.backfsize[0] : 1;
    field->nbackfy = field->nbacky>1 ? prefs.backfsize[1] : 1;
/*--  Set the back_type flag if absolute background is selected */
    if (((flags & DETECT_FIELD) && prefs.back_type[0]==BACK_ABSOLUTE)
	|| ((flags & MEASURE_FIELD) && prefs.back_type[1]==BACK_ABSOLUTE))
      field->back_type = BACK_ABSOLUTE;
    }

/* Add a comfortable margin for local background estimates */
  margin = (prefs.pback_type == LOCAL)? prefs.pback_size + prefs.mem_bufsize/4
					: 0;

  field->stripheight = prefs.mem_bufsize + margin;
  if (field->stripheight>field->height)
    field->stripheight = field->height;
/* Compute the image buffer size */
/* Basically, only one margin line is sufficient... */
  field->stripmargin = 1 + margin;
/* ...but : */
  if (prefs.filter_flag)
    {
/*-- If filtering is on, one should consider the height of the conv. mask */
/*-- + 1 line for detectinhg zero-weight neighbours */
    if (field->stripheight < thefilter->convh)
      field->stripheight = thefilter->convh;
    if (field->stripmargin < (margin = (thefilter->convh-1)/2+1))
      field->stripmargin = margin;
    }

  return field;
  }


/****** field_inherit ********************************************************
PROTO	fieldstruct *field_inherit(fieldstruct *infield, int flags)
PURPOSE	Make a copy of a field structure, e.g. for interpolation purposes.
INPUT	Pointer to input field structure,
	new field image flags (e.g. DETECT_FIELD, MEASURE_FIELD...).
OUTPUT	Pointer to a new malloc'ed field structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/12/2011
 ***/
fieldstruct	*field_inherit(fieldstruct *infield, int flags)

  {
   fieldstruct	*field;

/* First allocate memory for the new field (and nullify pointers) */
  QCALLOC(field, fieldstruct, 1);

/* Copy what is important and reset the remaining */
  *field = *infield;
  field->flags = flags;
  if (infield->wcs)
    field->wcs = copy_wcs(infield->wcs);
  field->interp_flag = 0;
  field->assoc = NULL;
  field->strip = NULL;
  field->fstrip = NULL;
  field->reffield = infield;

  return field;
  }


/****** field_end ************************************************************
PROTO	void field_end(fieldstruct *field)
PURPOSE	Free and close everything related to an image field structure.
INPUT	Pointer to field structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/12/2011
 ***/
void	field_end(fieldstruct *field)

  {
/* Free cat only if associated with an open file */
  if (field->cat &&  field->cat->file)
    free_cat(&field->cat, 1);
  free(field->strip);
  free(field->fstrip);
  if (field->wcs)
    end_wcs(field->wcs);
  if (field->interp_flag)
    end_interpolate(field);
  back_end(field);
  free(field);

  return;
  }

