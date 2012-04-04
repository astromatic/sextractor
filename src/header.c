/*
*				header.c
*
* Read external ASCII headers.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2010-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		26/03/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "field.h"
#include "header.h"
#include "prefs.h"

/****** header_readasc ******************************************************
PROTO	int header_readasc(char *filename, int frameno, tabstruct *tab)
PURPOSE	Read a ASCII header file and update the current field's tab
INPUT	Name of the ASCII file,
	Frame number (if extensions),
	Tab structure.
OUTPUT	RETURN_OK if the file was found, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/03/2012
 ***/
int     header_readasc(char *filename, int frameno, tabstruct *tab)
  {
   char         keyword[88],data[88],comment[88], str[88];
   FILE         *file;
   h_type       htype;
   t_type       ttype;
   int          i, flag;

  if ((file=fopen(filename, "r")))
    {
/*- Skip previous ENDs in multi-FITS extension headers */
    for (i=frameno; i--;)
      while (fgets(str, MAXCHAR, file)
		&& strncmp(str,"END ",4)
		&& strncmp(str,"END\n",4));
    memset(str, ' ', 80);
    flag = RETURN_ERROR;
    while (fgets(str, 81, file) && strncmp(str,"END ",4)
			&& strncmp(str,"END\n",4))
      {
      if (fitspick(str, keyword, data, &htype, &ttype, comment) != RETURN_OK)
        {
        memset(str, ' ', 80);
        continue;
        }
/*---- Block critical keywords */
      if (!wstrncmp(keyword, "SIMPLE  ", 8)
	||!wstrncmp(keyword, "BITPIX  ", 8)
	||!wstrncmp(keyword, "NAXIS   ", 8)
	||!wstrncmp(keyword, "BSCALE  ", 8)
	||!wstrncmp(keyword, "BZERO   ", 8))
        continue;
      addkeywordto_head(tab, keyword, comment);
      fitswrite(tab->headbuf, keyword, data, htype, ttype);
      memset(str, ' ', 80);
      flag = RETURN_OK;
      }
    fclose(file);
/*-- Update the tab data */
    readbasic_head(tab);
    return flag;
    }
  else
    return RETURN_ERROR;
  }


/****** header_readima ******************************************************
PROTO	void header_readima(fieldstruct *field)
PURPOSE	Read an image field header
INPUT	Pointer to field structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/03/2012
 ***/
void	header_readima(fieldstruct *field)
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


