 /*
 				header.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Read extern ASCII headers.
*
*	Last modify:	05/02/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "header.h"
#include "prefs.h"

/****** read_aschead ********************************************************
PROTO	int	read_aschead(char *filename, int frameno, tabstruct *tab)
PURPOSE	Read a ASCII header file and update the current field's tab
INPUT	Name of the ASCII file,
	Frame number (if extensions),
	Tab structure.
OUTPUT	RETURN_OK if the file was found, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	05/02/2010
 ***/
int     read_aschead(char *filename, int frameno, tabstruct *tab)
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

