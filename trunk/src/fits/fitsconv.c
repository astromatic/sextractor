/*
*				fitsconv.c
*
* Functions for converting LDAC FITS catalogs.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic FITS/LDAC library
*
*	Copyright:		(C) 1995-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		09/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"fitscat_defs.h"
#include	"fitscat.h"

/****** asc2bin_tab ***********************************************************
PROTO	tabstruct *asc2bin_tab(catstruct *catin, char *tabinname, 
        catstruct *catout, char *taboutname)
PURPOSE	Convert an ASCII table to a BINARY table
	which is then stored in a destination catalog.
INPUT	Pointer to the input catalog,
	Name of the input ASCII table,
	Pointer to the output catalog,
	Name of the output binary table.
OUTPUT	RETURN_OK if the ASCII table was transformed, and RETURN_ERROR
	otherwise.
NOTES	This function can be used to stick the binary translation of
        similar ASCII tables.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	25/09/2004
 ***/
tabstruct *asc2bin_tab(catstruct *catin, char *tabinname, catstruct *catout,
		char *taboutname)

  {
   catstruct	*tabcat;
   keystruct	*key;
   tabstruct	*tabin,*tabout;
   char		comment[82], keyword[16], ptr[82];
   h_type	htype;
   t_type	ttype;
   char		*buf, *lptr;
   int		i;

/*Convert the table name to a pointer*/
  if (!(tabin = name_to_tab(catin, tabinname, 0)))
    return NULL;

/*Get the original catalog*/
  tabcat = tabin->cat;

/*Create a new binary table*/
  tabout = new_tab(taboutname);

/*Alloc. mem. for the whole ASCII table at once (should not be very large)*/
  QMALLOC(buf, char, tabin->tabsize);
/*Now read all the elements from the original table*/
  if (open_cat(tabcat, READ_ONLY) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: Cannot access ", tabcat->filename);
  QFSEEK(tabcat->file, tabin->bodypos, SEEK_SET, tabcat->filename);
  QFREAD(buf, tabin->tabsize, tabcat->file, tabcat->filename);
  if (close_cat(tabcat) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: Problem while closing", tabcat->filename);
  lptr = buf;
  for (i=tabin->tabsize/80; i-- && strncmp(lptr, "END     ", 8);)
    {
/*Interprete the next FITS line */
    if (fitspick(lptr, keyword, ptr, &htype, &ttype, comment) != RETURN_OK) {
      char line[81];
      int  qflag=1;
      strncpy(line, lptr, 80);
      line[80] = '\0';
      QFPRINTF(OUTPUT, line);
      warning("*Warning*: incorrect FITS field will be ignored in ",
	tabcat->filename);
    }
    if (htype != H_COMMENT)
      {
/*----Create a new key and fill it with the right parameters*/
      key = new_key(keyword);
      strcpy(key->comment, comment+strspn(comment, " "));
      key->htype = htype;
      key->ttype = ttype;
      key->nbytes = t_size[ttype];
/*----!!Temporary (?)  solution for STRINGS*/
      if (htype==H_STRING)
        {
        key->naxis = 1;
        QMALLOC(key->naxisn, int, 1);
        key->naxisn[0] = 32;
        key->nbytes *= key->naxisn[0];
        }
      key->nobj = 1;
/*----Allocate memory and copy data in the same time*/
      QMEMCPY(ptr, key->ptr, char, key->nbytes);
      if (add_key(key, tabout, 0)==RETURN_ERROR)
        {
        sprintf(comment, "%s keyword found twice in ",
		keyword);
        warning(comment, tabcat->filename);
        }
      }
    lptr += 80;
    }

  free(buf);

  update_tab(tabout);
  return tabout;
  }


/****** ttypeconv ************************************************************
PROTO	void    ttypeconv(void *ptrin, void *ptrout,
		t_type ttypein, t_type ttypeout)
PURPOSE	Convert data from one type to another.
INPUT	Pointer to element to convert,
	destination pointer to the converted element,
	t_type of the element to convert,
	t_type of the converted element.
OUTPUT	-.
NOTES	ttypeconv does not yet handle arrays.
AUTHOR  E. Bertin (IAP)
VERSION 19/11/2009
 ***/

void	ttypeconv(void *ptrin, void *ptrout, t_type ttypein, t_type ttypeout)

  {
//   union	{char tbyte; short tshort; int tlong; LONGLONG tlonglong;
//		 float tfloat; double tdouble; char tstring;} ival;

#ifdef HAVE_LONG_LONG_INT
#define		OUTCONV(x, y)			\
      switch(y)					\
         {					\
         case T_BYTE:				\
         case T_STRING:				\
	   *((char *)ptrout) = (char)x;		\
           break;				\
         case T_SHORT:				\
           *((short *)ptrout) = (short)x;	\
           break;				\
         case T_LONG:				\
           *((int *)ptrout) = (int)x;		\
           break;				\
         case T_LONGLONG:			\
           *((LONGLONG *)ptrout) = (LONGLONG)x;	\
           break;				\
         case T_FLOAT:				\
           *((float *)ptrout) = (float)x;	\
           break;				\
         case T_DOUBLE:				\
           *((double *)ptrout) = (double)x;	\
           break;				\
         default:				\
	   break;				\
         }
#else
#define		OUTCONV(x, y)			\
      switch(y)					\
         {					\
         case T_BYTE:				\
         case T_STRING:				\
	   *((char *)ptrout) = (char)x;		\
           break;				\
         case T_SHORT:				\
           *((short *)ptrout) = (short)x;	\
           break;				\
         case T_LONG:				\
           *((int *)ptrout) = (int)x;		\
           break;				\
         case T_FLOAT:				\
           *((float *)ptrout) = (float)x;	\
           break;				\
         case T_DOUBLE:				\
           *((double *)ptrout) = (double)x;	\
           break;				\
         default:				\
	   break;				\
         }
#endif

  switch(ttypein)
    {
    case T_BYTE:
    case T_STRING:
      OUTCONV(*(char *)ptrin, ttypeout);
      break;
    case T_SHORT:
      OUTCONV(*(short *)ptrin, ttypeout);
      break;
    case T_LONG:
      OUTCONV(*(int *)ptrin, ttypeout);
      break;
#ifdef HAVE_LONG_LONG_INT
    case T_LONGLONG:
      OUTCONV(*(LONGLONG *)ptrin, ttypeout);
      break;
#endif
    case T_FLOAT:
      OUTCONV(*(float *)ptrin, ttypeout);
      break;
    case T_DOUBLE:
      OUTCONV(*(double *)ptrin, ttypeout);
      break;
    default:
      break;
    }

  return;
  }

