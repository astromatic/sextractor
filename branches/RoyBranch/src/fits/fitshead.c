/*
*				fitshead.c
*
* General functions for handling FITS file headers
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic FITS/LDAC library
*
*	Copyright:		(C) 1995-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		29/08/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef	HAVE_CONFIG_H
#include "config.h"
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"fitscat_defs.h"
#include	"fitscat.h"

extern	char	histokeys[][12];
const int	t_size[] = {1, 2, 4, 8, 4, 8, 1};/* size in bytes per t_type */

/******* get_head *************************************************************
PROTO	int get_head(tabstruct *tab)
PURPOSE	Read a FITS header.
INPUT	Table structure.
OUTPUT	RETURN_OK if a FITS header has been found and loaded, or RETURN_ERROR
	otherwise.
NOTES	The file must be opened, and the file pointer must be located at
	the beginning of a header.
	The headbuf pointer in the catstruct is reallocated.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	08/02/96
 ***/
int	get_head(tabstruct *tab)

  {
   catstruct	*cat;
   int		i;
   char		*buf;

  buf = tab->headbuf;
  if (!(cat = tab->cat))
    error(EXIT_FAILURE, "*Internal Error*: Table has no parent catalog","!");

  QFREE(buf);
  QMALLOC(buf, char, FBSIZE);

/*Read the first block and check that it is FITS */
  if (!fread(buf, FBSIZE, 1, cat->file))
    {
    QFREE(buf);
    return RETURN_ERROR;
    }

  if (strncmp(buf, "SIMPLE  ", 8) && strncmp(buf, "XTENSION", 8))
    {
    QFREE(buf);
    return RETURN_ERROR;
    }

/*Find the number of FITS blocks of the header while reading it */
  for (i=1; !fitsnfind(buf,"END     ", i); i++)
    {
    QREALLOC(buf, char, FBSIZE*(i+1));
    QFREAD(&buf[FBSIZE*i], FBSIZE, cat->file, cat->filename);
    }

  tab->headnblock = i;
  tab->headbuf = buf;

  return  RETURN_OK;
  }


/****** readbasic_head ********************************************************
PROTO	void readbasic_head(tabstruct *tab)
PURPOSE	Read the current FITS header basic keywords.
INPUT	pointer to catstruct.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/06/2012
 ***/
void	readbasic_head(tabstruct *tab)

  {
   char		str[88];
   char		key[12], name[16],
		*filename;
   int		i;
   KINGSIZE_T	tabsize;

  filename = (tab->cat? tab->cat->filename : strcpy(name, "internal header"));

// CFITSIO
  char NAXIS_KEYWORD[14];
  char BITPIX_KEYWORD[14];

  if (fitsread(tab->headbuf, "ZIMAGE  ", str, H_STRING, T_STRING) == RETURN_OK) {

          tab->isTileCompressed = 1;

          strcpy(NAXIS_KEYWORD, "ZNAXIS  ");
          strcpy(BITPIX_KEYWORD, "ZBITPIX ");
  }
  else {
          tab->isTileCompressed = 0;

          strcpy(NAXIS_KEYWORD, "NAXIS   ");
          strcpy(BITPIX_KEYWORD, "BITPIX  ");
  }

  if (fitsread(tab->headbuf, BITPIX_KEYWORD, &tab->bitpix, H_INT, T_LONG)
	==RETURN_ERROR)
    error(EXIT_FAILURE, "*Error*: Corrupted FITS header in ", filename);

  tab->bytepix = tab->bitpix>0?(tab->bitpix/8):(-tab->bitpix/8);

  if (fitsread(tab->headbuf, NAXIS_KEYWORD, &tab->naxis, H_INT, T_LONG)
                  ==RETURN_ERROR) {
    error(EXIT_FAILURE, "*Error*: Corrupted FITS header in ", filename);
  }


  tabsize = 0;
  if (tab->naxis>0)
    {
    QFREE(tab->naxisn);
    QMALLOC(tab->naxisn, int, tab->naxis);
/*--get the size of the array*/
    tabsize = 1;
    for (i=0; i<tab->naxis && i<999; i++)
      {
       // CFITSIO
       if (tab->isTileCompressed)
               sprintf(key,"ZNAXIS%-2d", (i+1));
       else
               sprintf(key,"NAXIS%-3d", (i+1));
      if (fitsread(tab->headbuf, key, &tab->naxisn[i], H_INT, T_LONG)
		==RETURN_ERROR)
        error(EXIT_FAILURE, "*Error*: incoherent FITS header in ", filename);
      tabsize *= tab->naxisn[i];
      }
    }

/*random groups parameters (optional)*/
  tab->pcount = 0;
  fitsread(tab->headbuf, "PCOUNT  ", &tab->pcount, H_INT, T_LONG);

  // CFITSIO TODO HACK
  if (tab->isTileCompressed) tab->pcount = 0;

  tab->gcount = 1;
  fitsread(tab->headbuf, "GCOUNT  ", &tab->gcount, H_INT, T_LONG);

/*number of fields (only for tables)*/
  tab->tfields = 0;
  fitsread(tab->headbuf, "TFIELDS ", &tab->tfields, H_INT, T_LONG);

/*in case of a non-primary header*/
  tab->xtension[0] = (char)'\0';
  fitsread(tab->headbuf, "XTENSION", tab->xtension, H_STRING, T_STRING);
  tab->extname[0] = (char)'\0';
  fitsread(tab->headbuf, "EXTNAME ", tab->extname, H_STRING, T_STRING);

  tab->tabsize = tab->bytepix*tab->gcount*((size_t)tab->pcount+tabsize);

/* Scaling parameters for basic FITS integer arrays */
  tab->bscale = 1.0;
  fitsread(tab->headbuf, "BSCALE ", &tab->bscale, H_FLOAT, T_DOUBLE);
  tab->bzero = 0.0;
  fitsread(tab->headbuf, "BZERO  ", &tab->bzero, H_FLOAT, T_DOUBLE);
  tab->blankflag =
    (fitsread(tab->headbuf,"BLANK   ",&tab->blank,H_INT,T_LONG) == RETURN_OK)?
	1 : 0;

/* Custom basic FITS parameters */
  tab->bitsgn = (tab->bitpix==BP_BYTE) ? 0 : 1;
  fitsread(tab->headbuf, "BITSGN  ", &tab->bitsgn, H_INT, T_LONG);

  if (fitsread(tab->headbuf, "IMAGECOD", str, H_STRING, T_STRING)==RETURN_OK)
    {
    if (!strcmp(str, "NONE"))
      tab->compress_type = COMPRESS_NONE;
    else if (!strcmp(str, "BASEBYTE"))
      tab->compress_type = COMPRESS_BASEBYTE;
    else if (!strcmp(str, "PREV_PIX"))
      tab->compress_type = COMPRESS_PREVPIX;
    else
      warning("Compression skipped: unknown IMAGECOD parameter:", str);
    }

/* Checksum */
  if (fitsread(tab->headbuf, "DATASUM ", str, H_STRING, T_STRING)==RETURN_OK)
    tab->bodysum = (unsigned int)atoi(str);

  return;
  }


/******* readbintabparam_head *************************************************
PROTO	int readbintabparam_head(tabstruct *tab)
PURPOSE	Read the current FITS header parameters concerning the binary-table.
INPUT	pointer to tabstruct.
OUTPUT	RETURN_OK if a binary table was found and mapped, RETURN_ERROR
	otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	20/07/2010
 ***/
int	readbintabparam_head(tabstruct *tab)

  {
   catstruct	*cat;
   keystruct	*key, *prevkey;
   char		strf[88], strk[16];
   char		*str;
   int		naxisn[32];
   int		i,j, larray, nfields,narray, pos;

  if (!(cat = tab->cat))
    error(EXIT_FAILURE, "*Internal Error*: Table has no parent catalog","!");

/*We are expecting a 2D binary-table, and nothing else*/
  if ((tab->naxis != 2)
	|| (tab->bitpix!=8)
	|| (tab->tfields == 0)
	|| strncmp(tab->xtension, "BINTABLE", 8))
    return RETURN_ERROR;

/*Size and number of lines in the binary table*/
  larray = tab->naxisn[0];
  nfields= tab->nkey = tab->tfields;
  narray = tab->naxisn[1];

  prevkey = NULL;
/*For each of the data fields...*/
  pos = 0;
  for (i=0; i<nfields; i++)
    {
/*--manage the chaining of keys*/
    QCALLOC(key, keystruct, 1);
    if (prevkey)
       {
       prevkey->nextkey = key;
       key->prevkey = prevkey;
       }
    else
       tab->key = key;
     prevkey = key;

/*--map binary-table fields*/

    sprintf(strk, "TTYPE%-3d", i+1);
    if (fitsread(tab->headbuf, strk, key->name, H_STRING, T_STRING)
	!= RETURN_OK) {
      error(EXIT_FAILURE,
	"*Error*: Incorrect FITS binary-table header in ", cat->filename); 
    }
    fitsread(tab->headbuf, strk, key->comment, H_HCOMMENT, T_STRING);

    sprintf(strk, "TUNIT%-3d", i+1);
    fitsread(tab->headbuf, strk, key->unit, H_STRING, T_STRING);
    sprintf(strk, "TDISP%-3d", i+1);
    fitsread(tab->headbuf, strk, key->printf, H_STRING, T_STRING);
    if (*key->printf)
      tdisptoprintf(key->printf, key->printf);

    sprintf(strk, "TFORM%-3d", i+1);
    if (fitsread(tab->headbuf, strk, strf, H_STRING, T_STRING) != RETURN_OK) {
      error(EXIT_FAILURE,
	"*Error*: Incorrect FITS binary-table header in ", cat->filename); 
    }
    key->pos = pos;
    pos += (key->nbytes = tsizeof(strf));
    key->ttype = ttypeof(strf);
    switch(key->ttype)
      {
      case T_BYTE:
      case T_SHORT:
      case T_LONG:
      case T_LONGLONG:
        key->htype = H_INT;
        break;
      case T_FLOAT:
      case T_DOUBLE:
        key->htype = H_EXPO;
        break;
      case T_STRING:
        key->htype = H_STRING;
        break;
      default:
        // CFITSIO TODO dodgy
    	key->ttype = T_FLOAT;
        //error(EXIT_FAILURE, "*Error*: Unknown TFORM in ", cat->filename);
      }

/*--handle the special case of multimensional arrays*/
    if ((naxisn[0] = key->nbytes/t_size[key->ttype]) > 1)
      {
      sprintf(strk, "TDIM%-3d", i+1);
      if (fitsread(tab->headbuf, strk, strf, H_STRING, T_STRING) == RETURN_OK)
        {
        str = strf;
        for (j=0; (naxisn[j]=(int)strtol(str+1, &str, 10)); j++);
        key->naxis = j;
        }
      else
        key->naxis = 1;
      QMALLOC(key->naxisn, int, key->naxis);
      for (j=0; j<key->naxis; j++)
        key->naxisn[j] = naxisn[j];
      }
    else
      key->naxis = 0;

    key->nobj = narray;
    key->tab = tab;
    }

  if (pos != larray)
    error(EXIT_FAILURE,
	"*Error*: Malformed FITS binary-table header in ", cat->filename); 

/*make both ends of the chain meet*/
  prevkey->nextkey = tab->key;
  tab->key->prevkey = prevkey;

  return RETURN_OK;
  }


/****** update_head ***********************************************************
PROTO	int update_head(tabstruct *tab)
PURPOSE	Update a FITS header according to what's in the table.
INPUT	Table structure.
OUTPUT	RETURN_OK if tab is a binary table, or RETURN_ERROR otherwise.
NOTES	The headbuf pointer in the tabstruct might be reallocated.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	11/06/2007
 ***/
int	update_head(tabstruct *tab)

  {
   keystruct	*key;
   tabstruct	*ctab;
   int		i,j,n,naxis1;
   char		strk[88], str[88];
   char		*buf;

/*Update EXTNAME, the table name */
  if (*tab->extname)
    {
    addkeywordto_head(tab, "EXTNAME ", "EXTENSION NAME");
    fitswrite(tab->headbuf, "EXTNAME ", tab->extname, H_STRING, T_STRING);
    }

/* If not a binary table, do only a few basic things */
  if ((tab->naxis != 2)
	|| (tab->bitpix!=8)
	|| (tab->tfields == 0)
	|| strncmp(tab->xtension, "BINTABLE", 8))
    {
    addkeywordto_head(tab, "BITPIX  ", "BITS PER PIXEL");
    fitswrite(tab->headbuf, "BITPIX  ", &tab->bitpix, H_INT, T_LONG);
    addkeywordto_head(tab, "NAXIS   ", "NUMBER OF AXES");
    fitswrite(tab->headbuf, "NAXIS   ", &tab->naxis, H_INT, T_LONG);
    for (i=0; i<tab->naxis; i++)
      {
      sprintf(strk, "NAXIS%-3d", i+1);
      addkeywordto_head(tab, strk, "NUMBER OF ELEMENTS ALONG THIS AXIS");
      fitswrite(tab->headbuf, strk, &tab->naxisn[i], H_INT, T_LONG);
      }
    return RETURN_ERROR;
    }

/*First, remove all existing TTYPE, TFORM, etc...*/
  removekeywordfrom_head(tab, "TTYPE???");
  removekeywordfrom_head(tab, "TFORM???");
  removekeywordfrom_head(tab, "TUNIT???");
  removekeywordfrom_head(tab, "TZERO???");
  removekeywordfrom_head(tab, "TSCAL???");
  removekeywordfrom_head(tab, "TDIM???");
  removekeywordfrom_head(tab, "TDISP???");


/*Change NAXIS1 in order to take into account changes in width*/
  naxis1 = 0;
  key = tab->key;
  if (tab->nkey>1000) {
     for (i=0; i<MIN(999,tab->nkey); i++) {
        naxis1 += key->nbytes;
        key = key->nextkey;
     }
     fitswrite(tab->headbuf, "NAXIS1  ", &naxis1, H_INT, T_LONG);
  } else {
     fitswrite(tab->headbuf, "NAXIS1  ", &tab->naxisn[0], H_INT, T_LONG);
  }

/*Change NAXIS1 in the number of fields */
  tab->tfields = MIN(999,tab->tfields);
  fitswrite(tab->headbuf, "TFIELDS ", &tab->tfields, H_INT, T_LONG);

/*Changes in the number of elements (look for possible segments)*/
  for (ctab = tab, n = ctab->naxisn[1];
	(ctab=ctab->nexttab) && !ctab->nseg;)
    n += ctab->naxisn[1];
  fitswrite(tab->headbuf, "NAXIS2  ", &n, H_INT, T_LONG);

  key = tab->key;
  if (!key)
    return RETURN_ERROR;

  if (tab->nkey>1000)
     warning("Too many output keys, trashing the ones bejond 999", "");
  for (i=0; i<MIN(999,tab->nkey); i++)
    {
    sprintf(strk, "TTYPE%-3d", i+1);
    addkeywordto_head(tab, strk, key->comment);
    fitswrite(tab->headbuf, strk, key->name, H_STRING, T_STRING);
    sprintf(strk, "TFORM%-3d", i+1);
    addkeywordto_head(tab, strk, "");
    tformof(str, key->ttype, key->nbytes/t_size[key->ttype]);
    fitswrite(tab->headbuf, strk, str, H_STRING, T_STRING);
    if (key->naxis>1)
      {
       char	*str2, *str2lim;

      sprintf(strk, "TDIM%-3d", i+1);
      addkeywordto_head(tab, strk, "");
      sprintf(str, "(");
      str2 = str+1;
      str2lim = str+70;	/* Prevent an excessively large string */
      for (n=0; n<key->naxis && str2<str2lim; n++)
        {
        sprintf(str2, n?", %d%n":"%d%n", key->naxisn[n],&j);
        str2 += j;
        }
      sprintf(str2, ")");
      fitswrite(tab->headbuf, strk, str, H_STRING, T_STRING);
      }
    if (*key->unit)
      {
      sprintf(strk, "TUNIT%-3d", i+1);
      addkeywordto_head(tab, strk, "");
      fitswrite(tab->headbuf, strk, key->unit, H_STRING, T_STRING);
      }
    if (*key->printf)
      {
      sprintf(strk, "TDISP%-3d", i+1);
      addkeywordto_head(tab, strk, "");
      fitswrite(tab->headbuf, strk, printftotdisp(key->printf, str),
		H_STRING, T_STRING);
      }
    key = key->nextkey;
    }

/*Finally re-compute CHECKSUM if present */
  if (fitsfind(tab->headbuf, "CHECKSUM")==RETURN_OK)
    {
    unsigned int	sum;

    if (tab->bodysum)
      {
      sprintf(str, "%u", tab->bodysum);
      fitswrite(tab->headbuf, "DATASUM ", str, H_STRING, T_STRING);
      }
    sum = tab->bodysum;
/*-- Now the header */
    buf = tab->headbuf;
    for (i=tab->headnblock; i--; buf+=FBSIZE)
      sum = compute_blocksum(buf, sum);
/*-- Complement to 1 */
    encode_checksum(~sum, str);
    fitswrite(tab->headbuf, "CHECKSUM", str, H_STRING, T_STRING);
    }

/*That may be enough for now; to be continued...*/

  return RETURN_OK;
  }


/****** prim_head *************************************************************
PROTO	int prim_head(tabstruct *tab)
PURPOSE	Update a FITS header to make it "primary" (not extension)
INPUT	Table structure.
OUTPUT	RETURN_OK if tab header was already primary, or RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP) C. Marmo (IAP)
VERSION	30/08/2011
 ***/
int	prim_head(tabstruct *tab)

  {
  if (!tab->headbuf)
    return RETURN_ERROR;
  if (!strncmp(tab->headbuf, "XTENSION",8))
      {
      strncpy(tab->headbuf, "SIMPLE  =                    T "
	"/ This is a FITS file                            ", 80);
      removekeywordfrom_head(tab, "PCOUNT");      
      removekeywordfrom_head(tab, "GCOUNT");      
      removekeywordfrom_head(tab, "TFIELDS");      
      removekeywordfrom_head(tab, "EXTNAME");      
      *tab->extname = '\0';
      return RETURN_ERROR;
      }

  return RETURN_OK;
  }


/****** ext_head *************************************************************
PROTO	int ext_head(tabstruct *tab)
PURPOSE	Update a FITS header to make it "extension" (not primary)
INPUT	Table structure.
OUTPUT	RETURN_OK if tab header was already extension, or RETURN_ERROR
	otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory) C. Marmo (IAP)
VERSION	20/06/2007
 ***/
int	ext_head(tabstruct *tab)

  {
  if (!tab->headbuf)
    return RETURN_ERROR;
  if (!strncmp(tab->headbuf, "SIMPLE  ",8))
      {
      strncpy(tab->headbuf, "XTENSION= 'IMAGE   '           "
		"/ Image extension                                ", 80);
/* fitsverify 4.13 (CFITSIO V3.002) return an error
   if EXTEND are in an extension header (20/06/2007)*/
      removekeywordfrom_head(tab, "EXTEND");      
/* fitsverify 4.13 (CFITSIO V3.002) return an error
   if PCOUNT and GCOUNT are not in the extension header (23/05/2007) */
      addkeywordto_head(tab, "PCOUNT  ", "required keyword; must = 0");      
      addkeywordto_head(tab, "GCOUNT  ", "required keyword; must = 1");
      fitswrite(tab->headbuf,"PCOUNT  ", &tab->pcount, H_INT, T_LONG);     
      fitswrite(tab->headbuf,"GCOUNT  ", &tab->gcount, H_INT, T_LONG);
      return RETURN_ERROR;
      }

  return RETURN_OK;
  }


/****** addkeyto_head *********************************************************
PROTO	int addkeyto_head(tabstruct *tab, keystruct *key)
PURPOSE	Add a keyword and its value to a table header.
INPUT	Table structure,
	Key containing the keyword and its value.
OUTPUT	Line position in the FITS header.
NOTES	The headbuf pointer in the tabstruct might be reallocated.
	Pre-existing keywords are overwritten (but not their comments).
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	11/05/2002
 ***/
int	addkeyto_head(tabstruct *tab, keystruct *key)

  {
   int	n;

  n = addkeywordto_head(tab, key->name, key->comment);
  fitswrite(tab->headbuf, key->name, key->ptr, key->htype, key->ttype);

  return n;
  }


/****** addkeywordto_head *****************************************************
PROTO	int addkeywordto_head(tabstruct *tab, char *keyword, char *comment)
PURPOSE	Add a keyword and a comment to a table header.
INPUT	Table structure,
	String containing the keyword,
	String containing the comment.
OUTPUT	Line position in the FITS header.
NOTES	The headbuf pointer in the tabstruct might be reallocated.
	Pre-existing keywords are overwritten (but not their comments).
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	21/04/2003
 ***/
int	addkeywordto_head(tabstruct *tab, char *keyword, char *comment)

  {
   int	n;

  if ((fitsfind(tab->headbuf, keyword) == RETURN_ERROR
	|| findkey(keyword, (char *)histokeys, 12)!=RETURN_ERROR)
	&& (fitsfind(tab->headbuf, "END     ")+1)*80 >= tab->headnblock*FBSIZE)
    {
    tab->headnblock++;
    QREALLOC(tab->headbuf, char, tab->headnblock*FBSIZE);
    memset(tab->headbuf + (tab->headnblock-1)*FBSIZE, ' ', FBSIZE);
    }

  n = fitsadd(tab->headbuf, keyword, comment);

  return n;
  }


/****** removekeywordfrom_head ************************************************
PROTO	int removekeywordfrom_head(tabstruct *tab, char *keyword)
PURPOSE	Remove a keyword from a table header.
INPUT	Table structure,
	String containing the keyword.
OUTPUT	RETURN_OK if the keyword was found, RETURN_ERROR otherwise..
NOTES	The headbuf pointer in the tabstruct might be reallocated.
        '?' wildcard allowed; Don't remove the ``END'' keyword with this!!!
AUTHOR	E. Bertin (IAP)
VERSION	11/06/2007
 ***/
int	removekeywordfrom_head(tabstruct *tab, char *keyword)

  {
   int	nb;

  if (fitsremove(tab->headbuf, keyword) == RETURN_OK)
    {
    if ((nb=fitsfind(tab->headbuf, "END     ")/(FBSIZE/80)+1) < tab->headnblock)
      {
      tab->headnblock = nb;
      QREALLOC(tab->headbuf, char, tab->headnblock*FBSIZE);
      }
    return RETURN_OK;
    }
  else
    return RETURN_ERROR;
  }


/****** tformof ***************************************************************
PROTO	int tformof(char *str, t_type ttype, int n)
PURPOSE	Return the ``TFORM'' string corresponding to a t_type
	and the number of elements.
INPUT	a char pointer (to be filled with the T_FORM string),
	t_type,
	Number of elements.
OUTPUT	RETURN_OK if everything went as expected, or RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	28/10/2009
 ***/
int	tformof(char *str, t_type ttype, int n)

  {
   char	t;

  switch (ttype)
    {
    case T_BYTE:	t = 'B';
			break;
    case T_SHORT:	t = 'I';
			break;
    case T_LONG:	t = 'J';
			break;
    case T_LONGLONG:	t = 'K';
			break;
    case T_FLOAT:	t = 'E';
			break;
    case T_DOUBLE:	t = 'D';
			break;
    case T_STRING:	t = 'A';
			break;
    default:		return	RETURN_ERROR;
    }

  sprintf(str, "%d%c", n, t);

  return RETURN_OK;
  }


/****** tsizeof ***************************************************************
PROTO	int tsizeof(char *str)
PURPOSE	Return the size of a binary-table field from its ``TFORM''.
INPUT	TFORM string (see the FITS documentation).
OUTPUT	size in bytes, or RETURN_ERROR if the TFORM is unknown.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	10/11/2010
 ***/
int	tsizeof(char *str)

  {
   int	n;
   char	*str2;

  str2 = str;
  n = strtol(str, &str2, 10);
  if (str2==str)
    n = 1;

  switch ((int)*str2)
    {
    case 'L': case 'B': case 'A':		return	n;
    case 'X':					return	(n-1)/8+1;
    case 'I':					return	2*n;
    case 'J': case 'E':				return	4*n;
    case 'C': case 'D': case 'K': case 'P':	return	8*n;
    case 'M':					return	16*n;
    default:					return	RETURN_ERROR;
    }

  }


/****** ttypeof ***************************************************************
PROTO	t_type ttypeof(char *str)
PURPOSE	Give the ``t_type'' of a binary-table field from its ``TFORM''.
INPUT	TFORM string (see the FITS documentation).
OUTPUT	size in bytes, or RETURN_ERROR if the TFORM is unknown.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	29/08/2012
 ***/
t_type	ttypeof(char *str)

  {
   char	*str2;

  str2 = str;
  strtol(str, &str2, 10);
  switch ((int)*str2)
    {
    case 'L': case 'B': case 'X':	return	T_BYTE;
    case 'I':				return	T_SHORT;
    case 'J':				return	T_LONG;
    case 'K':				return	T_LONGLONG;
    case 'E':				return	T_FLOAT;
    case 'D':				return	T_DOUBLE;
    case 'A':				return	T_STRING;
    default:				return	(t_type)RETURN_ERROR;
    }

  }


/****** tdisptoprintf *********************************************************
PROTO	char	*tdisptoprintf(char *tdisp, char *str)
PURPOSE	Convert the ``TDISP'' FITS format to the printf() format.
INPUT	TDISP format string (see the FITS documentation),
	output string (allocated pointer).
OUTPUT	printf() format string (see e.g.  K&R).
NOTES	The present conversion does not handle binary or engineer notations.
	A NULL vector is returned if the conversion was unsuccessful.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	25/09/2004
 ***/
char	*tdisptoprintf(char *tdisp, char *str)

  {
   char		control[4];
   int		w,d, n;

  w = d = 0;
  n = 0;
  n=sscanf(tdisp,"%[ALIBOZFENSGD]%d.%d", control, &w, &d)-1;
  if (!w)
    {
    warning("Strange TDISP format: ", tdisp);
    return NULL;
    }
  switch ((int)*control)
    {
    case 'A':
      sprintf(str, "%%%dc",w);
      break;
    case 'L':
      sprintf(str, "%%%dd",w);
      break;      
    case 'I':
      if (n>1)
        sprintf(str, "%%%d.%dd",w,d);
      else
        sprintf(str, "%%%dd",w);
      break;
    case 'B': case 'Z':
      if (n>1)
        sprintf(str, "%%%d.%dx",w,d);
      else
        sprintf(str, "%%%dx",w);
      break;
    case 'O':
      if (n>1)
        sprintf(str, "%%%d.%do",w,d);
      else
        sprintf(str, "%%%do",w);
      break;
    case 'F':
      if (n>1)
        sprintf(str, "%%%d.%df",w,d);
      else
        sprintf(str, "%%%df",w);
      break;
    case 'E': case 'D': 
      if (n>1)
        sprintf(str, "%%%d.%dE",w,d);
      else
        sprintf(str, "%%%dE",w);
      break;
    case 'G': 
      if (n>1)
        sprintf(str, "%%%d.%dG",w,d);
      else
        sprintf(str, "%%%dG",w);
      break;
    default:
      warning("Unknown TDISP format: ", tdisp);
      return NULL;
    }

  return str;
  }


/****** printftotdisp *********************************************************
PROTO	char	*printftotdisp(char *tdisp, char *str)
PURPOSE	Convert the printf() format to the ``TDISP'' FITS format.
INPUT	printf() format string (see e.g.  K&R),
	output string (allocated pointer).
OUTPUT	TDISP format string (see the FITS documentation).
NOTES	The handling of C string formatting does not include the precision.
	NULL is returned in case of unsucessful conversion.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	25/09/2004
 ***/
char	*printftotdisp(char *cprintf, char *str)

  {
   char		*control;
   int		w,d,n;

  *str = 0;
  w = d = 0;
  if (!(control = strpbrk(cprintf, "cdueERfFgGoOxXs")))
    {
    warning("Unknown printf() format: ", cprintf);
    return NULL;
    }

  n = sscanf(cprintf,"%%%d.%d", &w, &d);
  w = abs(w);
  if (!n)
    {
    warning("Unconvertible printf() format: ", cprintf);
    return NULL;
    }

  switch ((int)*control)
    {
    case 'c':
      sprintf(str, "A%d",w);
      break;
    case 's':
      sprintf(str, "A%d",w);
      break;
    case 'd': case 'u':
      if (n>1)
        sprintf(str, "I%d.%d",w,d);
      else
        sprintf(str, "I%d",w);
      break;
    case 'o': case 'O':
      if (n>1)
        sprintf(str, "O%d.%d",w,d);
      else
        sprintf(str, "O%d",w);
      break;
    case 'x': case 'X':
      if (n>1)
        sprintf(str, "Z%d.%d",w,d);
      else
        sprintf(str, "Z%d",w);
      break;
    case 'f': case 'F':
      if (n>1)
        sprintf(str, "F%d.%d",w,d);
      else
        sprintf(str, "F%d",w);
      break;
    case 'e': case 'E':
      if (n>1)
        sprintf(str, "E%d.%d",w,d);
      else
        sprintf(str, "E%d",w);
      break;
    case 'g': case 'G':
      if (n>1)
        sprintf(str, "G%d.%d",w,d);
      else
        sprintf(str, "G%d",w);
      break;
    default:
      warning("Unknown printf() format: ", cprintf);
      return NULL;
    }

  return str;
  }

