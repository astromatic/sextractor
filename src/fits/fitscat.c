/*
*				fitscat.c
*
* Low-level functions for handling FITS images and tables.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic FITS/LDAC library
*
*	Copyright:		(C) 1994,1997 ESO
*	          		(C) 1995,1996 Leiden Observatory 
*	          		(C) 1998-2021 IAP/CNRS/SorbonneU
*	          		(C) 2021-2023 CFHT/CNRS
*	          		(C) 2023-2025 CEA/AIM/UParisSaclay
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
*	Last modified:		21/03/2025
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<sys/types.h>
#include	<sys/stat.h>
#include	<fcntl.h>
#include	<time.h>

#include	"fitscat_defs.h"
#include	"fitscat.h"

extern	float *fits_rand_value;

char		fits_str[MAXCHAR];

/****** about_cat **************************************************************
PROTO	int about_cat(catstruct *cat, FILE *stream)
PURPOSE	Print some info about a catalog.
INPUT	Catalog structure,
	output stream.
OUTPUT	RETURN_OK if everything went as expected, RETURN_ERROR otherwise.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	19/03/2002
 ***/
int	about_cat(catstruct *cat, FILE *stream)

  {
   tabstruct	*tab;
   int		i;

  fprintf(stream,"\n");

/*General info about the catalog itself*/
  fprintf(stream,
	"------------------Catalog information----------------\n");
  fprintf(stream,
	"Filename:..............%s\n", cat->filename);
  fprintf(stream,
	"Number of segments:....%d\n", cat->ntab);
  fprintf(stream,"\n");

/*Now for each table*/
  tab = cat->tab;
  for (i=0; i<cat->ntab; i++)
    {
    fprintf(stream,
	"******	Table #%d\n", i+1);
    fprintf(stream,
	"	Extension type:.........%s\n",
	tab->xtension[0]? tab->xtension: "(Primary HDU)");
    fprintf(stream,
	"	Extension name:.........%s\n", tab->extname);
    if (tab->naxis)
      {
      fprintf(stream,
	"	Number of dimensions:...%d\n", tab->naxis);
      fprintf(stream,
	"	Number of elements:.....%d\n", tab->naxisn[1]);
      if (tab->tfields)
        fprintf(stream,
	"	Number of data fields...%d\n", tab->tfields);
      fprintf(stream,
	"	Body size:..............%ld bytes\n",
		(unsigned long)tab->tabsize);
      }
    fprintf(stream,"\n");
    while (!(tab=tab->nexttab)->nseg);
    }

  fprintf(stream,"\n");

  return RETURN_OK;
  }


/****** addhistoryto_cat *******************************************************
PROTO	int addhistoryto_cat(catstruct *cat, char *str)
PURPOSE	Add a HISTORY line to a FITS catalog.
INPUT	A pointer to catalog structure, and the character string to insert.
OUTPUT	RETURN_OK if everything went as expected, RETURN_ERROR otherwise.
NOTES	The pointer to the primary header might be reallocated if necessary.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	25/09/2004
 ***/
int	addhistoryto_cat(catstruct *cat, char *str)

  {
   time_t		thetime;
   char			str2[82];
   tabstruct		*tab;
   int			n, headpos;

  tab = cat->tab;
  n = tab->headnblock;
  headpos = fitsfind(tab->headbuf, "END     ");
  if (headpos >= n*(FBSIZE/80) - 1)
    {
    QREALLOC(tab->headbuf, char, (n+1)*FBSIZE);
    memset(&tab->headbuf[n*FBSIZE], ' ', FBSIZE);
    tab->headnblock++;
    }

  if (time(&thetime)==-1)
    warning("No time available for history","");

  if (!strftime(str2, 16, "%d/%m/%Y %H:%M", localtime(&thetime)))
    error(EXIT_FAILURE, "*Internal Error*: Time/date string too long in ",
	"addhistoryto_cat()");
  sprintf(str2, "%s %.65s", str2, str);
  fitsadd(tab->headbuf, "HISTORY ", str2);

  return RETURN_OK;
  }


/****** close_cat **************************************************************
PROTO	int close_cat(catstruct *cat)
PURPOSE	Close a FITS catalog.
INPUT	catalog structure.
OUTPUT	RETURN_OK if everything went as expected, RETURN_ERROR otherwise.
NOTES	the file structure member is set to NULL;
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	21/03/2025
 ***/
int	close_cat(catstruct *cat)

  {
   int status = RETURN_OK;
  
  if (cat->file && fclose(cat->file))
    status = RETURN_ERROR;   

  cat->file = NULL;

  if (cat->cfitsio_infptr && close_cfitsio(cat))
    status = RETURN_ERROR;

  cat->cfitsio_infptr = NULL;

  return status;
  }


#ifdef	HAVE_CFITSIO
/****** open_cfitsio ***********************************************************
PROTO	int close_cfitsio(catstruct *cat)
PURPOSE	Open a file with cfitsio 
INPUT	Pointer to the catalog.
    access type (can be WRITE_ONLY or READ_ONLY).
OUTPUT	RETURN_OK if everything went as expected (exit in error otherwise).
NOTES	-.
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	25/03/2025
***/
int	open_cfitsio(catstruct *cat, access_type_t at) {

   int status;

  if (cat->cfitsio_flag && !cat->cfitsio_infptr) {
    // Trigger CFITSIO file opening
    fits_open_file(
        &cat->cfitsio_infptr,
        cat->filename,
        at==WRITE_ONLY? 1 : 0,
        &status
    );
    if (status != 0) {
      fits_report_error(stderr, status);
      return RETURN_ERROR;
    }
  }

  return RETURN_OK;
}


/****** close_cfitsio **********************************************************
PROTO	int close_cfitsio(catstruct *cat)
PURPOSE	Close a file previously opened by cfitsio 
INPUT	Catalog structure.
OUTPUT	RETURN_OK if everything went as expected (exit in error otherwise).
NOTES	- .
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	21/03/2025
***/
int	close_cfitsio(catstruct *cat) {

  if ((cat) && (cat->cfitsio_infptr)) {
    int status = 0; fits_close_file(cat->cfitsio_infptr, &status);

    if (status != 0) {
      fits_report_error(stderr, status);
      error(EXIT_FAILURE, "Could not close FITS file with cfitsio: ",
      	cat->filename);
    } else {
      // Successfully closed FITS file with cfitsio
      cat->cfitsio_infptr = NULL;
      // Free random seed in CFITSIO.
      QFREE(fits_rand_value);
    }
  } else
    return RETURN_ERROR;
  return RETURN_OK;
}
#endif // HAVE_CFITSIO


/****** free_cat ***************************************************************
PROTO	void free_cat(catstruct **cat, int ncat)
PURPOSE	Free all structures allocated for one or several FITS catalog.
INPUT	Pointer to a catalog structure,
	Number of catalogs.
OUTPUT	-.
NOTES	Unallocated pointers should have been put to NULL.
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	21/03/2025
 ***/
void	free_cat(catstruct **cat, int ncat)

  {
   catstruct	**thecat;
   int		i;

/*--free memory allocated within each catalog */
  thecat = cat;
  for (i=ncat; i--;)
    {
    if ((*thecat)->file)
      close_cat(*thecat);
    remove_tabs(*thecat);
    free(*(thecat++));
    }


  return;
  }


/****** inherit_cat ************************************************************
PROTO	int inherit_cat(catstruct *catin, catstruct *catout)
PURPOSE	Copy the primary table, and all other informations from one catalog
	to another, except those related to the associated file itself
	(filename, etc...),
INPUT	A pointer to both catalog structures.
OUTPUT	RETURN_OK if at least one table was copied, RETURN_ERROR otherwise.
NOTES	The output catalog should be ``cleaned'' before call.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	13/06/2002
 ***/
int	inherit_cat(catstruct *catin, catstruct *catout)

  {
   tabstruct	*tabin, *tabout, *prevtabout;
   int		j;

  catout->ntab = 1;
  tabin = catin->tab;

/*copy only one table: well it could be simpler, but let's stay general!*/
  prevtabout = NULL;
  for (j=tabin->nseg; j--;)
    {
    QCALLOC(tabout, tabstruct, 1);
    *tabout = *tabin;
    if (tabin->naxis)
      QMEMCPY(tabin->naxisn, tabout->naxisn, int, (size_t)tabin->naxis);
    if (tabin->headbuf)
      QMEMCPY(tabin->headbuf, tabout->headbuf, char,
	tabin->headnblock*FBSIZE);
    if (tabin->bodybuf)
      QMEMCPY(tabin->bodybuf, tabout->bodybuf, char, (size_t)tabin->tabsize);
    if (prevtabout)
      {
      tabout->prevtab = prevtabout;
      prevtabout->nexttab = tabout;
      }
    else
      {
      catout->tab = tabout;
      }
    prevtabout = tabout;
    tabin = tabin->nexttab;
    }

  if (prevtabout)
    {
    prevtabout->nexttab = catout->tab;
    catout->tab->prevtab = prevtabout;
    }
  else
    return RETURN_ERROR;

  return RETURN_OK;
  }

/****** init_cat ***************************************************************
PROTO	int init_cat(catstruct *cat)
PURPOSE	Initialize a catalog, "cleaning" any content if present
	 and adding the primary header "table".
INPUT	A pointer to the catalog structure.
OUTPUT	RETURN_OK if everything went as expected, RETURN_ERROR otherwise.
NOTES	The output catalog should be ``cleaned'' before call.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	28/05/2001
 ***/
int	init_cat(catstruct *cat)

  {
   static char	bintabtemplate[][80] = {
"SIMPLE  =                    T / This is a FITS file",
"BITPIX  =                    8 / ",
"NAXIS   =                    0 / ",
"EXTEND  =                    T / This file may contain FITS extensions",
"END                            "};
   tabstruct	*tab;
   char		*buf;
   int		i;

/* Initialize the primary header itself */
  QCALLOC(tab, tabstruct, 1);
  tab->naxis = 0;
  tab->bitpix = 8;
  tab->bytepix = 1;
  tab->pcount = 0;
  tab->gcount = 1;
  tab->seg = 1;
  tab->nseg = 1;
/* Provide a new header*/
  QCALLOC(tab->headbuf, char, FBSIZE);
  memcpy(tab->headbuf, bintabtemplate, sizeof(bintabtemplate));
  for (buf = tab->headbuf, i=0; i<FBSIZE; i++, buf++)
    if (!*buf)
      *buf = ' ';
  tab->headnblock = 1;
/* Clean catalog and add the table to it */
  remove_tabs(cat);
  cat->tab = tab->prevtab = tab->nexttab = tab;
  cat->ntab = 1;

  return RETURN_OK;
  }


/****** map_cat ****************************************************************
PROTO	int map_cat(catstruct *cat)
PURPOSE	Explores the whole FITS file
	and gets information for each of the FITS tables it contains.
INPUT	catalog structure.
OUTPUT	RETURN_OK if at least one table was found, RETURN_ERROR otherwise.
NOTES	Memory space for the array of fits structures is reallocated.
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	21/03/2025
 ***/
int	map_cat(catstruct *cat)

  {
   tabstruct		*tab, *prevtab;
   int			ntab;
#ifdef	HAVE_CFITSIO
   int  hdunum = 1;
#endif // HAVE_CFITSIO

/*scan through the file until we reach the end*/
  prevtab = NULL;
  QCALLOC(tab, tabstruct, 1);
  tab->cat = cat;
  QFTELL(cat->file, tab->headpos, cat->filename);


  for (ntab=0; !get_head(tab); ntab++)
    {
    readbasic_head(tab);
    readbintabparam_head(tab);
    QFTELL(cat->file, tab->bodypos, cat->filename);
    tab->nseg = tab->seg = 1;

#ifdef	HAVE_CFITSIO
    if (tab->isTileCompressed && !cat->cfitsio_flag) {
      cat->cfitsio_flag = 1;
      open_cfitsio(cat, READ_ONLY);
   }
   tab->cfitsio_hdunum = hdunum++;
#endif // HAVE_CFITSIO
    if (tab->tabsize)
      QFSEEK(cat->file, PADTOTAL(tab->tabsize), SEEK_CUR, cat->filename);

    if (prevtab)
      {
      tab->prevtab = prevtab;
      prevtab->nexttab = tab;
      }
    else
      cat->tab = tab;
    prevtab = tab;
    QCALLOC(tab, tabstruct, 1);
    tab->cat = cat;
    QFTELL(cat->file, tab->headpos, cat->filename);
    }

  cat->ntab = ntab;
  free(tab);
  if (prevtab)
    {
    prevtab->nexttab = cat->tab;
    cat->tab->prevtab = prevtab;
    }
  else
    return RETURN_ERROR;

  return RETURN_OK;
  }


/****** new_cat ****************************************************************
PROTO	catstruct *new_cat(int ncat)
PURPOSE	Initialize a structure for a FITS catalog.
INPUT	Number of catalogs.
OUTPUT	A pointer to the catalog array.
NOTES	All fields are initialized to 0.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	20/03/96
 ***/
catstruct	*new_cat(int ncat)

  {
   catstruct	*cat;

  QCALLOC(cat, catstruct, ncat);

  cat->access_type = WRITE_ONLY;

  return cat;
  }


/****** open_cat ***************************************************************
PROTO	int open_cat(catstruct *cat, access_type_t at)
PURPOSE	Open a FITS catalog with name filename.
INPUT	catalog structure,
	access type (can be WRITE_ONLY or READ_ONLY).
OUTPUT	RETURN_OK if the cat is found, RETURN_ERROR otherwise.
NOTES	If the file was already opened by this catalog, nothing is done.
AUTHOR	E. Bertin (CEA/AIM/UParisSaclay)
VERSION	21/03/2025
 ***/
int	open_cat(catstruct *cat, access_type_t at)

  {
   int status = RETURN_OK;

  if  (cat->access_type == READ_ONLY && at == WRITE_ONLY)
    error(EXIT_FAILURE, "*Internal Error*: Trying to write to the "
	"READ_ONLY catalog ", cat->filename);

  if (!cat->file)
    {
    if ((cat->file = fopen(cat->filename, at==WRITE_ONLY?"wb":"rb")) == NULL)
      status = RETURN_ERROR;
    cat->access_type = at;
    }

#ifdef	HAVE_CFITSIO
  if (cat->cfitsio_flag && open_cfitsio(cat, at))
    status = RETURN_ERROR;
#endif // HAVE_CFITSIO

  return status;
  }


