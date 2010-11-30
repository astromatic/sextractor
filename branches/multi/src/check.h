/*
*				check.h
*
* Include file for check.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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

/*--------------------------------- structures ------------------------------*/
/* Check-image parameters */
typedef struct structcheck
  {
  char		filename[MAXCHAR];	/* ptr to check-image filename */
  FILE		*file;			/* ptr to check-image file structure */
  char		*fitshead;		/* ptr to check-image FITS header */
  int		fitsheadsize;		/* size of check-image FITS header */
  void		*pix;			/* ptr to check-image pixmap */
  int		width, height, depth;	/* size of check-image */
  size_t	npix;			/* number of pixels in check-image */
  int		y;			/* current line in check-image */
  PIXTYPE	overlay;		/* intensity of the overlayed plots */
  PIXTYPE	*line;			/* buffered image line */
  checkenum	type;			/* CHECKIMAGE_TYPE */
  }	checkstruct;

/*------------------------------- functions ---------------------------------*/

checkstruct	*initcheck(char *, checkenum, int next);

void		addcheck(checkstruct *, float *, int,int, int,int, float),
		blankcheck(checkstruct *, PIXTYPE *, int,int,int,int,PIXTYPE),
		endcheck(checkstruct *),
		reendcheck(picstruct *field, checkstruct *),
		reinitcheck(picstruct *, checkstruct *),
		writecheck(checkstruct *, PIXTYPE *, int);
