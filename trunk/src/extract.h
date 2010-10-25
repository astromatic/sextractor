/*
*				extract.h
*
* Include file for extract.c.
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

/*------------------------------ definitions --------------------------------*/

#define	NOBJ			256		/* starting number of obj. */
#define	UNKNOWN			-1		/* flag for LUTZ */

/*--------------------------------- typedefs --------------------------------*/

typedef	enum		{COMPLETE, INCOMPLETE, NONOBJECT, OBJECT}
				status;	/* Extraction status */

/*--------------------------------- variables -------------------------------*/
PIXTYPE		*dumscan;

/*------------------------------- structures --------------------------------*/
/* Temporary object parameters during extraction */
typedef struct structinfo
  {
  LONG		pixnb;			/* Number of pixels included */
  LONG		firstpix;		/* Pointer to first pixel of pixlist */
  LONG		lastpix;		/* Pointer to last pixel of pixlist */
  short		flag;			/* Extraction flag */
  }       infostruct;


/*------------------------------- functions ---------------------------------*/
void		lutzalloc(int, int),
		lutzfree(void),
		lutzsort(infostruct *, objliststruct *),
		sortit(picstruct *, picstruct *, picstruct *, picstruct *,
			infostruct *, objliststruct *, PIXTYPE *, PIXTYPE *),
		update(infostruct *, infostruct *, pliststruct *);

int		lutz(objliststruct *, int, objstruct *, objliststruct *); 
