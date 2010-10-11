/*
*				assoc.h
*
* Include file for assoc.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1998-2010 IAP/CNRS/UPMC
*				(C) 1997 ESO
*
*	Author:			Emmanuel Bertin (IAP)
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

#define		ASSOC_BUFINC	131072	/* Assoc buffer increment (bytes) */

/*--------------------------------- typedefs --------------------------------*/

typedef struct structassoc
  {
  float		*list;			/* Pointer to the list of data */
  int		nobj;			/* Number of data rows */
  int		ncol;			/* Total number of columns per row */
  int		ndata;			/* Number of retained cols per row */
  int		*hash;			/* Pointer to the hash table */
  float		*data;			/* Copy of current parameters */
  float		radius;			/* Radius of search for association */
  }             assocstruct;

/*------------------------------ Prototypes ---------------------------------*/

assocstruct	*load_assoc(char *filename);

int		do_assoc(picstruct *field, float x, float y);

void		init_assoc(picstruct *field),
		end_assoc(picstruct *field),
		sort_assoc(picstruct *field, assocstruct *assoc);
