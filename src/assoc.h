/*
*				assoc.h
*
* Include file for assoc.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1997-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		03/10/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSWCS_H_
#include        "fitswcs.h"
#endif

#define		ASSOC_BUFINC	131072	/* Assoc buffer increment (bytes) */

/*--------------------------------- typedefs --------------------------------*/

typedef struct structassoc
  {
  double	*list;			/* Pointer to the list of data */
  int		nobj;			/* Number of data rows */
  int		ncol;			/* Total number of columns per row */
  int		ndata;			/* Number of retained cols per row */
  int		*hash;			/* Pointer to the hash table */
  double	radius;			/* Radius of search for association */
  }             assocstruct;

/*------------------------------ Prototypes ---------------------------------*/

assocstruct	*load_assoc(char *filename, wcsstruct *wcs);

int		do_assoc(fieldstruct *field, double x, double y, double *data);

void		init_assoc(fieldstruct *field),
		end_assoc(fieldstruct *field),
		sort_assoc(fieldstruct *field, assocstruct *assoc);
