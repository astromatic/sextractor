/*
*				poly.h
*
* Include file for poly.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic WCS library
*
*	Copyright:		(C) 1998-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _POLY_H_
#define _POLY_H_

/*--------------------------------- constants -------------------------------*/

#define	POLY_MAXDIM		4	/* Max dimensionality of polynom */
#define POLY_MAXDEGREE		10	/* Max degree of the polynom */

/*---------------------------------- macros ---------------------------------*/

/*--------------------------- structure definitions -------------------------*/

typedef struct poly
  {
  double	*basis;		/* Current values of the basis functions */
  double	*coeff;		/* Polynom coefficients */
  int		ncoeff;		/* Number of coefficients */
  int		*group;		/* Groups */
  int		ndim;		/* dimensionality of the polynom */
  int		*degree;	/* Degree in each group */
  int		ngroup;		/* Number of different groups */
  }	polystruct;

/*---------------------------------- protos --------------------------------*/

extern polystruct	*poly_init(int *group,int ndim,int *degree,int ngroup);

extern double			poly_func(polystruct *poly, double *pos);

extern int		cholsolve(double *a, double *b, int n),
			*poly_powers(polystruct *poly);

extern void		poly_addcste(polystruct *poly, double *cste),
			poly_end(polystruct *poly),
			poly_fit(polystruct *poly, double *x, double *y,
				double *w, int ndata, double *extbasis),
			poly_solve(double *a, double *b, int n);
#endif

