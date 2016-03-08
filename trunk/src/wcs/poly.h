/*
*				poly.h
*
* Include file for poly.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1998-2011 IAP/CNRS/UPMC
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
*	Last modified:		20/11/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _POLY_H_
#define _POLY_H_

/*--------------------------------- constants -------------------------------*/

#define	POLY_MAXDIM		4	/* Max dimensionality of polynom */
#define POLY_MAXDEGREE		40	/* Max degree of the polynom */
#define	POLY_TINY		1e-30	/* A tiny number */

/*---------------------------------- macros ---------------------------------*/

/*--------------------------- structure definitions -------------------------*/

typedef struct poly
  {
  double	*basis;		/* Current values of the basis functions */
  double	*orthobasis;	/* Curr orthonormalized basis function values */
  double	*coeff;		/* Polynom coefficients */
  int		ncoeff;		/* Number of coefficients */
  int		*group;		/* Groups */
  int		ndim;		/* dimensionality of the polynom */
  int		*degree;	/* Degree in each group */
  int		ngroup;		/* Number of different groups */
  double	*orthomat;	/* Orthonormalization matrix */
  double	*deorthomat;	/* "Deorthonormalization" matrix */
  }	polystruct;

/*---------------------------------- protos --------------------------------*/

extern polystruct	*poly_copy(polystruct *poly),
			*poly_init(int *group,int ndim,int *degree,int ngroup);

extern double		poly_func(polystruct *poly, double *pos);

extern int		cholsolve(double *a, double *b, int n),
			poly_fit(polystruct *poly, double *x, double *y,
				double *w, int ndata, double *extbasis,
				double regul),
			*poly_powers(polystruct *poly),
			poly_solve(double *a, double *b, int n);

extern void		poly_addcste(polystruct *poly, double *cste),
			poly_end(polystruct *poly);

extern double		*poly_deortho(polystruct *poly, double *datain,
				double *dataout),
			*poly_ortho(polystruct *poly, double *datain,
				double *dataout);
extern void		poly_initortho(polystruct *poly, double *data,
				int ndata);

#endif

