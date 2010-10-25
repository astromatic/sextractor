/*
*				tnx.h
*
* Manage polynomials.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic WCS library
*
*	Copyright:		(C) 2000-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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

#ifndef _TNX_H_
#define _TNX_H_

/*-------------------------------- macros -----------------------------------*/

#define		TNX_MAXCHARS	2048	/* Maximum FITS "WAT" string length */

/* TNX permitted types of surfaces */
#define		TNX_CHEBYSHEV	1
#define		TNX_LEGENDRE	2
#define		TNX_POLYNOMIAL	3

/* TNX cross-terms flags */
#define		TNX_XNONE	0	/* no x-terms (old no) */
#define		TNX_XFULL	1	/* full x-terms (new yes) */
#define		TNX_XHALF	2	/* half x-terms (new) */

/*----------------------------- Internal constants --------------------------*/

/*------------------------------- structures --------------------------------*/

typedef struct tnxaxis
  {
  int		type;			/* Projection correction type */
  int		xorder,yorder;		/* Polynomial orders */
  int		xterms;			/* Well... */
  int		ncoeff;			/* Number of polynom coefficients */
  double	xrange,yrange;		/* Coordinate ranges */
  double	xmaxmin,ymaxmin;	/* Well... */
  double	*coeff;			/* Polynom coefficients */
  double	*xbasis,*ybasis;	/* Basis function values */
  }	tnxaxisstruct;

/*------------------------------- functions ---------------------------------*/

tnxaxisstruct	*copy_tnxaxis(tnxaxisstruct *axis),
		*read_tnxaxis(char *tnxstr);

double		raw_to_tnxaxis(tnxaxisstruct *axis, double x, double y);

void		free_tnxaxis(tnxaxisstruct *axis);

#endif

