/*
*				neurro.h
*
* Include file for neurro.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993,1998-2010 IAP/CNRS/UPMC
*				(C) 1994,1997 ESO
*				(C) 1995,1996 Sterrewacht Leiden
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

/*--------------------------- Neural Network parameters ---------------------*/
#define		LAYERS		3	/* max. number of hidden+i/o layers */
#define		CONNEX		LAYERS-1
#define		NEURONS		10	/* maximum number of neurons/layer */

/*------------------------------- structures --------------------------------*/
typedef	struct
	{
	int	layersnb;
	int	nn[LAYERS];
	double	inbias[NEURONS];
	double	inscale[NEURONS];
	double	outbias[NEURONS];
	double	outscale[NEURONS];
	double	ni[NEURONS];
	double	no[NEURONS];
	double	n[LAYERS][NEURONS];
	double	w[CONNEX][NEURONS][NEURONS];
	double	b[CONNEX][NEURONS];
	}	brainstruct;

/*------------------------------- globals ----------------------------------*/

extern double	f(double);
