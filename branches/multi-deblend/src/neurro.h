/*
*				neurro.h
*
* Include file for neurro.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		28/03/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*--------------------------- Neural Network parameters ---------------------*/
#define		NEUR_LAYERS	3	/* max. number of hidden+i/o layers */
#define		NEUR_CONNEX	(NEUR_LAYERS-1)
#define		NEUR_NODES	10	/* maximum number of neurons/layer */

/*------------------------------- structures --------------------------------*/
typedef	struct
	{
	int	layersnb;
	int	nn[NEUR_LAYERS];
	double	inbias[NEUR_NODES];
	double	inscale[NEUR_NODES];
	double	outbias[NEUR_NODES];
	double	outscale[NEUR_NODES];
	double	ni[NEUR_NODES];
	double	no[NEUR_NODES];
	double	n[NEUR_LAYERS][NEUR_NODES];
	double	w[NEUR_CONNEX][NEUR_NODES][NEUR_NODES];
	double	b[NEUR_CONNEX][NEUR_NODES];
	}	brainstruct;

/*------------------------------- functions ---------------------------------*/

extern void	neur_close(void),
		neur_getnnw(char *filename),
		neur_init(void),
		neur_resp(double *input, double *output);

