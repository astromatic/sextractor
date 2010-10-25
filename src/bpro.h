/*
*				bpro.h
*
* Include file for bpro.c.
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

/*------------------------------- definitions -------------------------------*/
#define	SIGMOID(u)	((u)<15.0?((u)>-15.0?1/(1+expf(-(u))):0.0):1.0)
				/* In-line activation function */

/*---------------------------------- types ----------------------------------*/
typedef	float	NFLOAT;		/* Floating point units for neural data */

/*------------------------------- structures --------------------------------*/
typedef	struct structbpann
	{
	int	nlayers;		/* Number of "active" layers */
	int	*nn;			/* Nb of neurons per "active" layer */
/*------ The ANN itself */
	NFLOAT	**neuron;		/* Neuron array (layer,pos in layer) */
	NFLOAT	**weight;		/* Weight array (layer,pos in layer) */
	int	linearoutflag;		/* Flag: 0 if outputs are non-linear */
	}	bpannstruct;


/*------------------------------ Prototypes ---------------------------------*/

bpannstruct	*loadtab_bpann(tabstruct *tab, char *filename);

void		free_bpann(bpannstruct *bpann),
		play_bpann(bpannstruct *bpann, NFLOAT *invec, NFLOAT *outvec);

