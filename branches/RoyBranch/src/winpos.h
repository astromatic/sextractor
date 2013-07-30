/*
*				winpos.h
*
* Include file for winpos.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2005-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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

/*----------------------------- Internal constants --------------------------*/

#define	WINPOS_NITERMAX	16	/* Maximum number of steps */
#define	WINPOS_NSIG	4	/* Measurement radius */
#define	WINPOS_OVERSAMP	11	/* oversampling in each dimension */
#define	WINPOS_STEPMIN	0.0001	/* Minimum change in position for continueing*/
#define	WINPOS_FAC	2.0	/* Centroid offset factor (2 for a Gaussian) */

/* NOTES:
One must have:
	WINPOS_NITERMAX >= 1
	WINPOS_OVERSAMP >= 1
*/

/*------------------------------- functions ---------------------------------*/
extern void	compute_winpos(picstruct *field, picstruct *wfield,
			       objstruct *obj);
