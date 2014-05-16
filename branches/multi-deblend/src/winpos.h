/*
*				winpos.h
*
* Include file for winpos.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2005-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		02/04/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*--------------------------- WINPOS flags ---------------------------------*/

#define	WINFLAG_SINGULAR	0x0001	/* Singularity in WINdowed light dist.*/
#define	WINFLAG_NEGMOMENT	0x0002	/* Negative WINdowed light 2nd moments*/
#define	WINFLAG_NEGFLUX		0x0004	/* Negative WINdowed integrated flux */
#define	WINFLAG_APERT_PB	0x0008	/* Window area incomplete or truncated*/
#define	WINFLAG_NOTCOVERED	0x0010	/* Field does not include object */

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
extern void	win_pos(fieldstruct **fields, fieldstruct **wfields,
				int nfield, obj2struct *obj2);
