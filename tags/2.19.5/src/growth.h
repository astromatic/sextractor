/*
*				growth.h
*
* Include file for growth.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1995-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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

#define	GROWTH_NSTEP	64	/* number of growth curve samples */
#define	GROWTH_OVERSAMP	5	/* pixel oversampling in each dimension */
#define	GROWTH_NSIG	3*MARGIN_SCALE	/* MAG_AUTO analysis range (nsigmas) */
#define	GROWTH_MINHLRAD	0.5	/* Minimum internal half-light radius (pixels)*/

/* NOTES:
One must have:	GROWTH_SAMP >= 1
		GROWTH_OVERSAMP >= 1
		GROWTH_OVERSAMPRADIUS >= 0
		GROWTH_NSIG > 0.0
*/

/*------------------------------- functions ---------------------------------*/
extern void	endgrowth(void),
		initgrowth(void),
		makeavergrowth(picstruct *field, picstruct *wfield,
			objstruct *obj);

