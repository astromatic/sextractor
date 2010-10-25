/*
*				photom.h
*
* Include file for photom.c.
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

/*----------------------------- Internal constants --------------------------*/

#define	APER_OVERSAMP	5	/* oversampling in each dimension (MAG_APER) */
#define	KRON_NSIG	3*MARGIN_SCALE	/* MAG_AUTO analysis range (number */
					/* of sigma) */
#define	PETRO_NSIG	3*MARGIN_SCALE	/* MAG_PETRO analysis range (number */
					/* of sigma) */
#define	CROWD_THRESHOLD	0.1	/* The OBJ_CROWDED flag is set if photometric*/
				/* contamination may exceed this fraction of */
				/* flux */

/* NOTES:
One must have:	APER_OVERSAMP >= 1
		KRON_NSIG > 0.0
		PETRO_NSIG > 0.0
		CROWD_THRESHOLD >= 0
*/

/*------------------------------- functions ---------------------------------*/
extern void	computeaperflux(picstruct *, picstruct *, objstruct *, int),
		computeautoflux(picstruct *, picstruct *, picstruct *,
			picstruct *, objstruct *),
		computeisocorflux(picstruct *, objstruct *),
		computemags(picstruct *, objstruct *),
		computepetroflux(picstruct *, picstruct *, picstruct *,
				picstruct *, objstruct *);
