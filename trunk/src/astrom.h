/*
*				astrom.h
*
* Include file for astrom.c.
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

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

/*----------------------------- Internal constants --------------------------*/

#define		DEG	(PI/180.0)	/* 1 deg in radians */
#define		ARCSEC	(DEG/3600.0)	/* 1 arcsec in radians */
#define		MJD2000	51544.50000	/* Modified Julian date for J2000.0 */
#define		MJD1950	33281.92346	/* Modified Julian date for B1950.0 */
#define		JU2TROP	1.0000214	/* 1 Julian century in tropical units*/
#define		MAMA_CORFLEX	3.3e-5	/* MAMA coordinate correction factor */

/*------------------------------- structures --------------------------------*/
/*------------------------------- functions ---------------------------------*/
extern void		astrom_errparam(picstruct *, objstruct *),
			astrom_peakpos(picstruct *, objstruct *),
			astrom_pos(picstruct *, objstruct *),
			astrom_proferrparam(picstruct *, objstruct *),
			astrom_profpos(picstruct *, objstruct *),
			astrom_profshapeparam(picstruct *, objstruct *),
			astrom_psferrparam(picstruct *, objstruct *),
			astrom_psfpos(picstruct *, objstruct *),
			astrom_shapeparam(picstruct *, objstruct *),
			astrom_winerrparam(picstruct *, objstruct *),
			astrom_winpos(picstruct *, objstruct *),
			astrom_winshapeparam(picstruct *, objstruct *),
			initastrom(picstruct *),
			j2b(double, double, double, double *, double *),
			precess(double,double,double,double,double *,double *);

extern double		*compute_wcs(picstruct *, double, double);
