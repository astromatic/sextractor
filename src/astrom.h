/*
*				astrom.h
*
* Include file for astrom.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		06/10/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

/*----------------------------- Internal constants --------------------------*/

#define	MAMA_CORFLEX	3.3e-5		/* MAMA coordinate correction factor */

/*------------------------------- structures --------------------------------*/
/*------------------------------- functions ---------------------------------*/
extern void		astrom_errparam(picstruct *field, obj2struct *obj2),
			astrom_peakpos(picstruct *field, obj2struct *obj2),
			astrom_pos(picstruct *field, obj2struct *obj2),
			astrom_proferrparam(picstruct *field, obj2struct *obj2),
			astrom_profpos(picstruct *field, obj2struct *obj2),
			astrom_profshapeparam(picstruct *field,
					obj2struct *obj2),
			astrom_psferrparam(picstruct *field, obj2struct *obj2),
			astrom_psfpos(picstruct *field, obj2struct *obj2),
			astrom_shapeparam(picstruct *field, obj2struct *obj2),
			astrom_winerrparam(picstruct *field, obj2struct *obj2),
			astrom_winpos(picstruct *field, obj2struct *obj2),
			astrom_winshapeparam(picstruct *field, obj2struct *obj2),
			initastrom(picstruct *field);

