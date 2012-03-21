/*
*				astrom.h
*
* Include file for astrom.c.
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
*	Last modified:		06/03/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

/*----------------------------- Internal constants --------------------------*/

#define	MAMA_CORFLEX	3.3e-5		/* MAMA coordinate correction factor */

/*------------------------------- structures --------------------------------*/
/*------------------------------- functions ---------------------------------*/
extern void		astrom_errparam(fieldstruct *field, obj2struct *obj2),
			astrom_init(fieldstruct *field, double pixel_scale),
			astrom_peakpos(fieldstruct *field, obj2struct *obj2),
			astrom_pos(fieldstruct **fields, int nfield,
				obj2struct *obj2),
			astrom_proferrparam(fieldstruct *field,
				obj2struct *obj2),
			astrom_profpos(fieldstruct *field, obj2struct *obj2),
			astrom_profshapeparam(fieldstruct *field,
				obj2struct *obj2),
			astrom_psferrparam(fieldstruct *field,
				obj2struct *obj2),
			astrom_psfpos(fieldstruct *field, obj2struct *obj2),
			astrom_shapeparam(fieldstruct *field, obj2struct *obj2),
			astrom_winerrparam(fieldstruct **fields, int nfield,
				obj2struct *obj2),
			astrom_winpos(fieldstruct **fields, int nfield,
				obj2struct *obj2),
			astrom_winshapeparam(fieldstruct **fields, int nfield,
				obj2struct *obj2);

