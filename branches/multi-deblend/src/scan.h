/*
*				scan.h
*
* Include file for scan.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2014 IAP/CNRS/UPMC
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
*	Last modified:		26/06/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _LUTZ_H_
#include "lutz.h"
#endif

#ifndef _PLIST_H_
#include "plist.h"
#endif

/*------------------------------ definitions --------------------------------*/

#define	NOBJ			256		/* starting number of obj. */
#define	UNKNOWN			-1		/* flag for LUTZ */

/*---------------------------- preanalyse flags -----------------------------*/

#define		ANALYSE_FAST		0
#define		ANALYSE_FULL		1
#define		ANALYSE_ROBUST		2

/*--------------------------------- typedefs --------------------------------*/

/*------------------------------- functions ---------------------------------*/
void		scan_extract(fieldstruct *dfield, fieldstruct *dwfield,
			fieldstruct **fields, fieldstruct **wfields, int nfield,
			fieldstruct **ffields, int nffield),
		scan_output(fieldstruct **fields, fieldstruct **wfields,
			int nfield, infostruct *info, pliststruct *plist,
			objliststruct *cleanobjlist),
		scan_preanalyse(objstruct *obj, pliststruct *plist,
			int analyse_type);

