/*
*				analyse.h
*
* Include file for analyse.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		07/10/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*----------------------------- Internal constants --------------------------*/

#define	ANALYSE_NMULTITER	5	/* number of multi-model iterations */

/*--------------------------------- typedefs --------------------------------*/

/*------------------------------ Prototypes ---------------------------------*/

obj2struct	*analyse_obj2obj2(picstruct *field, picstruct *dfield,
			picstruct *wfield, picstruct *dwfield,
			objstruct *obj, obj2liststruct *obj2list);
int		analyse_overlapness(objliststruct *objlist, int iobj);

void		analyse_full(picstruct *field, picstruct *dfield,
			picstruct *wfield, picstruct *dwfield,
			obj2struct *obj2),
		analyse_group(picstruct *field, picstruct *dfield,
			picstruct *wfield, picstruct *dwfield,
			obj2struct *obj2),
		analyse_iso(picstruct *field, picstruct *dfield,
			objliststruct *objlist, int n);
		

