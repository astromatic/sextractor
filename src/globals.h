/*
*				globals.h
*
* Global declarations and variables.
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
*	Last modified:		21/12/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include	"types.h"
#ifndef _FIELD_H_
#include        "field.h"
#endif

/*----------------------- miscellaneous variables ---------------------------*/

sexcatstruct		thecat;
fieldstruct		thefield1,thefield2, thewfield1,thewfield2;
extern objstruct	flagobj;
extern obj2struct	flagobj2;
char			gstr[MAXCHAR];

/*------------------------------- functions ---------------------------------*/
extern void	allocparcelout(void),
		blankit(char *, int),
                closecheck(void),
		copydata(fieldstruct *, int, int),
		flagcleancrowded(int, objliststruct *),
		freeparcelout(void),
		getnnw(void),
		makeit(void),
		mergeobject(objstruct *, objstruct *),
		neurinit(void),
		neurclose(void),
		neurresp(double *, double *),
		readdata(fieldstruct *, PIXTYPE *, int),
		readidata(fieldstruct *, FLAGTYPE *, int),
		readimagehead(fieldstruct *),
		readprefs(char *, char **, char **, int),
		useprefs(void),
		write_error(char *msg1, char *msg2);

extern int	addobj(int, objliststruct *, objliststruct *),
		belong(int, objliststruct *, int, objliststruct *),
		gatherup(objliststruct *, objliststruct *),
		parcelout(objliststruct *, objliststruct *);

extern void	*loadstrip(fieldstruct *, fieldstruct *);

extern char	*readfitshead(FILE *, char *, int *);


