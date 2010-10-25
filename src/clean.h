/*
*				clean.h
*
* Include file for clean.c.
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

/*------------------------------ definitions --------------------------------*/

#define		CLEAN_ZONE		10.0	/* zone (in sigma) to */
						/* consider for processing */

/*------------------------------- variables ---------------------------------*/

objliststruct	*cleanobjlist;		/* laconic, isn't it? */

/*------------------------------- functions ---------------------------------*/

extern void	addcleanobj(objstruct *),
		endclean(void),
		initclean(void),
		subcleanobj(int);

extern int	clean(picstruct *field, picstruct *dfield,
			int, objliststruct *);

