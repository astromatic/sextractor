/*
*				catout.h
*
* Include file for catout.c.
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
*	Last modified:		18/07/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*--------------------------------- typedefs --------------------------------*/

/*------------------------------ Prototypes ---------------------------------*/

obj2liststruct	*catout_readparams(char **paramlist, int nparam,
						int nobj2);

void		catout_allocparams(obj2liststruct *obj2list),
		catout_changeparamsize(char *keyword, int *axisn, int naxis),		
		catout_dumpparams(void),
		catout_end(char *error),
		catout_endext(void),
		catout_freeparams(obj2liststruct *obj2list),
		catout_init(void),
		catout_initext(picstruct *field),
		catout_updateparamflags(void),
		catout_writeobj(obj2struct *obj2),
		catout_writevofields(FILE *file);

