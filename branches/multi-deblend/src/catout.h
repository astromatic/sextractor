/*
*				catout.h
*
* Include file for catout.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2011-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		09/06/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*--------------------------------- typedefs --------------------------------*/

/*------------------------------ Prototypes ---------------------------------*/

int		catout_allocobjother(objstruct *obj, void *flagobj2elem,
				int nbytes),
		catout_freeobjother(objstruct *obj, void *flagobj2elem);

void		catout_allocobjparams(objstruct *obj),
		catout_changeparamsize(char *keyword, int *axisn, int naxis),		
		catout_dumpparams(void),
		catout_end(char *error),
		catout_endext(void),
		catout_freeobjparams(objstruct *obj),
		catout_init(void),
		catout_initext(fieldstruct *field),
		catout_readparams(char **paramlist, int nparam),
		catout_updateparamflags(void),
		catout_writeobj(objstruct *obj),
		catout_writevofields(FILE *file);

