/**
* @file         dgeo.h
* @brief        Include file for dgeo.c.
* @date         02/12/2015
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       This file part of:      SExtractor
*
*       Copyright:              (C) 2011-2015 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*       License:                GNU General Public License
*
*       SExtractor is free software: you can redistribute it and/or modify
*       it under the terms of the GNU General Public License as published by
*       the Free Software Foundation, either version 3 of the License, or
*       (at your option) any later version.
*       SExtractor is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*       GNU General Public License for more details.
*       You should have received a copy of the GNU General Public License
*       along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _DGEO_H_
#define _DGEO_H_

//----------------------------- Internal constants ----------------------------
//------------------------------- functions -----------------------------------
extern int	dgeo_copy(picstruct *dgeofield, PIXTYPE *destx, PIXTYPE *desty,
		int w,int h, int ix,int iy);

#endif
