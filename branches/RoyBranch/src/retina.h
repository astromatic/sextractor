/*
*				retina.h
*
* Include file for retina.c. 
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1995-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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

/*------------------------------- structures --------------------------------*/

typedef struct structreti
  {
/*---- convolution */
  float		*pix;		/* Pointer to the copy of the pixel array */
  int		width, height;	/* x,y size of the mask */
  int		npix;		/* Number of pixels in the retina */
  float		minnorm;	/* Minimum normalisation factor */
  struct structbpann	*bpann;	/* The neural network */
  }     retistruct;

retistruct	*theretina;

/*------------------------------- functions ---------------------------------*/

retistruct	*getretina(char *filename);
float		readretina(picstruct *, retistruct *, float, float);
void		endretina(retistruct *retina);

