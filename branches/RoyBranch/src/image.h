/*
*				image.h
*
* Include file for image.c.
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
*	Last modified:		19/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*----------------------------- Internal constants --------------------------*/

#define INTERPW		8	/* Interpolation function range */
#define	INTERPFAC	4.0	/* Interpolation envelope factor */

#define	INTERPF(x)	(x<1e-5 && x>-1e-5? 1.0 \
			:(x>INTERPFAC?0.0:(x<-INTERPFAC?0.0 \
			:sinf(PI*x)*sinf(PI/INTERPFAC*x)/(PI*PI/INTERPFAC*x*x))))
				/* Lanczos approximation */

/*--------------------------- structure definitions -------------------------*/


/*----------------------------- Global variables ----------------------------*/

/*------------------------------- functions ---------------------------------*/
extern void    	addimage(picstruct *field, float *psf,
			int w,int h, int ix,int iy, float amplitude),
		addfrombig(float *pixbig, int wbig,int hbig,
			float *pixsmall, int wsmall, int hsmall,			
			int ix,int iy, float amplitude),
		addtobig(float *pixsmall, int wsmall,int hsmall,
			float *pixbig, int wbig, int hbig,
			int ix,int iy, float amplitude),
		addimage_center(picstruct *field, float *psf,
			int w,int h, float x, float y, float amplitude),
		blankimage(picstruct *, PIXTYPE *, int,int, int,int, PIXTYPE),
		pasteimage(picstruct *, PIXTYPE *, int ,int, int, int);

extern int	copyimage(picstruct *, PIXTYPE *, int, int, int, int),
		copyimage_center(picstruct *, PIXTYPE *, int,int, float,float),
		vignet_resample(float *pix1, int w1, int h1, float *pix2,
			int w2, int h2, float dx, float dy, float step2);

