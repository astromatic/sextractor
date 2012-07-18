/*
*				fft.h
*
* Include file for fft.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2007-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		18/07/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef FFTW3_H
#include FFTW_H
#endif

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

/*---------------------------- Internal constants ---------------------------*/

/*------------------------------- Other Macros ------------------------------*/
#define	QFFTWF_MALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)fftwf_malloc((size_t)(nel)*sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}
#define	QFFTWF_FREE(ptr) \
		{fftwf_free(ptr); ptr=NULL;}

/*--------------------------- structure definitions -------------------------*/
typedef	struct fftscratch
	{
        fftwf_plan   	fplan, bplan;	/* Forward and backward FFT plans */
	fftwf_complex	*fdata;		/* Scratch (Fourier) space */
	int		size[2];	/* Scratch space dimensions */
	}	fftscratchstruct;

/*---------------------------------- protos --------------------------------*/
extern void	fft_conv(float *data1, float *fdata2, int *size,
			fftscratchstruct **fftscratch),
		fft_end(),
		fft_init(int nthreads),
		fft_reset(void),
		fft_scratchend(fftscratchstruct *fftscratch);

extern float	*fft_rtf(float *data, int *size);
