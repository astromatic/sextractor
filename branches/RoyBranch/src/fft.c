/*
*				fft.c
*
* Single precision FFT functions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2007-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		17/04/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_THREADS
#include <pthread.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef FFTW3_H
#include FFTW_H
#endif

#include "define.h"
#include "globals.h"
#include "fft.h"
#include "prefs.h"
#ifdef USE_THREADS
#include "threads.h"
#endif

 fftwf_plan	fplan, bplan;
 int    	firsttimeflag;

#ifdef USE_THREADS
 pthread_mutex_t	fftmutex;
#endif
 fftwf_complex 		*fdata1;

/****** fft_init ************************************************************
PROTO	void fft_init(void)
PURPOSE	Initialize the FFT routines
INPUT	-.
OUTPUT	-.
NOTES	Global preferences are used for multhreading.
AUTHOR	E. Bertin (IAP)
VERSION	26/06/2009
 ***/
void    fft_init(int nthreads)
 {
  if (!firsttimeflag)
    {
#ifdef USE_THREADS
    if (nthreads > 1)
      {
      if (!fftw_init_threads())
        error(EXIT_FAILURE, "*Error*: thread initialization failed in ", "FFTW");
      fftwf_plan_with_nthreads(prefs.nthreads);
      }
#endif
    firsttimeflag = 1;
    }

  return;
  }


/****** fft_end ************************************************************
PROTO	void fft_init(void)
PURPOSE	Clear up stuff set by FFT routines
INPUT	-.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/06/2009
 ***/
void    fft_end(void)
 {

  if (firsttimeflag)
    {
    firsttimeflag = 0;
    fftwf_cleanup();
    }

  return;
  }


/****** fft_reset ***********************************************************
PROTO	void fft_reset(void)
PURPOSE	Reset the FFT plans
INPUT	-.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/10/2009
 ***/
void    fft_reset(void)
 {
  if (fplan)
    {
    QFFTWF_FREE(fdata1);
    fftwf_destroy_plan(fplan);
    }
  if (bplan)
    fftwf_destroy_plan(bplan);
  fplan = bplan = NULL;

  return;
  }


/****** fft_conv ************************************************************
PROTO	void fft_conv(float *data1, float *fdata2, int *size)
PURPOSE	Optimized 2-dimensional FFT convolution using the FFTW library.
INPUT	ptr to the first image,
	ptr to the Fourier transform of the second image,
	image size vector.
OUTPUT	-.
NOTES	For data1 and fdata2, memory must be allocated for
	size[0]* ... * 2*(size[naxis-1]/2+1) floats (padding required).
AUTHOR	E. Bertin (IAP)
VERSION	29/03/2013
 ***/
void    fft_conv(float *data1, float *fdata2, int *size)
  {
   float		*fdata1p,*fdata2p,
			real,imag, fac;
   int			i, npix,npix2;

/* Convert axis indexing to that of FFTW */
  npix = size[0]*size[1];
  npix2 = ((size[0]/2) + 1) * size[1];

/* Forward FFT "in place" for data1 */
  if (!fplan)
    {
    QFFTWF_MALLOC(fdata1, fftwf_complex, npix2);
    fplan = fftwf_plan_dft_r2c_2d(size[1], size[0], data1,
        (fftwf_complex *)fdata1, FFTW_ESTIMATE);
    }

  fftwf_execute_dft_r2c(fplan, data1, fdata1);

/* Actual convolution (Fourier product) */
  fac = 1.0/npix;  
  fdata1p = (float *)fdata1;
  fdata2p = fdata2;
#pragma ivdep
  for (i=npix2; i--;)
    {
    real = *fdata1p **fdata2p - *(fdata1p+1)**(fdata2p+1);
    imag = *(fdata1p+1)**fdata2p + *fdata1p**(fdata2p+1);
    *(fdata1p) = fac*real;
    *(fdata1p+1) = fac*imag;
    fdata1p+=2;
    fdata2p+=2;
    }

/* Reverse FFT */
  if (!bplan)
    bplan = fftwf_plan_dft_c2r_2d(size[1], size[0], (fftwf_complex *)fdata1, 
        data1, FFTW_ESTIMATE);
  fftwf_execute_dft_c2r(bplan, fdata1, data1);

//  fftwf_execute(plan);


  return;
  }


/****** fft_rtf ************************************************************
PROTO	float *fft_rtf(float *data, int *size)
PURPOSE	Optimized 2-dimensional FFT "in place" using the FFTW library.
INPUT	ptr to the image,
	ptr to image size vector.
OUTPUT	Pointer to the compressed, memory-allocated Fourier transform.
NOTES	Input data may end up corrupted.
AUTHOR	E. Bertin (IAP)
VERSION	12/07/2012
 ***/
float	*fft_rtf(float *data, int *size)
  {
   fftwf_plan   	plan;
   fftwf_complex	*fdata;
   int			npix2;

/* Convert axis indexing to that of FFTW */
  npix2 = ((size[0]/2) + 1) * size[1];

/* Forward FFT "in place" for data1 */
  QFFTWF_MALLOC(fdata, fftwf_complex, npix2);
  plan = fftwf_plan_dft_r2c_2d(size[1], size[0], data, fdata, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  return (float *)fdata;
  }


