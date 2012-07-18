/*
*				fft.c
*
* Single precision FFT functions.
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

 int    firsttimeflag;

#ifdef USE_THREADS
 pthread_mutex_t	fftmutex;
#endif

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
    QPTHREAD_MUTEX_INIT(&fftmutex, NULL);
#ifdef HAVE_FFTW_MP
    if (nthreads > 1)
      {
      if (!fftw_init_threads())
        error(EXIT_FAILURE, "*Error*: thread initialization failed in ", "FFTW");
      fftwf_plan_with_nthreads(prefs.nthreads);
      }
#endif
#endif
    firsttimeflag = 1;
    }

  return;
  }


/****** fft_end *************************************************************
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
#ifdef USE_THREADS
#ifdef HAVE_FFTW_MP
    fftwf_cleanup_threads();
#endif
    QPTHREAD_MUTEX_DESTROY(&fftmutex);
#endif
    fftwf_cleanup();
    }

  return;
  }


/****** fft_scratchend ******************************************************
PROTO	void fft_scratchend(fftscratchstruct *fftscratch)
PURPOSE	Reset FFT buffers (including plans and scratch space).
INPUT	Pointer to the fftscratch structure to reset.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/07/2012
 ***/
void    fft_scratchend(fftscratchstruct *fftscratch)
  {
  if (!fftscratch)
    return;

  if (fftscratch->fplan)
    fftwf_destroy_plan(fftscratch->fplan);
   if (fftscratch->bplan)
    fftwf_destroy_plan(fftscratch->bplan);
  if (fftscratch->fdata)
    QFFTWF_FREE(fftscratch->fdata);
  free(fftscratch);

  return;
  }


/****** fft_conv ************************************************************
PROTO	fftscratchstruct	*fft_conv(float *data1,float *fdata2, int *size,
					fftscratchstruct	**pfftscratch)
PURPOSE	Optimized 2-dimensional FFT convolution using the FFTW library.
INPUT	ptr to the first image,
	ptr to the Fourier transform of the second image,
	image size vector,
	ptr to an FFT scratch structure pointer. If the structure pointer points
	to NULL, initialization is performed and the returned pointer points to
	the new FFT scratch structure.
OUTPUT	ptr to a new FFT scratch structure if initialization is required.
NOTES	For data1 and fdata2, memory must be allocated for
	size[0]* ... * 2*(size[naxis-1]/2+1) floats (padding required).
AUTHOR	E. Bertin (IAP)
VERSION	18/07/2012
 ***/
void	fft_conv(float *data1, float *fdata2, int *size,
			fftscratchstruct **pfftscratch)
  {
   fftscratchstruct	*fftscratch;
   float		*fdata1p,*fdata2p,
			real,imag, fac;
   int			i, npix,npix2;

/* Convert axis indexing to that of FFTW */
  npix = size[0]*size[1];
  npix2 = ((size[0]/2) + 1) * size[1];

/* Forward FFT "in place" for data1 */
  if (!pfftscratch || !(fftscratch=*pfftscratch)
	|| fftscratch->size[0] != size[0] || fftscratch->size[1] != size[1])
    {
    QCALLOC(fftscratch, fftscratchstruct, 1);
    if (pfftscratch)
      {
      if (*pfftscratch)
        fft_scratchend(*pfftscratch);
      *pfftscratch = fftscratch;
      }
#ifdef USE_THREADS
    QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
    QFFTWF_MALLOC(fftscratch->fdata, fftwf_complex, npix2);
    fftscratch->fplan = fftwf_plan_dft_r2c_2d(size[1], size[0], data1,
        fftscratch->fdata, FFTW_ESTIMATE);
#ifdef USE_THREADS
    QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif
    }

  fftwf_execute_dft_r2c(fftscratch->fplan, data1, fftscratch->fdata);

/* Actual convolution (Fourier product) */
  fac = 1.0/npix;  
  fdata1p = (float *)fftscratch->fdata;
  fdata2p = fdata2;
  for (i=npix2; i--; fdata2p+=2)
    {
    real = *fdata1p **fdata2p - *(fdata1p+1)**(fdata2p+1);
    imag = *(fdata1p+1)**fdata2p + *fdata1p**(fdata2p+1);
    *(fdata1p++) = fac*real;
    *(fdata1p++) = fac*imag;
    }

/* Reverse FFT */
  if (!fftscratch->bplan)
    {
#ifdef USE_THREADS
    QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
    fftscratch->bplan = fftwf_plan_dft_c2r_2d(size[1], size[0],
		fftscratch->fdata, data1, FFTW_ESTIMATE);
#ifdef USE_THREADS
    QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif
    }

  fftwf_execute_dft_c2r(fftscratch->bplan, fftscratch->fdata, data1);

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
VERSION	18/07/2012
 ***/
float	*fft_rtf(float *data, int *size)
  {
   fftwf_plan   	plan;
   fftwf_complex	*fdata;
   int			npix2;

/* Convert axis indexing to that of FFTW */
  npix2 = ((size[0]/2) + 1) * size[1];

/* Forward FFT "in place" for data1 */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  QFFTWF_MALLOC(fdata, fftwf_complex, npix2);
  plan = fftwf_plan_dft_r2c_2d(size[1], size[0], data, fdata, FFTW_ESTIMATE);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif
  fftwf_execute(plan);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftwf_destroy_plan(plan);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  return (float *)fdata;
  }


