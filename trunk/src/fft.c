/*
                                  fft.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        A program that uses FFTs
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Routines dealing with float precision FFT.
*
*       Last modify:    17/11/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include FFTW_H

#include "define.h"
#include "globals.h"
#include "fft.h"
#include "prefs.h"
#ifdef USE_THREADS
#include "threads.h"
#endif

 fftwf_plan	fplan, bplan;
 int    firsttimeflag;
#ifdef USE_THREADS
 pthread_mutex_t	fftmutex;
#endif
 fftwf_complex 	*fdata1;
#define SWAP(a,b)       tempr=(a);(a)=(b);(b)=tempr

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
    QFFTWFREE(fdata1);
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
VERSION	08/10/2009
 ***/
void    fft_conv(float *data1, float *fdata2, int *size)
  {
   fftwf_complex	*fdata0;
   float		*fdata1p,*fdata2p,*data0,
			real,imag, fac;
   int			i, npix,npix2;

/* Convert axis indexing to that of FFTW */
  npix = size[0]*size[1];
  npix2 = (((size[0]/2) + 1)*2) * size[1];

/* Forward FFT "in place" for data1 */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  if (!fplan)
    {
    QFFTWMALLOC(fdata1, fftwf_complex, npix2);
    fplan = fftwf_plan_dft_r2c_2d(size[1], size[0], data1,
        (fftwf_complex *)fdata1, FFTW_ESTIMATE);
    }
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif
  fftwf_execute_dft_r2c(fplan, data1, fdata1);

//  fftwf_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
//  fftwf_destroy_plan(plan);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

/* Actual convolution (Fourier product) */
  fac = 1.0/npix;  
  fdata1p = (float *)fdata1;
  fdata2p = fdata2;
  for (i=npix2/2; i--; fdata2p+=2)
    {
    real = *fdata1p **fdata2p - *(fdata1p+1)**(fdata2p+1);
    imag = *(fdata1p+1)**fdata2p + *fdata1p**(fdata2p+1);
    *(fdata1p++) = fac*real;
    *(fdata1p++) = fac*imag;
    }

/* Reverse FFT */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  if (!bplan)
    bplan = fftwf_plan_dft_c2r_2d(size[1], size[0], (fftwf_complex *)fdata1, 
        data1, FFTW_ESTIMATE);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif
  fftwf_execute_dft_c2r(bplan, fdata1, data1);

//  fftwf_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
//  fftwf_destroy_plan(plan);
/* Free the fdata1 scratch array */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

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
VERSION	08/10/2009
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
  QFFTWMALLOC(fdata, fftwf_complex, npix2);
  plan = fftwf_plan_dft_r2c_2d(size[1], size[0], data,
        fdata, FFTW_ESTIMATE);
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


