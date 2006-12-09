/*
                                  fft.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        A program that uses FFTs
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Routines dealing with double precision FFT.
*
*       Last modify:    08/12/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <fftw3.h>

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

#define SWAP(a,b)       tempr=(a);(a)=(b);(b)=tempr

/****** fft_init ************************************************************
PROTO	void fft_init(void)
PURPOSE	Initialize the FFT routines
INPUT	-.
OUTPUT	-.
NOTES	Global preferences are used for multhreading.
AUTHOR	E. Bertin (IAP)
VERSION	29/11/2006
 ***/
void    fft_init(void)
 {
  if (!firsttimeflag)
    {
#ifdef USE_THREADS
    if (!fftw_init_threads())
      error(EXIT_FAILURE, "*Error*: thread initialization failed in ", "FFTW");
    fftw_plan_with_nthreads(prefs.nthreads);
    QPTHREAD_MUTEX_INIT(&fftmutex, NULL);
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
VERSION	29/11/2006
 ***/
void    fft_end(void)
 {

  if (firsttimeflag)
    {
    firsttimeflag = 0;
#ifdef USE_THREADS
    fftw_cleanup_threads();
    QPTHREAD_MUTEX_DESTROY(&fftmutex);
#endif
    fftw_cleanup();
    }

  return;
  }


/****** fft_conv ************************************************************
PROTO	void fft_conv(double *data1, double *fdata2, int *size)
PURPOSE	Optimized 2-dimensional FFT convolution using the FFTW library.
INPUT	ptr to the first image,
	ptr to the Fourier transform of the second image,
	image size vector.
OUTPUT	-.
NOTES	For data1 and fdata2, memory must be allocated for
	size[0]* ... * 2*(size[naxis-1]/2+1) doubles (padding required).
AUTHOR	E. Bertin (IAP)
VERSION	07/12/2006
 ***/
void    fft_conv(double *data1, double *fdata2, int *size)
  {
   fftw_plan	plan;
   double	*fdata1,*fdata1p,*fdata2p,
		real,imag, fac;
   int		i, npix,npix2;

/* Convert axis indexing to that of FFTW */
  npix = size[0]*size[1];
  npix2 = (((size[0]>>1) + 1)<< 1) * size[1];

/* Forward FFT "in place" for data1 */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  QFFTWMALLOC(fdata1, double, npix2);
  plan = fftw_plan_dft_r2c_2d(size[0], size[1], data1,
        (fftw_complex *)fdata1, FFTW_ESTIMATE|FFTW_FORWARD|FFTW_DESTROY_INPUT);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif
  fftw_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftw_destroy_plan(plan);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

/* Actual convolution (Fourier product) */
  fac = 1.0/npix;  
  fdata1p = fdata1;
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
  plan = fftw_plan_dft_c2r_2d(size[0], size[1], (fftw_complex *)fdata1, 
        data1, FFTW_ESTIMATE|FFTW_BACKWARD|FFTW_DESTROY_INPUT);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  fftw_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftw_destroy_plan(plan);
/* Free the fdata1 scratch array */
  QFFTWFREE(fdata1);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  return;
  }


/****** fft_rtf ************************************************************
PROTO	double *fft_rtf(double *data, int *size)
PURPOSE	Optimized 2-dimensional FFT "in place" using the FFTW library.
INPUT	ptr to the image,
	ptr to image size vector.
OUTPUT	Pointer to the compressed, memory-allocated Fourier transform.
NOTES	Input data may end up corrupted.
AUTHOR	E. Bertin (IAP)
VERSION	29/11/2006
 ***/
double	*fft_rtf(double *data, int *size)
  {
   fftw_plan   plan;
   double	*fdata;
   int		npix2;

/* Convert axis indexing to that of FFTW */
  npix2 = (((size[0]>>1) + 1)<< 1) * size[1];

/* Forward FFT "in place" for data1 */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  QFFTWMALLOC(fdata, double, npix2);
  plan = fftw_plan_dft_r2c_2d(size[0], size[1], data,
        (fftw_complex *)fdata, FFTW_ESTIMATE|FFTW_FORWARD|FFTW_DESTROY_INPUT);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  fftw_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftw_destroy_plan(plan);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  return fdata;
  }


