 /*
 				fft.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	A program that uses FFTs
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for fft.c.
*
*	Last modify:	28/05/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

/*---------------------------- Internal constants ---------------------------*/

/*------------------------------- Other Macros ------------------------------*/
#define	QFFTWMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)fftw_malloc((size_t)(nel)*sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}
#define	QFFTWFREE(ptr)	fftw_free(ptr)

/*--------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern void	fft_conv(double *data1, double *fdata2, int *size),
		fft_end(),
		fft_init(int nthreads);

extern double	*fft_rtf(double *data, int *size);
