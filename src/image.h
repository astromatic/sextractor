 /*
 				image.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for image.c.
*
*	Last modify:	12/09/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
		addimage_center(picstruct *field, float *psf,
			int w,int h, float x, float y, float amplitude),
		blankimage(picstruct *, PIXTYPE *, int,int, int,int, PIXTYPE),
		pasteimage(picstruct *, PIXTYPE *, int ,int, int, int);

extern int	copyimage(picstruct *, PIXTYPE *, int, int, int, int),
		copyimage_center(picstruct *, PIXTYPE *, int,int, float,float),
		vignet_resample(float *pix1, int w1, int h1, float *pix2,
			int w2, int h2, float dx, float dy, float step2);

