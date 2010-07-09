/*
 				weight.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	Include file for weight.c.
*
*	Last modify:	01/10/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*------------------------------- definitions -------------------------------*/

#define	WTHRESH_CONVFAC		1e-4	/* Factor to apply to weights when */
					/* thresholding filtered weight-maps */

/*---------------------------------- protos --------------------------------*/

extern picstruct	*newweight(char *filename, picstruct *reffield,
				weightenum wtype, int nok);

void			weight_count(objstruct *obj, pliststruct *pixel),
			weight_to_var(picstruct *wfield, PIXTYPE *data,
				int npix);

