 /*
 				filter.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN, IAP/Leiden
*
*	Contents:	functions dealing with on-line filtering of the image
*			(for detection).
*
*	Last modify:	01/10/2009
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*------------------------------- definitions -------------------------------*/

#define	MAXMASK		1024	/* Maximum number of mask elements (=32x32) */

/*------------------------------- structures --------------------------------*/

typedef struct structfilter
  {
/*---- convolution */
  float		*conv;		/* pointer to the convolution mask */
  int		nconv;		/* total number of elements */
  int		convw, convh;	/* x,y size of the mask */
  float		varnorm;
/*---- neural filtering */
  struct structbpann	*bpann;
  }	filterstruct;

filterstruct	*thefilter;

/*------------------------------- functions ---------------------------------*/
void		convolve(picstruct *, PIXTYPE *, int y),
		convolve_image(picstruct *field, float *vig1,
				float *vig2, int width, int height),
		filter(picstruct *, PIXTYPE *, int y),
		neurfilter(picstruct *, PIXTYPE *, int y),
		endfilter(void),
		getfilter(char *filename);

int		getconv(char *filename),
		getneurfilter(char *filename);
