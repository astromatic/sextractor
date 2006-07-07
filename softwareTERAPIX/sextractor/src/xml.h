 /*
 				xml.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	XML logging.
*
*	Last modify:	07/07/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

/*----------------------------- Internal constants --------------------------*/
/*
#define	XSL_URL	"file:///home/bertin/sources/sex/xsl/sex.xsl"
*/
#ifndef XSL_URL
#define	XSL_URL	"."
#endif
/* Alternate XSLT file at TERAPIX: */
/* will not work with recent browsers because of security limitations */
/*
#define	XSL_URL_ALT	"http://terapix.iap.fr/cplt/xsl/sex.xsl"
*/
/*--------------------------------- typedefs --------------------------------*/
typedef struct
  {
  int		currext;
  int		ndetect;
  int		ntotal;
  char 		ext_date[16],ext_time[16];		/* date and time */
  double        ext_elapsed;				/* processing time */
  char		dident[MAXCHAR],ident[MAXCHAR];		/* identifiants */
  float		dbackmean, backmean;			/* mean background */
  float		dbacksig, backsig;			/* mean back stddev */
  float		dsigfac, sigfac;			/* mean weight scaling*/
  float		dthresh, thresh;			/* thresholds (ADU) */
  double	dpixscale, pixscale;			/* pixel scale (deg2) */
  double	depoch, epoch;				/* epoch of coords */
  }	xmlstruct;

/*------------------------------- functions ---------------------------------*/

extern int		init_xml(int next),
			update_xml(sexcatstruct *sexcat, picstruct *dfield, 
				picstruct *field, picstruct *dwfield,
				picstruct *wfield),
			write_xml(void);

extern void		write_xmlerror(char *msg1, char *msg2);
