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
*	Last modify:	14/07/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
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
  char		ident[2][MAXCHAR];			/* identifiants */
  float		backmean[2];				/* mean background */
  float		backsig[2];				/* mean back stddev */
  float		sigfac[2];				/* mean weight scaling*/
  float		thresh[2];				/* thresholds (ADU) */
  double	pixscale[2];				/* pixel scale (deg2) */
  double	epoch[2];				/* epoch of coords */
  }	xmlstruct;

/*------------------------------- functions ---------------------------------*/

extern int		end_xml(void),
			init_xml(int next),
			update_xml(sexcatstruct *sexcat, picstruct *dfield, 
				picstruct *field, picstruct *dwfield,
				picstruct *wfield),
			write_xml(char *filename),
			write_xml_header(FILE *file),
			write_xml_meta(FILE *file, char *error);

extern void		write_xmlerror(char *filename, char *error);
