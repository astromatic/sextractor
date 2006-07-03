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
*	Last modify:	03/07/2005
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
  sexcatstruct	sexcat;
  }	xmlstruct;

/*------------------------------- functions ---------------------------------*/

extern int		init_xml(int next),
			update_xml(void),
			write_xml(void);
