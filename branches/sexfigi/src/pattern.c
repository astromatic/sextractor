 /*
 				pattern.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Authors:	E.BERTIN (IAP)
*
*	Contents:	Generate and handle image patterns for image fitting.
*
*	Last modify:	09/09/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#define _GNU_SOURCE
#include <math.h>
#endif

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"pattern.h"

/*------------------------------- variables ---------------------------------*/


/****** pattern_init ***********************************************************
PROTO	patternstruct pattern_init(pattern_type ptype, int nvec)
PURPOSE	Allocate and initialize a new pattern structure.
INPUT	Pattern type,
	Number of vectors.
OUTPUT	Pointer to the new pattern structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/09/2008
 ***/
patternstruct	*pattern_init(pattypenum ptype, int nvec)
  {
   patternstruct	*pattern;
   double		stepx,stepy, dx,dy,dy2, dr2,dr2min, lr2, invlr0, mod,
			ang, cosang,sinang;
   float		*cpix,*spix;
   int			p, nx,ny, width,height, npix;

  QCALLOC(pattern, patternstruct, 1);
  pattern->type = ptype;
  width = pattern->size[0] = height = pattern->size[1] = PATTERN_SIZE;
  npix = width*height;
  stepx = 2.0*PATTERN_RADIUS/(width-1);
  stepy = 2.0*PATTERN_RADIUS/(height-1);
  invlr0 = 1.0/log(PATTERN_SCALE);
  if ((dr2min = stepx*stepx) < stepy*stepy)
    dr2min = stepy*stepy;
  pattern->size[2] = nvec*2;
  QCALLOC(pattern->pix, float, pattern->size[0]*pattern->size[1]
	*pattern->size[2]);
  cpix = pattern->pix;
  spix = pattern->pix + npix;
  for (p=0; p<nvec; p++, cpix+=npix, spix+=npix)
    {
    dy = -PATTERN_RADIUS;
    for (ny=height; ny--; dy+=stepy)
      {
      dy2 =dy*dy;
      dx = -PATTERN_RADIUS;
      for (nx=width; nx--; dx+=stepx)
        {
        dr2 = dy2+dx*dx;
        lr2 = 0.5*log(dr2 > dr2min ? dr2 : dr2min)*invlr0 - 1.0;
        mod = exp(-0.5*lr2*lr2);
        ang = 2.0*atan2(dy,dx);
#ifdef HAVE_SINCOS
        sincos(ang, &sinang, &cosang);
#else
        sinang = sin(ang);
        cosang = cos(ang);
#endif
        *(cpix++) = mod*cosang;
        *(spix++) = mod*sinang;
        }
      }
    }

  return pattern;
  }  


/****** pattern_end ***********************************************************
PROTO	void pattern_end(patternstruct *pattern)
PURPOSE	End (deallocate) a pattern structure.
INPUT	Pattern structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	09/09/2008
 ***/
void	pattern_end(patternstruct *pattern)
  {
  free(pattern->pix);
  free(pattern);

  return;
  }

