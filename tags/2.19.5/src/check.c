/*
*				check.c
*
* Manage "check-images".
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SExtractor is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SExtractor is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		23/09/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"fits/fitscat.h"
#include	"fitswcs.h"
#include	"check.h"

/********************************* addcheck **********************************/
/*
Add a PSF to a CHECK-image (with a multiplicative factor).
Outside boundaries are taken into account.
*/
void	addcheck(checkstruct *check, float *psf,
			int w,int h, int ix,int iy, float amplitude)
  {
   PIXTYPE	*pix;
   int		x,y, xmin,xmax,ymin,ymax,w2, dwpsf;

/* Set the image boundaries */
  w2 = w;
  ymin = iy-h/2;
  ymax = ymin + h;
  if (ymin<0)
    {
    psf -= ymin*w;
    ymin = 0;
    }
  if (ymax>check->height)
    ymax = check->height;

  xmin = ix-w/2;
  xmax = xmin + w;
  if (xmax>check->width)
    {
    w2 -= xmax-check->width;
    xmax = check->width;
    }
  if (xmin<0)
    {
    psf -= xmin;
    w2 += xmin;
    xmin = 0;
    }

  dwpsf = w-w2;
/* Subtract the right pixels to the destination */
  for (y=ymin; y<ymax; y++, psf += dwpsf)
    {
    pix = (float *)check->pix+y*check->width+xmin;
    for (x=w2; x--;)
      *(pix++) += amplitude**(psf++);
    }

  return;
  }


/****** addcheck_resample *****************************************************
PROTO	void addcheck_resample(checkstruct *check, float *thumb, int w, int h,
			int ix,int iy, float zoom, float amplitude)
PURPOSE	Add a resampled thumbnail to a check image (with a multiplicative
	factor).
INPUT	Pointer to the check-image,
	pointer to the thumbnail,
	thumbnail width,
	thumbnail height,
	thumbnail center x coordinate,
	thumbnail center y coordinate,
	inverse zoom factor,
	flux scaling factor.
OUTPUT	-.
NOTES	Outside boundaries are taken into account.
AUTHOR	E. Bertin (IAP)
VERSION	23/09/2013
 ***/
void	addcheck_resample(checkstruct *check, float *thumb, int w, int h,
			int ix, int iy, float step1, float amplitude)
  {
   float	interpm[CHECKINTERPW*CHECKINTERPW],
		*pix2, *mask,*maskt,*pix12, *pixin,*pixin0, *pixout,*pixout0,
		xs1,ys1, x1,y1, x,y, dxm,dym, val, norm, step2;
   int		i,j,k,n,t, *start,*startt, *nmask,*nmaskt,
		ixs2,iys2, ix2,iy2, dix2,diy2, nx2,ny2, iys1a, ny1, hmw,hmh,
		ixs,iys, ix1,iy1, w2,h2;

  step2 = 1.0/step1;
  pix2 = check->pix;
  w2 = check->width;
  h2 = check->height;

  xs1 = -0.5 - ((float)ix - (int)(w*step1)/2 - 0.5)*step2; /* Im1 start coord */
  if ((int)xs1 >= w)
    return;
  ixs2 = 0;			/* Int part of Im2 start x-coord */
  if (xs1<0.0)
    {
    dix2 = (int)(1-xs1*step1);
/*-- Simply leave here if the images do not overlap in x */
    if (dix2 >= w2)
      return;
    ixs2 += dix2;
    xs1 += dix2*step2;
    }
  nx2 = (int)((w-1-xs1)*step1+1);/* nb of interpolated Im2 pixels along x */
  if (nx2>(ix2=w2-ixs2))
    nx2 = ix2;
  if (nx2<=0)
    return;

  ys1 = -0.5 - ((float)iy - (int)(h*step1)/2 - 0.5)*step2; /* Im1 start coord */
  if ((int)ys1 >= h)
    return;
  iys2 = 0;			/* Int part of Im2 start y-coord */
  if (ys1<0.0)
    {
    diy2 = (int)(1-ys1*step1);
/*-- Simply leave here if the images do not overlap in y */
    if (diy2 >= h2)
      return;
    iys2 += diy2;
    ys1 += diy2*step2;
    }
  ny2 = (int)((h-1-ys1)*step1+1);/* nb of interpolated Im2 pixels along y */
  if (ny2>(iy2=h2-iys2))
    ny2 = iy2;
  if (ny2<=0)
    return;

/* Set the yrange for the x-resampling with some margin for interpolation */
  iys1a = (int)ys1;		/* Int part of Im1 start y-coord with margin */
  hmh = CHECKINTERPW/2 - 1;		/* Interpolant start */
  if (iys1a<0 || ((iys1a -= hmh)< 0))
    iys1a = 0;
  ny1 = (int)(ys1+ny2*step2)+CHECKINTERPW-hmh;	/* Interpolated Im1 y size */
  if (ny1>h)					/* with margin */
    ny1 = h;
/* Express everything relative to the effective Im1 start (with margin) */
  ny1 -= iys1a;
  ys1 -= (float)iys1a;

/* Allocate interpolant stuff for the x direction */
  QMALLOC(mask, float, nx2*CHECKINTERPW);	/* Interpolation masks */
  QMALLOC(nmask, int, nx2);			/* Interpolation mask sizes */
  QMALLOC(start, int, nx2);			/* Int part of Im1 conv starts */
/* Compute the local interpolant and data starting points in x */
  hmw = CHECKINTERPW/2 - 1;
  x1 = xs1;
  maskt = mask;
  nmaskt = nmask;
  startt = start;
  for (j=nx2; j--; x1+=step2)
    {
    ixs = (ix1=(int)x1) - hmw;
    dxm = ix1 - x1 - hmw;	/* starting point in the interpolation func */
    if (ixs < 0)
      {
      n = CHECKINTERPW+ixs;
      dxm -= (float)ixs;
      ixs = 0;
      }
    else
      n = CHECKINTERPW;
    if (n>(t=w-ixs))
      n=t;
    *(startt++) = ixs;
    *(nmaskt++) = n;
    norm = 0.0;
    for (x=dxm, i=n; i--; x+=1.0)
      norm += (*(maskt++) = CHECKINTERPF(x));
    norm = norm>0.0? 1.0/norm : 1.0;
    maskt -= n;
    for (i=n; i--;)
      *(maskt++) *= norm;
    }

  QCALLOC(pix12, float, nx2*ny1);	/* Intermediary frame-buffer */

/* Make the interpolation in x (this includes transposition) */
  pixin0 = thumb+iys1a*w;
  pixout0 = pix12;
  for (k=ny1; k--; pixin0+=w, pixout0++)
    {
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    pixout = pixout0;
    for (j=nx2; j--; pixout+=ny1)
      {
      pixin = pixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)**(pixin++);
      *pixout = val;
      }
    }

/* Reallocate interpolant stuff for the y direction */
  QREALLOC(mask, float, ny2*CHECKINTERPW);	/* Interpolation masks */
  QREALLOC(nmask, int, ny2);			/* Interpolation mask sizes */
  QREALLOC(start, int, ny2);			/* Int part of Im1 conv starts */

/* Compute the local interpolant and data starting points in y */
  hmh = CHECKINTERPW/2 - 1;
  y1 = ys1;
  maskt = mask;
  nmaskt = nmask;
  startt = start;
  for (j=ny2; j--; y1+=step2)
    {
    iys = (iy1=(int)y1) - hmh;
    dym = iy1 - y1 - hmh;	/* starting point in the interpolation func */
    if (iys < 0)
      {
      n = CHECKINTERPW+iys;
      dym -= (float)iys;
      iys = 0;
      }
    else
      n = CHECKINTERPW;
    if (n>(t=ny1-iys))
      n=t;
    *(startt++) = iys;
    *(nmaskt++) = n;
    norm = 0.0;
    for (y=dym, i=n; i--; y+=1.0)
      norm += (*(maskt++) = CHECKINTERPF(y));
    norm = norm>0.0? 1.0/norm : 1.0;
    maskt -= n;
    for (i=n; i--;)
      *(maskt++) *= norm;
    }

/* Make the interpolation in y  and transpose once again */
  pixin0 = pix12;
  pixout0 = pix2+ixs2+iys2*w2;
  for (k=nx2; k--; pixin0+=ny1, pixout0++)
    {
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    pixout = pixout0;
    for (j=ny2; j--; pixout+=w2)
      {
      pixin = pixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)**(pixin++);
      *pixout += amplitude*val;
      }
    }

/* Free memory */
  free(pix12);
  free(mask);
  free(nmask);
  free(start);

  return;
  }


/********************************* blankcheck *******************************/
/*
Blank a part of the CHECK-image according to a mask.
*/
void	blankcheck(checkstruct *check, PIXTYPE *mask, int w,int h,
		int xmin,int ymin, PIXTYPE val)
  {
   PIXTYPE	*pixt;
   int		x,y, xmax,ymax,w2,wc;

/* Don't go further if out of frame!! */
  if (xmin+w<0 || xmin>=check->width
	|| ymin+h<0 || ymin>=check->height)
    return;
 
/* Set the image boundaries */
  w2 = w;
  ymax = ymin + h;
  if (ymin<0)
    {
    mask -= ymin*w;
    ymin = 0;
    }
  if (ymax>check->height)
    ymax = check->height;

  xmax = xmin + w;
  if (xmax>check->width)
    {
    w2 -= xmax - check->width;
    xmax = check->width;
    }
  if (xmin<0)
    {
    mask += -xmin;
    w2 -= -xmin;
    xmin = 0;
    }

  w -= w2;
  wc = check->width;
  ymin = ymin*wc+xmin;
  ymax = ymax*wc+xmin;

/* Blank the right pixels in the image */
  for (y=ymin; y<ymax; y+=wc, mask += w)
    {
    pixt = (float *)check->pix + y;
    for (x=w2; x--; pixt++)
      if (*(mask++) > -BIG)
        *pixt = val;
    }

  return;
  }


/******************************** initcheck **********************************/
/*
initialize check-image.
*/
checkstruct	*initcheck(char *filename, checkenum check_type, int next)

  {
   catstruct	*cat;
   checkstruct	*check;

  QCALLOC(check, checkstruct, 1);
  check->type = check_type;
  check->next = next;
  cat = check->cat = new_cat(1);
  strcpy(cat->filename, filename);

  if (next>1)
/*-- Create a "pure" primary HDU */
    {
    init_cat(cat);
    addkeywordto_head(cat->tab, "NEXTEND ", "Number of extensions");
    fitswrite(cat->tab->headbuf, "NEXTEND ", &next, H_INT, T_LONG);
    if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE,"*Error*: cannot open for writing ", filename);
    save_head(cat, cat->tab);
    remove_tabs(cat);
    }
  else
    open_cat(cat, WRITE_ONLY);

  return check;
  }


/******************************** reinitcheck ********************************/
/*
initialize check-image (for subsequent writing).
*/
void	reinitcheck(picstruct *field, checkstruct *check)

  {
   catstruct	*cat;
   tabstruct	*tab;
   wcsstruct	*wcs;
   char		*fitshead;
   PIXTYPE	*ptrf;
   double	dval;
   int		i;

  cat = check->cat;
/* Inherit the field FITS header */
  remove_tabs(cat);
  copy_tab_fromptr(field->tab, cat, 0);
  tab = cat->tab;
  tab->cat = cat;
  if (check->next<=1)
    prim_head(tab);
  check->y = 0;
  fitshead = tab->headbuf;
/* Neutralize possible scaling factors */
  tab->bscale = 1.0;
  tab->bzero = 0.0;
  fitswrite(fitshead, "BSCALE  ", &tab->bscale, H_FLOAT, T_DOUBLE);
  fitswrite(fitshead, "BZERO   ", &tab->bzero, H_FLOAT, T_DOUBLE);
  fitswrite(fitshead, "BITSGN  ", &tab->bitsgn, H_INT, T_LONG);
  if (tab->compress_type != COMPRESS_NONE)
    {
    tab->compress_type = COMPRESS_NONE;
    fitswrite(fitshead, "IMAGECOD", "NONE", H_STRING, T_STRING);
    }
  fitswrite(fitshead, "ORIGIN  ", BANNER, H_STRING, T_STRING);

  switch(check->type)
    {
    case CHECK_IDENTICAL:
    case CHECK_BACKGROUND:
    case CHECK_FILTERED:
    case CHECK_SUBTRACTED:
      tab->bitpix = BP_FLOAT;
      tab->bytepix = 4;
      tab->bitsgn = 0;
      tab->naxisn[0] = check->width = field->width;
      tab->naxisn[1] = check->height = field->height;
      check->npix = field->npix;
      QMALLOC(ptrf, PIXTYPE, check->width);
      check->pix = (void *)ptrf;
      save_head(cat, cat->tab);
      break;

    case CHECK_BACKRMS:
    case CHECK_SUBOBJECTS:
      tab->bitpix = BP_FLOAT;
      tab->bytepix = 4;
      tab->bitsgn = 0;
      tab->naxisn[0] = check->width = field->width;
      tab->naxisn[1] = check->height = field->height;
      check->npix = field->npix;
      QMALLOC(check->pix, PIXTYPE, check->width);
      save_head(cat, cat->tab);
/*---- Allocate memory for replacing the blanked pixels by 0 */
      if (!check->line)
        QMALLOC(check->line, PIXTYPE, field->width);
      break;

    case CHECK_OBJECTS:
    case CHECK_APERTURES:
    case CHECK_PSFPROTOS:
    case CHECK_SUBPSFPROTOS:
    case CHECK_PCPROTOS:
    case CHECK_SUBPCPROTOS:
    case CHECK_PCOPROTOS:
    case CHECK_PROFILES:
    case CHECK_SUBPROFILES:
    case CHECK_SPHEROIDS:
    case CHECK_SUBSPHEROIDS:
    case CHECK_DISKS:
    case CHECK_SUBDISKS:
    case CHECK_PATTERNS:
      tab->bitpix = -32;
      tab->bytepix = 4;
      tab->bitsgn = 0;
      tab->naxisn[0] = check->width = field->width;
      tab->naxisn[1] = check->height = field->height;
      check->npix = field->npix;
      check->overlay = 30*field->backsig;
      QCALLOC(check->pix, PIXTYPE, check->npix);
      save_head(cat, cat->tab);
      break;

    case CHECK_SEGMENTATION:
      tab->bitpix = BP_LONG;
      tab->bytepix = 4;
      tab->bitsgn = 1;
      tab->naxisn[0] = check->width = field->width;
      tab->naxisn[1] = check->height = field->height;
      check->npix = field->npix;
      QCALLOC(check->pix, ULONG, check->npix);
      save_head(cat, cat->tab);
      break;

    case CHECK_MASK:
    case CHECK_SUBMASK:
      tab->bitpix = BP_BYTE;
      tab->bytepix = 1;
      tab->bitsgn = 1;
      tab->naxisn[0] = check->width = field->width;
      tab->naxisn[1] = check->height = field->height;
      check->npix = field->npix;
      save_head(cat, cat->tab);
/*---- Allocate memory */
      if (!check->line)
        QMALLOC(check->line, FLAGTYPE, check->width);
      break;

    case CHECK_ASSOC:
      tab->bitpix = BP_FLOAT;
      tab->bytepix = 4;
      tab->bitsgn = 0;
      tab->naxisn[0] = check->width = field->width;
      tab->naxisn[1] = check->height = field->height;
      check->npix = field->npix;
      QMALLOC(check->pix, PIXTYPE, check->npix);
/*---- Initialize the pixmap to IEEE NaN */
      memset(check->pix, 0xFF, check->npix*sizeof(LONG));
      save_head(cat, cat->tab);
      break;

    case CHECK_MINIBACKGROUND:
    case CHECK_MINIBACKRMS:
      tab->bitpix = BP_FLOAT;
      tab->bytepix = 4;
      tab->bitsgn = 0;
      tab->naxisn[0] = check->width = field->nbackx;
      tab->naxisn[1] = check->height = field->nbacky;
/*---- Scale the WCS information if present */
      if ((wcs=field->wcs))
        {
        dval = wcs->cdelt[0]*field->backw;
        fitswrite(fitshead, "CDELT1  ", &dval, H_EXPO, T_DOUBLE);
        dval = wcs->cdelt[1]*field->backh;
        fitswrite(fitshead, "CDELT2  ", &dval, H_EXPO, T_DOUBLE);
        dval = (wcs->crpix[0]-0.5)/field->backw + 0.5;
        fitswrite(fitshead, "CRPIX1  ", &dval, H_EXPO, T_DOUBLE);
        dval = (wcs->crpix[1]-0.5)/field->backh + 0.5;
        fitswrite(fitshead, "CRPIX2  ", &dval, H_EXPO, T_DOUBLE);

        dval = wcs->cd[0]*field->backw;
        fitswrite(fitshead, "CD1_1   ", &dval, H_EXPO, T_DOUBLE);
        dval = wcs->cd[1]*field->backh;
        fitswrite(fitshead, "CD1_2  ", &dval, H_EXPO, T_DOUBLE);
        dval = wcs->cd[wcs->naxis]*field->backw;
        fitswrite(fitshead, "CD2_1   ", &dval, H_EXPO, T_DOUBLE);
        dval = wcs->cd[wcs->naxis+1]*field->backh;
        fitswrite(fitshead, "CD2_2  ", &dval, H_EXPO, T_DOUBLE);
        }
      check->npix = check->width*check->height;
      QMALLOC(check->pix, PIXTYPE, check->npix);
      if (check->type==CHECK_MINIBACKRMS)
        memcpy(check->pix, field->sigma, check->npix*sizeof(float));
      else
        memcpy(check->pix, field->back, check->npix*sizeof(float));
      save_head(cat, cat->tab);
      write_body(cat->tab, check->pix, check->npix);
      free(check->pix);
      break;

    case CHECK_MAPSOM:
      tab->bitpix = BP_FLOAT;
      tab->bytepix = 4;
      tab->bitsgn = 0;
      tab->naxisn[0] = check->width = field->width;
      tab->naxisn[1] = check->height = field->height;
      check->npix = field->npix;
      QMALLOC(ptrf, PIXTYPE, check->npix);
      check->pix = (void *)ptrf;
      for (i=check->npix; i--;)
        *(ptrf++) = -10.0;
      save_head(cat, cat->tab);
      break;

    case CHECK_OTHER:
      tab->bitpix = BP_FLOAT;
      tab->bytepix = 4;
      tab->bitsgn = 0;
      tab->naxisn[0] = check->width;
      tab->naxisn[1] = check->height;
      check->npix = check->width*check->height;
      QCALLOC(check->pix, PIXTYPE, check->npix);
      save_head(cat, cat->tab);
      break;

    default:
      error(EXIT_FAILURE, "*Internal Error* in ", "reinitcheck()!");
    }

  return;
  }


/******************************** writecheck *********************************/
/*
Write ONE line of npix pixels of a check-image.
*/
void	writecheck(checkstruct *check, PIXTYPE *data, int w)

  {
  if (check->type == CHECK_APERTURES || check->type == CHECK_SUBPSFPROTOS
	|| check->type == CHECK_SUBPCPROTOS || check->type == CHECK_PCOPROTOS
	|| check->type == CHECK_SUBPROFILES || check->type == CHECK_SUBSPHEROIDS
	|| check->type == CHECK_SUBDISKS || check->type == CHECK_OTHER)
    {
    memcpy((PIXTYPE *)check->pix + w*(check->y++), data, w*sizeof(PIXTYPE));
    return;
    }
  else if (check->type == CHECK_SUBOBJECTS)
    {
     int	i;
     PIXTYPE	*pixt;

    pixt = (PIXTYPE *)check->line;
    for (i=w; i--; data++)
      *(pixt++) = (*data>-BIG)? *data:0.0;
    write_body(check->cat->tab, (PIXTYPE *)check->line, w);
   }
  else if (check->type == CHECK_MASK)
    {
     int		i;
     FLAGTYPE		*pixt;

    pixt = (FLAGTYPE *)check->line;
    for (i=w; i--;)
      *(pixt++) = (*(data++)>-BIG)?0:1;
    write_ibody(check->cat->tab, (FLAGTYPE *)check->line, w);
    }
  else if (check->type == CHECK_SUBMASK)
    {
     int		i;
     FLAGTYPE		*pixt;

    pixt = (FLAGTYPE *)check->line;
    for (i=w; i--;)
      *(pixt++) = (*(data++)>-BIG)?1:0;
    write_ibody(check->cat->tab, (FLAGTYPE *)check->line, w);
    }
  else
    write_body(check->cat->tab, data, w);

  return;
  }


/********************************* reendcheck ********************************/
/*
Finish current check-image.
*/
void	reendcheck(picstruct *field, checkstruct *check)
  {
   catstruct	*cat;

  cat = check->cat;
  switch(check->type)
    {
    case CHECK_MINIBACKGROUND:
    case CHECK_MINIBACKRMS:
      break;

    case CHECK_IDENTICAL:
    case CHECK_BACKGROUND:
    case CHECK_BACKRMS:
    case CHECK_FILTERED:
    case CHECK_SUBTRACTED:
      free(check->pix);
      free(check->line);
      check->line = NULL;
      break;

    case CHECK_OBJECTS:
    case CHECK_APERTURES:
    case CHECK_PSFPROTOS:
    case CHECK_SUBPSFPROTOS:
    case CHECK_PCPROTOS:
    case CHECK_SUBPCPROTOS:
    case CHECK_PCOPROTOS:
    case CHECK_PROFILES:
    case CHECK_SUBPROFILES:
    case CHECK_SPHEROIDS:
    case CHECK_SUBSPHEROIDS:
    case CHECK_DISKS:
    case CHECK_SUBDISKS:
    case CHECK_ASSOC:
    case CHECK_PATTERNS:
    case CHECK_MAPSOM:
    case CHECK_OTHER:
      write_body(cat->tab, check->pix, check->npix);
      free(check->pix);
      break;

    case CHECK_SEGMENTATION:
      write_ibody(cat->tab, check->pix, check->npix);
      free(check->pix);
      break;

    case CHECK_MASK:
    case CHECK_SUBMASK:
      {
       int	y;

      for (y=field->ymin; y<field->ymax; y++)
        writecheck(check, &PIX(field, 0, y), field->width);
      free(check->line);
      check->line = NULL;
      break;
      }

    case CHECK_SUBOBJECTS:
      {
       int	y;

      for (y=field->ymin; y<field->ymax; y++)
        writecheck(check, &PIX(field, 0, y), field->width);
      free(check->pix);
      free(check->line);
      check->line = NULL;
      break;
      }

    default:
      error(EXIT_FAILURE, "*Internal Error* in ", "endcheck()!");
    }

  pad_tab(cat, check->npix*sizeof(PIXTYPE));

  return;
  }

/********************************* endcheck **********************************/
/*
close check-image.
*/
void	endcheck(checkstruct *check)
  {
  free_cat(&check->cat,1);
  free(check);

  return;
  }

