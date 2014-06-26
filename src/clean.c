/*
*				clean.c
*
* Remove spurious detections from the catalogue.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		26/06/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"check.h"
#include	"clean.h"
#include	"flag.h"
#include	"image.h"
#include	"objlist.h"


/****** clean_init ***********************************************************
PROTO   objliststruct	*clean_init(void)
PURPOSE Initialize things for CLEANing.
INPUT   -.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 26/06/2014
 ***/
objliststruct	*clean_init(void)
  {
   objliststruct	*cleanobjlist;

  cleanobjlist = objlist_new();

  if (prefs.clean_flag)
    QMALLOC(cleanobjlist->objindex, int, prefs.clean_stacksize);

  return cleanobjlist;
  }


/****** clean_end ************************************************************
PROTO   void clean_end(objliststruct *cleanobjlist)
PURPOSE End things related to CLEANing.
INPUT   -.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 26/06/2014
 ***/
void	clean_end(objliststruct *cleanobjlist)
  {
  objlist_end(cleanobjlist);

  return;
  }


/****** clean_process *********************************************************
PROTO   int clean_process(objliststruct *cleanobjlist, fieldstruct *field,
			objstruct *objin)
PURPOSE Examine object in frame-buffer and put it in the "clean object list" if
	necessary.
INPUT   Pointer to image field,
        Object (source).
OUTPUT  0 if the object was CLEANed, 1 otherwise.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 26/06/2014
 ***/
int	clean_process(objliststruct *cleanobjlist, fieldstruct *field,
			objstruct *objin)
  {
   objstruct		*obj;
   subimagestruct	*subimage;
   double		amp,ampin,alpha,alphain, unitarea,unitareain,beta,val;
   float	       	dx,dy,rlim;
   int			*cleanvictim,
			i,j,k;

  beta = prefs.clean_param;
  unitareain = PI*objin->a*objin->b;
  ampin = objin->fdflux/(2*unitareain*objin->abcor);
  alphain = (pow(ampin/objin->dthresh, 1.0/beta)-1)*unitareain/objin->fdnpix;
  j=0;
  obj = cleanobjlist->obj;
  cleanvictim = cleanobjlist->objindex;
  for (i=0; i<cleanobjlist->nobj; i++, obj++)
    {
    dx = objin->mx - obj->mx;
    dy = objin->my - obj->my;
    rlim = objin->a+obj->a;
    rlim *= rlim;
    if (dx*dx+dy*dy<rlim*CLEAN_ZONE*CLEAN_ZONE)
      {
      if (obj->fdflux < objin->fdflux)
        {
        val = 1+alphain*(objin->cxx*dx*dx+objin->cyy*dy*dy+objin->cxy*dx*dy);
        if (val>1.0
	    && ((float)(val<1e10?ampin*pow(val,-beta) : 0.0) > obj->mthresh))
/*------- the newcomer puts that object in its menu! */
          cleanvictim[j++] = i;
        }
      else
        {
        unitarea = PI*obj->a*obj->b;
        amp = obj->fdflux/(2*unitarea*obj->abcor);
        alpha = (pow(amp/obj->dthresh, 1.0/beta) - 1)*unitarea/obj->fdnpix;
        val = 1+alpha*(obj->cxx*dx*dx+obj->cyy*dy*dy+obj->cxy*dx*dy);
        if (val>1.0
	    && ((float)(val<1e10?amp*pow(val,-beta) : 0.0) > objin->mthresh))
          {
/*------- the newcomer is eaten!! */
          clean_merge(objin, obj);
          if (prefs.blank_flag && (subimage = objin->isoimage))
            {
/*---------- Paste back ``CLEANed'' object pixels before forgetting them */
            pasteimage(field, subimage->image,
			subimage->size[0], subimage->size[1],
			subimage->ipos[0],subimage->ipos[1]);
            QFREE(objin->isoimage);
            }

          return 0;
          }
        }
      }
    }

/* the newcomer eats the content of the menu */
  for (i=j; i--;)
    {
    k = cleanvictim[i];
    obj = cleanobjlist->obj + k;
    clean_merge(obj, objin);
    if (prefs.blank_flag && (subimage = objin->isoimage))
      {
/*---- Paste back ``CLEANed'' object pixels before forgetting them */
      pasteimage(field, subimage->image,
			subimage->size[0], subimage->size[1],
			subimage->ipos[0],subimage->ipos[1]);
      QFREE(objin->isoimage);
      }
    clean_sub(cleanobjlist, k);
    }

  return 1;
  }


/****** clean_add ************************************************************
PROTO   void clean_add(objliststruct *cleanobjlist, objstruct *objin)
PURPOSE Add an object to the "clean object list".
INPUT   Pointer to object.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 25/06/2014
 ***/
void	clean_add(objliststruct *cleanobjlist, objstruct *objin)

  {
   int		margin, y;
   float	hh1,hh2;

/* Compute the max. vertical extent of the object: */
/* First from 2nd order moments, compute y-limit of the 3-sigma ellips... */
  hh1 = objin->cyy - objin->cxy*objin->cxy/(4.0*objin->cxx);
  hh1 = hh1 > 0.0 ? 1/sqrt(3*hh1) : 0.0;
/* ... then from the isophotal limit, which should not be TOO different... */
  hh2 = (objin->ymax-objin->ymin+1.0);
  margin = (int)((hh1>hh2?hh1:hh2)*MARGIN_SCALE+MARGIN_OFFSET+0.49999);
  objin->ycmax = objin->ymax+margin;
/* ... and finally compare with the predefined margin */
  if ((y=(int)(objin->my+0.49999)+prefs.cleanmargin)>objin->ycmax)
    objin->ycmax = y;
  objin->ycmin = objin->ymin-margin;
  if ((y=(int)(objin->my+0.49999)-prefs.cleanmargin)<objin->ycmin)
    objin->ycmin = y;

  if (objlist_addobj(cleanobjlist, objin, NULL) != RETURN_OK)
    error(EXIT_FAILURE, "Not enough memory for ", "CLEANing");

  return;
  }


/****** clean_merge **********************************************************
PROTO   void clean_merge(objstruct *objin, objstruct *objout)
PURPOSE Merge the content of two objects.
INPUT   Pointer to the input object,
	pointer to the output object.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 15/03/2012
 ***/
void	clean_merge(objstruct *objin, objstruct *objout)

  {
   checkstruct	*check;

  if ((check = prefs.check[CHECK_SEGMENTATION]))
    {
     ULONG	*pix;
     ULONG	colorin = objin->number,
		colorout = objout->number;
     int	dx,dx0,dy,dpix;

    pix = (ULONG *)check->pix+check->width*objin->ymin + objin->xmin;
    dx0 = objin->xmax-objin->xmin+1;
    dpix = check->width-dx0;
    for (dy=objin->ymax-objin->ymin+1; dy--; pix += dpix)
      for (dx=dx0; dx--; pix++)
        if (*pix==colorin)
          *pix = colorout;
    }

  objout->fdnpix += objin->fdnpix;
  objout->dnpix += objin->dnpix;
  objout->fdflux += objin->fdflux;
  objout->dflux += objin->dflux;
  objout->dfluxerr += objin->dfluxerr;

  if (objin->fdpeak>objout->fdpeak)
    {
    objout->fdpeak = objin->fdpeak;
    objout->dpeakx = objin->dpeakx;
    objout->dpeaky = objin->dpeaky;
    }
  if (objin->dpeak>objout->dpeak)
    objout->dpeak = objin->dpeak;

  if (objin->xmin<objout->xmin)
    objout->xmin = objin->xmin;
  if (objin->xmax>objout->xmax)
    objout->xmax = objin->xmax;

  if (objin->ymin<objout->ymin)
    objout->ymin = objin->ymin;
  if (objin->ymax>objout->ymax)
    objout->ymax = objin->ymax;

  objout->flag |= (objin->flag & (~(OBJ_MERGED|OBJ_CROWDED)));
  flag_merge(objout, objin);

  return;
  }


/****** clean_sub ************************************************************
PROTO   void clean_sub(objliststruct *cleanobjlist, int objindex)
PURPOSE Remove an object from the "clean object list".
INPUT   Object index.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 26/06/2014
 ***/
void	clean_sub(objliststruct *cleanobjlist, int objindex)
  {
   int	code;

  code = objlist_subobj(cleanobjlist, objindex);

/* Update the object list */
  if (code == RETURN_ERROR)
    error(EXIT_FAILURE, "*Internal Error*: no CLEAN object to remove ",
	"in clean_sub()");
  else if (code == RETURN_FATAL_ERROR)
    error(EXIT_FAILURE, "Not enough memory for ", "CLEANing");

  return;
  }

