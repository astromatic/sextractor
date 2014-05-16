/**
* @file		lutz.c
* @brief	Lutz (1980) algorithm to extract connected pixels from an image
		raster.
* @date		14/05/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2014 IAP/CNRS/UPMC
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
#include	"lutz.h"
#include	"plist.h"
#include	"scan.h"

/****** lutz_subextract **************************************************//**
PROTO	int lutz_subextract(subimagestruct *subimage, objstruct *objparent,
		objliststruct *objlist)
PURPOSE	C implementation of R.K LUTZ' algorithm for the extraction of
	8-connected pixels in a sub-image around a parent object
INPUT	Pointer to sub-image,
	pointer to parent object,
	pointer to the object list where new detections will be added.
OUTPUT	RETURN_OK if no memory allocation problem occured, RETURN_FATAL_ERROR
	otherwise.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
TODO    Check propagation of flags
VERSION	14/05/2014
 ***/
int	lutz_subextract(subimagestruct *subimage, objstruct *objparent,
		objliststruct *objlist)

  {
   infostruct		curpixinfo,initinfo,
			*info, *store;
   objstruct		*obj;
   pliststruct		*plist,*pixel;
   status		*psstack;

   char			*marker,
			newmarker;
   int			*start, *end,
			xmax,ymax, subw,subh,scansize,
			cn, co, luflag, objnb, pstop, xl,xl2,yl,
			out, minarea, stx,sty,enx,eny, step,
			nobjm = NSUBOBJ_START,
			inewsymbol;
   short		trunflag;
   PIXTYPE		*scan,*cscan,*dscan,
			thresh;
   status		cs, ps;


  stx = objparent->xmin;
  sty = objparent->ymin;
  enx = objparent->xmax + 1;
  eny = objparent->ymax + 1;
  xmax = subimage->field->width-1;
  ymax = subimage->field->height-1;
  subw = enx - stx;				// Sub-subimage width
  subh = eny - sty;				// Sub-subimage height
  scansize = subw + 1;				// Sub-subimage scan size
  QMALLOC(info, infostruct, scansize);
  QMALLOC(store, infostruct, scansize);
  QCALLOC(marker, char, scansize);		// Must be initialized to 0
  QMALLOC(psstack, status, scansize);
  QMALLOC(start, int, scansize);
  QMALLOC(end, int, scansize);
  QMALLOC(dscan, int, scansize);
  dscant = dscan;
  for (i=stacksize; i--;)
    *(dscant++) = -BIG;

  out = RETURN_OK;
  minarea = prefs.deb_maxarea;

  thresh = objlist->dthresh;
  initinfo.pixnb = 0;
  initinfo.flag = 0;
  initinfo.firstpix = initinfo.lastpix = -1;
  cn = 0;
  imsize = subimage->imsize[0];
  scan = subimage->image
	+ (sty-subimage->immin[1])*imsize
	+ (stx-subimage->immin[0]);
  cscan = (subimage->fimage? subimage->fimage : subimage->image)
	+ (sty-subimage->immin[1])*imsize
	+ (stx-subimage->immin[0]);
// As we only analyse a fraction of the subimage, a step occurs between lines
  step = imsize - subw;

// Allocate memory to store object data */
  free(objlist->obj);				// Free existing object if any
  if (!(obj=objlist->obj=(objstruct *)malloc(nobjm*sizeof(objstruct))))
    {
    out = RETURN_FATAL_ERROR;
    plist = NULL;				// Avoid gcc -Wall warnings
    goto exit_lutz;
    }

// Allocate memory for the pixel list */
  free(objlist->plist);
  if (!(objlist->plist = (pliststruct *)malloc(subw*subh*plistsize)))
    {
    out = RETURN_FATAL_ERROR;
    plist = NULL;				// Avoid gcc -Wall warnings
    goto exit_lutz;
    }

  pixel = plist = objlist->plist;
  objnb = objlist->nobj = 0;
  co = pstop = 0;
  curpixinfo.pixnb = 1;

  for (yl=sty; yl<=eny; yl++, cscan += step, scan += imsize)
    {
    ps = COMPLETE;
    cs = NONOBJECT;
    trunflag =  (yl==0 || yl==ymax) ? OBJ_TRUNC : 0;
    if (yl==eny)
      cscan = scan = dscan;

    for (xl=stx; xl<=enx; xl++)
      {
      newmarker = marker[xl];
      marker[xl] = 0;
      if ((cnewsymbol = (xl!=enx)?*(cscan++):-BIG) < 0)
        luflag = 0;
      else
        {
        curpixinfo.flag = trunflag;
        luflag =  > thresh?1:0);
        }
      if (luflag)
        {
        if (xl==0 || xl==xmax)
          curpixinfo.flag |= OBJ_TRUNC;
        PLIST(pixel, nextpix) = -1;
        PLIST(pixel, x) = xl;
        PLIST(pixel, y) = yl;
        PLIST(pixel, value) = scan[xl-stx];
        if (PLISTEXIST(cvalue))
          PLISTPIX(pixel, cvalue) = cnewsymbol;
        curpixinfo.lastpix = curpixinfo.firstpix = cn;
        cn += plistsize;
        pixel += plistsize;
        if (cs != OBJECT)
/*------------------------------- Start Segment -----------------------------*/
          {
          cs = OBJECT;
          if (ps == OBJECT)
              {
              if (start[co] == UNKNOWN)
                {
                marker[xl] = 'S';
                start[co] = xl;
                }
              else  marker[xl] = 's';
              }
          else
            {
            psstack[pstop++] = ps;
            marker[xl] = 'S';
            start[++co] = xl;
            ps = COMPLETE;
            info[co] = initinfo;
            }
          }
        }
/*---------------------------------------------------------------------------*/
      if (newmarker)
/*---------------------------- Process New Marker ---------------------------*/

        {
        if (newmarker == 'S')
          {
          psstack[pstop++] = ps;
          if (cs == NONOBJECT)
            {
            psstack[pstop++] = COMPLETE;
            info[++co] = store[xl];
            start[co] = UNKNOWN;
          }
          else
            lutz_update(&info[co],&store[xl], plist);
          ps = OBJECT;
          }
        else if (newmarker == 's')
          {
          if ((cs == OBJECT) && (ps == COMPLETE))
            {
            pstop--;
            xl2 = start[co];
            lutz_update(&info[co-1],&info[co], plist);
          if (start[--co] == UNKNOWN)
              start[co] = xl2;
            else
              marker[xl2] = 's';
            }
          ps = OBJECT;
          }
        else if (newmarker == 'f')
          ps = INCOMPLETE;
        else if (newmarker == 'F')
          {
          ps = psstack[--pstop];
          if ((cs == NONOBJECT) && (ps == COMPLETE))
            {
          if (start[co] == UNKNOWN)
              {
              if ((int)info[co].pixnb >= minarea)
                {
                if (objlist->nobj>=nobjm)
                  if (!(obj = objlist->obj = (objstruct *)
  			realloc(obj, (nobjm+=nobjm/2)* sizeof(objstruct))))
                    {
                    out = RETURN_FATAL_ERROR;
                    goto exit_lutz;
                    }
                lutz_output(&info[co], objlist);
                }
              }
            else
              {
              marker[end[co]] = 'F';
              store[start[co]] = info[co];
              }
            co--;
            ps = psstack[--pstop];
            }
          }
        }
  
/*---------------------------------------------------------------------------*/

      if (luflag)
        lutz_update(&info[co],&curpixinfo, plist);
    else
        {
        if (cs == OBJECT)
/*-------------------------------- End Segment ------------------------------*/
          {
          cs = NONOBJECT;
          if (ps != COMPLETE)
            {
            marker[xl] = 'f';
            end[co] = xl;
            }
          else
            {
            ps = psstack[--pstop];
            marker[xl] = 'F';
            store[start[co]] = info[co];
            co--;
            }
          }
        }
/*---------------------------------------------------------------------------*/
      }
    }

exit_lutz:

   if (objlist->nobj && out == RETURN_OK)
    {
    if (!(objlist->obj=(objstruct *)realloc(obj,
		objlist->nobj*sizeof(objstruct))))
      error(EXIT_FAILURE,"problem with mem. realloc. in lutz()","");
    }
  else
    {
    free(obj);
    objlist->obj = NULL;
    }

  if (cn && out == RETURN_OK)
    {
    if (!(objlist->plist=(pliststruct *)realloc(plist,cn)))
      error(EXIT_FAILURE,"problem with mem. realloc. in lutz()","");
    }
  else
    {
    free(objlist->plist);
    objlist->plist = NULL;
    }

  free(dscan);
  free(info);
  free(store);
  free(marker);
  free(psstack);
  free(start);
  free(end);

  return  out;
  }


/****** lutz_update **********************************************************
PROTO	void lutz_update(infostruct *infoptr1, infostruct *infoptr2,
		pliststruct *pixel)
PURPOSE	Update the properties of a detection each time one of its pixels is
	scanned.
INPUT	Pointer to detection info,
	pointer to pointer to new object info,
	pointer to pixel list.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	15/03/2012
 ***/
void	lutz_update(infostruct *infoptr1, infostruct *infoptr2,
		pliststruct *pixel)

  {
  infoptr1->pixnb += infoptr2->pixnb;
  infoptr1->flag |= infoptr2->flag;
  if (infoptr1->firstpix == -1)
    {
    infoptr1->firstpix = infoptr2->firstpix;
    infoptr1->lastpix = infoptr2->lastpix;
    }
  else if (infoptr2->lastpix != -1)
    {
    PLIST(pixel+infoptr1->lastpix, nextpix) = infoptr2->firstpix;
    infoptr1->lastpix = infoptr2->lastpix;
    }

  return;
  }


/****** lutz_output **********************************************************
PROTO	void lutz_output(infostruct *info, objliststruct *objlist)
PURPOSE	Convert an output detection to a new (sub-)object.
INPUT	-.
OUTPUT	-.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	15/03/2012
 ***/
/*
Build the object structure.
*/
void  lutz_output(infostruct *info, objliststruct *objlist)

  {
  objstruct  *obj = objlist->obj+objlist->nobj;

  memset(obj, 0, (size_t)sizeof(objstruct));
  obj->firstpix = info->firstpix;
  obj->lastpix = info->lastpix;
  obj->flag = info->flag;
  objlist->npix += info->pixnb;

  scan_preanalyse(objlist, objlist->nobj, ANALYSE_FAST);

  objlist->nobj++;

  return;
  }

