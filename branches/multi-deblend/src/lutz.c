/**
* @file		lutz.c
* @brief	Lutz (1980) algorithm to extract connected pixels from an image
		raster.
* @date		04/06/2014
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
PROTO	objliststruct	*lutz_subextract(subimagestruct *subimage,
			PIXTYPE thresh, int xmin, int xmax, int ymin, int ymax)
PURPOSE	C implementation of R.K LUTZ' algorithm for the extraction of
	8-connected pixels in a sub-image above a given threshold and within
	given limits
INPUT	Pointer to sub-image,
	detection threshold,
	minimum x pixel coordinate,
	maximum x pixel coordinate,
	minimum y pixel coordinate,
	maximumy pixel coordinate.
OUTPUT	Pointer to the objlist if no memory allocation problem occured,
	NULL otherwise.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
TODO    Check propagation of flags
VERSION	04/06/2014
 ***/
objliststruct	*lutz_subextract(subimagestruct *subimage, PIXTYPE thresh,
			int xmin, int xmax, int ymin, int ymax) {
   infostruct		curpixinfo,initinfo,
			*info, *store;
   pliststruct		*plist,*pixel;
   status		*psstack;

   char			*marker,
			newmarker;
   int			*start, *end,
			wminus1,hminus1, subw,subh,scansize,
			cn, co, luflag, objnb, pstop, xl,xl2,yl,
			out, minarea, step, nobjm,
			inewsymbol;
   short		trunflag;
   PIXTYPE		*scan,*cscan,*dscan;
   status		cs, ps;


  wminus1 = subimage->field->width-1;
  hminus1 = subimage->field->height-1;
  subw = xmax - xmin;				// Sub-subimage width
  subh = ymax - ymin;				// Sub-subimage height
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

  initinfo.pixnb = 0;
  initinfo.flag = 0;
  initinfo.firstpix = initinfo.lastpix = -1;
  cn = 0;
  imsize = subimage->size[0];
  scan = subimage->image
	+ (ymin - subimage->xmin[1])*imsize
	+ (xmin - subimage->xmin[0]);
  cscan = (subimage->fimage? subimage->fimage : subimage->image)
	+ (ymin - subimage->xmin[1])*imsize
	+ (xmin - subimage->xmin[0]);
// As we only analyse a fraction of the subimage, a step occurs between lines
  step = imsize - subw;

// Allocate memory to store object data */
  objlist = objlist_new();
  nobjm = NSUBOBJ_START;
  if (!(objlist->obj=(objstruct *)malloc(nobjm*sizeof(objstruct))))
    {
    out = RETURN_FATAL_ERROR;
    goto exit_lutz;
    }

// Allocate memory for the pixel list */
  if (!(objlist->plist = (pliststruct *)malloc(subw*subh*plistsize)))
    {
    out = RETURN_FATAL_ERROR;
    goto exit_lutz;
    }

  pixel = plist = objlist->plist;
  objnb = objlist->nobj = 0;
  co = pstop = 0;
  curpixinfo.pixnb = 1;

  for (yl=ymin; yl<=ymax; yl++, cscan += step, scan += imsize)
    {
    ps = COMPLETE;
    cs = NONOBJECT;
    trunflag =  (yl==0 || yl==hminus1) ? OBJ_TRUNC : 0;
    if (yl==ymax)
      cscan = scan = dscan;

    for (xl=xmin; xl<=xmax; xl++)
      {
      newmarker = marker[xl];
      marker[xl] = 0;
      if ((cnewsymbol = (xl!=xmax)?*(cscan++):-BIG) < 0)
        luflag = 0;
      else
        {
        curpixinfo.flag = trunflag;
        luflag =  > thresh?1:0);
        }
      if (luflag)
        {
        if (xl==0 || xl==wminus1)
          curpixinfo.flag |= OBJ_TRUNC;
        PLIST(pixel, nextpix) = -1;
        PLIST(pixel, x) = xl;
        PLIST(pixel, y) = yl;
        PLIST(pixel, value) = scan[xl - xmin];
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
                  if (!(objlist->obj = (objstruct *)
  			realloc(objlist->obj,
			(nobjm+=nobjm/2)* sizeof(objstruct))))
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

  free(dscan);
  free(info);
  free(store);
  free(marker);
  free(psstack);
  free(start);
  free(end);

  if (out == RETURN_OK) {
    if (objlist->nobj) {
      if (!(objlist->obj=(objstruct *)realloc(objlist->obj,
		objlist->nobj*sizeof(objstruct))))
        error(EXIT_FAILURE,"problem with mem. realloc. in lutz()","");
    } else {
      free(objlist->obj);
      objlist->obj = NULL;
      }
    if (cn) {
      if (!(objlist->plist=(pliststruct *)realloc(plist,cn)))
        error(EXIT_FAILURE,"problem with mem. realloc. in lutz()","");
    } else {
      free(objlist->plist);
      objlist->plist = NULL;
    }
    return objlist;
  } else {
    objlist_end(objlist);
    return NULL;
  }
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

