/**
* @file		lutz.c
* @brief	Lutz (1980) algorithm to extract connected pixels from an image
		raster.
* @date		05/01/2015
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2015 IAP/CNRS/UPMC
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
#include	"objlist.h"
#include	"plist.h"
#include	"scan.h"
#include	"weight.h"

/****** lutz_subextract **************************************************//**
PROTO	objliststruct	*lutz_subextract(subimagestruct *subimage,
			PIXTYPE thresh, int xmin, int xmax, int ymin, int ymax,
			int extflags)
PURPOSE	C implementation of R.K LUTZ' algorithm for the extraction of
	8-connected pixels in a sub-image above a given threshold and within
	given limits
INPUT	Pointer to sub-image,
	detection threshold: relative if weighted threshold, absolute otherwise,
	minimum x pixel coordinate,
	maximum x pixel coordinate (+1),
	minimum y pixel coordinate,
	maximum y pixel coordinate (+1),
	binary combination of extraction flags.
OUTPUT	Pointer to the objlist if no memory allocation problem occured,
	NULL otherwise.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
TODO    Check propagation of flags and filtering.
VERSION	05/01/2015
 ***/
objliststruct	*lutz_subextract(subimagestruct *subimage, PIXTYPE thresh,
			int xmin, int xmax, int ymin, int ymax, int extflags) {

   objliststruct	*objlist;
   infostruct		curpixinfo,initinfo,
			*info, *store;
   pliststruct		*plist,*pixel;
   status		*psstack,
			cs, ps;
   short		trunflag;
   char			*marker,
			newmarker;
   PIXTYPE		*scan, *cscan, *varscan, *cvarscan, *dumscan, *dumscant,
			varbadthresh, relthresh, cnewsymbol;
   int			*start, *end,
			i, xl,xl2,yl, wminus1,hminus1, subw,subh,
			scansize, imsize, pos,
			varflag, varthreshflag,
			cn, co, luflag, pstop,
			out, minarea, inewsymbol;

// Flag the presence of a variance map
  varflag = subimage->imvar || subimage->fimvar;
// Flag variable thresholding
  varthreshflag = (extflags&SUBEX_VARTHRESH) && varflag;
// Relative threshold will be used only for weighted thresholding
  relthresh = thresh;

  wminus1 = subimage->field->width - 1 - xmin;
  hminus1 = subimage->field->height - 1 - ymin;
  subw = xmax - xmin;				// Sub-subimage width
  subh = ymax - ymin;				// Sub-subimage height
  scansize = subw + 1;				// Sub-subimage scan size
  QMALLOC(info, infostruct, scansize);
  QMALLOC(store, infostruct, scansize);
  QCALLOC(marker, char, scansize);		// Must be initialized to 0
  QMALLOC(psstack, status, scansize);
  QMALLOC(start, int, scansize);
  QMALLOC(end, int, scansize);
  QMALLOC(dumscan, PIXTYPE, scansize);
  dumscant = dumscan;
  for (i=scansize; i--;)
    *(dumscant++) = -BIG;

  out = RETURN_OK;
  minarea = prefs.deb_maxarea;

// Variance threshold for bad pixels (= infinite variance)
  varbadthresh = varflag && subimage->wfield?
			subimage->wfield->weight_thresh : 0.0;
  if (varbadthresh>BIG*WTHRESH_CONVFAC)
    varbadthresh = BIG*WTHRESH_CONVFAC;

  initinfo.pixnb = 0;
  initinfo.flag = 0;
  initinfo.firstpix = initinfo.lastpix = -1;
  cn = 0;
  imsize = subimage->size[0];
  pos = (ymin - subimage->xmin[1])*imsize + (xmin - subimage->xmin[0]);
  scan = subimage->image + pos;
  cscan = subimage->fimage? subimage->fimage + pos : scan;
  varscan = subimage->imvar? subimage->imvar + pos : NULL;
  cvarscan = subimage->fimvar? subimage->fimvar + pos : varscan;

// Allocate memory to store object data
  objlist = objlist_new();

// Allocate memory for the pixel list
  if (!(objlist->plist = (pliststruct *)malloc(subw*subh*plistsize))) {
    out = RETURN_FATAL_ERROR;
    goto exit_lutz;
  }

  pixel = plist = objlist->plist;
  co = pstop = 0;
  curpixinfo.pixnb = 1;

  for (yl=0; yl<=subh; yl++) {
    ps = COMPLETE;
    cs = NONOBJECT;
    trunflag =  (yl==-ymin || yl==hminus1) ? OBJ_TRUNC : 0;
    if (yl==subh)
      cvarscan = varscan = cscan = scan = dumscan;
    for (xl=0; xl<=subw; xl++) {
      newmarker = marker[xl];
      marker[xl] = 0;
      cnewsymbol = (xl==subw)? -BIG-1 : cscan[xl];

      if (varthreshflag)
        thresh = relthresh*sqrt((xl==subw || yl==subh)? 0.0 : cvarscan[xl]);

      luflag =  (cnewsymbol > thresh);

      if (luflag) {
        curpixinfo.flag = trunflag;
        if (xl==-xmin || xl==wminus1)
          curpixinfo.flag |= OBJ_TRUNC;
        PLIST(pixel, x) = xl + xmin;
        PLIST(pixel, y) = yl + ymin;
        PLIST(pixel, value) = scan[xl];
        if (PLISTEXIST(cvalue))
          PLISTPIX(pixel, cvalue) = cnewsymbol;
        if (PLISTEXIST(var))
          PLISTPIX(pixel, var) = varscan[xl];
        PLIST(pixel, nextpix) = -1;
        curpixinfo.lastpix = curpixinfo.firstpix = cn;
        cn += plistsize;
        pixel += plistsize;
        if (cs != OBJECT) {
// ----- Start Segment
          cs = OBJECT;
          if (ps == OBJECT) {
              if (start[co] == UNKNOWN) {
                marker[xl] = 'S';
                start[co] = xl;
              }
              else  marker[xl] = 's';
            } else {
            psstack[pstop++] = ps;
            marker[xl] = 'S';
            start[++co] = xl;
            ps = COMPLETE;
            info[co] = initinfo;
          }
        }
      }

      if (newmarker) {
//---- Process New Marker
        if (newmarker == 'S') {
          psstack[pstop++] = ps;
          if (cs == NONOBJECT) {
            psstack[pstop++] = COMPLETE;
            info[++co] = store[xl];
            start[co] = UNKNOWN;
          }
          else
            lutz_update(&info[co],&store[xl], plist);
          ps = OBJECT;
        } else if (newmarker == 's') {
          if ((cs == OBJECT) && (ps == COMPLETE)) {
            pstop--;
            xl2 = start[co];
            lutz_update(&info[co-1],&info[co], plist);
            if (start[--co] == UNKNOWN)
              start[co] = xl2;
            else
              marker[xl2] = 's';
          }
          ps = OBJECT;
        } else if (newmarker == 'f')
          ps = INCOMPLETE;
        else if (newmarker == 'F') {
          ps = psstack[--pstop];
          if ((cs == NONOBJECT) && (ps == COMPLETE)) {
            if (start[co] == UNKNOWN) {
              if ((int)info[co].pixnb >= minarea
                && lutz_output(&info[co], objlist) != RETURN_OK) {
                out = RETURN_FATAL_ERROR;
                goto exit_lutz;
              }
            } else {
              marker[end[co]] = 'F';
              store[start[co]] = info[co];
            }
            co--;
            ps = psstack[--pstop];
          }
        }
      }
  
      if (luflag)
        lutz_update(&info[co],&curpixinfo, plist);
      else {
        if (cs == OBJECT) {
//------ End Segment
          cs = NONOBJECT;
          if (ps != COMPLETE) {
            marker[xl] = 'f';
            end[co] = xl;
          } else {
            ps = psstack[--pstop];
            marker[xl] = 'F';
            store[start[co]] = info[co];
            co--;
          }
        }
      }
    }

    scan += imsize;
    cscan += imsize;
    if (varscan)
      varscan += imsize;
    if (cvarscan)
      cvarscan += imsize;
  }

exit_lutz:

  free(dumscan);
  free(info);
  free(store);
  free(marker);
  free(psstack);
  free(start);
  free(end);

  if (out == RETURN_OK) {
    if (objlist->nobj) {
      objlist->nobjmax = objlist->nobj;
      if (!(objlist->obj=(objstruct *)realloc(objlist->obj,
		objlist->nobjmax*sizeof(objstruct))))
        error(EXIT_FAILURE,"problem with mem. realloc. in lutz()","");
    } else {
      free(objlist->obj);
      objlist->nobjmax = 0;
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


/****** lutz_update ******************************************************//**
Update basic detection properties based on new pixel information
@param[in] infoout	Pointer to the destination detection info
@param[in] infoin	Pointer to the source detection info
@param[in] pixel	Pointer to the new pixel
@author 		E. Bertin (IAP)
@date			05/01/2015
 ***/
void	lutz_update(infostruct *infoout, infostruct *infoin,
			pliststruct *pixel) {

  infoout->pixnb += infoin->pixnb;
  infoout->flag |= infoin->flag;
  if (infoout->firstpix == -1) {
    infoout->firstpix = infoin->firstpix;
    infoout->lastpix = infoin->lastpix;
  } else if (infoin->lastpix != -1) {
    PLIST(pixel + infoout->lastpix, nextpix) = infoin->firstpix;
    infoout->lastpix = infoin->lastpix;
  }

  return;
}


/****** lutz_output ******************************************************//**
Convert a basic detection to a new object in the provided objlist
@param[in]  info	Pointer to the input detection info
@param[in]  objlist	Pointer to the destination objlist
@param[out] RETURN_OK, or RETURN_ERROR / RETURN_FATAL_ERROR in case of a problem
@author 		E. Bertin (IAP)
@date			24/06/2014
 ***/
int  lutz_output(infostruct *info, objliststruct *objlist) {

   objstruct	obj;

  memset(&obj, 0, (size_t)sizeof(objstruct));
  obj.firstpix = info->firstpix;
  obj.lastpix = info->lastpix;
  obj.flag = info->flag;
  objlist->npix += info->pixnb;

  scan_preanalyse(&obj, objlist->plist, ANALYSE_FAST);

  return objlist_addobj(objlist, &obj, NULL);	// Pixels are already in plist
}
