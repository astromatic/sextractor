/*
*				deblend.c
*
* Deblend sources based on their pixel lists.
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
*	Last modified:		21/10/2014
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
#include	"deblend.h"
#include	"lutz.h"
#include	"plist.h"
#include	"subimage.h"
#include	"scan.h"

/****** deblend_parcelout ****************************************************
PROTO	deblend_parcelout(objstruct *objin, pliststruct *plist)
PURPOSE	Divide a list of isophotal detections in several parts (deblending).
INPUT	Pointer to input obj,
	pointer to subimage,
	pointer to pixel list
OUTPUT	RETURN_OK if success, RETURN_FATAL_ERROR otherwise (memory overflow).
NOTES	Even if the object is not deblended, the output objlist threshold is
	recomputed if a variable threshold is used.
AUTHOR	E. Bertin (IAP)
VERSION	21/10/2014
TODO	Checkout "out" variable at exit
 ***/
objliststruct	*deblend_parcelout(objstruct *objin, subimagestruct *subimage,
				pliststruct *plist, PIXTYPE thresh)

  {
   objliststruct	**objlist,
			*debobjlist, *debobjlist2, *objlistout;
   objstruct		*obj;
   double		dthresh, value0;
   short		*son, *ok;
   int			h,i,j,k,m, xn,
			nbm = DEBLEND_NBRANCH,
			out;

  out = RETURN_OK;

  xn = prefs.deblend_nthresh;

// Allocate memory for son and selection tables
  QMALLOC(son, short, xn * DEBLEND_NSONMAX*DEBLEND_NBRANCH);
  QMALLOC(ok, short, xn * DEBLEND_NSONMAX);

// Initialize object list tree
  QMALLOC(objlist, objliststruct *,  xn);
  for (k=0; k<xn; k++)
    objlist[k] = objlist_new();

// Initialize working object list
  debobjlist2 = objlist_new();
  objlistout = objlist_new();

  if ((out = objlist_addobj(objlist[0], objin, plist)) == RETURN_FATAL_ERROR)
    goto exit_parcelout;
  if ((out = objlist_addobj(debobjlist2, objin, plist)) == RETURN_FATAL_ERROR)
    goto exit_parcelout;

// Calculate flux threshold from DEBLEND_MINCONT
  value0 = objin->fdflux*prefs.deblend_mincont;

  ok[0] = (short)1;

  for (k=1; k<xn; k++) {
//-- Calculate threshold
    dthresh = objin->fdpeak;
    if (dthresh>0.0) {
      if (prefs.detector_type[0] == DETECTOR_PHOTO)
        dthresh = thresh + (dthresh-thresh) * (double)k/xn;
      else
        dthresh = thresh * pow(dthresh/thresh, (double)k/xn);
    } else
      dthresh = thresh;

//-- Build tree (bottom->up)
    if (objlist[k-1]->nobj >= DEBLEND_NSONMAX) {
      out = RETURN_FATAL_ERROR;
      goto exit_parcelout;
    }

    for (i=0; i<objlist[k-1]->nobj; i++) {
      obj = &objlist[k-1]->obj[i];
      if (!(debobjlist = lutz_subextract(subimage, (PIXTYPE)dthresh,
		obj->xmin, obj->xmax + 1, obj->ymin, obj->ymax + 1,
		SUBEX_NONE))) {
        out = RETURN_FATAL_ERROR;
        goto exit_parcelout;
      }

      for (j=h=0; j<debobjlist->nobj; j++)
        if (deblend_belong(j, debobjlist, i, objlist[k-1])) {
          debobjlist->obj[j].dthresh = dthresh;
          if ((m = objlist[k]->nobj) >= DEBLEND_NSONMAX
		|| objlist_addobj(objlist[k], debobjlist->obj + j,
			debobjlist->plist) == RETURN_FATAL_ERROR) {
            out = RETURN_FATAL_ERROR;
            goto exit_parcelout;
          }
          if (h >= (nbm - 1) && !(son = (short *)realloc(son,
		xn*DEBLEND_NSONMAX*(nbm+=16)*sizeof(short)))) {
            out = RETURN_FATAL_ERROR;
            goto exit_parcelout;
          }
          son[k-1+xn*(i+DEBLEND_NSONMAX*(h++))] = (short)m;
          ok[k+xn*m] = (short)1;
        }
      son[k-1+xn*(i+DEBLEND_NSONMAX*h)] = (short)-1;
      objlist_end(debobjlist);
      debobjlist = NULL;
    }
  }

// Cut the right branches (top->down)
  for (k = xn-2; k>=0; k--) {
    obj = objlist[k+1]->obj;
    for (i=0; i<objlist[k]->nobj; i++) {
      for (m=h=0; (j=(int)son[k+xn*(i+DEBLEND_NSONMAX*h)])!=-1; h++) {
        if (obj[j].fdflux - obj[j].dthresh*obj[j].fdnpix > value0)
          m++;
        ok[k+xn*i] &= ok[k+1+xn*j];
      }
      if (m > 1) {
        for (h=0; (j=(int)son[k+xn*(i+DEBLEND_NSONMAX*h)])!=-1; h++)
          if (ok[k+1+xn*j] && obj[j].fdflux-obj[j].dthresh*obj[j].fdnpix
			> value0) {
            objlist[k+1]->obj[j].flag |= OBJ_MERGED	// Merge flag on
		| ((OBJ_ISO_PB|OBJ_APERT_PB|OBJ_OVERFLOW)
			& debobjlist2->obj[0].flag);
            if ((out = objlist_addobj(debobjlist2, objlist[k+1]->obj + j, 
				objlist[k+1]->plist)) == RETURN_FATAL_ERROR)
              goto exit_parcelout;
          }
        ok[k+xn*i] = (short)0;
      }
    }
  }

  if (ok[0])
    out = objlist_addobj(objlistout, debobjlist2->obj, debobjlist2->plist);
  else {
/*-- Now we have passed the deblending section, reset thresholds */
    obj = debobjlist2->obj + 1;
    for (i=debobjlist2->nobj; --i; obj++) {
      obj->dthresh = objin->dthresh;
      obj->thresh = objin->thresh;
    }
    out = deblend_gatherup(debobjlist2, objlistout);
  }

exit_parcelout:

  for (k=0; k<xn; k++)
    objlist_end(objlist[k]);

  if (debobjlist)  
    objlist_end(debobjlist);
  if (debobjlist2)  
    objlist_end(debobjlist2);

  free(son);
  free(ok);

  return objlistout;
}



/****** deblend_belong ******************************************************
PROTO	int deblend_belong(int corenb, objliststruct *coreobjlist,
	       int shellnb, objliststruct *shellobjlist)
PURPOSE	Tell if an object is "included" in another.
INPUT	"Core" object index,
	"core" objlist,
	"shell" object index,
	"shell" objlist.
OUTPUT	1 if the "core" object is included in the "shell" object; 0 otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	27/03/2012
 ***/
int	deblend_belong(int corenb, objliststruct *coreobjlist,
	       int shellnb, objliststruct *shellobjlist)

  {
   objstruct	*cobj = &(coreobjlist->obj[corenb]),
		*sobj = &(shellobjlist->obj[shellnb]);
   pliststruct	*cpl = coreobjlist->plist, *spl = shellobjlist->plist, *pixt;

   int		xc=PLIST(cpl+cobj->firstpix,x), yc=PLIST(cpl+cobj->firstpix,y);

  for (pixt = spl+sobj->firstpix; pixt>=spl; pixt = spl+PLIST(pixt,nextpix))
    if ((PLIST(pixt,x) == xc) && (PLIST(pixt,y) == yc))
      return 1;

  return 0;
  }


/****** deblend_gatherup *****************************************************
PROTO	int deblend_gatherup(objliststruct *objlistin,objliststruct *objlistout)
PURPOSE	Collect faint remaining pixels and allocate them to their most probable
	progenitor.
INPUT	Input objlist,
	output objlist,
OUTPUT	RETURN_OK if success, RETURN_FATAL_ERROR otherwise (memory overflow).
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/06/2014
 ***/
int	deblend_gatherup(objliststruct *objlistin, objliststruct *objlistout)

  {
   char		*bmp;
   float	*amp, *p, dx,dy, drand, dist, distmin;
   objstruct	*objin = objlistin->obj, *objout, *objt;

   pliststruct	*pixelin = objlistin->plist, *pixelout, *pixt,*pixt2;

   int		i,k,l, *n, iclst, npix, bmwidth,
		nobj = objlistin->nobj, xs,ys, x,y, out;

  out = RETURN_OK;

  QMALLOC(amp, float, nobj);
  QMALLOC(p, float, nobj);
  QMALLOC(n, int, nobj);

  for (i=1; i<nobj; i++)
    scan_preanalyse(objlistin->obj + i, objlistin->plist, ANALYSE_FULL);

  p[0] = 0.0;
  bmwidth = objin->xmax - (xs=objin->xmin) + 1;
  npix = bmwidth * (objin->ymax - (ys=objin->ymin) + 1);
  if (!(bmp = (char *)calloc(1, npix*sizeof(char))))
    {
    bmp = 0;
    out = RETURN_FATAL_ERROR;
    goto exit_gatherup;
    }

  for (objt = objin+(i=1); i<nobj; i++, objt++)
    {
/* ------------	flag pixels which are already allocated */
    for (pixt=pixelin+objin[i].firstpix; pixt>=pixelin;
	pixt=pixelin+PLIST(pixt,nextpix))
      bmp[(PLIST(pixt,x)-xs) + (PLIST(pixt,y)-ys)*bmwidth] = '\1';

    if ((n[i] = objlist_addobj(objlistout, objlistin->obj + i,
			objlistin->plist)) == RETURN_FATAL_ERROR)
      {
      out = RETURN_FATAL_ERROR;
      goto exit_gatherup;
      }
    dist = objt->fdnpix/(2*PI*objt->abcor*objt->a*objt->b);
    amp[i] = dist<70.0? objt->dthresh*expf(dist) : 4.0*objt->fdpeak;

/* ------------ limitate expansion ! */
    if (amp[i]>4.0*objt->fdpeak)
      amp[i] = 4.0*objt->fdpeak;
    }

  objout = objlistout->obj;		/* DO NOT MOVE !!! */

  if (!(pixelout=(pliststruct *)realloc(objlistout->plist,
	(objlistout->npix + npix)*plistsize)))
    {
    out = RETURN_FATAL_ERROR;
    goto exit_gatherup;
    }

  objlistout->plist = pixelout;
  k = objlistout->npix;
  iclst = 0;				/* To avoid gcc -Wall warnings */
  for (pixt=pixelin+objin->firstpix; pixt>=pixelin;
	pixt=pixelin+PLIST(pixt,nextpix))
    {
    x = PLIST(pixt,x);
    y = PLIST(pixt,y);
    if (!bmp[(x-xs) + (y-ys)*bmwidth])
      {
      pixt2 = pixelout + (l=(k++*plistsize));
      memcpy(pixt2, pixt, (size_t)plistsize);
      PLIST(pixt2, nextpix) = -1;
      distmin = 1e+31;
      for (objt = objin+(i=1); i<nobj; i++, objt++)
        {
        dx = x - objt->mx;
        dy = y - objt->my;
        dist=0.5*(objt->cxx*dx*dx+objt->cyy*dy*dy+objt->cxy*dx*dy)/objt->abcor;
        p[i] = p[i-1] + (dist<70.0?amp[i]*expf(-dist) : 0.0);
        if (dist<distmin)
          {
          distmin = dist;
          iclst = i;
          }
        }			
      if (p[nobj-1] > 1.0e-31)
        {
        drand = p[nobj-1]*rand()/RAND_MAX;
        for (i=1; i<nobj && p[i]<drand; i++);
        if (i==nobj)
          i=iclst;
	}
      else
        i = iclst;
      objout[n[i]].lastpix=PLIST(pixelout+objout[n[i]].lastpix,nextpix)=l;
      }
    }

  objlistout->npix = k;
  if (!(objlistout->plist = (pliststruct *)realloc(pixelout,
	objlistout->npix*plistsize)))
    error (-1, "Not enough memory to update pixel list in ", "gatherup()");

exit_gatherup:

  free(bmp);
  free(amp);
  free(p);
  free(n);

  return out;
  }

