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
*	Last modified:		09/06/2014
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

static short		*son, *ok;

/****** deblend_parcelout ****************************************************
PROTO	deblend_parcelout(objliststruct *objlistin, objliststruct *objlistout)
PURPOSE	Divide a list of isophotal detections in several parts (deblending).
INPUT	Input objlist,
	output objlist.
OUTPUT	RETURN_OK if success, RETURN_FATAL_ERROR otherwise (memory overflow).
NOTES	Even if the object is not deblended, the output objlist threshold is
	recomputed if a variable threshold is used.
AUTHOR	E. Bertin (IAP)
VERSION	09/06/2014
TODO	Remove dependency towards global variables son and ok.
 ***/
int	deblend_parcelout(objliststruct *objlistin, objliststruct *objlistout)

  {
   objliststruct	**objlist,
			*debobjlist, *debobjlist2;
   objstruct		*obj;
   subimagestruct	*subimage;
   double		dthresh, dthresh0, value0;
   int			h,i,j,k,l,m, xn,
			nbm = DEBLEND_NBRANCH,
			out;

  out = RETURN_OK;

  xn = prefs.deblend_nthresh;

// Initialize object list tree
  QMALLOC(objlist, objliststruct *,  prefs.deblend_nthresh);
  for (k=0; k<prefs.deblend_nthresh; k++)
    objlist[k] = objlist_new();

// Initialize working object list
  debobjlist2 = objlist_new();
  objlistout->thresh = debobjlist2->thresh = objlistin->thresh;
  subimage = objlistin->subimage;

  for (l=0; l<objlistin->nobj && out==RETURN_OK; l++) {
    dthresh0 = objlistin->obj[l].dthresh;

    objlistout->dthresh = debobjlist2->dthresh = dthresh0;
    if ((out = deblend_addobj(l, objlistin, objlist[0])) == RETURN_FATAL_ERROR)
      goto exit_parcelout;
    if ((out = deblend_addobj(l, objlistin, debobjlist2)) == RETURN_FATAL_ERROR)
      goto exit_parcelout;
    value0 = objlist[0]->obj[0].fdflux*prefs.deblend_mincont;
    ok[0] = (short)1;
    for (k=1; k<xn; k++) {
//---- Calculate threshold
      dthresh = objlistin->obj[l].fdpeak;
      if (dthresh>0.0) {
        if (prefs.detector_type[0] == DETECTOR_PHOTO)
          dthresh = dthresh0 + (dthresh-dthresh0) * (double)k/xn;
        else
          dthresh = dthresh0 * pow(dthresh/dthresh0,(double)k/xn);
      } else
        dthresh = dthresh0;

//-- Build tree (bottom->up)
      if (objlist[k-1]->nobj >= DEBLEND_NSONMAX) {
        out = RETURN_FATAL_ERROR;
        goto exit_parcelout;
      }

      for (i=0; i<objlist[k-1]->nobj; i++) {
        obj = &objlist[k-1]->obj[i];
        if (!(debobjlist = lutz_subextract(subimage, (PIXTYPE)dthresh,
		obj->xmin, obj->xmax, obj->ymin, obj->ymax))) {
          out = RETURN_FATAL_ERROR;
          goto exit_parcelout;
        }

        for (j=h=0; j<debobjlist->nobj; j++)
          if (deblend_belong(j, debobjlist, i, objlist[k-1])) {
            debobjlist->obj[j].dthresh = dthresh;
            m = deblend_addobj(j, debobjlist, objlist[k]);
            if (m == RETURN_FATAL_ERROR || m >= DEBLEND_NSONMAX) {
              out = RETURN_FATAL_ERROR;
              goto exit_parcelout;
            }
            if ( h >= nbm - 1 && !(son = (short *)realloc(son,
			xn*DEBLEND_NSONMAX*(nbm+=16)*sizeof(short)))) {
              out = RETURN_FATAL_ERROR;
              goto exit_parcelout;
            }
            son[k-1+xn*(i+DEBLEND_NSONMAX*(h++))] = (short)m;
            ok[k+xn*m] = (short)1;
          }
        son[k-1+xn*(i+DEBLEND_NSONMAX*h)] = (short)-1;
        objlist_end(debobjlist);
      }
    }

//-- Cut the right branches (top->down)
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
              if ((out = deblend_addobj(j, objlist[k+1], debobjlist2))
			== RETURN_FATAL_ERROR)
                goto exit_parcelout;
            }
          ok[k+xn*i] = (short)0;
        }
      }
    }

    if (ok[0])
      out = deblend_addobj(0, debobjlist2, objlistout);
    else
      out = deblend_gatherup(debobjlist2, objlistout);
  }

exit_parcelout:
  for (k=0; k<xn; k++)
    objlist_end(objlist[k]);

  if (debobjlist)  
    objlist_end(debobjlist);
  if (debobjlist2)  
    objlist_end(debobjlist2);

  return out;
}


/****** deblend_addobj ******************************************************
PROTO	int deblend_addobj(int objnb, objliststruct *objl1,objliststruct *objl2)
PURPOSE	Add an object from an object list to another object list.
INPUT	Object index,
	input objlist,
	output objlist,
OUTPUT	New object index, or RETURN_FATAL_ERROR if a memory overflow occurs.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	27/03/2012
 ***/
int	deblend_addobj(int objnb, objliststruct *objl1, objliststruct *objl2)

  {
   objstruct	*objl2obj;
   pliststruct	*plist1 = objl1->plist, *plist2 = objl2->plist;
   int		fp, i, j, npx, objnb2;

  j = (fp = objl2->npix)*plistsize;
  objnb2 = objl2->nobj;

/* Update the object list */
  if (objl2->nobj)
    {
    if (!(objl2obj = (objstruct *)realloc(objl2->obj,
		(++objl2->nobj) * sizeof(objstruct))))
      goto exit_addobj;
    }
  else
    if (!(objl2obj = (objstruct *)malloc((++objl2->nobj)*sizeof(objstruct))))
      goto exit_addobj;

/* Update the pixel list */
  npx = objl1->obj[objnb].fdnpix;
  if (fp)
    {
    if (!(plist2 = (pliststruct *)realloc(plist2,
		(objl2->npix+=npx) * plistsize)))
      goto exit_addobj;
    }
  else
    if (!(plist2=(pliststruct *)malloc((objl2->npix=npx)*plistsize)))
      goto exit_addobj;

  objl2->obj = objl2obj;
  objl2->plist = plist2;

  plist2 += j;
  for(i=objl1->obj[objnb].firstpix; i!=-1; i=PLIST(plist1+i,nextpix))
    {
    memcpy(plist2, plist1+i, (size_t)plistsize);
    PLIST(plist2,nextpix) = (j+=plistsize);
    plist2 += plistsize;
    }

  PLIST(plist2-=plistsize, nextpix) = -1;

  objl2->obj[objnb2] = objl1->obj[objnb];
  objl2->obj[objnb2].firstpix = fp*plistsize;
  objl2->obj[objnb2].lastpix = j-plistsize;
  return	objnb2;

exit_addobj:

  objl2->nobj--;
  objl2->npix = fp;
  return RETURN_FATAL_ERROR;
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


/****** deblend_alloc ********************************************************
PROTO	deblend_alloc(void)
PURPOSE	Allocate the memory allocated for global pointers in refine.c
INPUT	-,
OUTPUT	-.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	07/03/2012
 ***/
void	deblend_alloc(void)
  {
  QMALLOC(son, short,  prefs.deblend_nthresh*DEBLEND_NSONMAX*DEBLEND_NBRANCH);
  QMALLOC(ok, short,  prefs.deblend_nthresh*DEBLEND_NSONMAX);

  return;
  }


/****** deblend_free *********************************************************
PROTO	deblend_free(void)
PURPOSE	Free the memory allocated for global pointers in refine.c
INPUT	-,
OUTPUT	-.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	07/03/2012
 ***/
void	deblend_free(void)
  {
  QFREE(son);
  QFREE(ok);
  return;
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
VERSION	07/03/2012
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

  objlistout->dthresh = objlistin->dthresh;
  objlistout->thresh = objlistin->thresh;

  QMALLOC(amp, float, nobj);
  QMALLOC(p, float, nobj);
  QMALLOC(n, int, nobj);

  for (i=1; i<nobj; i++)
    scan_preanalyse(objlistin, i, ANALYSE_FULL);

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
/*-- Now we have passed the deblending section, reset thresholds */
    objt->dthresh = objlistin->dthresh;
    objt->thresh = objlistin->thresh;

/* ------------	flag pixels which are already allocated */
    for (pixt=pixelin+objin[i].firstpix; pixt>=pixelin;
	pixt=pixelin+PLIST(pixt,nextpix))
      bmp[(PLIST(pixt,x)-xs) + (PLIST(pixt,y)-ys)*bmwidth] = '\1';

    if ((n[i] = deblend_addobj(i, objlistin, objlistout)) == RETURN_FATAL_ERROR)
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

