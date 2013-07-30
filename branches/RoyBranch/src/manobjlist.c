/*
*				manobjlist.c
*
* Manage object lists.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		11/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"plist.h"

/********************************* belong ************************************/
/*
say if an object is "included" in another.
*/
int	belong(int corenb, objliststruct *coreobjlist,
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


/********************************* addobj ************************************/
/*
Add an object to an objlist.
*/
int	addobj(int objnb, objliststruct *objl1, objliststruct *objl2)

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



