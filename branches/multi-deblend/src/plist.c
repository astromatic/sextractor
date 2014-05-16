/*
*				plist.c
*
* Manage pixel lists.
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
*	Last modified:		14/05/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<stdio.h>
#include	<stdlib.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"plist.h"


/******************************** createblank *******************************
PROTO   int createblank(int no, objliststruct *objlist)
PURPOSE Create pixel map for BLANKing.
INPUT   objlist number,
        objlist pointer,
OUTPUT  RETURN_OK if success, RETURN_FATAL_ERROR otherwise (memory overflow).
NOTES   -.
AUTHOR  E. Bertin (IAP & Leiden & ESO)
VERSION 12/01/2012
 ***/
int	createblank(objliststruct *objlist, int no)

  {
   objstruct	*obj;
   pliststruct	*pixel, *pixt;
   int		i, n, pos, xmin,ymin, w, dflag;
   PIXTYPE	*pix, *dpix, *pt;

  obj = objlist->obj+no;
  pixel = objlist->plist;
  dpix = NULL;			/* To avoid gcc -Wall warnings */

  obj->subx = xmin = obj->xmin;
  obj->suby = ymin = obj->ymin;
  obj->subw = w = obj->xmax - xmin + 1;
  obj->subh = obj->ymax - ymin + 1;

  n = w*obj->subh;
  if (!(obj->blank = pix = (PIXTYPE *)malloc(n*sizeof(PIXTYPE))))
    return RETURN_FATAL_ERROR;
  pt = pix;
  for (i=n; i--;)
    *(pt++) = -BIG;

  for (i=obj->firstpix; i!=-1; i=PLIST(pixt,nextpix))
    {
    pixt = pixel+i;
    pos = (PLIST(pixt,x)-xmin) + (PLIST(pixt,y)-ymin)*w;
    *(pix+pos) = PLIST(pixt, value);
    }

  return RETURN_OK;
  }


/****************************** init_plist ************************************
PROTO	pliststruct *init_plist(void)
PURPOSE	initialize a pixel-list and its components.
INPUT	-.
OUTPUT  -.
NOTES   The preparation of components relies on the preferences.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 21/02/2012
 ***/
void	init_plist(void)

  {
   pbliststruct	*pbdum = NULL;
   int		i;

  plistsize = sizeof(pbliststruct);
  plistoff_value = (char *)&pbdum->value - (char *)pbdum;

  if (prefs.filter_flag)
    {
    plistexist_cvalue = 1;
    plistoff_cvalue = plistsize;
    plistsize += sizeof(PIXTYPE);
    }
  else
    {
    plistexist_cvalue = 0;
    plistoff_cvalue = plistoff_value;
    }

  if (FLAG(obj2.imaflags))
    {
    plistexist_flag = 1;
    for (i=0; i<prefs.nfimage; i++)
      {
      plistoff_flag[i] = plistsize;
      plistsize += sizeof(FLAGTYPE);
      }
    }
  else
    plistexist_flag = 0;

  if (prefs.weight_flag[0] && FLAG(obj2.wflag))
    {
    plistexist_wflag = 1;
    plistoff_wflag = plistsize;
    plistsize += sizeof(FLAGTYPE);
    }
  else
    plistexist_wflag = 0;

  if (prefs.weight_flag[0])
    {
    plistexist_var = 1;
    plistoff_var = plistsize;
    plistsize += sizeof(PIXTYPE);
    }
  else
    plistexist_var = 0;

  if (prefs.weight_flag[0])
    {
    plistexist_dthresh = 1;
    plistoff_dthresh = plistsize;
    plistsize += sizeof(PIXTYPE);
    }
  else
    plistexist_dthresh = 0;

  return;
  }

