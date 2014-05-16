/*
*				flag.c
*
* Manage external flag maps.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1997-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		22/02/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<limits.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"plist.h"
#include	"flag.h"

/****** flag_get ************************************************************
PROTO   void	flag_get(objstruct *obj, pliststruct *pixel)
PURPOSE Return the composited flags extracted from the flag-maps.
INPUT   obj structure,
	pixel list.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 22/02/2012
 ***/
void	flag_get(objstruct *obj, pliststruct *pixel)
  {
   pliststruct	*pixt;
   FLAGTYPE	imaflags,cimaflags,
		*flagstack, *fs;
   int		i,n,nmax,nc,nflag,nflag0,
		*nflagstack, *nfs;

  for (i=0; i<prefs.nfimage; i++)
    {
    nmax = 0;
    imaflags = 0;
    switch(prefs.flag_type[i])
      {
      case FLAG_OR:
        for (pixt=pixel+obj->firstpix;pixt>=pixel;
		pixt=pixel+PLIST(pixt,nextpix))
          if ((cimaflags = PLISTFLAG(pixt,flag[i])))
            {
            imaflags |= cimaflags;
            nmax++;
            }
        break;
      case FLAG_AND:
        for (pixt=pixel+obj->firstpix;pixt>=pixel;
		pixt=pixel+PLIST(pixt,nextpix))
          if ((cimaflags = PLISTFLAG(pixt,flag[i])))
            {
            imaflags &= cimaflags;
            nmax++;
            }
        break;
      case FLAG_MIN:
        imaflags = UINT_MAX;
        for (pixt=pixel+obj->firstpix;pixt>=pixel;
		pixt=pixel+PLIST(pixt,nextpix))
          if ((cimaflags = PLISTFLAG(pixt,flag[i])))
            {
            if (cimaflags<imaflags)
              {
              imaflags = cimaflags;
              nmax = 1;
              }
            else if (cimaflags==imaflags)
              nmax++;
            }
        if (!nmax)
          imaflags = 0;
        break;
      case FLAG_MAX:
        imaflags = 0;
        for (pixt=pixel+obj->firstpix;pixt>=pixel;
		pixt=pixel+PLIST(pixt,nextpix))
          if ((cimaflags = PLISTFLAG(pixt,flag[i])))
            {
            if (cimaflags>imaflags)
              {
              imaflags = cimaflags;
              nmax = 1;
              }
            else if (cimaflags==imaflags)
              nmax++;
            }
        if (!nmax)
          imaflags = 0;
        break;
      case FLAG_MOST:
/*------ Allocate memory for the buffers */
        nflag = FLAG_BUFSIZE;
        QCALLOC(flagstack, FLAGTYPE, nflag);
        QCALLOC(nflagstack, int, nflag);
/*------ Count flag values */
        for (pixt=pixel+obj->firstpix;pixt>=pixel;
		pixt=pixel+PLIST(pixt,nextpix))
          if ((cimaflags = PLISTFLAG(pixt,flag[i])))
            {
            for (n=nflag, fs=flagstack, nfs=nflagstack; n-- && *nfs; nfs++)
              if (*(fs++) == cimaflags)
                {
                (*nfs)++;
                break;
                }
            if (n<0)
              {
              nflag0 = nflag;
              nflag += FLAG_BUFSIZE;
              QREALLOC(flagstack, FLAGTYPE, nflag)
              fs = flagstack + nflag0;
              memset(fs, 0, (size_t)FLAG_BUFSIZE*sizeof(FLAGTYPE));
              QREALLOC(nflagstack, int, nflag)
              nfs = nflagstack + nflag0;
              memset(nfs, 0, (size_t)FLAG_BUFSIZE*sizeof(int));
              }
            if (!*nfs)
              {
              *fs = cimaflags;
              *nfs = 1;
              }
            }

/*------ Explore the list we have built and search for most frequent flags */
        for (n=nflag, fs=flagstack, nfs=nflagstack; n-- && *nfs; fs++)
          if ((nc=*(nfs++))>nmax)
            {
            nmax = nc;
            imaflags = *fs;
            }

/*------ Free memory allocated for the buffers */
        free(flagstack);
        free(nflagstack);
        break;
      default:
        error(EXIT_FAILURE, "*Internal Error*: Unknown FLAG_TYPE","");
      }

    if (i<prefs.nfimage)
      {
      obj->imaflags[i] = imaflags;
      obj->imanflags[i] = nmax;
      }
    }

  return;
  }


/****** flag_merge **********************************************************
PROTO   void	flag_merge(objstruct *objmaster, objstruct *objslave)
PURPOSE Composite flag extracted from the flag-maps.
INPUT   Destination object,
	source object.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 21/02/2012
 ***/
void	flag_merge(objstruct *objmaster, objstruct *objslave)
  {
   FLAGTYPE	*masterflag,*slaveflag;
   int		i, *masternflag,*slavenflag;

  masterflag = objmaster->imaflags;
  masternflag = objmaster->imanflags;
  slaveflag = objslave->imaflags;
  slavenflag = objslave->imanflags;
  for (i=0; i<prefs.nfimage; i++,
	masterflag++,masternflag++,slaveflag++,slavenflag++)
    switch(prefs.flag_type[i])
      {
      case FLAG_OR:
        *masterflag |= *slaveflag;
        *masternflag += *slavenflag;
        break;
      case FLAG_AND:
        *masterflag &= *slaveflag;
        *masternflag += *slavenflag;
        break;
      case FLAG_MIN:
        if (*slaveflag == *masterflag)
          *masternflag += *slavenflag;
        else if (*slaveflag<*masterflag)
          {
          *masterflag = *slaveflag;
          *masternflag = *slavenflag;
          }
        break;
      case FLAG_MAX:
        if (*slaveflag == *masterflag)
          *masternflag += *slavenflag;
        else if (*slaveflag>*masterflag)
          {
          *masterflag = *slaveflag;
          *masternflag = *slavenflag;
          }
        break;
      case FLAG_MOST:
        if (*slavenflag>*masternflag)
          {
          *masterflag = *slaveflag;
          *masternflag = *slavenflag;
          }
        break;
      default:
        error(EXIT_FAILURE, "*Internal Error*: Unknown FLAG_TYPE","");
      }

  return;
  }

