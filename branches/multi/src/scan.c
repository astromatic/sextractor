/*
*				scan.c
*
* Main image scanning routines.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		02/08/2012
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
#include	"prefs.h"
#include	"analyse.h"
#include	"back.h"
#include	"check.h"
#include	"clean.h"
#include	"deblend.h"
#include	"filter.h"
#include	"image.h"
#include	"lutz.h"
#include	"plist.h"
#include	"readimage.h"
#include	"scan.h"
#include	"weight.h"

#ifdef USE_THREADS
#include	"threads.h"
#endif

static void	scan_initmarkers(fieldstruct *field),
		scan_updatemarkers(fieldstruct *field, int yl);

/****** scan_extract *********************************************************
PROTO	void scan_extract(fieldstruct *dfield, fieldstruct *dwfield,
			fieldstruct **fields, fieldstruct **wfields, int nfield,
			fieldstruct **ffields, int nffield)
PURPOSE	Scan of the large pixmap(s). Main loop and heart of the program.
INPUT	Pointer to the detection image field,
	pointer to the detection weight field,
	pointer to an array of image field pointers,
	pointer to an array of weight-map field pointers,
	number of images,
	pointer to an array of flag map field pointers,
	number of flag maps.
OUTPUT	-.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	02/08/2012
 ***/
void	scan_extract(fieldstruct *dfield, fieldstruct *dwfield,
			fieldstruct **fields, fieldstruct **wfields, int nfield,
			fieldstruct **ffields, int nffield)

  {
   static infostruct	curpixinfo, *info, *store,
			initinfo, freeinfo, *victim;
   fieldstruct		*field,*ffield;
   checkstruct		*check;
   objliststruct       	objlist;
   objstruct		*cleanobj;
   pliststruct		*pixel, *pixt;

   char			*marker, newmarker, *blankpad, *bpt,*bpt0;
   int			co, i,j, flag, luflag,pstop, xl,xl2,yl, cn,
			nposize, stacksize, w, h, blankh, maxpixnb,
			varthreshflag, ontotal;
   short	       	trunflag;
   PIXTYPE		thresh, relthresh, cnewsymbol,cwthresh,wthresh,
			*scan,*cscan, *wscan,*wscanp,*wscann,
			*cwscan,*cwscanp,*cwscann, *wscand,*scant;
   FLAGTYPE		*fscan[MAXFLAG];
   status		cs, ps, *psstack;
   int			*start, *end, ymax;

/* Avoid gcc -Wall warnings */
  scan = cscan = cwscan = cwscann = cwscanp = wscan = wscann = wscanp = NULL;
  victim = NULL;			/* Avoid gcc -Wall warnings */
  blankh = 0;				/* Avoid gcc -Wall warnings */
/*----- Beginning of the main loop: Initialisations  */
  thecat.ntotal = thecat.ndetect = 0;
  
  wthresh = dwfield ? dwfield->weight_thresh : 0.0;
  if (wthresh>BIG*WTHRESH_CONVFAC);
    wthresh = BIG*WTHRESH_CONVFAC;

/* If WEIGHTing and no absolute thresholding, activate threshold scaling */
  varthreshflag = (dwfield && prefs.thresh_type[0]!=THRESH_ABSOLUTE);
  relthresh = varthreshflag ? prefs.dthresh[0] : 0.0;/* To avoid gcc warnings*/
  w = dfield->width;
  h = dfield->height;
  objlist.dthresh = dfield->dthresh;
  objlist.thresh = dfield->thresh;
  scan_initmarkers(dfield);
  scan_initmarkers(dwfield);

  for (i=1; i<nfield; i++)
    if (!(fields[i]->flags&MULTIGRID_FIELD))
      {
      scan_initmarkers(fields[i]);
      scan_initmarkers(wfields[i]);
      }

  if (nffield)
    for (i=0; i<nffield; i++)
      scan_initmarkers(ffields[i]);

/*Allocate memory for buffers */
  stacksize = w+1;
  QMALLOC(info, infostruct, stacksize);
  QCALLOC(store, infostruct, stacksize);
  QMALLOC(marker, char, stacksize);
  QMALLOC(dumscan, PIXTYPE, stacksize);
  QMALLOC(psstack, status, stacksize);
  QCALLOC(start, int, stacksize);
  QMALLOC(end, int, stacksize);
  blankpad = bpt = NULL;
  lutz_alloc(w,h);
  deblend_alloc();

/* Some initializations */

  thresh = objlist.dthresh;
  initinfo.pixnb = 0;
  initinfo.flag = 0;
  initinfo.firstpix = initinfo.lastpix = -1;

  for (xl=0; xl<stacksize; xl++)
    {
    marker[xl]  = 0 ;
    dumscan[xl] = -BIG ;
    }

  co = pstop = 0;
  objlist.nobj = 1;
  curpixinfo.pixnb = 1;

/* Init cleaning procedure */
  clean_init();

/*----- Allocate memory for the pixel list */
  init_plist();

  if (!(pixel = objlist.plist = malloc(nposize=prefs.mem_pixstack*plistsize)))
    error(EXIT_FAILURE, "Not enough memory to store the pixel stack:\n",
        "           Try to decrease MEMORY_PIXSTACK");

/*----- at the beginning, "free" object fills the whole pixel list */
  freeinfo.firstpix = 0;
  freeinfo.lastpix = nposize-plistsize;
  pixt = pixel;
  for (i=plistsize; i<nposize; i += plistsize, pixt += plistsize)
    PLIST(pixt, nextpix) = i;
  PLIST(pixt, nextpix) = -1;

/* Allocate memory for other buffers */
  if (prefs.filter_flag)
    {
    QMALLOC(cscan, PIXTYPE, stacksize);
    if (dwfield)
      {
      QCALLOC(cwscan, PIXTYPE, stacksize);
      if (PLISTEXIST(wflag))
        {
        QCALLOC(cwscanp, PIXTYPE, stacksize);
        QCALLOC(cwscann, PIXTYPE, stacksize);
        }
      }
/*-- One needs a buffer to protect filtering if source-blanking applies */
    if (prefs.blank_flag)
      {
      blankh = thefilter->convh/2+1;
      QMALLOC(blankpad, char, w*blankh);
      dfield->yblank -= blankh;
      bpt = blankpad;
      }
    }

#ifdef USE_THREADS
/*Setup measurement threads as we meet the 1st object; leave 1 for extraction */
  if (prefs.nthreads>1)
    pthread_init_obj2(fields, wfields, nfield, prefs.nthreads);
#endif

/*----- Here we go */
  for (yl=0; yl<=h;)
    {
    ps = COMPLETE;
    cs = NONOBJECT;
    if (yl==h)
      {
/*---- Need an empty line for Lutz' algorithm to end gracely */
      if (prefs.filter_flag)
        {
        free(cscan);
        if (dwfield)
          {
          if (PLISTEXIST(wflag))
            {
            free(cwscanp);
            free(cwscann);
            cwscanp = cwscan;
            }
          else
            free(cwscan);
          }
        }
      cwscan = cwscann = cscan = dumscan;
      }
    else
      {
      if (dwfield)
        {
/*------ Copy the previous weight line to track bad pixel limits */
        wscan = (dwfield->stripy==dwfield->stripysclim)?
		  (PIXTYPE *)readimage_loadstrip(dwfield, (fieldstruct *)NULL)
		: &dwfield->strip[dwfield->stripy*dwfield->width];
        if (PLISTEXIST(wflag))
          {
          if (yl>0)
            wscanp = &dwfield->strip[((yl-1)%dwfield->stripheight)
			*dwfield->width];
          if (yl<h-1)
            wscann = &dwfield->strip[((yl+1)%dwfield->stripheight)
			*dwfield->width];
            }
        }
      if (dfield->stripy==dfield->stripysclim)
        {
        scan = (PIXTYPE *)readimage_loadstrip(dfield, dwfield);
        for (i=1; i<nfield; i++)
          if (!(fields[i]->flags&MULTIGRID_FIELD))
            {
            readimage_loadstrip(fields[i], wfields[i]);
            if (wfields[i])
              readimage_loadstrip(wfields[i], (fieldstruct *)NULL);
            }
        }
      else
        scan = &dfield->strip[dfield->stripy*dfield->width];

      if (nffield)
        for (i=0; i<nffield; i++)
          {
          ffield = ffields[i];
          fscan[i] = (ffield->stripy==ffield->stripysclim)?
		  (FLAGTYPE *)readimage_loadstrip(ffield, (fieldstruct *)NULL)
		: &ffield->fstrip[ffield->stripy*ffield->width];
          }

      if (prefs.filter_flag)
        {
        filter(dfield, cscan, dfield->y);
        if (dwfield)
          {
          if (PLISTEXIST(wflag))
            {
            if (yl==0)
              filter(dwfield, cwscann, yl);
            wscand = cwscanp;
            cwscanp = cwscan;
            cwscan = cwscann;
            cwscann = wscand;
            if (yl < h-1)
              filter(dwfield, cwscann, yl + 1);
            }
          else
            filter(dwfield, cwscan, yl);
          }
        }
      else
        {
        cscan = scan;
        cwscan = wscan;
        if (PLISTEXIST(wflag))
          {
          cwscanp = wscanp;
          cwscann = wscann;
          }
        }

      if ((check=prefs.check[CHECK_FILTERED]))
        check_write(check, cscan, w);
      }

    trunflag = (yl==0 || yl==h-1)? OBJ_TRUNC:0;

    for (xl=0; xl<=w; xl++)
      {
      if (xl == w)
        cnewsymbol = -BIG;
      else
        cnewsymbol = cscan[xl];

      newmarker = marker[xl];
      marker[xl] = 0;

      curpixinfo.flag = trunflag;
      if (varthreshflag)
        thresh = relthresh*sqrt((xl==w || yl==h)? 0.0:cwscan[xl]);
      luflag = cnewsymbol > thresh?1:0;

      if (luflag)
        {
        if (xl==0 || xl==w-1)
          curpixinfo.flag |= OBJ_TRUNC;
        pixt = pixel + (cn=freeinfo.firstpix);
        freeinfo.firstpix = PLIST(pixt, nextpix);

/*------- Running out of pixels, the largest object becomes a "victim" ------*/

        if (freeinfo.firstpix==freeinfo.lastpix)
          {
          sprintf(gstr, "%d,%d", xl+1, yl+1);
          warning("Pixel stack overflow at position ", gstr);
          maxpixnb = 0;
          for (i=0; i<=w; i++)
            if (store[i].pixnb>maxpixnb)
              if (marker[i]=='S' || (newmarker=='S' && i==xl))
                {
                flag = 0;
                if (i<xl)
                  for (j=0; j<=co; j++)
                    flag |= (start[j]==i);
                if (!flag)
                  maxpixnb = (victim = &store[i])->pixnb;
                }
          for (j=1; j<=co; j++)
            if (info[j].pixnb>maxpixnb)
              maxpixnb = (victim = &info[j])->pixnb;

          if (!maxpixnb)
            error(EXIT_FAILURE, "*Fatal Error*: something is badly bugged in ",
		"scanimage()!");
          if (maxpixnb <= 1)
            error(EXIT_FAILURE, "Pixel stack overflow in ", "scanimage()");
          freeinfo.firstpix = PLIST(pixel+victim->firstpix, nextpix);
          PLIST(pixel+victim->lastpix, nextpix) = freeinfo.lastpix;
          PLIST(pixel+(victim->lastpix=victim->firstpix), nextpix) = -1;
          victim->pixnb = 1;
          victim->flag |= OBJ_OVERFLOW;
          }

/*---------------------------------------------------------------------------*/
        curpixinfo.lastpix = curpixinfo.firstpix = cn;
        PLIST(pixt, nextpix) = -1;
        PLIST(pixt, x) = xl;
        PLIST(pixt, y) = yl;
        PLIST(pixt, value) = scan[xl];
        if (PLISTEXIST(cvalue))
          PLISTPIX(pixt, cvalue) = cnewsymbol;
        if (PLISTEXIST(flag))
          for (i=0; i<nffield; i++)
            PLISTFLAG(pixt, flag[i]) = fscan[i][xl];
/*--------------------- Detect pixels with a low weight ---------------------*/
        if (PLISTEXIST(wflag) && wscan)
          {
	  PLISTFLAG(pixt, wflag) = 0;
          if (wscan[xl] >= wthresh)
            PLISTFLAG(pixt, wflag) |= OBJ_LOWWEIGHT;
          if (cwscan[xl] >= cwthresh)
            PLISTFLAG(pixt, wflag) |= OBJ_LOWDWEIGHT;

          if (yl>0)
            {
            if (cwscanp[xl] >= cwthresh)
              PLISTFLAG(pixt, wflag) |= OBJ_LOWDWEIGHT;
            if (xl>0 && cwscanp[xl-1]>=cwthresh)
              PLISTFLAG(pixt, wflag) |= OBJ_LOWDWEIGHT;
            if (xl<w-1 && cwscanp[xl+1]>=cwthresh)
              PLISTFLAG(pixt, wflag) |= OBJ_LOWDWEIGHT;
            }
          if (xl>0 && cwscan[xl-1]>=cwthresh)
              PLISTFLAG(pixt, wflag) |= OBJ_LOWDWEIGHT;
          if (xl<w-1 && cwscan[xl+1]>=cwthresh)
            PLISTFLAG(pixt, wflag) |= OBJ_LOWDWEIGHT;
          if (yl<h-1)
            {
            if (cwscann[xl] >= cwthresh)
              PLISTFLAG(pixt, wflag) |= OBJ_LOWDWEIGHT;
            if (xl>0 && cwscann[xl-1]>=cwthresh)
              PLISTFLAG(pixt, wflag) |= OBJ_LOWDWEIGHT;
            if (xl<w-1 && cwscann[xl+1]>=cwthresh)
              PLISTFLAG(pixt, wflag) |= OBJ_LOWDWEIGHT;
            }
          }
        if (PLISTEXIST(dthresh))
          PLISTPIX(pixt, dthresh) = thresh;
        if (PLISTEXIST(var))
          PLISTPIX(pixt, var) = wscan[xl];

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
            else
              marker[xl] = 's';
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

/*---------------------------------------------------------------------------*/
        }

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
            lutz_update(&info[co],&store[xl], pixel);
          ps = OBJECT;
          }
        else if (newmarker == 's')
          {
          if ((cs == OBJECT) && (ps == COMPLETE))
            {
            pstop--;
            xl2 = start[co];
            lutz_update(&info[co-1],&info[co], pixel);
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
              if ((int)info[co].pixnb >= prefs.ext_minarea)
                {
                scan_output(fields, wfields, nfield, &info[co], &objlist);
                }
/* ------------------------------------ free the chain-list */

              PLIST(pixel+info[co].lastpix, nextpix) = freeinfo.firstpix;
              freeinfo.firstpix = info[co].firstpix;
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
        lutz_update(&info[co],&curpixinfo, pixel);
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

      if (prefs.blank_flag && xl<w)
        {
        if (prefs.filter_flag)
	  *(bpt++) = (luflag)?1:0;
        else if (luflag)
          scan[xl] = -BIG;
        }
/*--------------------- End of the loop over the x's -----------------------*/
      }

/* Detected pixel removal at the end of each line */
    if (prefs.blank_flag && yl<h)
      {
      if (prefs.filter_flag)
        {
        bpt = bpt0 = blankpad + w*((yl+1)%blankh);
        if (dfield->yblank >= 0)
          {
          scant = &PIX(dfield, 0, dfield->yblank);
          for (i=w; i--; scant++)
            if (*(bpt++))
              *scant = -BIG;
          bpt = bpt0;
          }
        }
      dfield->yblank++;
      }

/*-- Prepare markers for the next line */
    yl++;
    scan_updatemarkers(dfield, yl);
    scan_updatemarkers(dwfield, yl);
    for (i=1; i<nfield; i++)
      if (!(fields[i]->flags&MULTIGRID_FIELD))
        {
        scan_updatemarkers(fields[i], yl);
        scan_updatemarkers(wfields[i], yl);
        }
    if (nffield)
      for (i=0; i<nffield; i++)
        scan_updatemarkers(ffields[i], yl);

/*-- Remove objects close to the ymin limit if ymin is ready to increase */
    if (dfield->stripy==dfield->stripysclim)
      {
      i = cleanobjlist->nobj;
      while (i--)
        {
        if (i>=cleanobjlist->nobj)
          i = cleanobjlist->nobj - 1;
        cleanobj = cleanobjlist->obj+i;
        if (cleanobj->ycmin <= dfield->ymin)
          analyse_final(fields, wfields, nfield, cleanobjlist, i);
        }
      }

#ifdef USE_THREADS
    if (prefs.nthreads>1)
      {
      QPTHREAD_MUTEX_LOCK(&pthread_countobj2mutex);
      thecat.nline = yl>h? h:yl;
      if ((prefs.prof_flag && !(thecat.ntotal%(10*prefs.nthreads))) ||
		!(thecat.nline%50))
        NPRINTF(OUTPUT, "\33[1M> Line:%5d  "
		"Objects: %8d detected / %8d sextracted\n\33[1A",
		thecat.nline, thecat.ndetect, thecat.ntotal);
      QPTHREAD_MUTEX_UNLOCK(&pthread_countobj2mutex);
      }
    else
#endif
      {
      thecat.nline = yl>h? h:yl;
      if ((prefs.prof_flag && !(thecat.ntotal%10)) || !(thecat.nline%50))
        NPRINTF(OUTPUT, "\33[1M> Line:%5d  "
		"Objects: %8d detected / %8d sextracted\n\33[1A",
		thecat.nline, thecat.ndetect, thecat.ntotal);
      }

/*--------------------- End of the loop over the y's -----------------------*/
    }

/* Removal or the remaining pixels */
  if (prefs.blank_flag && prefs.filter_flag && (dfield->yblank >= 0))
    for (j=blankh-1; j--; yl++)
      {
      bpt = bpt0 = blankpad + w*(yl%blankh);
      scant = &PIX(dfield, 0, dfield->yblank);
      for (i=w; i--; scant++)
        if (*(bpt++))
          *scant = -BIG;
      dfield->yblank++;
      }

/* Now that all "detected" pixels have been removed, analyse detections */
  while (cleanobjlist->nobj)
    analyse_final(fields, wfields, nfield, cleanobjlist,
	cleanobjlist->nobj-1);

  clean_end();

/*Free memory */
  if (prefs.filter_flag && dwfield && PLISTEXIST(wflag))
    free(cwscanp);
  deblend_free();
  free(pixel);
  lutz_free();
  free(info);
  free(store);
  free(marker);
  free(dumscan);
  free(psstack);
  free(start);
  free(end);
  if (prefs.blank_flag && prefs.filter_flag)
    free(blankpad);

#ifdef	USE_THREADS
  if (prefs.nthreads>1)
    pthread_end_obj2();
#endif

  return;
  }


/*i**** scan_initmarkers **************************************************
PROTO	static void scan_initmarkers(fieldstruct *field)
PURPOSE	Initialize scan markers of an image field for the first line.
INPUT	Pointer to the image field.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/05/2012
 ***/
static void	scan_initmarkers(fieldstruct *field)

  {
  if (field)
    field->y = field->stripy = field->ymin = field->stripylim
	= field->stripysclim = 0;

  return;
  }


/*i**** scan_updatemarkers **************************************************
PROTO	static void scan_updatemarkers(fieldstruct *field, int yl)
PURPOSE	Prepare scan markers of an image field for the next line.
INPUT	Pointer to the image field,
	scan line y coordinate.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/05/2012
 ***/
static void	scan_updatemarkers(fieldstruct *field, int yl)

  {
  if (field)
    field->stripy = (field->y=yl)%field->stripheight;

  return;
  }


/****** scan_output **********************************************************
PROTO	void scan_output(fieldstruct **fields, fieldstructs **wfields,
		int nfield, infostruct *info, objliststruct *objlist)
PURPOSE	Manage detection after primary extraction (deblending, cleaning,
	measurements), and add it to an object list.
INPUT	Pointer to an array of image field pointers,
	pointer to an array of weight-map field pointers,
	number of images,
	pointer to detection info,
	pointer to the output object list.
OUTPUT	-.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	07/03/2012
 ***/
void	scan_output(fieldstruct **fields, fieldstruct **wfields, int nfield,
		infostruct *info, objliststruct *objlist)
  {
   fieldstruct		*field;
   objliststruct	objlistd, *objlistout;
   obj2liststruct	*obj2list;
   static objstruct	obj;
   objstruct		*cobj, *vobj;
   obj2struct		*obj2, *firstobj2, *prevobj2;
   pliststruct		*pixel;
   int 			i,j,n,o;

  field = fields[0];

  pixel = objlist->plist;
  objlistd.obj = NULL;
  objlistd.plist = NULL;
  objlistd.nobj = objlistd.npix = 0;

/*----- Allocate memory to store object data */

  objlist->obj = &obj;
  objlist->nobj = 1;

  memset(&obj, 0, (size_t)sizeof(objstruct));
  objlist->npix = info->pixnb;
  obj.firstpix = info->firstpix;
  obj.lastpix = info->lastpix;
  obj.flag = info->flag;
  obj.dthresh = objlist->dthresh;
  obj.thresh = objlist->thresh;

  scan_preanalyse(objlist, 0, ANALYSE_FAST);

/*----- Check if the current strip contains the lower isophote... */
  if ((int)obj.ymin < field->ymin)
    obj.flag |= OBJ_ISO_PB;

  if (!(obj.flag & OBJ_OVERFLOW) && (createsubmap(objlist, 0) == RETURN_OK))
    {
    if (deblend_parcelout(objlist, &objlistd) == RETURN_OK)
      objlistout = &objlistd;
    else
      {
      objlistout = objlist;
      for (o=0; o<objlistout->nobj; o++)
        objlistout->obj[o].flag |= OBJ_DOVERFLOW;
      sprintf(gstr, "%.0f,%.0f", obj.mx+1, obj.my+1);
      warning("Deblending overflow for detection at ", gstr);
      }
    free(obj.submap);
    }
  else
    objlistout = objlist;

  ++thecat.nblend;			/* Parent blend index */
  for (o=0; o<objlistout->nobj; o++)
    {
/*-- Basic measurements */
    scan_preanalyse(objlistout, o, ANALYSE_FULL|ANALYSE_ROBUST);
    if (prefs.ext_maxarea && objlistout->obj[o].fdnpix > prefs.ext_maxarea)
      continue; 
    cobj = objlistout->obj + o;
    cobj->number = ++thecat.ndetect;
    cobj->blend = thecat.nblend;
/*--- Isophotal measurements */
    analyse_iso(fields, wfields, nfield, objlistout, o);
    if (prefs.blank_flag)
      {
      if (createblank(objlistout, o) != RETURN_OK)
        {
/*------ Not enough mem. for the BLANK vignet: flag the object now */
        cobj->flag |= OBJ_OVERFLOW;
        cobj->blank = NULL;
        sprintf(gstr, "%.0f,%.0f", cobj->mx+1, cobj->my+1);
        warning("Memory overflow during masking for detection at ", gstr);
        }
      }

    if ((n=cleanobjlist->nobj) >= prefs.clean_stacksize)
      {
       objstruct	*cleanobj;
       int		ymin, ymax, victim=0;

      ymin = 2000000000;	/* No image is expected to be that tall ! */
      cleanobj = cleanobjlist->obj;
      for (j=0; j<n; j++, cleanobj++)
        if (cleanobj->ycmax < ymin)
          {
          victim = j;
          ymin = cleanobj->ycmax;
          }

      cleanobj = &cleanobjlist->obj[victim];
/*---- Warn if there is a possibility for any aperture to be truncated */
      if (field->ymax < field->height)
        {
        if ((ymax=cleanobj->ycmax) > field->ymax)
          {
          sprintf(gstr, "Object at position %.0f,%.0f ",
		cleanobj->mx+1, cleanobj->my+1);
          QWARNING(gstr, "may have some apertures truncated:\n"
		"          You might want to increase MEMORY_OBJSTACK");
          }
        else if (ymax>field->yblank && prefs.blank_flag)
          {
          sprintf(gstr, "Object at position %.0f,%.0f ",
		cleanobj->mx+1, cleanobj->my+1);
          QWARNING(gstr, "may have some unBLANKed neighbours\n"
		"          You might want to increase MEMORY_OBJSTACK");
          }
        }

      analyse_final(fields, wfields, nfield, cleanobjlist, victim);
      }

/*-- Add the object only if it is not "swallowed" by cleaning */
    if (!prefs.clean_flag || clean_process(field, cobj))
      clean_add(cobj);
    }

  free(objlistd.plist);
  free(objlistd.obj);

  return;
  }


/****** scan_preanalyse ******************************************************
PROTO   void scan_preanalyse(objliststruct *objlist, int no, int analyse_type)
PURPOSE Compute basic image parameters from the pixel-list for each detection.
INPUT   Objlist pointer,
	objlist number,
        analysis switch flag.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 15/02/2012
 ***/
void  scan_preanalyse(objliststruct *objlist, int no, int analyse_type)

  {
   objstruct	*obj = &objlist->obj[no];
   pliststruct	*pixel = objlist->plist, *pixt;
   PIXTYPE	peak, cpeak, val, cval, minthresh, thresht;
   double	thresh,thresh2, t1t2,darea,
		mx,my, mx2,my2,mxy, rv, tv,
		xm,ym, xm2,ym2,xym,
		temp,temp2, theta,pmx2,pmy2;
   int		x, y, xmin,xmax, ymin,ymax,area2, fdnpix, dnpix;
  

/*-----  initialize stacks and bounds */
  thresh = obj->dthresh;
  if (PLISTEXIST(dthresh))
    minthresh = BIG;
  else
    minthresh = 0.0;
  fdnpix = dnpix = 0;
  rv = 0.0;
  peak = cpeak = -BIG;
  ymin = xmin = 2*MAXPICSIZE;    /* to be really sure!! */
  ymax = xmax = 0;

/*-----  integrate results */
  for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
    {
    x = PLIST(pixt, x);
    y = PLIST(pixt, y);
    val=PLISTPIX(pixt, value);
    if (cpeak < (cval=PLISTPIX(pixt, cvalue)))
      cpeak = cval;
    if (PLISTEXIST(dthresh) && (thresht=PLISTPIX(pixt, dthresh))<minthresh)
      minthresh = thresht;
    if (peak < val)
      peak = val;
    rv += cval;
    if (xmin > x)
      xmin = x;
    if (xmax < x)
      xmax = x;
    if (ymin > y)
      ymin = y;
    if (ymax < y)
      ymax = y;
    fdnpix++;
    }    

  if (PLISTEXIST(dthresh))
    obj->dthresh = thresh = minthresh;

/* copy some data to "obj" structure */

  obj->fdnpix = (LONG)fdnpix;
  obj->fdflux = (float)rv;
  obj->fdpeak = cpeak;
  obj->xmin = xmin;
  obj->xmax = xmax;
  obj->ymin = ymin;
  obj->ymax = ymax;

  if (analyse_type & ANALYSE_FULL)
    {
    mx = my = tv = 0.0;
    mx2 = my2 = mxy = 0.0;
    thresh2 = (thresh + peak)/2.0;
    area2 = 0;
    for (pixt=pixel+obj->firstpix; pixt>=pixel; pixt=pixel+PLIST(pixt,nextpix))
      {
      x = PLIST(pixt,x)-xmin;	/* avoid roundoff errors on big images */
      y = PLIST(pixt,y)-ymin;	/* avoid roundoff errors on big images */
      cval = PLISTPIX(pixt, cvalue);
      tv += (val = PLISTPIX(pixt, value));
      if (val>thresh)
        dnpix++;
      if (val > thresh2)
        area2++;
      mx += cval * x;
      my += cval * y;
      mx2 += cval * x*x;
      my2 += cval * y*y;
      mxy += cval * x*y;
      }

/*----- compute object's properties */
    xm = mx / rv;			/* mean x */
    ym = my / rv;			/* mean y */

/*-- In case of blending, use previous barycenters */
    if ((analyse_type&ANALYSE_ROBUST) && (obj->flag&OBJ_MERGED))
      {
       double	xn,yn;

      xn = obj->mx-xmin;
      yn = obj->my-ymin;
      xm2 = mx2 / rv + xn*xn - 2*xm*xn;
      ym2 = my2 / rv + yn*yn - 2*ym*yn;
      xym = mxy / rv + xn*yn - xm*yn - xn*ym;
      xm = xn;
      ym = yn;
      }
    else
      {
      xm2 = mx2 / rv - xm * xm;	/* variance of x */
      ym2 = my2 / rv - ym * ym;	/* variance of y */
      xym = mxy / rv - xm * ym;	/* covariance */
      }

/* Handle fully correlated x/y (which cause a singularity...) */
    if ((temp2=xm2*ym2-xym*xym)<0.00694)
      {
      xm2 += 0.0833333;
      ym2 += 0.0833333;
      temp2 = xm2*ym2-xym*xym;
      obj->singuflag = 1;
      }
    else
      obj->singuflag = 0;

    if ((fabs(temp=xm2-ym2)) > 0.0)
      theta = atan2(2.0 * xym,temp) / 2.0;
    else
      theta = PI/4.0;

    temp = sqrt(0.25*temp*temp+xym*xym);
    pmy2 = pmx2 = 0.5*(xm2+ym2);
    pmx2+=temp;
    pmy2-=temp;

    obj->dnpix = (obj->flag & OBJ_OVERFLOW)? obj->fdnpix:(LONG)dnpix;
    obj->mx = xm+xmin;	/* add back xmin */
    obj->my = ym+ymin;	/* add back ymin */
    obj->mx2 = xm2;
    obj->my2 = ym2;
    obj->mxy = xym;
    obj->a = (float)sqrt(pmx2);
    obj->b = (float)sqrt(pmy2);
    obj->theta = theta*180.0/PI;

    obj->cxx = (float)(ym2/temp2);
    obj->cyy = (float)(xm2/temp2);
    obj->cxy = (float)(-2*xym/temp2);

    darea = (double)area2 - dnpix;
    t1t2 = thresh/thresh2;
    if (t1t2>0.0)
      {
      obj->abcor = (darea<0.0?darea:-1.0)/(2*PI*log(t1t2<1.0?t1t2:0.99)
	*obj->a*obj->b);
      if (obj->abcor>1.0)
        obj->abcor = 1.0;
      }
    else
      obj->abcor = 1.0;
    }

  return;
  }

