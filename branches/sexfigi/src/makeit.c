/*
 				makeit.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN, IAP & Leiden observatory
*
*	Contents:	main program.
*
*	Last modify:	14/07/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include	"fits/fitscat.h"
#include	"assoc.h"
#include	"back.h"
#include	"check.h"
#include	"field.h"
#include	"filter.h"
#include	"growth.h"
#include	"interpolate.h"
#include	"psf.h"
#include	"som.h"
#include	"weight.h"
#include	"xml.h"

time_t	thetimet, thetimet2;

/******************************** makeit *************************************/
/*
Manage the whole stuff.
*/
void	makeit()

  {
   checkstruct		*check;
   picstruct		*dfield, *field,*pffield[MAXFLAG], *wfield,*dwfield;
   catstruct		*imacat;
   tabstruct		*imatab;
   static time_t        thetime1, thetime2;
   struct tm		*tm;
   int			i, nok, ntab, next;

/* Install error logging */
  error_installfunc(write_error);

/* Processing start date and time */
  thetimet = time(NULL);
  tm = localtime(&thetimet);
  sprintf(prefs.sdate_start,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_start,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT, "----- %s %s started on %s at %s with %d thread%s\n\n",
		BANNER,
		MYVERSION,
		prefs.sdate_start,
		prefs.stime_start,
		prefs.nthreads,
		prefs.nthreads>1? "s":"");

/* Initialize globals variables */
  initglob();

  NFPRINTF(OUTPUT, "Setting catalog parameters");
  readcatparams(prefs.param_name);
  useprefs();			/* update things accor. to prefs parameters */

  if (prefs.psf_flag)
    {
    NFPRINTF(OUTPUT, "Reading PSF information");
    thepsf = psf_load(prefs.psf_name[0]); 
    if (prefs.dpsf_flag)
      ppsf = psf_load(prefs.psf_name[1]);
 /*-- Need to check things up because of PSF context parameters */
    updateparamflags();
    useprefs();
    }

  if (prefs.filter_flag)
    {
    NFPRINTF(OUTPUT, "Reading detection filter");
    getfilter(prefs.filter_name);	/* get the detection filter */
    }

  if (FLAG(obj2.sprob))
    {
    NFPRINTF(OUTPUT, "Initializing Neural Network");
    neurinit();
    NFPRINTF(OUTPUT, "Reading Neural Network Weights");
    getnnw(); 
    }

  if (prefs.somfit_flag)
    {
     int	margin;

    thesom = som_load(prefs.som_name);
    if ((margin=(thesom->inputsize[1]+1)/2) > prefs.cleanmargin)
      prefs.cleanmargin = margin;
    if (prefs.somfit_vectorsize>thesom->neurdim)
      {
      prefs.somfit_vectorsize = thesom->neurdim;
      sprintf(gstr,"%d", prefs.somfit_vectorsize);
      warning("Dimensionality of the SOM-fit vector limited to ", gstr);
      }
    }

/* Prepare growth-curve buffer */
  if (prefs.growth_flag)
    initgrowth();

/* Compute the number of valid input extensions */
  if (!(imacat = read_cat(prefs.image_name[0])))
    error(EXIT_FAILURE, "*Error*: cannot open ", prefs.image_name[0]);
  close_cat(imacat);
  imatab = imacat->tab;
  next = 0;
  for (ntab = 0 ; ntab<imacat->ntab; ntab++, imatab = imatab->nexttab)
    {
/*--  Check for the next valid image extension */
    if ((imatab->naxis < 2)
	|| !strncmp(imatab->xtension, "BINTABLE", 8)
	|| !strncmp(imatab->xtension, "ASCTABLE", 8))
      continue;
    next++;
    }
  thecat.next = next;

/*-- Init the CHECK-images */
  if (prefs.check_flag)
    {
     checkenum	c;

    NFPRINTF(OUTPUT, "Initializing check-image(s)");
    for (i=0; i<prefs.ncheck_type; i++)
    if ((c=prefs.check_type[i]) != CHECK_NONE)
      {
      if (prefs.check[c])
         error(EXIT_FAILURE,"*Error*: 2 CHECK_IMAGEs cannot have the same ",
			" CHECK_IMAGE_TYPE");
      prefs.check[c] = initcheck(prefs.check_name[i], prefs.check_type[i],
			next);
      free(prefs.check_name[i]);
      }
    }

  NFPRINTF(OUTPUT, "Initializing catalog");
  initcat();

/* Initialize XML data */
  if (prefs.xml_flag || prefs.cat_type==ASCII_VO)
    init_xml(next);

/* Go through all images */
  nok = -1;
  for (ntab = 0 ; ntab<imacat->ntab; ntab++, imatab = imatab->nexttab)
    {
/*--  Check for the next valid image extension */
    if ((imatab->naxis < 2)
	|| !strncmp(imatab->xtension, "BINTABLE", 8)
	|| !strncmp(imatab->xtension, "ASCTABLE", 8))
      continue;
    nok++;

/*-- Initial time measurement*/
    time(&thetime1);
    thecat.currext = nok+1;

    dfield = field = wfield = dwfield = NULL;

    if (prefs.dimage_flag)
      {
/*---- Init the Detection and Measurement-images */
      dfield = newfield(prefs.image_name[0], DETECT_FIELD, nok);
      field = newfield(prefs.image_name[1], MEASURE_FIELD, nok);
      if ((field->width!=dfield->width) || (field->height!=dfield->height))
        error(EXIT_FAILURE, "*Error*: Frames have different sizes","");
/*---- Prepare interpolation */
      if (prefs.dweight_flag && prefs.interp_type[0] == INTERP_ALL)
        init_interpolate(dfield, -1, -1);
      if (prefs.interp_type[1] == INTERP_ALL)
        init_interpolate(field, -1, -1);
      }
    else
      {
      field = newfield(prefs.image_name[0], DETECT_FIELD | MEASURE_FIELD, nok);
/*-- Prepare interpolation */
      if ((prefs.dweight_flag || prefs.weight_flag)
	&& prefs.interp_type[0] == INTERP_ALL)
      init_interpolate(field, -1, -1);       /* 0.0 or anything else */
      }

/*-- Init the WEIGHT-images */
    if (prefs.dweight_flag || prefs.weight_flag) 
      {
       weightenum	wtype;
       PIXTYPE	interpthresh;

      if (prefs.nweight_type>1)
        {
/*------ Double-weight-map mode */
        if (prefs.weight_type[1] != WEIGHT_NONE)
          {
/*-------- First: the "measurement" weights */
          wfield = newweight(prefs.wimage_name[1],field,prefs.weight_type[1],
		nok);
          wtype = prefs.weight_type[1];
          interpthresh = prefs.weight_thresh[1];
/*-------- Convert the interpolation threshold to variance units */
          weight_to_var(wfield, &interpthresh, 1);
          wfield->weight_thresh = interpthresh;
          if (prefs.interp_type[1] != INTERP_NONE)
            init_interpolate(wfield,
		prefs.interp_xtimeout[1], prefs.interp_ytimeout[1]);
          }
/*------ The "detection" weights */
        if (prefs.weight_type[0] != WEIGHT_NONE)
          {
          interpthresh = prefs.weight_thresh[0];
          if (prefs.weight_type[0] == WEIGHT_FROMINTERP)
            {
            dwfield=newweight(prefs.wimage_name[0],wfield,prefs.weight_type[0],
		nok);
            weight_to_var(wfield, &interpthresh, 1);
            }
          else
            {
            dwfield = newweight(prefs.wimage_name[0], dfield?dfield:field,
		prefs.weight_type[0], nok);
            weight_to_var(dwfield, &interpthresh, 1);
            }
          dwfield->weight_thresh = interpthresh;
          if (prefs.interp_type[0] != INTERP_NONE)
            init_interpolate(dwfield,
		prefs.interp_xtimeout[0], prefs.interp_ytimeout[0]);
          }
        }
      else
        {
/*------ Single-weight-map mode */
        wfield = newweight(prefs.wimage_name[0], dfield?dfield:field,
			prefs.weight_type[0], nok);
        wtype = prefs.weight_type[0];
        interpthresh = prefs.weight_thresh[0];
/*------ Convert the interpolation threshold to variance units */
        weight_to_var(wfield, &interpthresh, 1);
        wfield->weight_thresh = interpthresh;
        if (prefs.interp_type[0] != INTERP_NONE)
          init_interpolate(wfield,
		prefs.interp_xtimeout[0], prefs.interp_ytimeout[0]);
        }
      }

/*-- Init the FLAG-images */
    for (i=0; i<prefs.nimaflag; i++)
      {
      pffield[i] = newfield(prefs.fimage_name[i], FLAG_FIELD, nok);
      if ((pffield[i]->width!=field->width)
	|| (pffield[i]->height!=field->height))
        error(EXIT_FAILURE,
	"*Error*: Incompatible FLAG-map size in ", prefs.fimage_name[i]);
      }

/*-- Compute background maps for `standard' fields */
    QPRINTF(OUTPUT, dfield? "Measurement image:"
			: "Detection+Measurement image: ");
    makeback(field, wfield);
    QPRINTF(OUTPUT, (dfield || (dwfield&&dwfield->flags^INTERP_FIELD))? "(M)   "
		"Background: %-10g RMS: %-10g / Threshold: %-10g \n"
		: "(M+D) "
		"Background: %-10g RMS: %-10g / Threshold: %-10g \n",
	field->backmean, field->backsig, (field->flags & DETECT_FIELD)?
	field->dthresh: field->thresh);
    if (dfield)
    {
      QPRINTF(OUTPUT, "Detection image: ");
      makeback(dfield, dwfield? dwfield
			: (prefs.weight_type[0] == WEIGHT_NONE?NULL:wfield));
      QPRINTF(OUTPUT, "(D)   "
		"Background: %-10g RMS: %-10g / Threshold: %-10g \n",
	dfield->backmean, dfield->backsig, dfield->dthresh);
      }
    else if (dwfield && dwfield->flags^INTERP_FIELD)
      {
      makeback(field, dwfield);
      QPRINTF(OUTPUT, "(D)   "
		"Background: %-10g RMS: %-10g / Threshold: %-10g \n",
	field->backmean, field->backsig, field->dthresh);
      }

/*-- For interpolated weight-maps, copy the background structure */
    if (dwfield && dwfield->flags&(INTERP_FIELD|BACKRMS_FIELD))
      copyback(dwfield->reffield, dwfield);
    if (wfield && wfield->flags&(INTERP_FIELD|BACKRMS_FIELD))
      copyback(wfield->reffield, wfield);

/*-- Prepare learn and/or associations */
    if (prefs.assoc_flag)
      init_assoc(field);                  /* initialize assoc tasks */

/*-- Update the CHECK-images */
    if (prefs.check_flag)
      for (i=0; i<MAXCHECK; i++)
        if ((check=prefs.check[i]))
          reinitcheck(field, check);

/*-- Initialize PSF contexts and workspace */
    if (prefs.psf_flag)
      {
      psf_readcontext(thepsf, field);
      psf_init(thepsf);
      if (prefs.dpsf_flag)
        {
        psf_readcontext(thepsf, dfield);
        psf_init(thepsf); /*?*/
        }
      }

/*-- Copy field structures to static ones (for catalog info) */
    if (dfield)
      {
      thefield1 = *field;
      thefield2 = *dfield;
      }
    else
      thefield1 = thefield2 = *field;

    if (wfield)
      {
      thewfield1 = *wfield;
      thewfield2 = dwfield? *dwfield: *wfield;
      }
    else if (dwfield)
      thewfield2 = *dwfield;

    reinitcat(field);

/*-- Start the extraction pipeline */
    NFPRINTF(OUTPUT, "Scanning image");
    scanimage(field, dfield, pffield, prefs.nimaflag, wfield, dwfield);

/*-- Finish the current CHECK-image processing */
    if (prefs.check_flag)
      for (i=0; i<MAXCHECK; i++)
        if ((check=prefs.check[i]))
          reendcheck(field, check);

/*-- Final time measurements*/
    if (time(&thetime2)!=-1)
      {
      if (!strftime(thecat.ext_date, 12, "%d/%m/%Y", localtime(&thetime2)))
        error(EXIT_FAILURE, "*Internal Error*: Date string too long ","");
      if (!strftime(thecat.ext_time, 10, "%H:%M:%S", localtime(&thetime2)))
        error(EXIT_FAILURE, "*Internal Error*: Time/date string too long ","");
      thecat.ext_elapsed = difftime(thetime2, thetime1);
      }

    reendcat();

/* Update XML data */
  if (prefs.xml_flag || prefs.cat_type==ASCII_VO)
    update_xml(&thecat, dfield? dfield:field, field,
	dwfield? dwfield:wfield, wfield);


/*-- Close ASSOC routines */
    end_assoc(field);

    for (i=0; i<prefs.nimaflag; i++)
      endfield(pffield[i]);
    endfield(field);
    if (dfield)
      endfield(dfield);
    if (wfield)
      endfield(wfield);
    if (dwfield)
      endfield(dwfield);

    QPRINTF(OUTPUT, "Objects: detected %-8d / sextracted %-8d               \n",
	thecat.ndetect, thecat.ntotal);
    }

  if (nok<0)
    error(EXIT_FAILURE, "Not enough valid FITS image extensions in ",
	prefs.image_name[0]);
  free_cat(&imacat, 1);

  NFPRINTF(OUTPUT, "Closing files");

/* End CHECK-image processing */
  if (prefs.check_flag)
    for (i=0; i<MAXCHECK; i++)
      {
      if ((check=prefs.check[i]))
        endcheck(check);
      prefs.check[i] = NULL;
      }

  if (prefs.filter_flag)
    endfilter();

  if (prefs.somfit_flag)
    som_end(thesom);

  if (prefs.growth_flag)
    endgrowth();

  if (prefs.psf_flag)
    psf_end(thepsf,thepsfit); /*?*/

  if (prefs.dpsf_flag)
    psf_end(ppsf,ppsfit);

  if (FLAG(obj2.sprob))
    neurclose();

/* Processing end date and time */
  thetimet2 = time(NULL);
  tm = localtime(&thetimet2);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);
  prefs.time_diff = difftime(thetimet2, thetimet);

/* Write XML */
  if (prefs.xml_flag)
    write_xml(prefs.xml_name);

  endcat((char *)NULL);

  if (prefs.xml_flag || prefs.cat_type==ASCII_VO)
    end_xml();

  return;
  }


/******************************** initglob ***********************************/
/*
Initialize a few global variables
*/
void	initglob()
  {
   int	i;

  for (i=0; i<37; i++)
    {
    ctg[i] = cos(i*PI/18);
    stg[i] = sin(i*PI/18);
    }


  return;
  }


/****** write_error ********************************************************
PROTO	int	write_error(char *msg1, char *msg2)
PURPOSE	Manage files in case of a catched error
INPUT	a character string,
	another character string
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	14/07/2006
 ***/
void	write_error(char *msg1, char *msg2)
  {
   char			error[MAXCHAR];

  sprintf(error, "%s%s", msg1,msg2);
  if (prefs.xml_flag)
    write_xmlerror(prefs.xml_name, error);

/* Also close existing catalog */
  endcat(error);

  end_xml();

  return;
  }

