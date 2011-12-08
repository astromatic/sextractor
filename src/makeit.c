/*
*				makeit.c
*
* Main program.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		08/12/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
#include	"catout.h"
#include	"check.h"
#include	"fft.h"
#include	"field.h"
#include	"filter.h"
#include	"growth.h"
#include	"interpolate.h"
#include	"pattern.h"
#include	"psf.h"
#include	"profit.h"
#include	"som.h"
#include	"weight.h"
#include	"xml.h"

static int		selectext(char *filename);
time_t			thetimet, thetimet2;
extern char		profname[][32];
double			dtime;

/****** makeit ***************************************************************
PROTO	void main(void)
PURPOSE	Manage the whole processing.
INPUT	-.
OUTPUT	-.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	06/12/2011
 ***/
void	makeit(void)
  {
   profitstruct		*profit;
   checkstruct		*check;
   fieldstruct		**fields,**wfields,**ffields;
   catstruct		*imacat;
   tabstruct		*imatab;
   patternstruct	*pattern;
   psfstruct		**psfs;
   somstruct		*som;
   static time_t        thetime1, thetime2;
   struct tm		*tm;
   PIXTYPE		interpthresh;
   int			ext_flag[MAXFLAG],ext_image[MAXIMAGE],
			ext_wimage[MAXIMAGE], nparam2[2],
			i,t, nok, ntab, next, ntabmax, forcextflag,
			npat, nmodels;

/* Install error logging */
  error_installfunc(write_error);

/* Processing start date and time */
  dtime = counter_seconds();
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

  NFPRINTF(OUTPUT, "Setting catalog parameters");
  thecat.obj2list = catout_readparams(prefs.param, prefs.nparam,
					prefs.obj2_stacksize);
  useprefs();			/* update things accor. to prefs parameters */

  if (prefs.psf_flag)
    {
    NFPRINTF(OUTPUT, "Reading PSF information");
    QCALLOC(psfs, psfstruct *, prefs.nband);
    for (i=0; i<prefs.nband; i++)
      psfs[i] = psf_load(prefs.psf_name[i]);
 /*-- Need to check things up because of PSF context parameters */
    catout_updateparamflags();
    useprefs();
    }

  if (prefs.prof_flag)
    {
#ifdef USE_MODEL
    fft_init(prefs.nthreads);
/* Create profiles at full resolution */
    NFPRINTF(OUTPUT, "Preparing profile models");
    prefs.prof_modelflags = 
	 (FLAG(obj2.prof_offset_flux)? MODEL_BACK : MODEL_NONE)
	|(FLAG(obj2.prof_dirac_flux)? MODEL_DIRAC : MODEL_NONE)
	|(FLAG(obj2.prof_spheroid_flux)?
		(FLAG(obj2.prof_spheroid_sersicn)?
			MODEL_SERSIC : MODEL_DEVAUCOULEURS) : MODEL_NONE)
	|(FLAG(obj2.prof_disk_flux)? MODEL_EXPONENTIAL : MODEL_NONE)
	|(FLAG(obj2.prof_bar_flux)? MODEL_BAR : MODEL_NONE)
	|(FLAG(obj2.prof_arms_flux)? MODEL_ARMS : MODEL_NONE);

/*-- Setup a minimum profit structure for component names and number of params*/
    QCALLOC(profit, profitstruct, 1);
    QMALLOC(profit->prof, profstruct *, MODEL_NMAX);
    nmodels = 0;
    QPRINTF(OUTPUT, "Fitting model: ");
    for (t=1; t<(1<<MODEL_NMAX); t<<=1)
      if (prefs.prof_modelflags&t)
        {
        profit->prof[nmodels] = prof_init(profit, t);
        if (nmodels)
          QPRINTF(OUTPUT, " + ");
        QPRINTF(OUTPUT, "%s", profit->prof[nmodels]->name);
        nmodels++;
        }
    QPRINTF(OUTPUT, "\n");
    catout_changeparamsize("VECTOR_MODEL", &profit->nparam, 1);
    catout_changeparamsize("VECTOR_MODELERR", &profit->nparam, 1);
    nparam2[0] = nparam2[1] = profit->nparam;
    catout_changeparamsize("MATRIX_MODELERR", nparam2, 2);
    if (prefs.pattern_flag)
      {
      npat = prefs.prof_disk_patternvectorsize;
      if (npat<prefs.prof_disk_patternmodvectorsize)
        npat = prefs.prof_disk_patternmodvectorsize;
      if (npat<prefs.prof_disk_patternargvectorsize)
        npat = prefs.prof_disk_patternargvectorsize;
/*---- Do a copy of the original number of pattern components */
      prefs.prof_disk_patternncomp = npat;
      pattern = pattern_init(profit, prefs.pattern_type, npat);
      if (FLAG(obj2.prof_disk_patternvector))
        {
        npat = pattern->size[2];
        catout_changeparamsize("DISK_PATTERN_VECTOR", &npat, 1);
        }
      if (FLAG(obj2.prof_disk_patternmodvector))
        {
        npat = pattern->ncomp*pattern->nfreq;
        catout_changeparamsize("DISK_PATTERNMOD_VECTOR", &npat, 1);
        }
      if (FLAG(obj2.prof_disk_patternargvector))
        {
        npat = pattern->ncomp*pattern->nfreq;
        catout_changeparamsize("DISK_PATTERNARG_VECTOR", &npat, 1);
        }
      pattern_end(pattern);
      }
    profit_end(profit);
#else
    error(EXIT_FAILURE,
		"*Error*: model-fitting is not supported in this build.\n",
			" Please check your configure options");
#endif
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

    som = som_load(prefs.som_name);
    if ((margin=(som->inputsize[1]+1)/2) > prefs.cleanmargin)
      prefs.cleanmargin = margin;
    if (prefs.somfit_vector_size[0]>som->neurdim)
      {
      prefs.somfit_vector_size[0] = som->neurdim;
      sprintf(gstr,"%d", prefs.somfit_vector_size[0]);
      warning("Dimensionality of the SOM-fit vector limited to ", gstr);
      }
    }

/* Allocate memory for multidimensional catalog parameter arrays */
  catout_allocparams(thecat.obj2list);
  useprefs();

/* Read extension numbers and remove them from the filenames if present */
 for (i=0; i<prefs.nimage_name; i++)
   ext_image[i] = selectext(prefs.image_name[i]);
 for (i=0; i<prefs.nwimage_name; i++)
   ext_wimage[i] = selectext(prefs.wimage_name[i]);
 for (i=0; i<prefs.nfimage_name; i++)
   ext_flag[i] = selectext(prefs.fimage_name[i]);

/* Use the first image to probe the number of extensions to analyse */
  if (ext_image[0] != RETURN_ERROR)
    {
    forcextflag = 1;
    ntabmax = next = 1;
    }
  else
    forcextflag = 0;

  if (!(imacat = read_cat(prefs.image_name[0])))
    error(EXIT_FAILURE, "*Error*: cannot open ", prefs.image_name[0]);
  close_cat(imacat);
  imatab = imacat->tab;
  if (!forcextflag)
    {
    ntabmax = imacat->ntab;
/*-- Compute the number of valid input extensions */
    next = 0;
    for (ntab = 0 ; ntab<imacat->ntab; ntab++, imatab = imatab->nexttab)
      {
/*---- Check for the next valid image extension */
      if ((imatab->naxis < 2)
	|| !strncmp(imatab->xtension, "BINTABLE", 8)
	|| !strncmp(imatab->xtension, "ASCTABLE", 8))
        continue;
      next++;
      }
    }

  thecat.next = next;

/* Initialize the CHECK-images */
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
        prefs.check[c] = check_init(prefs.check_name[i], prefs.check_type[i],
			next, prefs.nband);
        }
    }

/* Initialize catalog output */
  NFPRINTF(OUTPUT, "Initializing catalog");
  catout_init();

/* Initialize XML data */
  if (prefs.xml_flag || prefs.cat_type==ASCII_VO)
    init_xml(next);

/* Go through all images */
  nok = -1;
  for (ntab = 0 ; ntab<ntabmax; ntab++, imatab = imatab->nexttab)
    {
/*--  Check for the next valid image extension */
    if (!forcextflag && ((imatab->naxis < 2)
	|| !strncmp(imatab->xtension, "BINTABLE", 8)
	|| !strncmp(imatab->xtension, "ASCTABLE", 8)))
      continue;
    nok++;

/*-- Initial time measurement*/
    time(&thetime1);
    thecat.currext = nok+1;

/*-- Allocate memory to store field (image) information */
    QCALLOC(fields, fieldstruct *, prefs.nband);
    QCALLOC(wfields, fieldstruct *, prefs.nband);

    for (i=0; i<prefs.nband; i++)
      {
      if (!i)
/*------ Load detection field (image) information */
        fields[i] = field_init(prefs.image_name[i], DETECT_FIELD | MEASURE_FIELD,
			ext_image[i]<0? nok:ext_image[0]);
      else
        {
/*------ Load measurement field (image) information */
        fields[i] = field_init(prefs.image_name[i], MEASURE_FIELD,
			ext_image[i]<0? nok:ext_image[i]);
        if ((fields[i]->width!=fields[0]->width)
		|| (fields[i]->height!=fields[i]->height))
          error(EXIT_FAILURE, "*Error*: Frames have different sizes","");
        }
/*---- Prepare image interpolation */
      if (prefs.weight_flag[i] && prefs.interp_type[i] == INTERP_ALL)
        init_interpolate(fields[i], -1, -1);
      if (prefs.weight_flag[i]) 
        {
/*------ Load weight image information */
        wfields[i] = weight_init(prefs.wimage_name[i], fields[i],
			prefs.weight_type[i],
			ext_wimage[i]<0? nok:ext_wimage[i]);
        interpthresh = prefs.weight_thresh[i];
/*------ Convert the interpolation threshold to variance units */
        weight_to_var(wfields[i], &interpthresh, 1);
        wfields[i]->weight_thresh = interpthresh;
          if (prefs.interp_type[i] != INTERP_NONE)
            init_interpolate(wfields[i],
			prefs.interp_xtimeout[i], prefs.interp_ytimeout[i]);
        }
      }

/*-- Init the FLAG-images */
    QCALLOC(ffields, fieldstruct *, prefs.nimaflag);
    for (i=0; i<prefs.nimaflag; i++)
      {
      ffields[i] = field_init(prefs.fimage_name[i], FLAG_FIELD,
		ext_flag[i]<0? nok:ext_flag[i]);
      if ((ffields[i]->width!=fields[0]->width)
	|| (ffields[i]->height!=fields[0]->height))
        error(EXIT_FAILURE,
	"*Error*: Incompatible FLAG-map size in ", prefs.fimage_name[i]);
      }

/*-- Compute background maps for `standard' fields */
    for (i=0; i<prefs.nband; i++)
      {
      QPRINTF(OUTPUT, i? "Measurement image:":"Detection+Measurement image: ");
      back_map(fields[i], wfields[i], prefs.wscale_flag[i]);
      if (i)
        QPRINTF(OUTPUT,
		"Background: %-10g RMS: %-10g / Analysis threshold: %-10g \n",
		fields[i]->backmean, fields[i]->backsig, fields[i]->thresh);
      else
        QPRINTF(OUTPUT,
		"Background: %-10g RMS: %-10g / Detection threshold: %-10g"
		" / Analysis threshold: %-10g \n",
		fields[i]->backmean, fields[i]->backsig,
		fields[i]->dthresh, fields[i]->thresh);

/*---- For interpolated weight-maps, copy the background structure */
      if (wfields[i] && wfields[i]->flags&(INTERP_FIELD|BACKRMS_FIELD))
        back_copy(wfields[i]->reffield, wfields[i]);
      }

/*-- Prepare learn and/or associations */
    if (prefs.assoc_flag)
      init_assoc(fields[0]);

/*-- Update CHECK-images */
    if (prefs.check_flag)
      for (i=0; i<MAXCHECK; i++)
        if ((check=prefs.check[i]))
          check_reinit(fields[0], check);	/* FIX */

/*-- Initialize PSF contexts and workspace */
    if (prefs.psf_flag)
      for (i=0; i<prefs.nband; i++)
        if (psfs[i])
          {
          psf_readcontext(psfs[i], fields[i]);
          psf_init(psfs[i]);
          fields[i]->psf = &psf[i];
          }

    catout_initext(fields[0]);

/*-- Start the extraction pipeline */
    NFPRINTF(OUTPUT, "Scanning image");
    scanimage(fields, wfields, prefs.nband, ffields, prefs.nimaflag);

/*-- Finish the current CHECK-image processing */
    if (prefs.check_flag)
      for (i=0; i<MAXCHECK; i++)
        if ((check=prefs.check[i]))
          check_reend(fields, check);

/*-- Final time measurements*/
    if (time(&thetime2)!=-1)
      {
      if (!strftime(thecat.ext_date, 12, "%d/%m/%Y", localtime(&thetime2)))
        error(EXIT_FAILURE, "*Internal Error*: Date string too long ","");
      if (!strftime(thecat.ext_time, 10, "%H:%M:%S", localtime(&thetime2)))
        error(EXIT_FAILURE, "*Internal Error*: Time/date string too long ","");
      thecat.ext_elapsed = difftime(thetime2, thetime1);
      }

    catout_endext();

/* --Update XML data */
    if (prefs.xml_flag || prefs.cat_type==ASCII_VO)
      update_xml(&thecat, fields, wfields);


/*-- Close ASSOC routines */
    end_assoc(fields[0]);

/*-- End flag-images */
    for (i=0; i<prefs.nimaflag; i++)
      field_end(ffields[i]);
/*-- End science images and weight maps */
    for (i=0; i<prefs.nband; i++)
      {
      field_end(fields[i]);
      if (wfields[i])
        field_end(wfields[i]);
      }
    QPRINTF(OUTPUT,"      Objects: detected %-8d / sextracted %-8d        \n\n",
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
        check_end(check);
      prefs.check[i] = NULL;
      }

/* End detection filter */
  if (prefs.filter_flag)
    endfilter();

/* End som-fitting */
  if (prefs.somfit_flag)
    som_end(som);

#ifdef USE_MODEL
  if (prefs.prof_flag)
    fft_end();
#endif

/* End PSFs */
  if (prefs.psf_flag)
    for (i=0; i<prefs.nband; i++)
      if (psfs[i])
        psf_end(psfs[i]);
  free(psfs);

/* End classification neural network */
  if (FLAG(obj2.sprob))
    neurclose();

/* Processing end date and time */
  thetimet2 = time(NULL);
  tm = localtime(&thetimet2);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);
  prefs.time_diff = counter_seconds() - dtime;

/* Write XML */
  if (prefs.xml_flag)
    write_xml(prefs.xml_name);

/* End catalog */
  catout_end((char *)NULL);

  if (prefs.xml_flag || prefs.cat_type==ASCII_VO)
    end_xml();

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
VERSION	18/07/2011
 ***/
void	write_error(char *msg1, char *msg2)
  {
   char			error[MAXCHAR];

  sprintf(error, "%s%s", msg1,msg2);
  if (prefs.xml_flag)
    write_xmlerror(prefs.xml_name, error);

/* Also close existing catalog */
  catout_end(error);

  end_xml();

  return;
  }


