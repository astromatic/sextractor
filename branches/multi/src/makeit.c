/*
*				makeit.c
*
* Main program.
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
*	Last modified:		06/05/2012
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
#include	"misc.h"
#include	"neurro.h"
#include	"pattern.h"
#include	"photom.h"
#include	"psf.h"
#include	"profit.h"
#include	"scan.h"
#include	"som.h"
#include	"weight.h"
#include	"xml.h"

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
VERSION	06/05/2012
 ***/
void	makeit(void)
  {
   profitstruct		*profit;
   checkstruct		*check;
   fieldstruct		**fields,**wfields,**ffields,
			**efields,**wefields,**fefields,
			*field,*wfield,*ffield, *dfield,*dwfield;
   catstruct		*imacat;
   tabstruct		*imatab;
   patternstruct	*pattern;
   psfstruct		**psfs;
   somstruct		*som;
   static time_t        thetime1, thetime2;
   struct tm		*tm;
   PIXTYPE		interpthresh;
   char			str[82];
   int			ext_fimage[MAXFLAG],ext_image[MAXIMAGE],
			ext_wimage[MAXIMAGE], nparam2[2],
			*file_index,*wfile_index,*ffile_index,*file_next,
			*ffile_next,
			e,e2,i,j,t, nok, ntab, previndex, ne,ext,
			next,next0,nfext, ntabmax, forcextflag, npat, nmodels,
			nfield,nffield, nfieldmax,nffieldmax, nimage,nfimage;

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
  prefs_use();	/* update things a first time according to prefs parameters */

  if (prefs.filter_flag)
    {
    NFPRINTF(OUTPUT, "Reading detection filter");
    getfilter(prefs.filter_name);	/* get the detection filter */
    }


/* Read extension numbers and remove them from the filenames if present */
 for (i=0; i<prefs.nimage; i++)
   ext_image[i] = selectext(prefs.image_name[i]);
 for (i=0; i<prefs.nwimage; i++)
   ext_wimage[i] = selectext(prefs.wimage_name[i]);
 for (i=0; i<prefs.nfimage; i++)
   ext_fimage[i] = selectext(prefs.fimage_name[i]);

/* Use the first image to probe the number of extensions to analyse */
  if (ext_image[0] != RETURN_ERROR)
    {
    forcextflag = 1;
    ntabmax = next0 = 1;
    }
  else
    forcextflag = 0;

  if (!(imacat = read_cat(prefs.image_name[0])))
    error(EXIT_FAILURE, "*Error*: cannot open ", prefs.image_name[0]);
  close_cat(imacat);
  imatab = imacat->tab;

  thecat.next = next0;


/* Examine all extensions of all input images */
  nfieldmax = nimage = prefs.nimage;
  QMALLOC(fields, fieldstruct *, nfieldmax);
  QMALLOC(wfields, fieldstruct *, nfieldmax);
  QMALLOC(file_index, int, nimage);
  QMALLOC(wfile_index, int, nimage);
  QMALLOC(file_next, int, nimage);
  ext = 0;
  previndex = -1;
  if ((prefs.multigrids_flag))
    QPRINTF(OUTPUT, "Multi-grid mode activated\n");

  for (i=0; i<nimage; i++, ext+=next)
    {
    if (!(imacat = read_cat(prefs.image_name[i])))
      error(EXIT_FAILURE, "*Error*: cannot open ", prefs.image_name[i]);
    if (!forcextflag)
      {
/*---- Compute the number of valid input extensions */
      next = 0;
      for (ntab = 0 ; ntab<imacat->ntab; ntab++, imatab = imatab->nexttab)
        {
/*------ Check for the next valid image extension */
        if ((imatab->naxis < 2)
		|| !strncmp(imatab->xtension, "BINTABLE", 8)
		|| !strncmp(imatab->xtension, "ASCTABLE", 8))
          continue;
        next++;
        }
      }
    if (ext+next>nfieldmax)
      {
      nfieldmax = (ext+next)*2;
      QREALLOC(fields, fieldstruct *, nfieldmax);
      QREALLOC(wfields, fieldstruct *, nfieldmax);
      }
    for (e=0; e<next; e++)
      {
      ne = ext+e;
      if (next>1)
        sprintf(str, "Examining %s %.40s[%d/%d]...",
		i? "measurement image":"detection image",
		prefs.image_name[i], e+1,next);
      else
        sprintf(str, "Examining %s %.40s...",
		i? "measurement image":"detection image",
		prefs.image_name[i]);
      NFPRINTF(OUTPUT, str);
      field = fields[ne] = field_init(prefs.image_name[i], i,
		forcextflag? ext_image[i]:e,
		i? MEASURE_FIELD|(prefs.multigrid_flag[i]*MULTIGRID_FIELD)
		: DETECT_FIELD|MEASURE_FIELD);
        
      if (previndex>=0 && (field->width != fields[previndex+e]->width
		|| field->height != fields[previndex+e]->height))
        error(EXIT_FAILURE, "*Error*: unexpected image size in ",
		prefs.image_name[i]);
/*---- Prepare image interpolation */
      if (prefs.weight_flag[i])
        {
        if (prefs.interp_type[i] == INTERP_ALL)
          init_interpolate(fields[ne], -1, -1);
        if (next>1)
          sprintf(str, "Examining %s %.40s[%d/%d]...",
		"weight image",
		prefs.wimage_name[i], e+1,next);
        else
          sprintf(str, "Examining %s %.40s...",
		"weight image",
		prefs.wimage_name[i]);
        NFPRINTF(OUTPUT, str);
        wfield = wfields[ne] = weight_init(prefs.wimage_name[i], field, i,
		forcextflag? ext_wimage[i]:e,
		prefs.weight_type[i]);
        interpthresh = prefs.weight_thresh[i];
/*------ Convert the interpolation threshold to variance units */
        weight_to_var(wfield, &interpthresh, 1);
        wfield->weight_thresh = interpthresh;
          if (prefs.interp_type[i] != INTERP_NONE)
            init_interpolate(wfield, prefs.interp_xtimeout[i],
				prefs.interp_ytimeout[i]);
        }
      else
        wfields[ne] = NULL;
      }

    previndex = wfile_index[i] = file_index[i] = ext;
    file_next[i] = next;
    }

  next0 = file_next[0];
  QMALLOC(efields, fieldstruct *, nfieldmax);
  QMALLOC(wefields, fieldstruct *, nfieldmax);

  NFPRINTF(OUTPUT,"");

  nfield = ext;
  photom_printinstruinfo();

/* Init the FLAG-images */
  nfimage = 0;
  if ((prefs.imaflags_flag))
    {
    nffieldmax = nfimage = prefs.nfimage;
    nfext = next0;
    if ((nfimage))
      {
      QMALLOC(fefields, fieldstruct *, nffieldmax)
      QMALLOC(ffields, fieldstruct *, nfimage);
      QMALLOC(ffile_index, int, nfimage);
      QMALLOC(ffile_next, int, nfimage);
      }
    else
      fefields = ffields = NULL;
    ext = 0;
    previndex = -1;
    for (i=0; i<nfimage; i++, ext+=nfext)
      {
      if (!(imacat = read_cat(prefs.fimage_name[i])))
        error(EXIT_FAILURE, "*Error*: cannot open ", prefs.fimage_name[i]);
      if (!forcextflag)
        {
/*------ Compute the number of valid input extensions */
        nfext = 0;
        for (ntab = 0 ; ntab<imacat->ntab; ntab++, imatab = imatab->nexttab)
          {
/*-------- Check for the next valid image extension */
          if ((imatab->naxis < 2)
		|| !strncmp(imatab->xtension, "BINTABLE", 8)
		|| !strncmp(imatab->xtension, "ASCTABLE", 8))
            continue;
          nfext++;
          }
        }
      if (ext+nfext>nffieldmax)
        {
        nffieldmax *= 2;
        QREALLOC(ffields, fieldstruct *, nffieldmax);
        }
      for (e=0; e<nfext; e++)
        {
        ne = ext+e;
        if (nfext>1)
          sprintf(str, "Examining %s %.40s[%d/%d]...",
		"flag image",
		prefs.fimage_name[i], e+1,nfext);
        else
          sprintf(str, "Examining %s %.40s...",
		"flag image",
		prefs.fimage_name[i]);
        NFPRINTF(OUTPUT, str);
        ffield = ffields[ne] = field_init(prefs.fimage_name[i], i,
		forcextflag? ext_fimage[i]:e, FLAG_FIELD);
        if ((previndex>=0) && (ffield->width != ffields[previndex+e]->width
		|| ffield->height != ffields[previndex+e]->height))
          {
          sprintf(str, "*Error*: %s does not match detection image size\n",
			prefs.fimage_name[i]);
          error(EXIT_FAILURE, str,
		"         You might want to set MULTIGRID to Y");
          }
        }
      previndex = ffile_index[i] = ext;
      ffile_next[i] = nfext;
      }

    nffield = ext;
    }

/*-- Initialize the CHECK-images */
  if (prefs.check_flag)
    {
     checkenum	c;

    NFPRINTF(OUTPUT, "Initializing check-image(s)...");
    for (i=0; i<prefs.ncheck_type; i++)
      if ((c=prefs.check_type[i]) != CHECK_NONE)
        {
        if (prefs.check[c])
           error(EXIT_FAILURE,"*Error*: 2 CHECK_IMAGEs cannot have the same ",
			" CHECK_IMAGE_TYPE");
        prefs.check[c] = check_init(prefs.check_name[i], prefs.check_type[i],
			next0, nimage);
        }
    }

  if (prefs.psf_flag)
    {
    NFPRINTF(OUTPUT, "Reading PSF information...");
    QCALLOC(psfs, psfstruct *, nfield);
    for (i=0; i<nfield; i++)
      psfs[i] = psf_load(prefs.psf_name[i]);
 /*-- Need to check things up because of PSF context parameters */
    catout_updateparamflags();
    prefs_use();
    }

  if (prefs.prof_flag)
    {
#ifdef USE_MODEL
    fft_init(prefs.nthreads);
/* Create profiles at full resolution */
    NFPRINTF(OUTPUT, "Preparing profile models...");
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
    QCALLOC(profit->subprofit, subprofitstruct, 1);
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

  if (FLAG(obj2.sprob))
    {
    NFPRINTF(OUTPUT, "Initializing classification neural network...");
    neur_init();
    NFPRINTF(OUTPUT, "Reading neural network weights...");
    neur_getnnw(prefs.nnw_name); 
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
/* Allocate memory for other arrays (not catalogue measurements) */
/* Sky background */
  catout_allocother(thecat.obj2list, &flagobj2.dbkg, nimage*sizeof(float));
  catout_allocother(thecat.obj2list, &flagobj2.sigbkg, nimage*sizeof(float));
/* Flux combination */
  catout_allocother(thecat.obj2list, &flagobj2.cflux,
		prefs.nphotinstru*sizeof(double));
  catout_allocother(thecat.obj2list, &flagobj2.cfluxw,
		prefs.nphotinstru*sizeof(double));
/* Position combination */
  catout_allocother(thecat.obj2list, &flagobj2.cposx,
		prefs.nphotinstru*sizeof(double));
  catout_allocother(thecat.obj2list, &flagobj2.cposy,
		prefs.nphotinstru*sizeof(double));
  catout_allocother(thecat.obj2list, &flagobj2.cposw,
		prefs.nphotinstru*sizeof(double));
  prefs_use();


/* Initialize catalog output */
  NFPRINTF(OUTPUT, "Initializing catalogue...");
  catout_init();

/* Initialize XML data */
  NFPRINTF(OUTPUT, "Initializing XML output...");
  if (prefs.xml_flag || prefs.cat_type==ASCII_VO)
    xml_init(next0);

/* Initial time measurement*/
  time(&thetime1);

/* Pre-process background subtraction in multigrid mode */
  if ((prefs.multigrids_flag))
    for (i=0; i<nimage; i++)
      {
      next = file_next[i];
      field = fields[file_index[i]];
      wfield = wfields[wfile_index[i]];
      for (e=0; e<next; e++, field++, wfield++)
        {
        if (next>1)
          sprintf(str, "[%d/%d]", e+1, next);
        NFPRINTF(OUTPUT, "");
        QPRINTF(OUTPUT, " \n");
        QPRINTF(OUTPUT, "----- Image %.60s%s:\n",
		field->rfilename, next>1? str:"");
        field_printinfo(field, wfield);
        back_map(field, wfield, prefs.wscale_flag[i]);
        if ((i))
          {
          QPRINTF(OUTPUT,
		"    Background: %.6g   RMS: %.5g"
		"   Analysis threshold: %.5g \n",
		field->backmean, field->backsig, field->thresh);
          }
        else
          {
          QPRINTF(OUTPUT,
		"    Background: %.6g   RMS: %.5g   Analysis threshold: %.5g"
		"   Detection threshold: %.5g \n",
		field->backmean, field->backsig,
		field->thresh, field->dthresh);
          }
        }
      }

/* Process one extension at a time */
  for (e=0; e<next0; e++)
    {
    dfield = fields[e];
    dwfield = wfields[e];
    if (!(prefs.multigrids_flag))
      {
/*---- Multigrid is off: select only one extension from every file */
      nfield = nimage;
      for (i=0; i<nimage; i++)
        {
        field = efields[i] = fields[i*next0+e];
        wfield = wefields[i] = wfields[i*next0+e];
        if (next0>1)
          sprintf(str, "[%d/%d]", e+1, next0);
        NFPRINTF(OUTPUT, "");
        QPRINTF(OUTPUT, " \n");
        QPRINTF(OUTPUT, "----- Image %.60s%s:\n",
		field->rfilename, next0>1? str:"");
        field_printinfo(field, wfield);
        back_map(field, wfield, prefs.wscale_flag[i]);
        if ((i))
          {
          QPRINTF(OUTPUT,
		"    Background: %.6g   RMS: %.5g"
		"   Analysis threshold: %.5g \n",
		field->backmean, field->backsig, field->thresh);
          }
        else
          {
          QPRINTF(OUTPUT,
		"    Background: %.6g   RMS: %.5g   Analysis threshold: %.5g"
		"   Detection threshold: %.5g \n",
		field->backmean, field->backsig,
		field->thresh, field->dthresh);
          }

/*------ For interpolated weight-maps, copy the background structure */
        if (wfield && wfield->flags&(INTERP_FIELD|BACKRMS_FIELD))
          back_copy(wfield->reffield, wfield);

/*------ Initialize PSF contexts and workspace */
        if (prefs.psf_flag)
          if (psfs[i])
            {
            psf_readcontext(psfs[i], field);
            psf_init(psfs[i]);
            fields[i]->psf = psfs[i];
            }
        }
      }
    else
      {
/*---- Multigrid is on: hybrid selection of fields */
      j=0;
      for (i=0; i<nimage; i++)
        {
        if (i && prefs.multigrid_flag[i])
          {
/*-------- Identify extensions overlapping current detection extension */
          field = fields[file_index[i]];
          wfield = wfields[file_index[i]];
          for (e2=file_next[i]; e2--; field++, wfield++)
            {
            if ((frame_wcs(field->wcs, dfield->wcs)))
              {
              efields[j] = field;
              wefields[j++] = wfield;
              }
            }
          }
        else
          {
/*-------- Basic (unique) extension selection */
          efields[j] = fields[i*next0+e];
          wefields[j++] = wfields[i*next0+e];
          }
        }
      nfield = j;
      if (next0>1)
        sprintf(str, "[%d/%d]", e+1, next0);
      NFPRINTF(OUTPUT, "");
      QPRINTF(OUTPUT, " \n");
      QPRINTF(OUTPUT, "----- Image %.60s%s:\n",
	fields[e]->rfilename, next0>1? str:"");
      }

/*-- Flag maps */
    if ((prefs.imaflags_flag))
      for (i=0; i<nfimage; i++)
        fefields[i]= ffields[i*next0+e];

/*-- Prepare learn and/or associations */
    if (prefs.assoc_flag)
      init_assoc(dfield);

/*-- Update CHECK-images */
    if (prefs.check_flag)
      for (i=0; i<MAXCHECK; i++)
        if ((check=prefs.check[i]))
          check_reinit(dfield, check);		/* FIX */

    thecat.currext++;
    catout_initext(dfield);

/*-- Start the extraction pipeline */
    NFPRINTF(OUTPUT, "Scanning image");
    scan_extract(dfield,dwfield, efields,wefields,nfield, fefields,nfimage);

    thecat.ntotalsum += thecat.ntotal;
    thecat.nlinesum += dfield->height;

/*-- Finish the current CHECK-image processing */
    if (prefs.check_flag)
      for (i=0; i<MAXCHECK; i++)
        if ((check=prefs.check[i]))
          check_reend(dfield, check);

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
      xml_update(&thecat, efields, wefields);


/*-- Close ASSOC routines */
    end_assoc(dfield);

/*-- End flag-images */
    if ((nfimage))
      for (i=0; i<nfimage; i++)
        field_end(fefields[i]);
/*-- End science images and weight maps */
    for (i=0; i<nimage; i++)
      {
      field_end(efields[i]);
      if (wefields[i])
        field_end(wefields[i]);
      }
    QPRINTF(OUTPUT,"      Objects: detected %-8d / sextracted %-8d        \n\n",
	thecat.ndetect, thecat.ntotal);
    }

  free(efields);
  free(wefields);
  free(fields);
  free(wfields);
  free(file_index);
  free(wfile_index);
  free(file_next);
  if ((nfimage))
    {
    free(fefields);
    free(ffile_next);
    }

  if (thecat.currext==0)
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
    {
    for (i=0; i<nfield; i++)
      if (psfs[i])
        psf_end(psfs[i], NULL);
    free(psfs);
    }

/* End classification neural network */
  if (FLAG(obj2.sprob))
    neur_end();

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
    xml_write(prefs.xml_name);

/* Free memory allocated for arrays that are not catalogue measurements */
  catout_freeother(thecat.obj2list, &flagobj2.dbkg);
  catout_freeother(thecat.obj2list, &flagobj2.sigbkg);
  catout_freeother(thecat.obj2list, &flagobj2.cflux);
  catout_freeother(thecat.obj2list, &flagobj2.cfluxw);
  catout_freeother(thecat.obj2list, &flagobj2.cposx);
  catout_freeother(thecat.obj2list, &flagobj2.cposy);
  catout_freeother(thecat.obj2list, &flagobj2.cposw);

/* End catalog */
  catout_end((char *)NULL);


  if (prefs.xml_flag || prefs.cat_type==ASCII_VO)
    xml_end();

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
VERSION	24/02/2012
 ***/
void	write_error(char *msg1, char *msg2)
  {
   char			error[MAXCHAR];

  sprintf(error, "%s%s", msg1,msg2);
  if (prefs.xml_flag)
    xml_write_error(prefs.xml_name, error);

/* Also close existing catalog */
  catout_end(error);

  xml_end();

  return;
  }


