 /*
				xml.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	XML logging.
*
*	Last modify:	06/07/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "field.h"
#include "key.h"
#include "prefs.h"
#include "xml.h"

extern pkeystruct	key[];	
extern char		keylist[][32];
xmlstruct		*xmlstack = NULL;
int			nxml, nxmlmax;

/****** init_xml ************************************************************
PROTO	int	init_xml(void)
PURPOSE	Initialize a set of meta-data kept in memory before being written to the
	XML file
INPUT	Number of image extensions.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/07/2006
 ***/
int	init_xml(int next)
  {
  QMALLOC(xmlstack, xmlstruct, next);
  nxml = 0;
  nxmlmax = next;

  return EXIT_SUCCESS;
  }

/****** update_xml ***********************************************************
PROTO	int	update_xml(void)
PURPOSE	Update a set of meta-data kept in memory before being written to the
	XML file
INPUT	-.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/07/2006
 ***/
int	update_xml(sexcatstruct *sexcat, picstruct *dfield, picstruct *field,
		picstruct *dwfield, picstruct *wfield)
  {
   xmlstruct	*x;

  if (nxml >= nxmlmax)
    error(EXIT_FAILURE, "*Internal Error*: too many extensions in XML stack",
			"");
  x = &xmlstack[nxml++];
  x->currext = sexcat->currext;
  x->ndetect = sexcat->ndetect;
  x->ntotal = sexcat->ntotal;
  strcpy(x->ext_date, sexcat->ext_date);
  strcpy(x->ext_time, sexcat->ext_time);
  x->ext_elapsed = sexcat->ext_elapsed;
  strcpy(x->dident, dfield->ident); 
  strcpy(x->ident, field->ident); 
  x->dbackmean = dfield->backmean;
  x->backmean = field->backmean;
  x->dbacksig = dfield->backsig;
  x->backsig = field->backsig;
  x->dsigfac = dfield->sigfac;
  x->sigfac = field->sigfac;
  x->dthresh = dfield->dthresh;
  x->thresh = field->thresh;
  x->dpixscale = dfield->pixscale;
  x->pixscale = field->pixscale;
  x->depoch = dfield->epoch;
  x->epoch = field->epoch;

  return EXIT_SUCCESS;
  }

/****** write_xml ************************************************************
PROTO	int	write_xml(void)
PURPOSE	Save meta-data to XML files
INPUT	-.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/07/2006
 ***/
int	write_xml(void)
  {
   FILE			*file;
   char			*pspath,*psuser, *pshost;
   int			n;

/* Username */
  psuser = pspath = pshost = NULL;
#ifdef HAVE_GETENV
  if (!(psuser=getenv("USERNAME")))	/* Cygwin,... */
    psuser = getenv("LOGNAME");		/* Linux,... */
  pspath = getenv("PWD");
  pshost = getenv("HOSTNAME");
#endif

  if ((file = fopen(prefs.xml_name, "wb")) == NULL)
    return RETURN_ERROR;
  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n",
	prefs.xsl_name);
/*
  fprintf(file, "<!DOCTYPE INFO SYSTEM \"instfile.dtd\">\n");
*/
  fprintf(file, "<SOURCE_EXTRACTION>\n");
  fprintf(file, " <software type=\"char\">%s</software>\n", BANNER);
  fprintf(file, " <version type=\"char\">%s</version>\n", MYVERSION);
  fprintf(file, " <date type=\"char\">%s</date>\n", prefs.sdate_end);
  fprintf(file, " <time type=\"char\">%s</time>\n", prefs.stime_end);
  fprintf(file, " <duration type=\"int\" unit=\"s\">%.0f</duration>\n",
    	prefs.time_diff);
  fprintf(file, " <nthreads type=\"int\">%d</nthreads>\n",
    	prefs.nthreads);
  fprintf(file, " <user type=\"char\">%s</user>\n", psuser);
  fprintf(file, " <host type=\"char\">%s</host>\n", pshost);
  fprintf(file, " <path type=\"char\">%s</path>\n", pspath);
  fprintf(file, " <DETECTION_IMAGE>\n");
  fprintf(file, "  <image_name type=\"char\">%s</image_name>\n",
    	prefs.image_name[0]);
  fprintf(file, " </DETECTION_IMAGE>\n");
  fprintf(file, " <MEASUREMENT_IMAGE>\n");
  fprintf(file, "  <image_name type=\"char\">%s</image_name>\n",
    	prefs.image_name[1]? prefs.image_name[1] : prefs.image_name[0]);
  fprintf(file, " </MEASUREMENT_IMAGE>\n");
  fprintf(file, " <CONFIG>\n");
  fprintf(file, "  <catalog_name type=\"char\">%s</catalog_name>\n",
    	prefs.cat_name);
  fprintf(file, "  <catalog_type type=\"char\">%s</catalog_type>\n",
    	key[findkeys("CATALOG_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.cat_type]);
  fprintf(file, "  <parameters_name type=\"char\">%s</parameters_name>\n",
    	prefs.param_name);
  fprintf(file, "  <detect_type type=\"char\">%s</detect_type>\n",
    key[findkeys("DETECT_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.detect_type]);
  fprintf(file, "  <detect_minarea type=\"int\">%d</detect_minarea>\n",
	prefs.ext_minarea);
  fprintf(file, "  <deblend_nthresh type=\"int\">%d</deblend_nthresh>\n",
	prefs.deblend_nthresh);
  fprintf(file, "  <deblend_mincont type=\"float\">%g</deblend_mincont>\n",
	prefs.deblend_mincont);
  fprintf(file, "  <filter_flag type=\"bool\">%c</filter_flag>\n",
    	prefs.filter_flag? 'T':'F');
  fprintf(file, "  <filter_name type=\"char\">%s</filter_name>\n",
    	prefs.filter_name);
  fprintf(file, "  <clean_flag type=\"bool\">%c</clean_flag>\n",
    	prefs.clean_flag? 'T':'F');
  fprintf(file, "  <clean_param type=\"float\">%g</clean_param>\n",
	prefs.clean_param);
  fprintf(file, "  <mask_type type=\"char\">%s</mask_type>\n",
	key[findkeys("MASK_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.mask_type]);
  fprintf(file, "  <weight_gain type=\"bool\">%c</weight_gain>\n",
    	prefs.weightgain_flag? 'T':'F');
  for (n=0; n<prefs.naper; n++)
    fprintf(file, "  <phot_apertures%d type=\"float\">%g</phot_apertures%d>\n",
	n+1, prefs.apert[n], n+1);
  for (n=0; n<prefs.nautoparam; n++)
    fprintf(file, "  <phot_autoparams%d type=\"float\">%g</phot_autoparams%d>\n",
	n+1, prefs.autoparam[n], n+1);
  for (n=0; n<prefs.npetroparam; n++)
    fprintf(file, "  <phot_petroparams%d type=\"float\">%g</phot_petroparams%d>\n",
	n+1, prefs.petroparam[n], n+1);
  for (n=0; n<prefs.nautoaper; n++)
    fprintf(file, "  <phot_autoapers%d type=\"float\">%g</phot_autoapers%d>\n",
	n+1, prefs.autoaper[n], n+1);
  for (n=0; n<prefs.nflux_frac; n++)
    fprintf(file, "  <phot_fluxfrac%d type=\"float\">%g</phot_fluxfrac%d>\n",
	n+1, prefs.flux_frac[n], n+1);
  fprintf(file, "  <satur_level type=\"float\">%g</satur_level>\n",
	prefs.satur_level);
  fprintf(file, "  <mag_zeropoint type=\"float\">%g</mag_zeropoint>\n",
	prefs.mag_zeropoint);
  fprintf(file, "  <mag_gamma type=\"float\">%g</mag_gamma>\n",
	prefs.mag_gamma);
  fprintf(file, "  <gain type=\"float\">%g</gain>\n",
	prefs.gain);
  fprintf(file, "  <pixel_scale type=\"float\">%g</pixel_scale>\n",
	prefs.pixel_scale);
  fprintf(file, "  <seeing_fwhm type=\"float\">%g</seeing_fwhm>\n",
	prefs.seeing_fwhm);
  fprintf(file, "  <starnnw_name type=\"char\">%s</starnnw_name>\n",
    	prefs.nnw_name);
  for (n=0; n<prefs.nbacksize; n++)
    fprintf(file, "  <back_size%d type=\"int\">%d</back_size%d>\n",
	n+1, prefs.backsize[n], n+1);
  for (n=0; n<prefs.nbackfsize; n++)
    fprintf(file, "  <back_filtersize%d type=\"int\">%d</back_filtersize%d>\n",
	n+1, prefs.backfsize[n], n+1);
  fprintf(file, "  <backphoto_type type=\"char\">%s</backphoto_type>\n",
    key[findkeys("BACKPHOTO_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.pback_type]);
  fprintf(file, "  <backphoto_thick type=\"int\">%d</backphoto_thick>\n",
	prefs.pback_size);
  fprintf(file, "  <back_filtthresh type=\"float\">%g</back_filtthresh>\n",
	prefs.backfthresh);
  for (n=0; n<prefs.ncheck_type; n++)
    {
    fprintf(file, "  <checkimage_type%d type=\"char\">%s</checkimage_type%d>\n",
	n+1, key[findkeys("CHECKIMAGE_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.check_type[n]], n+1);
    if (prefs.check_type[n] != CHECK_NONE)
      fprintf(file, " <checkimage_name%d type=\"char\">%s</checkimage_name%d>\n",
    	n+1, prefs.check_name[n], n+1);
    }
  fprintf(file, "  <memory_objstack type=\"int\">%d</memory_objstack>\n",
	prefs.clean_stacksize);
  fprintf(file, "  <memory_pixstack type=\"int\">%d</memory_pixstack>\n",
	prefs.mem_pixstack);
  fprintf(file, "  <memory_bufsize type=\"int\">%d</memory_bufsize>\n",
	prefs.mem_bufsize);
  fprintf(file, "  <assoc_name type=\"char\">%s</assoc_name>\n",
    	prefs.assoc_name);
  for (n=0; n<prefs.nassoc_data; n++)
    fprintf(file, "  <assoc_data%d type=\"int\">%d</assoc_data%d>\n",
	n+1, prefs.assoc_data[n], n+1);
  for (n=0; n<prefs.nassoc_param; n++)
    fprintf(file, "  <assoc_params%d type=\"int\">%d</assoc_params%d>\n",
	n+1, prefs.assoc_param[n], n+1);
  fprintf(file, "  <assoc_radius type=\"float\">%g</assoc_radius>\n",
	prefs.assoc_radius);
  fprintf(file, "  <assoc_type type=\"char\">%s</assoc_type>\n",
    key[findkeys("ASSOC_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.assoc_type]);
  fprintf(file, "  <assocselec_type type=\"char\">%s</assocselec_type>\n",
    key[findkeys("ASSOCSELEC_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.assocselec_type]);
  fprintf(file, "  <verbose_type type=\"char\">%s</verbose_type>\n",
    key[findkeys("VERBOSE_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.verbose_type]);
  fprintf(file, "  <fits_unsigned type=\"bool\">%c</fits_unsigned>\n",
    	prefs.fitsunsigned_flag? 'T':'F');
  fprintf(file, "  <psf_name type=\"char\">%s</psf_name>\n",
    	prefs.psf_name[0]);
  fprintf(file, "  <psf_nmax type=\"int\">%d</psf_nmax>\n",
	prefs.psf_npsfmax);
  fprintf(file, "  <psfdisplay_type type=\"char\">%s</psfdisplay_type>\n",
    key[findkeys("PSFDISPLAY_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.psfdisplay_type]);
  fprintf(file, "  <som_name type=\"char\">%s</som_name>\n",
    	prefs.som_name);
  fprintf(file, "  <DETECTION_IMAGE>\n");
  fprintf(file, "   <thresh_type type=\"char\">%s</thresh_type>\n",
	key[findkeys("THRESH_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.thresh_type[0]]);
  fprintf(file, "   <weight_type type=\"char\">%s</weight_type>\n",
	key[findkeys("WEIGHT_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.weight_type[0]]);
  if (prefs.weight_type[0] != WEIGHT_NONE
	&& prefs.weight_type[0] != WEIGHT_FROMBACK)
    fprintf(file, "   <weight_image type=\"char\">%s</weight_image>\n",
    	prefs.wimage_name[0]);
  fprintf(file, "   <weight_thresh type=\"float\">%g</weight_thresh>\n",
	prefs.weight_thresh[0]);
  fprintf(file, "   <back_type type=\"char\">%s</back_type>\n",
	key[findkeys("BACK_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.back_type[0]]);
  fprintf(file, "   <back_value type=\"float\">%g</back_value>\n",
	prefs.back_val[0]);
  fprintf(file, "   <interp_maxxlag>%d</interp_maxxlag>\n",
	prefs.interp_xtimeout[0]);
  fprintf(file, "   <interp_maxylag>%d</interp_maxylag>\n",
	prefs.interp_ytimeout[0]);
  fprintf(file, "   <interp_type type=\"char\">%s</interp_type>\n",
    key[findkeys("INTERP_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.interp_type[0]]);
  fprintf(file, "  </DETECTION_IMAGE>\n");
  fprintf(file, "  <MEASUREMENT_IMAGE>\n");
  fprintf(file, "   <thresh_type type=\"char\">%s</thresh_type>\n",
	key[findkeys("THRESH_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.thresh_type[1]]);
  fprintf(file, "   <weight_type type=\"char\">%s</weight_type>\n",
	key[findkeys("WEIGHT_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.weight_type[1]]);
  if (prefs.weight_type[1] != WEIGHT_NONE
	&& prefs.weight_type[1] != WEIGHT_FROMBACK)
    fprintf(file, "   <weight_image type=\"char\">%s</weight_image>\n",
    	prefs.wimage_name[1]);
  fprintf(file, "   <weight_thresh type=\"float\">%g</weight_thresh>\n",
	prefs.weight_thresh[1]);
  fprintf(file, "   <back_type type=\"char\">%s</back_type>\n",
	key[findkeys("BACK_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.back_type[1]]);
  fprintf(file, "   <back_value type=\"float\">%g</back_value>\n",
	prefs.back_val[1]);
  fprintf(file, "   <interp_maxxlag>%d</interp_maxxlag>\n",
	prefs.interp_xtimeout[1]);
  fprintf(file, "   <interp_maxylag>%d</interp_maxylag>\n",
	prefs.interp_ytimeout[1]);
  fprintf(file, "   <interp_type type=\"char\">%s</interp_type>\n",
    key[findkeys("INTERP_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.interp_type[1]]);
  fprintf(file, "  </MEASUREMENT_IMAGE>\n");
  for (n=0; n<prefs.nimaflag; n++)
    {
    fprintf(file, " <FLAG_IMAGE>\n");
    fprintf(file, "  <image_name type=\"char\">%s</image_name>\n",
    	prefs.fimage_name[n]);
    fprintf(file, "  <flag_type type=\"char\">%s</flag_type>\n",
	key[findkeys("FLAG_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.flag_type[n]]);
    fprintf(file, " </FLAG_IMAGE>\n");
    }
  fprintf(file, " </CONFIG>\n");
  fprintf(file, " <nextens type=\"int\">%d</nextens>\n", nxmlmax);
  for (n=0; n<nxmlmax; n++)
    {
    fprintf(file, " <EXT_PROPS>\n");
    fprintf(file, "  <extension type=\"int\">%d</extension>\n",
    	xmlstack[n].currext);
    fprintf(file, "  <duration type=\"int\" unit=\"s\">%.0f</duration>\n",
    	xmlstack[n].ext_elapsed);
    fprintf(file, "  <ndetect type=\"int\">%d</ndetect>\n",
    	xmlstack[n].ndetect);
    fprintf(file, "  <nsextracted type=\"int\">%d</nsextracted>\n",
    	xmlstack[n].ntotal);
    fprintf(file, "  <DETECTION_IMAGE>\n");
    fprintf(file, "   <image_ident type=\"char\">%s</image_ident>\n",
    	xmlstack[n].dident);
    fprintf(file, "   <background_mean type=\"float\">%g</background_mean>\n",
	xmlstack[n].dbackmean);
    fprintf(file, "   <background_stddev type=\"float\">%g</background_stddev>\n",
	xmlstack[n].dbacksig);
    fprintf(file, "   <threshold type=\"float\">%g</threshold>\n",
	xmlstack[n].dthresh);
    fprintf(file, "   <weight_scaling type=\"float\">%g</weight_scaling>\n",
	xmlstack[n].dsigfac);
    fprintf(file, "   <pixel_scale type=\"float\" unit=\"deg\">%g</pixel_scale>\n",
	xmlstack[n].dpixscale/3600.0);
    fprintf(file, "   <epoch type=\"float\">%g</epoch>\n",
	xmlstack[n].depoch);
    fprintf(file, "  </DETECTION_IMAGE>\n");
    fprintf(file, "  <MEASUREMENT_IMAGE>\n");
    fprintf(file, "   <image_ident type=\"char\">%s</image_ident>\n",
    	xmlstack[n].ident);
    fprintf(file, "   <background_mean type=\"float\">%g</background_mean>\n",
	xmlstack[n].backmean);
    fprintf(file, "   <background_stddev type=\"float\">%g</background_stddev>\n",
	xmlstack[n].backsig);
    fprintf(file, "   <threshold type=\"float\">%g</threshold>\n",
	xmlstack[n].thresh);
    fprintf(file, "   <weight_scaling type=\"float\">%g</weight_scaling>\n",
	xmlstack[n].sigfac);
    fprintf(file, "   <pixel_scale type=\"float\" unit=\"deg\">%g</pixel_scale>\n",
	xmlstack[n].pixscale/3600.0);
    fprintf(file, "   <epoch type=\"float\">%g</epoch>\n",
	xmlstack[n].epoch);
    fprintf(file, "  </MEASUREMENT_IMAGE>\n");
    fprintf(file, " </EXT_PROPS>\n");
    }
  fprintf(file, "</SOURCE_EXTRACTION>\n");
  fclose(file);

  free(xmlstack);

  return RETURN_OK;
  }


/****** write_xmlerror ********************************************************
PROTO	int	write_xmlerror(char *msg1, char *msg2)
PURPOSE	Save meta-data to a simplified XML file in case of a catched error
INPUT	a character string,
	another character string
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/07/2006
 ***/
void	write_xmlerror(char *msg1, char *msg2)
  {
   FILE			*file;
   time_t		thetime;
   struct tm		*tm;
   char			*pspath,*psuser, *pshost;

  if (!prefs.xml_flag)
    return;
/* Processing date and time */
  thetime = time(NULL);
  tm = localtime(&thetime);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);

/* Username */
  psuser = pspath = pshost = NULL;
#ifdef HAVE_GETENV
  if (!(psuser=getenv("USERNAME")))	/* Cygwin,... */
    psuser = getenv("LOGNAME");		/* Linux,... */
  pspath = getenv("PWD");
  pshost = getenv("HOSTNAME");
#endif

  if ((file = fopen(prefs.xml_name, "wb")) == NULL)
    return;
  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n",
	prefs.xsl_name);
/*
  fprintf(file, "<!DOCTYPE INFO SYSTEM \"instfile.dtd\">\n");
*/
  fprintf(file, "<SOURCE_EXTRACTION>\n");
  fprintf(file, " <software type=\"char\">%s</software>\n", BANNER);
  fprintf(file, " <version type=\"char\">%s</version>\n", MYVERSION);
  fprintf(file, " <date type=\"char\">%s</date>\n", prefs.sdate_end);
  fprintf(file, " <time type=\"char\">%s</time>\n", prefs.stime_end);
  fprintf(file, " <nthreads type=\"int\">%d</nthreads>\n",
    	prefs.nthreads);
  fprintf(file, " <user type=\"char\">%s</user>\n", psuser);
  fprintf(file, " <host type=\"char\">%s</host>\n", pshost);
  fprintf(file, " <path type=\"char\">%s</path>\n", pspath);
  fprintf(file, " <DETECTION_IMAGE>\n");
  fprintf(file, "  <image_name type=\"char\">%s</image_name>\n",
    	prefs.image_name[0]);
  fprintf(file, " </DETECTION_IMAGE>\n");
  fprintf(file, " <MEASUREMENT_IMAGE>\n");
  fprintf(file, "  <image_name type=\"char\">%s</image_name>\n",
    	prefs.image_name[1]? prefs.image_name[1] : prefs.image_name[0]);
  fprintf(file, " </MEASUREMENT_IMAGE>\n");
  fprintf(file, " <ERROR_MSG type=\"char\">%s%s</ERROR_MSG>\n", msg1,msg2);
  fprintf(file, "</SOURCE_EXTRACTION>\n");
  fclose(file);

  free(xmlstack);

  return;
  }

