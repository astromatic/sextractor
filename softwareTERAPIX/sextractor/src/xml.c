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
*	Last modify:	04/07/2006
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
xmlstruct		*xmlstack;
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
PROTO	int	write_xml(void);
PURPOSE	Save meta-data to XML files
INPUT	Pointer to the array of fields,
	number of fields,
	pointer to the array of field groups,
	number of field groups,
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/07/2006
 ***/
int	write_xml(void)
  {
   FILE			*file;
   time_t		thetime;
   char			*pspath,*psuser, *pshost;
   int			n;

/* Processing date and time */
  thetime = time(NULL);

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
  fprintf(file, " <catalog_name type=\"char\">%s</catalog_name>\n",
    	prefs.cat_name);
  fprintf(file, " <catalog_type type=\"char\">%s</catalog_type>\n",
    	key[findkeys("CATALOG_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.cat_type]);
  fprintf(file, " <parameters_name type=\"char\">%s</parameters_name>\n",
    	prefs.param_name);
  fprintf(file, " <detect_type type=\"char\">%s</detect_type>\n",
    key[findkeys("DETECT_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.detect_type]);
  fprintf(file, " <detect_minarea type=\"int\">%d</detect_minarea>\n",
	prefs.ext_minarea);
  fprintf(file, " <deblend_nthresh type=\"int\">%d</deblend_nthresh>\n",
	prefs.deblend_nthresh);
  fprintf(file, " <deblend_mincont type=\"float\">%g</deblend_mincont>\n",
	prefs.deblend_mincont);
  fprintf(file, " <filter_flag type=\"bool\">%c</filter_flag>\n",
    	prefs.filter_flag? 'T':'F');
  fprintf(file, " <filter_name type=\"char\">%s</filter_name>\n",
    	prefs.filter_name);
  fprintf(file, " <clean_flag type=\"bool\">%c</clean_flag>\n",
    	prefs.clean_flag? 'T':'F');
  fprintf(file, " <clean_param type=\"float\">%g</clean_param>\n",
	prefs.clean_param);
  fprintf(file, " <mask_type type=\"char\">%s</mask_type>\n",
	key[findkeys("MASK_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.mask_type]);
  fprintf(file, " <weight_gain type=\"bool\">%c</weight_gain>\n",
    	prefs.weightgain_flag? 'T':'F');
  fprintf(file, " <nextens type=\"int\">%d</nextens>\n", nxmlmax);
  fprintf(file, " <DETECTION_IMAGE>\n");
  fprintf(file, "  <image_name type=\"char\">%s</image_name>\n",
    	prefs.image_name[0]);
  fprintf(file, "  <thresh_type type=\"char\">%s</thresh_type>\n",
	key[findkeys("THRESH_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.thresh_type[0]]);
  fprintf(file, "  <weight_type type=\"char\">%s</weight_type>\n",
	key[findkeys("WEIGHT_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.weight_type[0]]);
  fprintf(file, "  <weight_image type=\"char\">%s</weight_image>\n",
    	prefs.wimage_name[0]);
  fprintf(file, "  <weight_thresh type=\"float\">%g</weight_thresh>\n",
	prefs.weight_thresh[0]);
  fprintf(file, " </DETECTION_IMAGE>\n");
  fprintf(file, " <MEASUREMENT_IMAGE>\n");
  fprintf(file, "  <image_name type=\"char\">%s</image_name>\n",
    	prefs.image_name[1]? prefs.image_name[1] : prefs.image_name[0]);
  fprintf(file, "  <thresh_type type=\"char\">%s</thresh_type>\n",
	key[findkeys("THRESH_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.thresh_type[1]]);
  fprintf(file, "  <weight_type type=\"char\">%s</weight_type>\n",
	key[findkeys("WEIGHT_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.weight_type[1]]);
  fprintf(file, "  <weight_image type=\"char\">%s</weight_image>\n",
    	prefs.wimage_name[1]? prefs.wimage_name[1] : prefs.wimage_name[0]);
  fprintf(file, "  <weight_thresh type=\"float\">%g</weight_thresh>\n",
	prefs.weight_thresh[1]);
  fprintf(file, " </MEASUREMENT_IMAGE>\n");
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


