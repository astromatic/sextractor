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
*	Last modify:	11/07/2006
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

extern time_t		thetimet,thetimet2;	/* from makeit.c */
extern pkeystruct	key[];			/* from preflist.h */
extern char		keylist[][32];		/* from preflist.h */
xmlstruct		*xmlstack = NULL;
int			nxml=0, nxmlmax=0;

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
VERSION	11/07/2006
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
  strcpy(x->ident[0], dfield->ident); 
  strcpy(x->ident[1], field->ident); 
  x->backmean[0] = dfield->backmean;
  x->backmean[1] = field->backmean;
  x->backsig[0] = dfield->backsig;
  x->backsig[1] = field->backsig;
  x->sigfac[0] = dfield->sigfac;
  x->sigfac[1] = field->sigfac;
  x->thresh[0] = dfield->dthresh;
  x->thresh[1] = field->thresh;
  x->pixscale[0] = dfield->pixscale;
  x->pixscale[1] = field->pixscale;
  x->epoch[0] = dfield->epoch;
  x->epoch[1] = field->epoch;

  return EXIT_SUCCESS;
  }

/****** write_xml ************************************************************
PROTO	int	write_xml(void)
PURPOSE	Save meta-data to XML files
INPUT	-.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/07/2006
 ***/
int	write_xml(void)
  {
   FILE			*file;
   char			*pspath,*psuser, *pshost, *str;
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
/*
  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n",
	prefs.xsl_name);
  fprintf(file, "<!DOCTYPE INFO SYSTEM \"instfile.dtd\">\n");
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
  fprintf(file, "  <prefs_name type=\"char\">%s</prefs_name>\n",
    	prefs.prefs_name);
  fprintf(file, "  <command_line type=\"char\">%s", prefs.command_line[0]);
  for (n=1; n<prefs.ncommand_line; n++)
    fprintf(file, " %s", prefs.command_line[n]);
  fprintf(file, "</command_line>\n");
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
    fprintf(file, "  <date type=\"char\">%s</date>\n", xmlstack[n].ext_date);
    fprintf(file, "  <time type=\"char\">%s</time>\n", xmlstack[n].ext_time);
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

  for (str = warning_history(); *str; str = warning_history())
    fprintf(file, " <WARNING>%s</WARNING>\n", str);

  fprintf(file, "</SOURCE_EXTRACTION>\n");
*/

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n",
	prefs.xsl_name);
  fprintf(file, "<VOTABLE "
	"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
	"xsi:noNamespaceSchemaLocation="
	"\"xmlns=http://www.ivoa.net/xml/VOTable/v1.1\">\n");
  fprintf(file, " <RESOURCE name=\"source_extraction\"/>\n");
  fprintf(file, "  <DESCRIPTION>%s meta-data</DESCRIPTION>\n", BANNER);
  fprintf(file, "  <PARAM name=\"software\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.title;meta.software\" value=\"%s\"/>\n",
	BANNER);
  fprintf(file, "  <PARAM name=\"version\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.version;meta.software\" value=\"%s\"/>\n",
	MYVERSION);
  fprintf(file, "  <PARAM name=\"soft_url\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.ref.url;meta.software\" value=\"%s\"/>\n",
	WEBSITE);
  fprintf(file, "  <PARAM name=\"soft_auth\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.author;meta.software\" value=\"%s\"/>\n",
	"Emmanuel Bertin");
  fprintf(file, "  <PARAM name=\"soft_ref\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.bibcode;meta.software\" value=\"%s\"/>\n",
	"1996A&AS..117..393B");
  fprintf(file, "  <PARAM name=\"nthreads\" datatype=\"int\""
	" ucd=\"meta.number;meta.software\" value=\"%d\"/>\n",
    	prefs.nthreads);
  fprintf(file, "  <PARAM name=\"date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.sdate_end);
  fprintf(file, "  <PARAM name=\"time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.stime_end);
  fprintf(file, "  <PARAM name=\"duration\" datatype=\"float\""
	" ucd=\"time.event;meta.software\" value=\"%.0f\" unit=\"s\"/>\n",
	prefs.time_diff);

  fprintf(file, "  <PARAM name=\"user\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	psuser);
  fprintf(file, "  <PARAM name=\"host\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	pshost);
  fprintf(file, "  <PARAM name=\"path\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset\" value=\"%s\"/>\n",
	pspath);

  fprintf(file, "  <RESOURCE name=\"config\">\n");
  fprintf(file, "   <DESCRIPTION>%s configuration</DESCRIPTION>\n", BANNER);
  fprintf(file,
	"   <PARAM name=\"command_line\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param\" value=\"%s",
	prefs.command_line[0]);
  for (n=1; n<prefs.ncommand_line; n++)
    fprintf(file, " %s", prefs.command_line[n]);
  fprintf(file, "\"/>\n");
  fprintf(file,
	"   <PARAM name=\"prefs_name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param;meta.file\" value=\"%s\"/>\n",
	prefs.prefs_name);
  fprintf(file,
	"   <PARAM name=\"catalog_type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\" value=\"%s\"/>\n",
    	key[findkeys("CATALOG_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.cat_type]);
  fprintf(file,
	"   <PARAM name=\"catalog_name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file\" value=\"%s\"/>\n",
	prefs.cat_name);
  fprintf(file,
	"   <PARAM name=\"parameters_name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param;meta.file\" value=\"%s\"/>\n",
	prefs.param_name);
  fprintf(file,
	"   <PARAM name=\"detect_type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;instr.det;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("DETECT_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.detect_type]);
  fprintf(file, "   <PARAM name=\"detect_minarea\" datatype=\"int\""
	" ucd=\"phys.area;obs.param\" value=\"%d\" unit=\"pix2\"/>\n",
    	prefs.ext_minarea);

  fprintf(file,
	"   <PARAM name=\"thresh_type\" datatype=\"char\" arraysize=\"*s,\""
	" ucd=\"meta.code;instr.sensitivity;obs.param\" value=\"%s",
    	key[findkeys("THRESH_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.thresh_type[0]]);
  if (prefs.nthresh_type>1)
    fprintf(file, ",%s", key[findkeys("THRESH_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.thresh_type[1]]);
  fprintf(file, "\"/>\n");

  fprintf(file, "   <PARAM name=\"detect_thresh\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.param\" value=\"%g",
	prefs.ndthresh, prefs.dthresh[0]);
  if (prefs.ndthresh>1)
    fprintf(file, " %g", prefs.dthresh[1]);
  fprintf(file, "\"/>\n");

  fprintf(file, "   <PARAM name=\"analysis_thresh\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.param\" value=\"%g",
	prefs.nthresh, prefs.thresh[0]);
  if (prefs.nthresh>1)
    fprintf(file, " %g", prefs.thresh[1]);
  fprintf(file, "\"/>\n");

  fprintf(file,
	"   <PARAM name=\"filter\" datatype=\"bool\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.filter_flag? 'T':'F');
  fprintf(file,
	"   <PARAM name=\"filter_name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file;obs.param\" value=\"%s\"/>\n",
	prefs.filter_name);

  if (prefs.nfilter_thresh)
    {
    fprintf(file, "   <PARAM name=\"filter_thresh\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.param\" value=\"%g",
	prefs.nfilter_thresh, prefs.filter_thresh[0]);
    if (prefs.nfilter_thresh>1)
      fprintf(file, " %g", prefs.filter_thresh[1]);
    fprintf(file, "\"/>\n");
    }

  fprintf(file, "   <PARAM name=\"deblend_nthresh\" datatype=\"int\""
	" ucd=\"meta.number;obs.param\" value=\"%d\"/>\n",
    	prefs.deblend_nthresh);
  fprintf(file, "   <PARAM name=\"deblend_mincont\" datatype=\"float\""
	" ucd=\"obs.param;arith.ratio\" value=\"%g\"/>\n",
    	prefs.deblend_mincont);
  fprintf(file,
	"   <PARAM name=\"clean\" datatype=\"bool\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.clean_flag? 'T':'F');
  fprintf(file, "   <PARAM name=\"clean_param\" datatype=\"float\""
	" ucd=\"meta\" value=\"%g\"/>\n",
    	prefs.clean_param);
  fprintf(file,
	"   <PARAM name=\"mask_type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param;\" value=\"%s\"/>\n",
    	key[findkeys("MASK_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.mask_type]);

  fprintf(file,
	"   <PARAM name=\"weight_type\" datatype=\"char\" arraysize=\"*s,\""
	" ucd=\"meta.code;obs.param\" value=\"%s",
    	key[findkeys("WEIGHT_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.weight_type[0]]);
  if (prefs.nweight_type>1)
    fprintf(file, ",%s", key[findkeys("WEIGHT_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.weight_type[1]]);
  fprintf(file, "\"/>\n");

  fprintf(file, "   <PARAM name=\"weight_thresh\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.param\" value=\"%g",
	prefs.nweight_thresh, prefs.weight_thresh[0]);
  if (prefs.nweight_thresh>1)
    fprintf(file, " %g", prefs.weight_thresh[1]);
  fprintf(file, "\"/>\n");

  if ((prefs.weight_type[0] != WEIGHT_NONE
		&& prefs.weight_type[0] != WEIGHT_FROMBACK)
	|| (prefs.weight_type[1] != WEIGHT_NONE
		&& prefs.weight_type[1] != WEIGHT_FROMBACK))
    {
    fprintf(file,
	"   <PARAM name=\"weight_image\" datatype=\"char\" arraysize=\"*s,\""
	" ucd=\"obs.image;meta.fits;obs.param\" value=\"%s",
    	(prefs.weight_type[0] != WEIGHT_NONE
	&& prefs.weight_type[0] != WEIGHT_FROMBACK) ?
		prefs.wimage_name[0] : NULL);
    if (prefs.weight_type[1] != WEIGHT_NONE
		&& prefs.weight_type[1] != WEIGHT_FROMBACK)
      fprintf(file, ",%s", prefs.wimage_name[1]);
    fprintf(file, "\"/>\n");
    }

  fprintf(file,
	"   <PARAM name=\"weight_gain\" datatype=\"bool\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.weightgain_flag? 'T':'F');

  if (prefs.nimaflag)
    {
    fprintf(file,
	"   <PARAM name=\"flag_image\" datatype=\"char\" arraysize=\"*s,\""
	" ucd=\"obs.image;meta.fits\" value=\"%s",
    	prefs.fimage_name[0]);
    for (n=1; n<prefs.nimaflag; n++)
      fprintf(file, ",%s", prefs.fimage_name[n]);
    fprintf(file, "\"/>\n");
    fprintf(file,
	"   <PARAM name=\"flag_type\" datatype=\"char\" arraysize=\"*s,\""
	" ucd=\"meta.code\" value=\"%s",
    	key[findkeys("FLAG_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.flag_type[0]]);
    for (n=1; n<prefs.nimaflag; n++)
      fprintf(file, ",%s", key[findkeys("FLAG_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.flag_type[n]]);
    fprintf(file, "\"/>\n");
    }

  fprintf(file, "   <PARAM name=\"phot_apertures\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"obs.param\" value=\"%g",
	prefs.naper, prefs.apert[0]);
  for (n=1; n<prefs.naper; n++)
    fprintf(file, " %g", prefs.apert[n]);
  fprintf(file, "\" unit=\"pix\"/>\n");

  fprintf(file, "   <PARAM name=\"phot_autoparams\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"obs.param;phot\" value=\"%g",
	prefs.nautoparam, prefs.autoparam[0]);
  for (n=1; n<prefs.nautoparam; n++)
    fprintf(file, " %g", prefs.autoparam[n]);
  fprintf(file, "\"/>\n");

  fprintf(file, "   <PARAM name=\"phot_petroparams\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"obs.param;phot\" value=\"%g",
	prefs.npetroparam, prefs.petroparam[0]);
  for (n=1; n<prefs.npetroparam; n++)
    fprintf(file, " %g", prefs.petroparam[n]);
  fprintf(file, "\"/>\n");

  fprintf(file, "   <PARAM name=\"phot_autoapers\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"obs.param;phot\" value=\"%g",
	prefs.nautoaper, prefs.autoaper[0]);
  for (n=1; n<prefs.nautoaper; n++)
    fprintf(file, " %g", prefs.autoaper[n]);
  fprintf(file, "\"/>\n");

  fprintf(file, "   <PARAM name=\"phot_fluxfrac\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"arith.factor;obs.param;phot\" value=\"%g",
	prefs.nflux_frac, prefs.flux_frac[0]);
  for (n=1; n<prefs.nflux_frac; n++)
    fprintf(file, " %g", prefs.flux_frac[n]);
  fprintf(file, "\"/>\n");

  fprintf(file, "   <PARAM name=\"satur_level\" datatype=\"float\""
	" ucd=\"instr.saturation;phot.count;obs.param\" value=\"%g\""
	" unit=\"ct\"/>\n", prefs.satur_level);
  fprintf(file, "   <PARAM name=\"mag_zeropoint\" datatype=\"float\""
	" ucd=\"phot.calib;phot.mag;obs.param\" value=\"%g\" unit=\"mag\"/>\n",
    	prefs.mag_zeropoint);
  fprintf(file, "   <PARAM name=\"mag_gamma\" datatype=\"float\""
	" ucd=\"phot.calib;obs.param\" value=\"%g\"/>\n",
    	prefs.mag_gamma);
  fprintf(file, "   <PARAM name=\"gain\" datatype=\"float\""
	" ucd=\"instr.param;obs.param\" value=\"%g\"/>\n",
    	prefs.gain);
  fprintf(file, "   <PARAM name=\"pixel_scale\" datatype=\"float\""
	" ucd=\"instr.scale;obs.param\" value=\"%g\" unit=\"arcsec\"/>\n",
    	prefs.pixel_scale);
  fprintf(file, "   <PARAM name=\"seeing_fwhm\" datatype=\"float\""
	" ucd=\"instr.det.psf;stat.mean;obs.param\" value=\"%g\""
	" unit=\"pix\"/>\n", prefs.seeing_fwhm);
  fprintf(file,
	"   <PARAM name=\"starnnw_name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file;obs.param\" value=\"%s\"/>\n",
	prefs.nnw_name);

  fprintf(file, "   <PARAM name=\"back_size\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"obs.param\" value=\"%d",
	prefs.nbacksize, prefs.backsize[0]);
  for (n=1; n<prefs.nbacksize; n++)
    fprintf(file, " %d", prefs.backsize[n]);
  fprintf(file, "\" unit=\"pix\"/>\n");

  fprintf(file, "   <PARAM name=\"back_filtersize\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"obs.param\" value=\"%d",
	prefs.nbackfsize, prefs.backfsize[0]);
  for (n=1; n<prefs.nbackfsize; n++)
    fprintf(file, " %d", prefs.backfsize[n]);
  fprintf(file, "\"/>\n");

  fprintf(file,
	"   <PARAM name=\"backphoto_type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param;\" value=\"%s\"/>\n",
    	key[findkeys("BACKPHOTO_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.pback_type]);

  fprintf(file, "   <PARAM name=\"backphoto_thick\" datatype=\"int\""
	" ucd=\"obs.param\" value=\"%d\" unit=\"pix\"/>\n",
    	prefs.pback_size);

  fprintf(file, "   <PARAM name=\"back_filtthresh\" datatype=\"float\""
	" ucd=\"phot.count;arith.ratio;obs.param\" value=\"%g\"/>\n",
    	prefs.backfthresh);

  fprintf(file,
	"   <PARAM name=\"checkimage_type\" datatype=\"char\" arraysize=\"*s,\""
	" ucd=\"meta.code\" value=\"%s",
    	key[findkeys("CHECKIMAGE_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.check_type[0]]);
  for (n=1; n<prefs.ncheck_type; n++)
    fprintf(file,
	",%s",
    	key[findkeys("CHECKIMAGE_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.check_type[n]]);
  fprintf(file, "\"/>\n");

  fprintf(file,
	"   <PARAM name=\"checkimage_name\" datatype=\"char\" arraysize=\"*s,\""
	" ucd=\"meta.file\" value=\"%s",
    	prefs.check_name[0]);
  for (n=1; n<prefs.ncheck_type; n++)
    if (prefs.check_type[n] != CHECK_NONE)
      fprintf(file,
	",%s", prefs.check_name[n]);
  fprintf(file, "\"/>\n");

  fprintf(file, "   <PARAM name=\"memory_objstack\" datatype=\"int\""
	" ucd=\"meta.number;src;obs.param\" value=\"%d\"/>\n",
    	prefs.clean_stacksize);
  fprintf(file, "   <PARAM name=\"memory_pixstack\" datatype=\"int\""
	" ucd=\"meta.number;obs.param\" value=\"%d\"/>\n",
    	prefs.mem_pixstack);
  fprintf(file, "   <PARAM name=\"memory_bufsize\" datatype=\"int\""
	" ucd=\"meta.number;obs.param\" value=\"%d\"/>\n",
    	prefs.mem_bufsize);

  fprintf(file,
	"   <PARAM name=\"assoc_name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file\" value=\"%s\"/>\n",
	prefs.assoc_name);
  if (prefs.nassoc_data)
    {
    fprintf(file, "   <PARAM name=\"assoc_data\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.code;obs.param\" value=\"%d",
	prefs.nassoc_data, prefs.assoc_data[0]);
    for (n=1; n<prefs.nassoc_data; n++)
      fprintf(file, " %d", prefs.assoc_data[n]);
    fprintf(file, "\"/>\n");
    }
  if (prefs.nassoc_param)
    {
    fprintf(file, "   <PARAM name=\"assoc_data\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.code;obs.param\" value=\"%d",
	prefs.nassoc_data, prefs.assoc_data[0]);
    for (n=1; n<prefs.nassoc_data; n++)
      fprintf(file, " %d", prefs.assoc_data[n]);
    fprintf(file, "\"/>\n");
    }
  fprintf(file, "   <PARAM name=\"assoc_radius\" datatype=\"float\""
	" ucd=\"phys.size.radius;obs.param\" value=\"%g\" unit=\"pix\"/>\n",
    	prefs.assoc_radius);
  fprintf(file,
	"   <PARAM name=\"assoc_type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("ASSOC_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.assoc_type]);
  fprintf(file,
	"   <PARAM name=\"assocselec_type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("ASSOCSELEC_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.assocselec_type]);

  fprintf(file,
	"   <PARAM name=\"verbose_type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code\" value=\"%s\"/>\n",
    	key[findkeys("VERBOSE_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.verbose_type]);

  fprintf(file,
	"   <PARAM name=\"fits_unsigned\" datatype=\"bool\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.fitsunsigned_flag? 'T':'F');

  fprintf(file,
	"   <PARAM name=\"psf_name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file;obs.param\" value=\"%s\"/>\n",
	prefs.psf_name[0]);
  fprintf(file, "   <PARAM name=\"psf_nmax\" datatype=\"int\""
	" ucd=\"meta.number;obs.param\" value=\"%d\"/>\n",
    	prefs.psf_npsfmax);
  fprintf(file,
	"   <PARAM name=\"psfdisplay_type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("PSFDISPLAY_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.psfdisplay_type]);

  fprintf(file,
	"   <PARAM name=\"som_name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file;obs.param\" value=\"%s\"/>\n",
	prefs.som_name);

  fprintf(file, "  </RESOURCE>\n");

/* Meta-data for each extension */
  fprintf(file, "  <TABLE name=\"extension_data\">\n");
  fprintf(file, "   <DESCRIPTION>Data gathered by %s for every FITS"
	" extension</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <PARAM name=\"nextensions\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	nxmlmax);
  fprintf(file, "   <FIELD name=\"extension\" datatype=\"int\""
	" ucd=\"meta.record\"/>\n");
  fprintf(file, "   <FIELD name=\"date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.record;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.record;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"duration\" datatype=\"float\""
	" ucd=\"meta.record;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"ndetect\" datatype=\"int\""
	" ucd=\"meta.number;src.sample\"/>\n");
  fprintf(file, "   <FIELD name=\"nsextracted\" datatype=\"int\""
	" ucd=\"meta.number;src.sample\"/>\n");
  fprintf(file, "   <FIELD name=\"image_ident\" datatype=\"char\""
	" arraysize=\"*s,\" ucd=\"meta.id;obs\"/>\n");
  fprintf(file, "   <FIELD name=\"background_mean\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.skyLevel;obs.image;stat.median\""
	" unit=\"ct\"/>\n", prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"background_stdev\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.stdev;obs.image;stat.median\""
	" unit=\"ct\"/>\n", prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"threshold\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.image;stat.median\""
	" unit=\"ct\"/>\n", prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"weight_scaling\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"arith.factor;obs.image;stat.median\"/>\n",
	prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"pixel_scale\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.scale;obs.image;stat.mean\""
	" unit=\"arcsec\"/>\n", prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"epoch\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"time.epoch;obs\" unit=\"yr\"/>\n",
	prefs.nimage_name);
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (n=0; n<nxmlmax; n++)
    if (prefs.nimage_name>1)
      fprintf(file, "    <TR>\n"
	"     <TD>%d</TD><TD>%s</TD><TD>%s</TD><TD>%.0f</TD>"
	"<TD>%d</TD><TD>%d</TD>\n"
	"     <TD>%s,%s</TD>\n"
	"     <TD>%g %g</TD><TD>%g %g</TD><TD>%g %g</TD>"
	"<TD>%g %g</TD><TD>%g %g</TD>\n"
	"    </TR>\n",
	xmlstack[n].currext,
	xmlstack[n].ext_date,
	xmlstack[n].ext_time,
	xmlstack[n].ext_elapsed,
	xmlstack[n].ndetect,
	xmlstack[n].ntotal,
	xmlstack[n].ident[0], xmlstack[n].ident[1],
	xmlstack[n].backmean[0], xmlstack[n].backmean[1],
	xmlstack[n].backsig[0], xmlstack[n].backsig[1],
	xmlstack[n].sigfac[0], xmlstack[n].sigfac[1],
	xmlstack[n].thresh[0], xmlstack[n].thresh[1],
	xmlstack[n].pixscale[0], xmlstack[n].pixscale[1]);
    else
      fprintf(file, "    <TR>\n"
	"     <TD>%d</TD><TD>%s</TD><TD>%s</TD><TD>%.0f</TD>"
	"<TD>%d</TD><TD>%d</TD>\n"
	"     <TD>%s</TD>\n"
	"     <TD>%g</TD><TD>%g</TD><TD>%g</TD><TD>%g</TD><TD>%g</TD>\n"
	"    </TR>\n",
	xmlstack[n].currext,
	xmlstack[n].ext_date,
	xmlstack[n].ext_time,
	xmlstack[n].ext_elapsed,
	xmlstack[n].ndetect,
	xmlstack[n].ntotal,
	xmlstack[n].ident[0],
	xmlstack[n].backmean[0],
	xmlstack[n].backsig[0],
	xmlstack[n].sigfac[0],
	xmlstack[n].thresh[0],
	xmlstack[n].pixscale[0]);
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* Warnings */
  fprintf(file, "  <TABLE name=\"warnings\">\n");
  fprintf(file,
	"   <DESCRIPTION>%s warnings (limited to the last %d)</DESCRIPTION>\n",
	BANNER, WARNING_NMAX);
  fprintf(file, "   <FIELD name=\"date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"msg\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\"/>\n");
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (str = warning_history(); *str; str = warning_history())
    fprintf(file, "    <TR><TD>%10.10s</TD><TD>%8.8s</TD><TD>%s</TD></TR>\n",
	str, str+11, str+22);
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");
  fprintf(file, " </RESOURCE>\n");
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
VERSION	10/07/2006
 ***/
void	write_xmlerror(char *msg1, char *msg2)
  {
   FILE			*file;
   struct tm		*tm;
   char			*pspath,*psuser, *pshost, *str;

  if (!prefs.xml_flag)
    return;
/* Processing date and time */
  thetimet2 = time(NULL);
  tm = localtime(&thetimet2);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);
  prefs.time_diff = difftime(thetimet2, thetimet);

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
  fprintf(file, " <nextens type=\"int\">%d</nextens>\n", nxmlmax);
  fprintf(file, " <currextens type=\"int\">%d</currextens>\n",
	nxml<nxmlmax? nxml+1 : nxml);
  for (str = warning_history(); *str; str = warning_history())
    fprintf(file, " <WARNING>%s</WARNING>\n", str);
  fprintf(file, " <ERROR_MSG type=\"char\">%s%s</ERROR_MSG>\n", msg1,msg2);
  fprintf(file, "</SOURCE_EXTRACTION>\n");
  fclose(file);

  free(xmlstack);

  return;
  }

