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
*	Last modify:	14/07/2006
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


/****** end_xml ************************************************************
PROTO	int	end_xml(void)
PURPOSE	Free the set of meta-data kept in memory.
INPUT	-.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	12/07/2006
 ***/
int	end_xml(void)
  {
  free(xmlstack);

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
PROTO	int	write_xml(char *filename)
PURPOSE	Save meta-data to an XML file/stream.
INPUT	XML file name.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	14/07/2006
 ***/
int	write_xml(char *filename)
  {
   FILE		*file;

  if (!(file = fopen(prefs.xml_name, "w")))
    return RETURN_ERROR;

  write_xml_header(file);
  write_vo_fields(file);

  fprintf(file, "   <DATA>\n");
  if (prefs.cat_type == FITS_LDAC || prefs.cat_type == FITS_TPX
	|| prefs.cat_type == FITS_10)
    fprintf(file,
	"   <FITS extnum=\"%d\"><STREAM href=\"%s%s\" /> </FITS>",
	prefs.cat_type == FITS_10? 1:2,
	prefs.cat_name[0] == '/'? "file://" : "file:",
	prefs.cat_name);
  fprintf(file, "   </DATA>\n");
  fprintf(file, "  </TABLE>\n");

  write_xml_meta(file, (char *)NULL);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

  fclose(file);

  return RETURN_OK;
  }


/****** write_xml_header ******************************************************
PROTO	int	write_xml_header(FILE *file)
PURPOSE	Save an XML-VOtable header to an XML file/stream
INPUT	file or stream pointer.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	14/07/2006
 ***/
int	write_xml_header(FILE *file)
  {
   char		*filename, *rfilename;

/* A short, "relative" version of the filename */
  filename = prefs.image_name[prefs.nimage_name>1? 1:0];
  if (!(rfilename = strrchr(filename, '/')))
    rfilename = filename;
  else
    rfilename++;

  fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(file, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n",
	prefs.xsl_name);
  fprintf(file, "<VOTABLE "
	"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
	"xsi:noNamespaceSchemaLocation="
	"\"http://www.ivoa.net/xml/VOTable/v1.1\">\n");
  fprintf(file, "<DESCRIPTION>produced by %s</DESCRIPTION>\n", BANNER);
  fprintf(file, "<!-- VOTable description at "
	"http://www.ivoa.net/Documents/latest/VOT.html -->\n");
  fprintf(file, "<RESOURCE ID=\"%s\" name=\"%s\">\n", BANNER, rfilename);
  fprintf(file, " <DESCRIPTION>Catalog of sources extracted with %s"
	"</DESCRIPTION>\n", BANNER);
  fprintf(file, " <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n");
  fprintf(file, " <COOSYS ID=\"J2000\" equinox=\"J2000\""
	" epoch=\"J%.10g\" system=\"%s\"/>\n", prefs.epoch, prefs.coosys);
  fprintf(file, " <TABLE ID=\"Source_List\" name=\"%s/out\">\n", rfilename);
  fprintf(file,
	"  <DESCRIPTION>Table of sources detected in image</DESCRIPTION>\n");
  fprintf(file,
	"  <!-- Now comes the definition of each %s parameter -->\n", BANNER);

  return RETURN_OK;
  }


/****** write_xml_meta ********************************************************
PROTO	int	write_xml_meta(FILE *file, char *error)
PURPOSE	Save meta-data to an XML-VOTable file or stream
INPUT	Pointer to the output file (or stream),
	Pointer to an error msg (or NULL).
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	14/07/2006
 ***/
int	write_xml_meta(FILE *file, char *error)
  {
   char			*pspath,*psuser, *pshost, *str;
   struct tm		*tm;
   int			n;

/* Processing date and time if msg error present */
  if (error)
    {
    thetimet2 = time(NULL);
    tm = localtime(&thetimet2);
    sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);
    prefs.time_diff = difftime(thetimet2, thetimet);
    }

/* Username */
  psuser = pspath = pshost = NULL;
#ifdef HAVE_GETENV
  if (!(psuser=getenv("USERNAME")))	/* Cygwin,... */
    psuser = getenv("LOGNAME");		/* Linux,... */
  pspath = getenv("PWD");
  pshost = getenv("HOSTNAME");
#endif

  fprintf(file, " <RESOURCE ID=\"MetaData\" name=\"MetaData\">\n");
  fprintf(file, "  <DESCRIPTION>%s meta-data</DESCRIPTION>\n", BANNER);
  fprintf(file, "  <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n");
  fprintf(file, "  <PARAM name=\"Software\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.title;meta.software\" value=\"%s\"/>\n",
	BANNER);
  fprintf(file, "  <PARAM name=\"Version\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.version;meta.software\" value=\"%s\"/>\n",
	MYVERSION);
  fprintf(file, "  <PARAM name=\"Soft_URL\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.ref.url;meta.software\" value=\"%s\"/>\n",
	WEBSITE);
  fprintf(file, "  <PARAM name=\"Soft_Auth\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.author;meta.software\" value=\"%s\"/>\n",
	"Emmanuel Bertin");
  fprintf(file, "  <PARAM name=\"Soft_Ref\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.bibcode;meta.software\" value=\"%s\"/>\n",
	"1996A&amp;AS..117..393B");
  fprintf(file, "  <PARAM name=\"NThreads\" datatype=\"int\""
	" ucd=\"meta.number;meta.software\" value=\"%d\"/>\n",
    	prefs.nthreads);
  fprintf(file, "  <PARAM name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.sdate_end);
  fprintf(file, "  <PARAM name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.stime_end);
  fprintf(file, "  <PARAM name=\"Duration\" datatype=\"float\""
	" ucd=\"time.event;meta.software\" value=\"%.0f\" unit=\"s\"/>\n",
	prefs.time_diff);

  fprintf(file, "  <PARAM name=\"User\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	psuser);
  fprintf(file, "  <PARAM name=\"Host\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	pshost);
  fprintf(file, "  <PARAM name=\"Path\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset\" value=\"%s\"/>\n",
	pspath);

  fprintf(file,
	"  <PARAM name=\"Image_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.image;meta.fits\" value=\"%s", prefs.image_name[0]);
  if (prefs.nimage_name>1)
    fprintf(file, ",%s", prefs.image_name[1]);
  fprintf(file, "\"/>\n");

  if (error)
    {
    fprintf(file, "\n  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!! an Error occured"
	" !!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file,"  <PARAM name=\"Error_Msg\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\" value=\"%s\"/>\n", error);
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n\n");
    }

/* Meta-data for each extension */
  fprintf(file, "  <TABLE ID=\"Extension_Data\" name=\"Extension_Data\">\n");
  fprintf(file, "   <DESCRIPTION>Data gathered by %s for every FITS"
	" extension</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <!-- NExtensions may be 0"
	" if an error occurred early in the processing -->\n");
  fprintf(file, "   <PARAM name=\"NExtensions\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	nxmlmax);
  fprintf(file, "   <!-- CurrExtension may differ fromq n_extensions"
	" if an error occurred -->\n");
  fprintf(file, "   <PARAM name=\"CurrExtension\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	nxml);
  fprintf(file, "   <FIELD name=\"Extension\" datatype=\"int\""
	" ucd=\"meta.record\"/>\n");
  fprintf(file, "   <FIELD name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.record;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.record;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Duration\" datatype=\"float\""
	" ucd=\"meta.record;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"NDetect\" datatype=\"int\""
	" ucd=\"meta.number;src.sample\"/>\n");
  fprintf(file, "   <FIELD name=\"NSextracted\" datatype=\"int\""
	" ucd=\"meta.number;src.sample\"/>\n");
  fprintf(file, "   <FIELD name=\"Image_Ident\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;obs\"/>\n");
  fprintf(file, "   <FIELD name=\"Background_Mean\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.skyLevel;obs.image;stat.median\""
	" unit=\"ct\"/>\n", prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"Background_StDev\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.stdev;obs.image;stat.median\""
	" unit=\"ct\"/>\n", prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"Threshold\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.image;stat.median\""
	" unit=\"ct\"/>\n", prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"Weight_Scaling\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"arith.factor;obs.image;stat.median\"/>\n",
	prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"Pixel_Scale\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.scale;obs.image;stat.mean\""
	" unit=\"arcsec\"/>\n", prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"Epoch\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"time.epoch;obs\" unit=\"yr\"/>\n",
	prefs.nimage_name);
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (n=0; n<nxml; n++)
    if (prefs.nimage_name>1)
      fprintf(file, "     <TR>\n"
	"      <TD>%d</TD><TD>%s</TD><TD>%s</TD><TD>%.0f</TD>"
	"<TD>%d</TD><TD>%d</TD>\n"
	"      <TD>%s,%s</TD>\n"
	"      <TD>%g %g</TD><TD>%g %g</TD><TD>%g %g</TD>"
	"<TD>%g %g</TD><TD>%g %g</TD>\n"
	"     </TR>\n",
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
  fprintf(file, "  <TABLE ID=\"Warnings\" name=\"Warnings\">\n");
  fprintf(file,
	"   <DESCRIPTION>%s warnings (limited to the last %d)</DESCRIPTION>\n",
	BANNER, WARNING_NMAX);
  fprintf(file, "   <FIELD name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Msg\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\"/>\n");
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (str = warning_history(); *str; str = warning_history())
    fprintf(file, "    <TR><TD>%10.10s</TD><TD>%8.8s</TD><TD>%s</TD></TR>\n",
	str, str+11, str+22);
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* Configuration file */
  fprintf(file, "  <RESOURCE ID=\"Config\" name=\"Config\">\n");
  fprintf(file, "   <DESCRIPTION>%s configuration</DESCRIPTION>\n", BANNER);
  fprintf(file,
	"   <PARAM name=\"Command_Line\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param\" value=\"%s",
	prefs.command_line[0]);
  for (n=1; n<prefs.ncommand_line; n++)
    fprintf(file, " %s", prefs.command_line[n]);
  fprintf(file, "\"/>\n");
  fprintf(file,
	"   <PARAM name=\"Prefs_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param;meta.file\" value=\"%s\"/>\n",
	prefs.prefs_name);

  if (!error)
    {
    fprintf(file,
	"   <PARAM name=\"Catalog_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\" value=\"%s\"/>\n",
    	key[findkeys("CATALOG_TYPE",keylist,
			FIND_STRICT)].keylist[prefs.cat_type]);
    fprintf(file,
	"   <PARAM name=\"Catalog_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file\" value=\"%s\"/>\n",
	prefs.cat_name);
    fprintf(file,
	"   <PARAM name=\"Parameters_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param;meta.file\" value=\"%s\"/>\n",
	prefs.param_name);
    fprintf(file,
	"   <PARAM name=\"Detect_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;instr.det;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("DETECT_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.detect_type]);
    fprintf(file, "   <PARAM name=\"Detect_MinArea\" datatype=\"int\""
	" ucd=\"phys.area;obs.param\" value=\"%d\" unit=\"pix2\"/>\n",
    	prefs.ext_minarea);

    fprintf(file,
	"   <PARAM name=\"Thresh_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;instr.sensitivity;obs.param\" value=\"%s",
    	key[findkeys("THRESH_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.thresh_type[0]]);
    if (prefs.nthresh_type>1)
      fprintf(file, ",%s", key[findkeys("THRESH_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.thresh_type[1]]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Detect_Thresh\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.param\" value=\"%g",
	prefs.ndthresh, prefs.dthresh[0]);
    if (prefs.ndthresh>1)
      fprintf(file, " %g", prefs.dthresh[1]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Analysis_Thresh\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.param\" value=\"%g",
	prefs.nthresh, prefs.thresh[0]);
    if (prefs.nthresh>1)
      fprintf(file, " %g", prefs.thresh[1]);
    fprintf(file, "\"/>\n");

    fprintf(file,
	"   <PARAM name=\"Filter\" datatype=\"boolean\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.filter_flag? 'T':'F');
    fprintf(file,
	"   <PARAM name=\"Filter_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file;obs.param\" value=\"%s\"/>\n",
	prefs.filter_name);

    if (prefs.nfilter_thresh)
      {
      fprintf(file, "   <PARAM name=\"Filter_Thresh\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.sensitivity;obs.param\" value=\"%g",
	prefs.nfilter_thresh, prefs.filter_thresh[0]);
      if (prefs.nfilter_thresh>1)
        fprintf(file, " %g", prefs.filter_thresh[1]);
      fprintf(file, "\"/>\n");
      }

    fprintf(file, "   <PARAM name=\"Deblend_NThresh\" datatype=\"int\""
	" ucd=\"meta.number;obs.param\" value=\"%d\"/>\n",
    	prefs.deblend_nthresh);
    fprintf(file, "   <PARAM name=\"Deblend_MinCont\" datatype=\"float\""
	" ucd=\"obs.param;arith.ratio\" value=\"%g\"/>\n",
    	prefs.deblend_mincont);
    fprintf(file,
	"   <PARAM name=\"Clean\" datatype=\"boolean\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.clean_flag? 'T':'F');
    fprintf(file, "   <PARAM name=\"Clean_Param\" datatype=\"float\""
	" ucd=\"meta\" value=\"%g\"/>\n",
    	prefs.clean_param);
    fprintf(file,
	"   <PARAM name=\"Mask_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param;\" value=\"%s\"/>\n",
    	key[findkeys("MASK_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.mask_type]);

    fprintf(file,
	"   <PARAM name=\"Weight_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\" value=\"%s",
    	key[findkeys("WEIGHT_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.weight_type[0]]);
    if (prefs.nweight_type>1)
      fprintf(file, ",%s", key[findkeys("WEIGHT_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.weight_type[1]]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Weight_Thresh\" datatype=\"float\""
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
	"   <PARAM name=\"Weight_Image\" datatype=\"char\" arraysize=\"*\""
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
	"   <PARAM name=\"Weight_Gain\" datatype=\"boolean\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.weightgain_flag? 'T':'F');

    if (prefs.nimaflag)
      {
      fprintf(file,
	"   <PARAM name=\"Flag_Image\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.image;meta.fits\" value=\"%s",
    	prefs.fimage_name[0]);
      for (n=1; n<prefs.nimaflag; n++)
        fprintf(file, ",%s", prefs.fimage_name[n]);
      fprintf(file, "\"/>\n");
      fprintf(file,
	"   <PARAM name=\"Flag_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code\" value=\"%s",
    	key[findkeys("FLAG_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.flag_type[0]]);
      for (n=1; n<prefs.nimaflag; n++)
        fprintf(file, ",%s", key[findkeys("FLAG_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.flag_type[n]]);
      fprintf(file, "\"/>\n");
      }

    fprintf(file, "   <PARAM name=\"Phot_Apertures\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"obs.param\" value=\"%g",
	prefs.naper, prefs.apert[0]);
    for (n=1; n<prefs.naper; n++)
      fprintf(file, " %g", prefs.apert[n]);
    fprintf(file, "\" unit=\"pix\"/>\n");

    fprintf(file, "   <PARAM name=\"Phot_AutoParams\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"obs.param;phot\" value=\"%g",
	prefs.nautoparam, prefs.autoparam[0]);
    for (n=1; n<prefs.nautoparam; n++)
      fprintf(file, " %g", prefs.autoparam[n]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Phot_PetroParams\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"obs.param;phot\" value=\"%g",
	prefs.npetroparam, prefs.petroparam[0]);
    for (n=1; n<prefs.npetroparam; n++)
      fprintf(file, " %g", prefs.petroparam[n]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Phot_AutoApers\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"obs.param;phot\" value=\"%g",
	prefs.nautoaper, prefs.autoaper[0]);
    for (n=1; n<prefs.nautoaper; n++)
      fprintf(file, " %g", prefs.autoaper[n]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Phot_FluxFrac\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"arith.factor;obs.param;phot\" value=\"%g",
	prefs.nflux_frac, prefs.flux_frac[0]);
    for (n=1; n<prefs.nflux_frac; n++)
      fprintf(file, " %g", prefs.flux_frac[n]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Satur_Level\" datatype=\"float\""
	" ucd=\"instr.saturation;phot.count;obs.param\" value=\"%g\""
	" unit=\"ct\"/>\n", prefs.satur_level);
    fprintf(file, "   <PARAM name=\"Mag_ZeroPoint\" datatype=\"float\""
	" ucd=\"phot.calib;phot.mag;obs.param\" value=\"%g\" unit=\"mag\"/>\n",
    	prefs.mag_zeropoint);
    fprintf(file, "   <PARAM name=\"Mag_Gamma\" datatype=\"float\""
	" ucd=\"phot.calib;obs.param\" value=\"%g\"/>\n",
    	prefs.mag_gamma);
    fprintf(file, "   <PARAM name=\"Gain\" datatype=\"float\""
	" ucd=\"instr.param;obs.param\" value=\"%g\"/>\n",
    	prefs.gain);
    fprintf(file, "   <PARAM name=\"Pixel_Scale\" datatype=\"float\""
	" ucd=\"instr.scale;obs.param\" value=\"%g\" unit=\"arcsec\"/>\n",
    	prefs.pixel_scale);
    fprintf(file, "   <PARAM name=\"Seeing_FWHM\" datatype=\"float\""
	" ucd=\"instr.det.psf;stat.mean;obs.param\" value=\"%g\""
	" unit=\"pix\"/>\n", prefs.seeing_fwhm);
    fprintf(file,
	"   <PARAM name=\"StarNNW_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file;obs.param\" value=\"%s\"/>\n",
	prefs.nnw_name);

    fprintf(file, "   <PARAM name=\"Back_Size\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"obs.param\" value=\"%d",
	prefs.nbacksize, prefs.backsize[0]);
    for (n=1; n<prefs.nbacksize; n++)
      fprintf(file, " %d", prefs.backsize[n]);
    fprintf(file, "\" unit=\"pix\"/>\n");

    fprintf(file, "   <PARAM name=\"Back_FilterSize\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"obs.param\" value=\"%d",
	prefs.nbackfsize, prefs.backfsize[0]);
    for (n=1; n<prefs.nbackfsize; n++)
      fprintf(file, " %d", prefs.backfsize[n]);
    fprintf(file, "\"/>\n");

    fprintf(file,
	"   <PARAM name=\"BackPhoto_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param;\" value=\"%s\"/>\n",
    	key[findkeys("BACKPHOTO_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.pback_type]);

    fprintf(file, "   <PARAM name=\"BackPhoto_Thick\" datatype=\"int\""
	" ucd=\"obs.param\" value=\"%d\" unit=\"pix\"/>\n",
    	prefs.pback_size);

    fprintf(file, "   <PARAM name=\"Back_FiltThresh\" datatype=\"float\""
	" ucd=\"phot.count;arith.ratio;obs.param\" value=\"%g\"/>\n",
    	prefs.backfthresh);

    fprintf(file,
	"   <PARAM name=\"CheckImage_Type\" datatype=\"char\" arraysize=\"*\""
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
	"   <PARAM name=\"CheckImage_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.file\" value=\"%s",
    	prefs.check_name[0]);
    for (n=1; n<prefs.ncheck_type; n++)
      if (prefs.check_type[n] != CHECK_NONE)
        fprintf(file, ",%s", prefs.check_name[n]);
    fprintf(file, "\"/>\n");

    fprintf(file, "   <PARAM name=\"Memory_ObjStack\" datatype=\"int\""
	" ucd=\"meta.number;src;obs.param\" value=\"%d\"/>\n",
    	prefs.clean_stacksize);
    fprintf(file, "   <PARAM name=\"Memory_PixStack\" datatype=\"int\""
	" ucd=\"meta.number;obs.param\" value=\"%d\"/>\n",
    	prefs.mem_pixstack);
    fprintf(file, "   <PARAM name=\"Memory_BufSize\" datatype=\"int\""
	" ucd=\"meta.number;obs.param\" value=\"%d\"/>\n",
    	prefs.mem_bufsize);

    fprintf(file,
	"   <PARAM name=\"Assoc_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file\" value=\"%s\"/>\n",
	prefs.assoc_name);
    if (prefs.nassoc_data)
      {
      fprintf(file, "   <PARAM name=\"Assoc_Data\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.code;obs.param\" value=\"%d",
	prefs.nassoc_data, prefs.assoc_data[0]);
      for (n=1; n<prefs.nassoc_data; n++)
        fprintf(file, " %d", prefs.assoc_data[n]);
      fprintf(file, "\"/>\n");
      }
    if (prefs.nassoc_param)
      {
      fprintf(file, "   <PARAM name=\"Assoc_Params\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.code;obs.param\" value=\"%d",
	prefs.nassoc_param, prefs.assoc_param[0]);
      for (n=1; n<prefs.nassoc_param; n++)
        fprintf(file, " %d", prefs.assoc_param[n]);
      fprintf(file, "\"/>\n");
      }
    fprintf(file, "   <PARAM name=\"Assoc_Radius\" datatype=\"float\""
	" ucd=\"phys.size.radius;obs.param\" value=\"%g\" unit=\"pix\"/>\n",
    	prefs.assoc_radius);
    fprintf(file,
	"   <PARAM name=\"Assoc_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("ASSOC_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.assoc_type]);
    fprintf(file,
	"   <PARAM name=\"AssocSelec_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("ASSOCSELEC_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.assocselec_type]);

    fprintf(file,
	"   <PARAM name=\"Verbose_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code\" value=\"%s\"/>\n",
    	key[findkeys("VERBOSE_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.verbose_type]);

    fprintf(file,
	"   <PARAM name=\"FITS_Unsigned\" datatype=\"boolean\""
	" ucd=\"meta.code;obs.param\" value=\"%c\"/>\n",
    	prefs.fitsunsigned_flag? 'T':'F');

    fprintf(file,
	"   <PARAM name=\"PSF_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file;obs.param\" value=\"%s\"/>\n",
	prefs.psf_name[0]);
    fprintf(file, "   <PARAM name=\"PSF_NMax\" datatype=\"int\""
	" ucd=\"meta.number;obs.param\" value=\"%d\"/>\n",
    	prefs.psf_npsfmax);
    fprintf(file,
	"   <PARAM name=\"PSFDisplay_Type\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.code;obs.param\" value=\"%s\"/>\n",
    	key[findkeys("PSFDISPLAY_TYPE", keylist,
			FIND_STRICT)].keylist[prefs.psfdisplay_type]);

    fprintf(file,
	"   <PARAM name=\"SOM_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset;meta.file;obs.param\" value=\"%s\"/>\n",
	prefs.som_name);
    }

  fprintf(file, "  </RESOURCE>\n");
  fprintf(file, " </RESOURCE>\n");

  return RETURN_OK;
  }




/****** write_xmlerror ******************************************************
PROTO	int	write_xmlerror(char *error)
PURPOSE	Save meta-data to a simplified XML file in case of a catched error
INPUT	a character string.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	14/07/2006
 ***/
void	write_xmlerror(char *filename, char *error)
  {
   FILE			*file;

  if (!(file = fopen(filename, "w")))
    return;

  write_xml_header(file);

  fprintf(file, " </TABLE>\n");

  write_xml_meta(file, error);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

  fclose(file);

  return;
  }


