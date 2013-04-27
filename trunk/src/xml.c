/*
*				xml.c
*
* Manage XML metadata.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2006-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		27/04/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
VERSION	05/02/2010
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
  x->headflag[0] = dfield->headflag;
  x->headflag[1] = field->headflag;
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
  x->gain[0] = dfield->gain;
  x->gain[1] = field->gain;
  x->satur_level[0] = dfield->satur_level;
  x->satur_level[1] = field->satur_level;

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
VERSION	27/04/2013
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
  fprintf(file, "<VOTABLE version=\"1.1\" "
//        "xmlns=\"http://www.ivoa.net/xml/VOTable/v1.1\" "
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
VERSION	25/04/2013
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
  fprintf(file, "   <!-- CurrExtension may differ from Nextensions"
	" if an error occurred -->\n");
  fprintf(file, "   <PARAM name=\"CurrExtension\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	nxml);
  fprintf(file, "   <FIELD name=\"Extension\" datatype=\"int\""
	" ucd=\"meta.record\"/>\n");
  fprintf(file, "   <FIELD name=\"External_Header\" datatype=\"boolean\""
	" arraysize=\"%d\" ucd=\"meta.code\"/>\n",
	prefs.nimage_name);
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
  fprintf(file, "   <FIELD name=\"Gain\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.param;obs.param\"/>\n",
	prefs.nimage_name);
  fprintf(file, "   <FIELD name=\"Satur_Level\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"instr.saturation;phot.count\" unit=\"ct\"/>\n",
	prefs.nimage_name);
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (n=0; n<nxml; n++)
    if (prefs.nimage_name>1)
      fprintf(file, "     <TR>\n"
	"      <TD>%d</TD><TD>%c %c</TD><TD>%s</TD><TD>%s</TD><TD>%.0f</TD>"
	"<TD>%d</TD><TD>%d</TD>\n"
	"      <TD>%s,%s</TD><TD>%g %g</TD>\n"
	"      <TD>%g %g</TD><TD>%g %g</TD><TD>%g %g</TD>"
	"<TD>%g %g</TD><TD>%f %f</TD>\n"
	"      <TD>%g %g</TD><TD>%g %g</TD>\n"
	"     </TR>\n",
	xmlstack[n].currext,
	xmlstack[n].headflag[0]?'T':'F',xmlstack[n].headflag[1]?'T':'F',
	xmlstack[n].ext_date,
	xmlstack[n].ext_time,
	xmlstack[n].ext_elapsed,
	xmlstack[n].ndetect,
	xmlstack[n].ntotal,
	xmlstack[n].ident[0], xmlstack[n].ident[1],
	xmlstack[n].backmean[0], xmlstack[n].backmean[1],
	xmlstack[n].backsig[0], xmlstack[n].backsig[1],
	xmlstack[n].thresh[0], xmlstack[n].thresh[1],
	xmlstack[n].sigfac[0], xmlstack[n].sigfac[1],
	xmlstack[n].pixscale[0], xmlstack[n].pixscale[1],
	xmlstack[n].epoch[0], xmlstack[n].epoch[1],
	xmlstack[n].gain[0], xmlstack[n].gain[1],
	xmlstack[n].satur_level[0], xmlstack[n].satur_level[1]);
    else
      fprintf(file, "    <TR>\n"
	"     <TD>%d</TD><TD>%c</TD><TD>%s</TD><TD>%s</TD><TD>%.0f</TD>"
	"<TD>%d</TD><TD>%d</TD>\n"
	"     <TD>%s</TD><TD>%g</TD>\n"
	"     <TD>%g</TD><TD>%g</TD><TD>%g</TD><TD>%g</TD><TD>%f</TD>\n"
	"     <TD>%g</TD><TD>%g</TD>\n"
	"    </TR>\n",
	xmlstack[n].currext,
	xmlstack[n].headflag[0]?'T':'F',
	xmlstack[n].ext_date,
	xmlstack[n].ext_time,
	xmlstack[n].ext_elapsed,
	xmlstack[n].ndetect,
	xmlstack[n].ntotal,
	xmlstack[n].ident[0],
	xmlstack[n].backmean[0],
	xmlstack[n].backsig[0],
	xmlstack[n].thresh[0],
	xmlstack[n].sigfac[0],
	xmlstack[n].pixscale[0],
        xmlstack[n].epoch[0],
	xmlstack[n].gain[0],
	xmlstack[n].satur_level[0]);
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
    write_xmlconfigparam(file, "Catalog_Type", "", "meta.code;meta.table", "");
    write_xmlconfigparam(file, "Catalog_Name", "",
	"meta;meta.table;meta.file", "");
    write_xmlconfigparam(file, "Parameters_Name", "", "meta;obs.param", "");

    write_xmlconfigparam(file, "Detect_Type", "",
	"meta.code;obs.param;instr.det", "");
    write_xmlconfigparam(file, "Detect_MinArea", "pix",
	"obs.param;phys.area;stat.min", "%d");
    write_xmlconfigparam(file, "Detect_MaxArea", "pix",
	"obs.param;phys.area;stat.max", "%d");
    write_xmlconfigparam(file, "Thresh_Type", "",
	"meta.code;obs.param;instr.sensitivity", "");
    write_xmlconfigparam(file, "Detect_Thresh", "",
	"obs.param;instr.sensitivity", "%g");
    write_xmlconfigparam(file, "Analysis_Thresh", "",
	"obs.param;instr.sensitivity", "%g");
    write_xmlconfigparam(file, "Filter", "", "meta.code;obs.param", "");
    write_xmlconfigparam(file, "Filter_Name","", "obs.param;meta.file", "");
    write_xmlconfigparam(file, "Filter_Thresh", "", "obs.param", "%g");
    write_xmlconfigparam(file, "Deblend_NThresh", "",
	"meta.number;obs.param", "%d");
    write_xmlconfigparam(file, "Deblend_MinCont", "",
	"obs.param;arith.ratio", "%g");
    write_xmlconfigparam(file, "Clean", "", "meta.code;obs.param", "");
    write_xmlconfigparam(file, "Clean_Param","","obs.param;arith.ratio", "%g");
    write_xmlconfigparam(file, "Mask_Type", "", "meta.code;obs.param", "");

    write_xmlconfigparam(file, "Weight_Type", "", "meta.code;obs.param", "");
    write_xmlconfigparam(file, "Rescale_Weights","","meta.code;obs.param", "");
    write_xmlconfigparam(file, "Weight_Suffix", "",
	"meta;meta.file;meta.fits", "");
    write_xmlconfigparam(file, "Weight_Thresh", "", "obs.param", "%g");
    write_xmlconfigparam(file, "Weight_Image", "",
	"stat.weight;meta.file;meta.fits", "");
    write_xmlconfigparam(file, "Weight_Gain", "", "meta.code;obs.param", "");

    write_xmlconfigparam(file, "Flag_Image", "",
	"meta.code;meta.file;meta.fits", "");
    write_xmlconfigparam(file, "Flag_Type", "", "meta.code;obs.param", "");

    write_xmlconfigparam(file, "Phot_Apertures", "pix",
	"obs.param;phys.size.diameter", "%g");
    write_xmlconfigparam(file, "Phot_AutoParams", "",
	"meta.code;obs.param", "%g");
    write_xmlconfigparam(file, "Phot_PetroParams", "",
	"meta.code;obs.param", "%g");
    write_xmlconfigparam(file, "Phot_AutoApers", "pix",
	"obs.param;phys.size.diameter", "%g");
    write_xmlconfigparam(file, "Phot_FluxFrac", "",
	"obs.param;arith.ratio", "%g");
    write_xmlconfigparam(file, "Satur_Level", "ct",
	"obs.param;instr.saturation", "%g");
    write_xmlconfigparam(file, "Satur_Key", "",
	"meta.code;obs.param;instr.saturation", "");
    write_xmlconfigparam(file, "Mag_ZeroPoint", "mag",
	"obs.param;phot.calib;phot.mag", "%.4f");
    write_xmlconfigparam(file, "Mag_Gamma", "",
	"obs.param;phot.calib;instr.plate.emulsion", "%.2f");
    write_xmlconfigparam(file, "Gain", "/ct", "obs.param;instr.param", "%g");
    write_xmlconfigparam(file, "Gain_Key", "",
	"meta.code;obs.param;instr.param", "");

    write_xmlconfigparam(file, "Pixel_Scale", "arcsec/pix",
	"obs.param;instr.scale", "%g");
    write_xmlconfigparam(file, "Seeing_FWHM", "arcsec",
	"obs.param;instr.obsty.seeing", "%g");
    write_xmlconfigparam(file, "StarNNW_Name", "",
	"obs.param;meta.dataset;meta.file", "");

    write_xmlconfigparam(file, "Back_Type", "", "meta.code;obs.param", "");
    write_xmlconfigparam(file, "Back_Value", "ct", "obs.param", "%g");
    write_xmlconfigparam(file, "Back_Size", "pix", "obs.param", "%d");
    write_xmlconfigparam(file, "Back_FilterSize", "", "obs.param", "%d");
    write_xmlconfigparam(file, "BackPhoto_Type", "","meta.code;obs.param", "");
    write_xmlconfigparam(file, "BackPhoto_Thick", "pix", "obs.param", "%d");
    write_xmlconfigparam(file, "Back_FiltThresh", "ct", "obs.param;", "%g");

    write_xmlconfigparam(file, "CheckImage_Type","","meta.code;obs.param", "");
    write_xmlconfigparam(file, "CheckImage_Name", "",
	"meta.id;meta.file;meta.fits", "");

    write_xmlconfigparam(file, "Memory_ObjStack", "",
	"meta.number;src;obs.param", "%d");
    write_xmlconfigparam(file, "Memory_PixStack", "pix",
	"meta.number;obs.param", "%d");
    write_xmlconfigparam(file, "Memory_BufSize", "pix",
	"meta.number;obs.param", "%d");

    write_xmlconfigparam(file, "Assoc_Name", "",
	"meta;meta.table;meta.file", "");
    write_xmlconfigparam(file, "Assoc_Data", "", "meta;meta.table", "%d");
    write_xmlconfigparam(file, "Assoc_Params", "", "meta;meta.table", "%d");
    write_xmlconfigparam(file, "AssocCoord_Type", "", "meta.code", "");
    write_xmlconfigparam(file, "Assoc_Radius", "pix",
	"meta;phys.size.radius", "");
    write_xmlconfigparam(file, "Assoc_Type", "", "meta.code", "");
    write_xmlconfigparam(file, "AssocSelec_Type", "", "meta.code", "");

    write_xmlconfigparam(file, "Verbose_Type", "", "meta.code", "");
    write_xmlconfigparam(file, "Header_Suffix", "", "meta.id;meta.file", "");
    write_xmlconfigparam(file, "Write_XML", "", "meta.code", "");
    write_xmlconfigparam(file, "XML_Name", "", "meta;meta.file", "");
    write_xmlconfigparam(file, "XSL_URL", "", "meta.ref.url;meta.file", "");
    write_xmlconfigparam(file, "NThreads", "", "meta.number", "%d");
    write_xmlconfigparam(file, "FITS_Unsigned", "", "meta.code;obs.param", "");

    write_xmlconfigparam(file, "PSF_Name", "",
	"instr.det.psf;meta.file;meta.fits", "");
    }

  fprintf(file, "  </RESOURCE>\n");
  fprintf(file, " </RESOURCE>\n");

  return RETURN_OK;
  }


/****** write_xmlconfigparam **************************************************
PROTO	int write_xmlconfigparam(FILE *file, char *name, char *unit,
		char *ucd, char *format)
PURPOSE	Write to a VO-table the configuration parameters.
INPUT	Output stream (file) pointer,
	Name of the parameter keyword,
	unit,
	UCD string,
	printf() format to use in "value".
OUTPUT	RETURN_OK if the keyword exists, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/04/2013
 ***/
int	write_xmlconfigparam(FILE *file, char *name, char *unit,
		 char *ucd, char *format)
  {
   char		value[MAXCHAR], uunit[MAXCHAR];
   int		i,j,n;

  for (i=0; key[i].name[0] && cistrcmp(name, key[i].name, FIND_STRICT); i++);
  if (!key[i].name[0])
    return RETURN_ERROR;

  if (*unit)
    sprintf(uunit, " unit=\"%s\"", unit);
  else
    *uunit = '\0';
  switch(key[i].type)
    {
    case P_FLOAT:
      sprintf(value, format, *((double *)key[i].ptr));
      fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, uunit, ucd, value);
      break;
    case P_FLOATLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, format, ((double *)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
		" arraysize=\"%d\" ucd=\"%s\" value=\"%s",
		name, uunit, n, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, format, ((double *)key[i].ptr)[j]);
          fprintf(file, " %s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
		" ucd=\"%s\" value=\"\"/>\n",
		name, uunit, ucd);
      break;
    case P_INT:
      sprintf(value, format, *((int *)key[i].ptr));
      fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"int\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, uunit, ucd, value);
      break;
    case P_INTLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, format, ((int *)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"int\""
		" arraysize=\"%d\" ucd=\"%s\" value=\"%s",
		name, uunit, n, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, format, ((int *)key[i].ptr)[j]);
          fprintf(file, " %s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
		" ucd=\"%s\" value=\"\"/>\n",
		name, uunit, ucd);
      break;
    case P_BOOL:
      sprintf(value, "%c", *((int *)key[i].ptr)? 'T':'F');
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"boolean\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, value);
      break;
    case P_BOOLLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, "%c", ((int *)key[i].ptr)[0]? 'T':'F');
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"boolean\""
		" arraysize=\"%d\" ucd=\"%s\" value=\"%s",
		name, n, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, "%c", ((int *)key[i].ptr)[j]? 'T':'F');
          fprintf(file, " %s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"boolean\""
		" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    case P_STRING:
      sprintf(value, (char *)key[i].ptr);
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, *value? value: " ");
      break;
    case P_STRINGLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, ((char **)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"%s",
		name, ucd, *value? value: " ");
        for (j=1; j<n; j++)
          {
          sprintf(value, ((char **)key[i].ptr)[j]);
          fprintf(file, ",%s", *value? value: " ");
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    case P_KEY:
      sprintf(value, key[i].keylist[*((int *)key[i].ptr)]);
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, value);
      break;
    case P_KEYLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, key[i].keylist[((int *)key[i].ptr)[0]]);
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"%s",
		name, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, key[i].keylist[((int *)key[i].ptr)[j]]);
          fprintf(file, ",%s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    default:
        error(EXIT_FAILURE, "*Internal Error*: Type Unknown",
		" in write_xmlconfigparam()");
    }

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


