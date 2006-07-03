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
*	Last modify:	03/07/2006
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
#include "prefs.h"
#include "xml.h"

xmlstruct	*xmlstack;
int		nxml, nxmlmax;

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
VERSION	03/07/2006
 ***/
int	update_xml(void)
  {
  if (nxml < nxmlmax)
    xmlstack[nxml++].sexcat = thecat;
  else
    error(EXIT_FAILURE, "*Internal Error*: too many extensions in XML stack",
			"");

  return EXIT_SUCCESS;
  }

/****** write_xml ************************************************************
PROTO	int	write_xml(fieldstruct **fields, int nfield)
			fgroupstruct **fgroups, int ngroup);
PURPOSE	Save meta-data to XML files
INPUT	Pointer to the array of fields,
	number of fields,
	pointer to the array of field groups,
	number of field groups,
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	10/03/2006
 ***/
int	write_xml(void)
  {
   FILE			*file;
   time_t		thetime;
   char			*pspath,*psuser, *pshost;
   int			f;

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
  fprintf(file, "<SEXTRACTION>\n");
  fprintf(file, " <software type=\"char\">%s</software>\n", BANNER);
  fprintf(file, " <version type=\"char\">%s</version>\n", MYVERSION);
  fprintf(file, " <date type=\"char\">%s</date>\n", prefs.sdate_end);
  fprintf(file, " <time type=\"char\">%s</time>\n", prefs.stime_end);
  fprintf(file, " <duration type=\"int\" unit=\"s\">%d</duration>\n",
    	prefs.time_diff);
  fprintf(file, " <nthreads type=\"int\">%d</nthreads>\n",
    	prefs.nthreads);
  fprintf(file, "<user type=\"char\">%s</user>\n", psuser);
  fprintf(file, "<host type=\"char\">%s</host>\n", pshost);
  fprintf(file, "<path type=\"char\">%s</path>\n", pspath);
  fprintf(file, " <filter_flag type=\"bool\">%c</filter_flag>\n",
    	prefs.filter_flag? 'T':'F');
  fprintf(file, "  <nextens type=\"int\">%d</nextens>\n", nxmlmax);
  for (f=0; f<nxmlmax; f++)
    {
    fprintf(file, " <IMAGE_PROPS>\n");
    fprintf(file, " </IMAGE_PROPS>\n");
    }
  fprintf(file, "</SEXTRACTION>\n");
  fclose(file);

  free(xmlstack);

  return RETURN_OK;
  }


