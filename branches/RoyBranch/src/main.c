/*
*				main.c
*
* Command line parsing.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		04/06/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<ctype.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"prefs.h"
#include "pattern.h"
#define		SYNTAX \
EXECUTABLE " <image> [<image2>][-c <configuration_file>][-<keyword> <value>]\n" \
"> to dump a default configuration file:          " EXECUTABLE " -d \n" \
"> to dump a default extended configuration file: " EXECUTABLE " -dd \n" \
"> to dump a full list of measurement parameters: " EXECUTABLE " -dp \n"

extern const char       notokstr[];
extern keystruct	objkey[];

/********************************** main ************************************/
int	main(int argc, char *argv[])

  {
   double	tdiff, lines, dets;
   int		a, narg, nim, opt, opt2;
   char		str[MAXCHARL],
		**argkey, **argval,
		*pstr;

setlinebuf(stdout);
 if (argc<2)
    {
    fprintf(OUTPUT, "\n         %s  version %s (%s)\n", BANNER,MYVERSION,DATE);
    fprintf(OUTPUT, "\nWritten by %s\n", AUTHORS);
    fprintf(OUTPUT, "Copyright %s\n", COPYRIGHT);
    fprintf(OUTPUT, "\nvisit %s\n", WEBSITE);
    fprintf(OUTPUT, "\n%s\n", DISCLAIMER);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }
  QMALLOC(argkey, char *, argc);
  QMALLOC(argval, char *, argc);

/*default parameters */
  prefs.command_line = argv;
  prefs.ncommand_line = argc;
  prefs.pipe_flag = 0;
  prefs.nimage_name = 1;
  prefs.image_name[0] = "image";
  strcpy(prefs.prefs_name, "default.sex");
  narg = nim = 0;

  for (a=1; a<argc; a++)
    {
    if (*(argv[a]) == '-')
      {
      opt = (int)argv[a][1];
      if (strlen(argv[a])<4 || opt == '-')
        {
        opt2 = (int)tolower((int)argv[a][2]);
        if (opt == '-')
          {
          opt = opt2;
          opt2 = (int)tolower((int)argv[a][3]);
          }
        switch(opt)
          {
          case 'c':
            if (a<(argc-1))
              strcpy(prefs.prefs_name, argv[++a]);
            break;
          case 'd':
            if (opt2=='d')
              dumpprefs(1);
            else if (opt2=='p')
              dumpparams();
            else
              dumpprefs(0);
            exit(EXIT_SUCCESS);
            break;
          case 'v':
            printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
            exit(0);
            break;
          case 'h':
          default:
            error(EXIT_SUCCESS,"SYNTAX: ", SYNTAX);
          }
        }
      else
        {
/*------ Config parameters */
        argkey[narg] = &argv[a][1];
        argval[narg++] = argv[++a];
        }       
      }
    else
      {
/*---- The input image filename(s) */
      for(; (a<argc) && (*argv[a]!='-'); a++)
        {
        strncpy(str, argv[a], MAXCHARL-1);
        for (pstr=NULL;(pstr=strtok(pstr?NULL:str, notokstr)); nim++)
          if (nim<MAXIMAGE)
            prefs.image_name[nim] = pstr;
          else
            error(EXIT_FAILURE, "*Error*: Too many input images: ", pstr);
        }
      prefs.nimage_name = nim;
      a--;
      }
    }

  readprefs(prefs.prefs_name, argkey, argval, narg);
  preprefs();

  free(argkey);
  free(argval);

  makeit();

  endprefs();
  NFPRINTF(OUTPUT, "");
  tdiff = prefs.time_diff>0.0? prefs.time_diff : 0.001;
  lines = (double)thefield1.height/tdiff;
  dets = (double)thecat.ntotal/tdiff;
  NPRINTF(OUTPUT,
	"> All done (in %.1f s: %.1f line%s/s , %.1f detection%s/s)\n",
	prefs.time_diff, lines, lines>1.0? "s":"", dets, dets>1.0? "s":"");

  return EXIT_SUCCESS;
  }

