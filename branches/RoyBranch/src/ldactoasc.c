/*
*				ldactoasc.c
*
* Stand-along tool that converts FITS-LDAC binary files to ASCII.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2007-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		11/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ldactoasc.h"
#include "fits/fitscat.h"

#define		SYNTAX  "ldactoasc catalog\n"
extern const char	notokstr[];

/********************************** main ************************************/
int	main(int argc, char *argv[])

  {
   catstruct		*cat;
   tabstruct		*tab;
   unsigned short	ashort=1;
   char			catname[MAXCHAR];
   int			a, t, opt,opt2, flag;

  if (argc<2)
    {
    fprintf(OUTPUT, "\n                %s  Version %s (%s)\n",
		BANNER, MYVERSION, DATE);
    fprintf(OUTPUT, "\nFor information, please contact: %s\n", COPYRIGHT);
    error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
    }

/* Test if byteswapping will be needed */
  bswapflag = *((char *)&ashort);

/* Default parameters */
  for (a=1; a<argc; a++)
    {
    if (*(argv[a]) == '-')
      {
      opt = (int)argv[a][1];
      {
      opt2 = (int)tolower((int)argv[a][2]);
      if (opt == '-')
        {
        opt = opt2;
        opt2 = (int)tolower((int)argv[a][3]);
        }
      switch(opt)
        {
        case 'v':
          printf("%s version %s (%s)\n", BANNER,MYVERSION,DATE);
          exit(EXIT_SUCCESS);
          break;
        case 'h':
          fprintf(OUTPUT, "\nSYNTAX: %s", SYNTAX);
          exit(EXIT_SUCCESS);
          break;
        default:
          error(EXIT_SUCCESS, "SYNTAX: ", SYNTAX);
          }
        }
      }
    else
      strcpy(catname, argv[a]);
    }

  flag = 1;	/* display banner of first extension */
  if ((cat = read_cat(catname)))
    {
    tab = cat->tab;
    for (t=cat->ntab; t--; tab=tab->nexttab)
      if (!strcmp("LDAC_OBJECTS", tab->extname)
	  || !strcmp("OBJECTS", tab->extname))
        {
        show_keys(tab, NULL, NULL, 0, NULL, stdout, 1, flag, 0, SHOW_ASCII);
        flag = 0;
        }
    free_cat(&cat, 1);
    }
  else
    error(EXIT_FAILURE,"Cannot open ",catname);

  return EXIT_SUCCESS;
  }

