/*
 				ldactoasc.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	LDACtoASC
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Convert LDAC binary format to ASCII.
*
*	Last modify:	04/06/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
   int			a, t, opt,opt2;

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

  cat = read_cat(catname);
  tab = cat->tab;
  for (t=cat->ntab; t--; tab=tab->nexttab)
    if (!strcmp("LDAC_OBJECTS", tab->extname)
	|| !strcmp("OBJECTS", tab->extname))
      show_keys(tab, NULL, NULL, 0, NULL, stdout, 1, 1, 0, SHOW_ASCII);
  free_cat(&cat, 1);

  exit(EXIT_SUCCESS);
  }

