 /*
 				fitsmisc.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	The LDAC Tools
*
*	Author:		E.BERTIN, DeNIS/LDAC
*
*	Contents:	miscellaneous functions.
*
*	Last modify:	18/07/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<ctype.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<time.h>

#include	"fitscat_defs.h"
#include	"fitscat.h"

static void	(*errorfunc)(char *msg1, char *msg2) = NULL;
static char	warning_historystr[WARNING_NMAX][192]={""};
static int	nwarning = 0, nwarning_history = 0, nerror = 0;

/********************************* error ************************************/
/*
I hope it will never be used!
*/
void	error(int num, char *msg1, char *msg2)
  {
  fprintf(stderr, "\n> %s%s\n\n",msg1,msg2);
  if (num && errorfunc && !nerror)
    {
    nerror = 1;
    errorfunc(msg1, msg2);
    }
  exit(num);
  }


/**************************** error_installfunc *****************************/
/*
I hope it will never be used!
*/
void	error_installfunc(void (*func)(char *msg1, char *msg2))
  {
  if (func)
    errorfunc = func;

  return;
  }


/********************************* warning **********************************/
/*
Print a warning message on screen.
*/
void    warning(char *msg1, char *msg2)
  {
   time_t	warntime;
   struct tm	*tm;

  warntime = time(NULL);
  tm = localtime(&warntime);
 
  fprintf(stderr, "\n> WARNING: %s%s\n\n",msg1,msg2);
  sprintf(warning_historystr[(nwarning++)%WARNING_NMAX],
	"%04d-%02d-%02d %02d:%02d:%02d : %.80s%.80s",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday,
	tm->tm_hour, tm->tm_min, tm->tm_sec,
	msg1, msg2);


  return;
  }


/****************************** warning_history ******************************/
/*
Return warning.
*/
char    *warning_history(void)
  {
   char		*str;

  if (nwarning_history >= WARNING_NMAX)
    {
    nwarning_history = 0;	/* So it can be accessed later on */
    return "";
    }

  str = warning_historystr[((nwarning>WARNING_NMAX? (nwarning%WARNING_NMAX):0)
	+ nwarning_history++)%WARNING_NMAX];
  if (!*str)
    nwarning_history = 0;	/* So it can be accessed later on */

  return str;
  }


/******************************* swapbytes **********************************/
/*
Swap bytes for doubles, longs and shorts (for DEC machines or PC for inst.).
*/
void    swapbytes(void *ptr, int nb, int n)
  {
   char	*cp,
	c;
   int	j;

  cp = (char *)ptr;

  if (nb&4)
    {
    for (j=n; j--; cp+=4)
      {
      c = cp[3];
      cp[3] = cp[0];
      cp[0] = c;
      c = cp[2];
      cp[2] = cp[1];
      cp[1] = c;
      }
    return;
    }

  if (nb&2)
    {
    for (j=n; j--; cp+=2)
      {
      c = cp[1];
      cp[1] = cp[0];
      cp[0] = c;
      }
    return;
    }

  if (nb&1)
    return;

  if (nb&8)
    {
    for (j=n; j--; cp+=8)
      {
      c = cp[7];
      cp[7] = cp[0];
      cp[0] = c;
      c = cp[6];
      cp[6] = cp[1];
      cp[1] = c;
      c = cp[5];
      cp[5] = cp[2];
      cp[2] = c;
      c = cp[4];
      cp[4] = cp[3];
      cp[3] = c;
      }
    return;
    }

  error(EXIT_FAILURE, "*Internal Error*: Unknown size in ", "swapbytes()");

  return;
  }


/****** wstrncmp ***************************************************************
PROTO	int wstrncmp(char *cs, char *ct, int n)
PURPOSE	simple wildcard strcmp.
INPUT	character string 1,
	character string 2,
	maximum number of characters to be compared.
OUTPUT	comparison integer (same meaning as strcmp).
NOTES	-.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	15/02/96
 ***/
int	wstrncmp(char *cs, char *ct, int n)

  {
   int	diff,i;

  i = n;
  diff = 0;
  do
    {
    diff = ((*cs=='?'&&*ct)||(*ct=='?'&&*cs))?0:*cs-*ct;
    } while (!diff && --i && *(cs++) && *(ct++));

  return diff;
  }


/****** findkey ****************************************************************
PROTO	int findkey(char *str, char *key, int size)
PURPOSE	Find an item within a list of keywords.
INPUT	character string,
	an array of character strings containing the list of keywords,
	offset (in char) between each keyword.
OUTPUT	position in the list (0 = first) if keyword matched,
	RETURN_ERROR otherwise.
NOTES	the matching is case-sensitive.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	15/02/96
 ***/
int	findkey(char *str, char *key, int size)

  {
  int i;

  for (i=0; key[0]; i++, key += size)
    if (!strcmp(str, key))
      return i;

  return RETURN_ERROR;
  }


/********************************* findnkey **********************************
PROTO	int findnkey(char *str, char *key, int size, int nkey)
PURPOSE	Find an item within a list of nkey keywords.
INPUT	character string,
	an array of character strings containing the list of keywords,
	offset (in char) between each keyword.
	number of keywords.
OUTPUT	position in the list (0 = first) if keyword matched,
	RETURN_ERROR otherwise.
NOTES	the matching is case-sensitive.
AUTHOR	E. Bertin (IAP & Leiden observatory)
VERSION	15/02/96
 ***/
int	findnkey(char *str, char *key, int size, int nkey)

  {
  int i;

  for (i=0; i<nkey; i++, key += size)
    if (!strcmp(str, key))
      return i;

  return RETURN_ERROR;
  }


