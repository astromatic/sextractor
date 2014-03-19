/*
*				prefs.c
*
* Functions related to run-time configurations.
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
*	Last modified:		18/06/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<ctype.h>
#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include        <unistd.h>

#if defined(USE_THREADS) \
&& (defined(__APPLE__) || defined(FREEBSD) || defined(NETBSD))	/* BSD, Apple */
 #include	<sys/types.h>
 #include	<sys/sysctl.h>
#elif defined(USE_THREADS) && defined(HAVE_MPCTL)		/* HP/UX */
 #include	<sys/mpctl.h>
#endif

#ifdef HAVE_MKL
 #include MKL_H
#endif

#include	"define.h"
#include	"globals.h"
#include	"back.h"
#include	"prefs.h"
#include	"preflist.h"
#include	"fits/fitscat.h"


/********************************* dumpprefs ********************************/
/*
Print the default preference parameters.
*/
void    dumpprefs(int state)
  {
   char **dp;

  dp = default_prefs;
  while (**dp)
    if (**dp != '*')
      printf("%s\n",*(dp++));
    else if (state)
      printf("%s\n",*(dp++)+1);
    else
      dp++;
  return;
  }


/********************************* readprefs ********************************/
/*
Read a configuration file in ``standard'' format (see the SExtractor
documentation)
*/
void    readprefs(char *filename, char **argkey, char **argval, int narg)

  {
   FILE          *infile;
   char          *cp, str[MAXCHARL], *keyword, *value, **dp;
   int           i, ival, nkey, warn, argi, flagc, flagd, flage, flagz;
   float         dval;
#ifndef	NO_ENVVAR
   static char	value2[MAXCHARL],envname[MAXCHAR];
   char		*dolpos;
#endif


  if ((infile = fopen(filename,"r")) == NULL)
    {
    flage = 1;
    warning(filename, " not found, using internal defaults");
    }
  else
    flage = 0;

/*Build the keyword-list from pkeystruct-array */

  for (i=0; key[i].name[0]; i++)
    strcpy(keylist[i], key[i].name);
  keylist[i][0] = '\0';


/*Scan the configuration file*/

  argi=0;
  flagc = 0;
  flagd = 1;
  dp = default_prefs;
  for (warn=0;;)
    {
    if (flagd)
      {
      if (**dp)
        {
        if (**dp=='*')
          strcpy(str, *(dp++)+1);
        else
          strcpy(str, *(dp++));
        }
      else
        flagd = 0;
      }
    if (!flagc && !flagd)
      if (flage || !fgets(str, MAXCHARL, infile))
        flagc=1;

    if (flagc)
      {
      if (argi<narg)
        {
        sprintf(str, "%s %s", argkey[argi], argval[argi]);
        argi++;
        }
      else
        break;
      }

    keyword = strtok(str, notokstr);
    if (keyword && keyword[0]!=0 && keyword[0]!=(char)'#')
      {
     if (warn>=10)
        error(EXIT_FAILURE, "*Error*: No valid keyword found in ", filename);
      nkey = findkeys(keyword, keylist, FIND_STRICT);
      if (nkey!=RETURN_ERROR)
        {
        value = strtok((char *)NULL, notokstr);
#ifndef	NO_ENVVAR
/*------ Expansion of environment variables (preceded by '$') */
        if (value && (dolpos=strchr(value, '$')))
          {
           int	nc;
           char	*valuet,*value2t, *envval;

          value2t = value2;
          valuet = value;
          while (dolpos)
            {
            while (valuet<dolpos)
              *(value2t++) = *(valuet++);	/* verbatim copy before '$' */
            if (*(++valuet) == (char)'{')
              valuet++;
            strncpy(envname, valuet, nc=strcspn(valuet,"}/:\"\'\\"));
            *(envname+nc) = (char)'\0';
            if (*(valuet+=nc) == (char)'}')
              valuet++;
            if (!(envval=getenv(envname)))
              error(EXIT_FAILURE, "Environment variable not found: ",
				envname);
            while(*envval)			/* Copy the ENV content */
              *(value2t++) = *(envval++);
            while(*valuet && *valuet!=(char)'$')/* Continue verbatim copy */
              *(value2t++) = *(valuet++);
            if (*valuet)
              dolpos = valuet;
            else
              {
              dolpos = NULL;
              *value2t = (char)'\0';
              }
	    }

          value = value2;
          }
#endif
        switch(key[nkey].type)
          {
          case P_FLOAT:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            dval = atof(value);
            if (dval>=key[nkey].dmin && dval<=key[nkey].dmax)
              *(double *)(key[nkey].ptr) = dval;
            else
              error(EXIT_FAILURE, keyword," keyword out of range");
            break;

          case P_INT:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            ival = atoi(value);
            if (ival>=key[nkey].imin && ival<=key[nkey].imax)
              *(int *)(key[nkey].ptr) = ival;
            else
              error(EXIT_FAILURE, keyword, " keyword out of range");
            break;

          case P_STRING:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," string is empty!");
            strcpy((char *)key[nkey].ptr, value);
            break;

          case P_BOOL:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            if ((cp = strchr("yYnN", (int)value[0])))
              *(int *)(key[nkey].ptr) = (tolower((int)*cp)=='y')?1:0;
            else
              error(EXIT_FAILURE, keyword, " value must be Y or N");
            break;

          case P_KEY:
            if (!value || value[0]==(char)'#')
              error(EXIT_FAILURE, keyword," keyword has no value!");
            if ((ival = findkeys(value, key[nkey].keylist,FIND_STRICT))
			!= RETURN_ERROR)
              *(int *)(key[nkey].ptr) = ival;
            else
              error(EXIT_FAILURE, keyword, " set to an unknown keyword");
            break;

          case P_BOOLLIST:
            for (i=0; i<MAXLIST&&value&&value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              if ((cp = strchr("yYnN", (int)value[0])))
                ((int *)(key[nkey].ptr))[i] = (tolower((int)*cp)=='y')?1:0;
              else
                error(EXIT_FAILURE, keyword, " value must be Y or N");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_INTLIST:
            for (i=0; i<MAXLIST&&value&&value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              ival = strtol(value, NULL, 0);
              if (ival>=key[nkey].imin && ival<=key[nkey].imax)
                ((int *)key[nkey].ptr)[i] = ival;
              else
                error(EXIT_FAILURE, keyword, " keyword out of range");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_FLOATLIST:
            for (i=0; i<MAXLIST&&value&&value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              dval = atof(value);
              if (dval>=key[nkey].dmin && dval<=key[nkey].dmax)
                ((double *)key[nkey].ptr)[i] = dval;
              else
                error(EXIT_FAILURE, keyword, " keyword out of range");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_KEYLIST:
            for (i=0; i<MAXLIST && value && value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
	      if ((ival = findkeys(value, key[nkey].keylist, FIND_STRICT))
			!= RETURN_ERROR)
                ((int *)(key[nkey].ptr))[i] = ival;
              else
                error(EXIT_FAILURE, keyword, " set to an unknown keyword");
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = i;
            break;

          case P_STRINGLIST:
            if (!value || value[0]==(char)'#')
              {
              value = "";
              flagz = 1;
              }
            else
              flagz = 0;
            for (i=0; i<MAXLIST && value && value[0]!=(char)'#'; i++)
              {
              if (i>=key[nkey].nlistmax)
                error(EXIT_FAILURE, keyword, " has too many members");
              free(((char **)key[nkey].ptr)[i]);
              QMALLOC(((char **)key[nkey].ptr)[i], char, MAXCHAR);
              strcpy(((char **)key[nkey].ptr)[i], value);
              value = strtok((char *)NULL, notokstr);
              }
            if (i<key[nkey].nlistmin)
              error(EXIT_FAILURE, keyword, " list has not enough members");
            *(key[nkey].nlistptr) = flagz?0:i;
            break;

          default:
            error(EXIT_FAILURE, "*Internal ERROR*: Type Unknown",
				" in readprefs()");
            break;
          }
        key[nkey].flag = 1;
        }
      else
        {
        warning(keyword, " keyword unknown");
        warn++;
        }
      }
    }

  for (i=0; key[i].name[0]; i++)
    if (!key[i].flag)
      error(EXIT_FAILURE, key[i].name, " configuration keyword missing");
  if (!flage)
    fclose(infile);

  return;
  }


/********************************** findkeys **********************************/
/*
 find an item within a list of keywords, SExtractor version.
*/
int	findkeys(char *str, char keyw[][32], int mode)

  {
  int i;

  for (i=0; keyw[i][0]; i++)
    if (!cistrcmp(str, keyw[i], mode))
      return i;

  return RETURN_ERROR;
  }


/******************************* cistrcmp ***********************************/
/*
case-insensitive strcmp.
*/
int     cistrcmp(char *cs, char *ct, int mode)

  {
   int  i, diff;

  if (mode)
    {
    for (i=0; cs[i]&&ct[i]; i++)
      if ((diff=tolower((int)cs[i])-tolower((int)ct[i])))
        return diff;
    }
  else
    {
    for (i=0; cs[i]||ct[i]; i++)
      if ((diff=tolower((int)cs[i])-tolower((int)ct[i])))
        return diff;
    }

  return 0;
  }


/********************************* preprefs **********************************/
/*
Set number of threads and endianity.
*/
void	preprefs()

  {
   char			str[80];
   unsigned short	ashort=1;
#ifdef USE_THREADS
   int			nproc;
#endif

/* Test if byteswapping will be needed */
  bswapflag = *((char *)&ashort);

/* Multithreading */
#ifdef USE_THREADS
  if (prefs.nthreads <= 0)
    {
/*-- Get the number of processors for parallel builds */
/*-- See, e.g. http://ndevilla.free.fr/threads */
    nproc = -1;
#if defined(_SC_NPROCESSORS_ONLN)		/* AIX, Solaris, Linux */
    nproc = (int)sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROCESSORS_CONF)
    nproc = (int)sysconf(_SC_NPROCESSORS_CONF);
#elif defined(__APPLE__) || defined(FREEBSD) || defined(NETBSD)	/* BSD, Apple */
    {
     int        mib[2];
     size_t     len;

     mib[0] = CTL_HW;
     mib[1] = HW_NCPU;
     len = sizeof(nproc);
     sysctl(mib, 2, &nproc, &len, NULL, 0);
     }
#elif defined (_SC_NPROC_ONLN)			/* SGI IRIX */
    nproc = sysconf(_SC_NPROC_ONLN);
#elif defined(HAVE_MPCTL)			/* HP/UX */
    nproc =  mpctl(MPC_GETNUMSPUS_SYS, 0, 0);
#endif

    if (nproc>0)
      prefs.nthreads = ((prefs.nthreads) && nproc>(-prefs.nthreads))?
		-prefs.nthreads : nproc;
    else
      {
      prefs.nthreads = prefs.nthreads? -prefs.nthreads : 2;
      sprintf(str, "NTHREADS defaulted to %d", prefs.nthreads);
      warning("Cannot find the number of CPUs on this system:", str);
      }
    }
#ifndef HAVE_ATLAS_MP
   if (prefs.nthreads>1)
     warning("This executable has been compiled using a version of the ATLAS "
	"library without support for multithreading. ",
	"Performance will be degraded.");
#endif

#else
  if (prefs.nthreads != 1)
    {
    prefs.nthreads = 1;
    warning("NTHREADS != 1 ignored: ",
	"this build of " BANNER " is single-threaded");
    }
#endif

/* We temporarily allow threads for non-multithreaded code (defaulted to 1)*/
#ifdef HAVE_MKL
  mkl_set_num_threads(prefs.nthreads);
#endif

/* Override INTEL CPU detection routine to help performance on 3rd-party CPUs */
#if defined(__INTEL_COMPILER) && defined (USE_CPUREDISPATCH)
  __get_cpuid(1, &eax, &ebx, &ecx, &edx);
  if (ecx&bit_AVX)
    __intel_cpu_indicator = 0x20000;
  else if (ecx&bit_SSE4_2)
    __intel_cpu_indicator = 0x8000;
  else if (ecx&bit_SSE4_1)
    __intel_cpu_indicator = 0x2000;
  else if (ecx&bit_SSSE3)
    __intel_cpu_indicator = 0x1000;
  else if (ecx&bit_SSE3)
    __intel_cpu_indicator = 0x0800;
  else if (edx&bit_SSE2)
    __intel_cpu_indicator = 0x0200;
  else if (edx&bit_SSE)
    __intel_cpu_indicator = 0x0080;
  else if (edx&bit_MMX)
    __intel_cpu_indicator = 0x0008;
  else
    __intel_cpu_indicator = 0x0001;
#endif

  }


/********************************* useprefs **********************************/
/*
Update various structures according to the prefs.
*/
void	useprefs()

  {
   int			i, margin, naper;
   char			*str;

/*-------------------------------- Images ----------------------------------*/
  prefs.dimage_flag = (prefs.nimage_name>1);

/*--------------------------------- ASSOC ----------------------------------*/
  prefs.assoc_flag = FLAG(obj2.assoc) || FLAG(obj2.assoc_number);

/*-------------------------------- Extracting ------------------------------*/
  if (prefs.nthresh_type<2)
    prefs.thresh_type[1] = prefs.thresh_type[0];

/*-------------------------------- Deblending ------------------------------*/
  prefs.deb_maxarea = (prefs.ext_minarea<MAXDEBAREA ?
		prefs.ext_minarea:MAXDEBAREA);

/*-------------------------------- Astrometry ------------------------------*/
  prefs.world_flag = FLAG(obj2.mxw) || FLAG(obj2.mamaposx)
		|| FLAG(obj2.peakxw) || FLAG(obj2.winpos_xw)
		|| FLAG(obj2.mx2w) || FLAG(obj2.win_mx2w)
		|| FLAG(obj2.xw_prof) || FLAG(obj2.poserrmx2w_prof)
		|| FLAG(obj2.poserr_mx2w) || FLAG(obj2.winposerr_mx2w)
		|| FLAG(obj2.area_flagw) || FLAG(obj2.prof_flagw)
		|| FLAG(obj2.fwhmw_psf)
		|| (FLAG(obj2.sprob) && prefs.pixel_scale == 0.0);

/* Default astrometric settings */
  strcpy(prefs.coosys, "ICRS");
  prefs.epoch = 2000.0;

/*-------------------------------- Photometry ------------------------------*/

/* Find the largest APERture-photometry vector */
  if (FLAG(obj2.flux_aper))
    {
    naper = prefs.flux_apersize;
    if (prefs.fluxerr_apersize>naper)
      naper = prefs.fluxerr_apersize;
    if (prefs.mag_apersize>naper)
      naper = prefs.mag_apersize;
    if (prefs.magerr_apersize>naper)
      naper = prefs.magerr_apersize;
    if (naper>prefs.naper)
      {
      warning("Not enough apertures provided in config.:\n",
	"         some APER photometric values will remain blank ");
      naper = prefs.naper;
      }
    else
      prefs.naper = naper;
    }
  else
    naper = 0;			/* To avoid gcc -Wall warnings */

/* Find the largest "minimum margin" necessary for apertures */
  prefs.cleanmargin = 0;
  if (FLAG(obj2.vignet)
	&& (margin=(prefs.vignetsize[1]+1)/2) > prefs.cleanmargin)
    prefs.cleanmargin = margin;
  if (FLAG(obj2.vigshift)
	&& (margin=(prefs.vigshiftsize[1]+1)/2+3)>prefs.cleanmargin)
    prefs.cleanmargin = margin;
  if (FLAG(obj2.flux_aper))
    for (i=0; i<naper; i++)
      if ((margin=(int)((prefs.apert[i]+1)/2)+1) > prefs.cleanmargin)
        prefs.cleanmargin = margin;

/* Growth-curve flag */
  if (FLAG(obj2.flux_growth)
	|| FLAG(obj2.mag_growth)
	|| FLAG(obj2.flux_radius)
	|| FLAG(obj2.hl_radius)
	|| FLAG(obj2.flux_growthstep)
	|| FLAG(obj2.mag_growthstep))
    prefs.growth_flag = 1;

  if (FLAG(obj2.flux_radius) && prefs.flux_radiussize)
    if (prefs.nflux_frac>prefs.flux_radiussize)
      prefs.nflux_frac = prefs.flux_radiussize;

/*------------------------------- MASKing ----------------------------------*/
  prefs.blank_flag = (prefs.mask_type!=MASK_NONE);

/*--------------------------- SOM-fitting ----------------------------------*/
  prefs.somfit_flag = FLAG(obj2.flux_somfit);

/*------------------------------ Background --------------------------------*/
  if (prefs.nback_type<2)
    prefs.back_type[1] = prefs.back_type[0];
  if (prefs.nback_val<2)
    prefs.back_val[1] = prefs.back_val[0];

  if (prefs.nbacksize<2)
    prefs.backsize[1] = prefs.backsize[0];
  if (prefs.nbackfsize<2)
    prefs.backfsize[1] = prefs.backfsize[0];

/*------------------------------ FLAG-images -------------------------------*/
  prefs.nimaisoflag = (prefs.imaflag_size > prefs.imanflag_size) ?
			prefs.imaflag_size : prefs.imanflag_size;
  prefs.nimaflag = (prefs.nimaisoflag < prefs.nfimage_name) ?
		prefs.nimaisoflag : prefs.nfimage_name;

/*----------------------------- CHECK-images -------------------------------*/
  prefs.check_flag = 0;
  for (i=0; i<prefs.ncheck_type; i++)
    if (prefs.check_type[i] != CHECK_NONE)	/* at least 1 is not NONE */
      prefs.check_flag = 1;

  if (prefs.check_flag && prefs.ncheck_name!=prefs.ncheck_type)
    error(EXIT_FAILURE, "*Error*: CHECKIMAGE_NAME(s) and CHECKIMAGE_TYPE(s)",
		" are not in equal number");

/*---------------------------- PSF-fitting ---------------------------------*/
  if (FLAG(obj2.flux_psf) )
    {
    prefs.psffit_flag = 1;
/* We deactivate double-PSF fits for now */
/*
    prefs.dpsffit_flag = (prefs.npsf_name>1);
*/
    }
  if (prefs.check_flag)
    for (i=0; i<prefs.ncheck_type; i++)
      if (prefs.check_type[i] == CHECK_SUBPSFPROTOS
		|| prefs.check_type[i] == CHECK_PSFPROTOS)
        prefs.psffit_flag = 1;
  if (prefs.psffit_flag)
    prefs.psf_flag = 1;
  if (prefs.dpsffit_flag)
    prefs.dpsf_flag = 1;

/*---------------------------- PC-fitting ----------------------------------*/
/* PC-fitting is possible only if a PSF file is loaded */
  if (prefs.psf_flag)
    {
    prefs.pc_flag = FLAG(obj2.mx2_pc);
    if (prefs.check_flag)
      for (i=0; i<prefs.ncheck_type; i++)
        if (prefs.check_type[i] == CHECK_SUBPCPROTOS
		|| prefs.check_type[i] == CHECK_PCPROTOS
		|| prefs.check_type[i] == CHECK_PCOPROTOS)
          prefs.pc_flag = 1;
    }

/*----------------------------- Model-fitting -------------------------------*/
  if (prefs.check_flag)
    for (i=0; i<prefs.ncheck_type; i++)
      if (prefs.check_type[i] == CHECK_PROFILES
	|| prefs.check_type[i] == CHECK_SUBPROFILES
	|| prefs.check_type[i] == CHECK_SPHEROIDS
	|| prefs.check_type[i] == CHECK_SUBSPHEROIDS
	|| prefs.check_type[i] == CHECK_DISKS
	|| prefs.check_type[i] == CHECK_SUBDISKS)
        prefs.prof_flag = 1;
  if (prefs.prof_flag)
    prefs.psf_flag = 1;
  if (prefs.dprof_flag)
    prefs.dpsf_flag = 1;

/*-------------------------- Tracking the PSF FWHM --------------------------*/
  if (prefs.seeing_fwhm == 0 && (FLAG(obj2.sprob) || FLAG(obj2.fwhm_psf)))
    prefs.psf_flag = 1;

/*-------------------------- Pattern-fitting -------------------------------*/
/* Profile-fitting is possible only if a PSF file is loaded */
  if (prefs.check_flag)
    for (i=0; i<prefs.ncheck_type; i++)
      if (prefs.check_type[i] == CHECK_PATTERNS)
        prefs.pattern_flag = 1;

/*----------------------------- WEIGHT-images ------------------------------*/
  if (prefs.nweight_type<2)
    prefs.weight_type[1] = prefs.weight_type[0];

  prefs.dweight_flag = (prefs.weight_type[0]!= WEIGHT_NONE);
  prefs.weight_flag = (prefs.weight_type[1]!= WEIGHT_NONE);

  if (prefs.dweight_flag || prefs.weight_flag)
    {
/*-- Handle the default weight-threshold values */
    if (prefs.nweight_thresh<2)
      for (i=2; --i >= prefs.nweight_thresh;)
        prefs.weight_thresh[i] = (prefs.weight_type[i]==WEIGHT_FROMWEIGHTMAP)?
					0.0 : BIG;
/*-- Weight rescaling flags */
    if (prefs.nwscale_flag<2)
       prefs.wscale_flag[1] = (prefs.weight_type[1]==WEIGHT_FROMBACK)?
					BACK_WSCALE : prefs.wscale_flag[0];

/*-- Check WEIGHT_IMAGE parameter(s) */
    if ((!prefs.nwimage_name
	&& ((prefs.weight_type[0]!=WEIGHT_FROMBACK
		&& prefs.weight_type[0]!=WEIGHT_NONE)
	|| (prefs.weight_type[1]!=WEIGHT_FROMBACK
		&& prefs.weight_type[1]!=WEIGHT_NONE)))
	|| (prefs.nwimage_name<2
	&& prefs.weight_type[0]!=WEIGHT_FROMBACK
	&& prefs.weight_type[0]!=WEIGHT_NONE
	&& prefs.weight_type[1]!=WEIGHT_FROMBACK
	&& prefs.weight_type[1]!=WEIGHT_NONE
	&& prefs.weight_type[0]!=prefs.weight_type[1]))
      error(EXIT_FAILURE, "*Error*: WEIGHT_IMAGE missing","");

    if (prefs.nwimage_name && prefs.nwimage_name<2)
      prefs.wimage_name[1] = prefs.wimage_name[0];
    if (prefs.nwimage_name==2 && prefs.nweight_type==1)
      prefs.nweight_type = 2;

/*-- If detection-only interpolation is needed with 1 Weight image... */
/*-- ...pretend we're using 2, with only one being interpolated */
    if (prefs.nweight_type==1
	&& prefs.nwimage_name && prefs.wimage_name[1]==prefs.wimage_name[0]
	&& prefs.interp_type[0]==INTERP_VARONLY )
      {
      prefs.nweight_type = 2;
      prefs.weight_type[1] = prefs.weight_type[0];
      prefs.weight_type[0] = WEIGHT_FROMINTERP;
      prefs.wimage_name[1] = prefs.wimage_name[0];
      prefs.interp_type[1] = INTERP_NONE;
      prefs.dweight_flag = 1;
      if (prefs.nweight_thresh<2)
        {
        prefs.nweight_thresh = 2;
        prefs.weight_thresh[1] = prefs.weight_thresh[0];
        }
      }
    }

/*------------------------------ Catalogue ---------------------------------*/

  if (!strcmp(prefs.cat_name, "STDOUT"))
    prefs.pipe_flag = 1;

  if ((str=strrchr(prefs.filter_name, '/')))
    strcpy(thecat.filter_name, str+1);
  else
    strcpy(thecat.filter_name, prefs.filter_name);

  if ((str=strrchr(prefs.prefs_name, '/')))
    strcpy(thecat.prefs_name, str+1);
  else
    strcpy(thecat.prefs_name, prefs.prefs_name);

  if ((str=strrchr(prefs.nnw_name, '/')))
    strcpy(thecat.nnw_name, str+1);
  else
    strcpy(thecat.nnw_name, prefs.nnw_name);

  if ((str=strrchr(prefs.image_name[prefs.nimage_name-1], '/')))
    strcpy(thecat.image_name, str+1);
  else
    strcpy(thecat.image_name, prefs.image_name[prefs.nimage_name-1]);

  sprintf(thecat.soft_name, "%s %s", BANNER, VERSION);

  return;
  }


/********************************* endprefs *********************************/
/*
Mostly free memory allocate for static arrays.
*/
void	endprefs(void)

  {
    int i;

  for (i=0; i<prefs.nfimage_name; i++)
      free(prefs.fimage_name[i]);
  for (i=0; i<prefs.nwimage_name; i++)
      free(prefs.wimage_name[i]);
  for (i=0; i<prefs.npsf_name; i++)
      free(prefs.psf_name[i]);
  for (i=0; i<prefs.ncheck_name; i++)
      free(prefs.check_name[i]);

  return;
  }
