/*
*				fitswcs.c
*
* Manage World Coordinate System data.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1993-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		13/07/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"fits/fitscat_defs.h"
#include	"fits/fitscat.h"
#include	"fitswcs.h"
#include	"wcscelsys.h"
#include	"wcs/wcs.h"
#include	"wcs/lin.h"
#include	"wcs/tnx.h"
#include	"wcs/poly.h"

/******* copy_wcs ************************************************************
PROTO	wcsstruct *copy_wcs(wcsstruct *wcsin)
PURPOSE	Copy a WCS (World Coordinate System) structure.
INPUT	WCS structure to be copied.
OUTPUT	pointer to a copy of the input structure.
NOTES	Actually, only FITS parameters are copied. Lower-level structures
	such as those created by the WCS or TNX libraries are generated.
AUTHOR	E. Bertin (IAP)
VERSION	31/08/2002
 ***/
wcsstruct	*copy_wcs(wcsstruct *wcsin)

  {
   wcsstruct	*wcs;

/* Copy the basic stuff */
  QMEMCPY(wcsin, wcs, wcsstruct, 1);
/* The PROJP WCS parameters */
  QMEMCPY(wcsin->projp, wcs->projp, double, wcs->naxis*100);

/* Set other structure pointers to NULL (they'll have to be reallocated) */
  wcs->wcsprm = NULL;
  wcs->lin = NULL;
  wcs->cel = NULL;
  wcs->prj = NULL;
  wcs->tnx_lngcor = copy_tnxaxis(wcsin->tnx_lngcor);
  wcs->tnx_latcor = copy_tnxaxis(wcsin->tnx_latcor);
  wcs->inv_x = wcs->inv_y = NULL;

  QCALLOC(wcs->wcsprm, struct wcsprm, 1);
/* Test if the WCS is recognized and a celestial pair is found */
  wcsset(wcs->naxis,(const char(*)[9])wcs->ctype, wcs->wcsprm);

/* Initialize other WCS structures */
  init_wcs(wcs);
/* Invert projection corrections */
  invert_wcs(wcs);
/* Find the range of coordinates */
  range_wcs(wcs);  

  return wcs;
  }


/******* create_wcs ***********************************************************
PROTO	wcsstruct *create_wcs(char **ctype, double *crval, double *crpix,
			double *cdelt, int *naxisn, int naxis)
PURPOSE	Generate a simple WCS (World Coordinate System) structure.
INPUT	Pointer to an array of char strings with WCS projection on each axis,
	pointer to an array of center coordinates (double),
	pointer to an array of device coordinates (double),
	pointer to an array of pixel scales (double),
	pointer to an array of image dimensions (int),
	number of dimensions.
OUTPUT	pointer to a WCS structure.
NOTES	If a pointer is set to null, the corresponding variables are set to
	default values.
AUTHOR	E. Bertin (IAP)
VERSION	09/08/2006
 ***/
wcsstruct	*create_wcs(char **ctype, double *crval, double *crpix,
			double *cdelt, int *naxisn, int naxis)

  {
   wcsstruct	*wcs;
   int		l;

  QCALLOC(wcs, wcsstruct, 1);
  wcs->naxis = naxis;
  QCALLOC(wcs->projp, double, naxis*100);
  wcs->nprojp = 0;

  wcs->longpole = wcs->latpole = 999.0;
  for (l=0; l<naxis; l++)
    {
    wcs->naxisn[l] = naxisn? naxisn[l] : 360.0;
/*-- The default WCS projection system is an all-sky Aitoff projection */
    if (ctype)
      strncpy(wcs->ctype[l], ctype[l], 8);
    else if (l==0)
      strncpy(wcs->ctype[l], "RA---AIT", 8);
    else if (l==1)
      strncpy(wcs->ctype[l], "DEC--AIT", 8);
    wcs->crval[l] = crval? crval[l]: 0.0;
    wcs->crpix[l] = crpix? crpix[l]: 0.0;
    wcs->cdelt[l] = 1.0;
    wcs->cd[l*(naxis+1)] = cdelt? cdelt[l] : 1.0;
    }

  wcs->epoch = wcs->equinox = 2000.0;
  QCALLOC(wcs->wcsprm, struct wcsprm, 1);

/* Test if the WCS is recognized and a celestial pair is found */
  wcsset(wcs->naxis,(const char(*)[9])wcs->ctype, wcs->wcsprm);

/* Initialize other WCS structures */
  init_wcs(wcs);
/* Invert projection corrections */
  invert_wcs(wcs);
/* Find the range of coordinates */
  range_wcs(wcs);  

  return wcs;
  }


/******* init_wcs ************************************************************
PROTO	void init_wcs(wcsstruct *wcs)
PURPOSE	Initialize astrometry and WCS (World Coordinate System) structures.
INPUT	WCS structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	17/05/2007
 ***/
void	init_wcs(wcsstruct *wcs)

  {
   int		l,n,lng,lat,naxis;

  naxis = wcs->naxis;
  if (wcs->lin)
    {
    free(wcs->lin->cdelt);
    free(wcs->lin->crpix);
    free(wcs->lin->pc);
    free(wcs->lin->piximg);
    free(wcs->lin->imgpix);
    free(wcs->lin);
    }
  QCALLOC(wcs->lin, struct linprm, 1);
  QCALLOC(wcs->lin->cdelt, double, naxis);
  QCALLOC(wcs->lin->crpix, double, naxis);
  QCALLOC(wcs->lin->pc, double, naxis*naxis);


  if (wcs->cel)
    free(wcs->cel);
  QCALLOC(wcs->cel, struct celprm, 1);

  if (wcs->prj)
    free(wcs->prj);
  QCALLOC(wcs->prj, struct prjprm, 1);

  if (wcs->inv_x)
    {
    poly_end(wcs->inv_x);
    wcs->inv_x = NULL;
    }
  if (wcs->inv_y)
    {
    poly_end(wcs->inv_y);
    wcs->inv_y = NULL;
    }

/* Set WCS flags to 0: structures will be reinitialized by the WCS library */
  wcs->lin->flag = wcs->cel->flag = wcs->prj->flag = 0;
  wcs->lin->naxis = naxis;

/* wcsprm structure */
  lng = wcs->lng = wcs->wcsprm->lng;
  lat = wcs->lat = wcs->wcsprm->lat;

/* linprm structure */
  for (l=0; l<naxis; l++)
    {
    wcs->lin->crpix[l] = wcs->crpix[l];
    wcs->lin->cdelt[l] = 1.0;
    }

  for (l=0; l<naxis*naxis; l++)
    wcs->lin->pc[l] = wcs->cd[l];

/* celprm structure */
  if (lng>=0)
    {
    wcs->cel->ref[0] = wcs->crval[lng];
    wcs->cel->ref[1] = wcs->crval[lat];
    }
  else
    {
    wcs->cel->ref[0] = wcs->crval[0];
    wcs->cel->ref[1] = wcs->crval[1];
    }
  wcs->cel->ref[2] = wcs->longpole;
  wcs->cel->ref[3] = wcs->latpole;

/* prjprm structure */
  wcs->prj->r0 = wcs->r0;
  wcs->prj->tnx_lngcor = wcs->tnx_lngcor;
  wcs->prj->tnx_latcor = wcs->tnx_latcor;
  if (lng>=0)
    {
    n = 0;
    for (l=100; l--;)
      {
      wcs->prj->p[l] = wcs->projp[l+lat*100];	/* lat comes first for ... */
      wcs->prj->p[l+100] = wcs->projp[l+lng*100];/* ... compatibility reasons */
      if (!n && (wcs->prj->p[l] || wcs->prj->p[l+100]))
        n = l+1;
      }
    wcs->nprojp = n;
    }

/* Check-out chirality */
  wcs->chirality = wcs_chirality(wcs);

/* Initialize Equatorial <=> Celestial coordinate system transforms */
  init_wcscelsys(wcs);

  return;
  }


/******* init_wcscelsys *******************************************************
PROTO	void init_wcscelsys(wcsstruct *wcs)
PURPOSE	Initialize Equatorial <=> Celestial coordinate system transforms.
INPUT	WCS structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/07/2006
 ***/
void	init_wcscelsys(wcsstruct *wcs)

  {
  double	*mat,
		a0,d0,ap,dp,ap2,y;
  int		s,lng,lat;

  lng = wcs->wcsprm->lng;
  lat = wcs->wcsprm->lat;
/* Is it a celestial system? If not, exit! */
  if (lng==lat)
    {
    wcs->celsysconvflag = 0;
    return;
    }
/* Find the celestial system */
  for (s=0; *celsysname[s][0] && strncmp(wcs->ctype[lng], celsysname[s][0], 4);
	s++);
/* Is it a known, non-equatorial system? If not, exit! */
  if (!s || !*celsysname[s][0])
    {
    wcs->celsysconvflag = 0;
    return;
    }
  wcs->celsys = (celsysenum)s;
/* Some shortcuts */
  a0 = celsysorig[s][0]*DEG;
  d0 = celsysorig[s][1]*DEG;
  ap = celsyspole[s][0]*DEG;
  dp = celsyspole[s][1]*DEG;
/* First compute in the output referential the longitude of the south pole */
  y = sin(ap - a0);
/*
  x = cos(d0)*(cos(d0)*sin(dp)*cos(ap-a0)-sin(d0)*cos(dp));
  ap2 = atan2(y,x);
*/
  ap2 = asin(cos(d0)*y) ;
/* Equatorial <=> Celestial System transformation parameters */
  mat = wcs->celsysmat;
  mat[0] = ap;
  mat[1] = ap2;
  mat[2] = cos(dp);
  mat[3] = sin(dp);

  wcs->celsysconvflag = 1;
  return;
  }


/******* read_wcs *************************************************************
PROTO	wcsstruct *read_wcs(tabstruct *tab)
PURPOSE	Read WCS (World Coordinate System) info in the FITS header.
INPUT	tab structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/06/2012
 ***/
wcsstruct	*read_wcs(tabstruct *tab)

  {
#define	FITSREADF(buf, k, val, def) \
		{if (fitsread(buf,k, &val, H_FLOAT,T_DOUBLE) != RETURN_OK) \
		   val = def; \
		}

#define	FITSREADI(buf, k, val, def) \
		{if (fitsread(buf,k, &val, H_INT,T_LONG) != RETURN_OK) \
		   val = def; \
		}

#define	FITSREADS(buf, k, str, def) \
		{if (fitsread(buf,k,str, H_STRING,T_STRING) != RETURN_OK) \
		   strcpy(str, (def)); \
		}
   char		str[MAXCHARS];
   char		wstr1[TNX_MAXCHARS], wstr2[TNX_MAXCHARS];

   wcsstruct	*wcs;
   double	drota;
   int		j, l, naxis;
   char		name[16],
		*buf, *filename, *ptr;

  buf = tab->headbuf;
  filename = (tab->cat? tab->cat->filename : strcpy(name, "internal header"));

  FITSREADS(buf, "OBJECT  ", str, "Unnamed");

  QCALLOC(wcs, wcsstruct, 1);
  if (tab->naxis > NAXIS)
    {
    warning("Maximum number of dimensions supported by this version of the ",
	"software exceeded\n");
    tab->naxis = 2;
    }

  wcs->naxis = naxis = tab->naxis;
  QCALLOC(wcs->projp, double, naxis*100);

  for (l=0; l<naxis; l++)
    {
    wcs->naxisn[l] = tab->naxisn[l];
    sprintf(str, "CTYPE%-3d", l+1);
    FITSREADS(buf, str, str, "");
    strncpy(wcs->ctype[l], str, 8);
    sprintf(str, "CUNIT%-3d", l+1);
    FITSREADS(buf, str, str, "deg");
    strncpy(wcs->cunit[l], str, 32);
    sprintf(str, "CRVAL%-3d", l+1);
    FITSREADF(buf, str, wcs->crval[l], 0.0);
    sprintf(str, "CRPIX%-3d", l+1);
    FITSREADF(buf, str, wcs->crpix[l], 1.0);
    sprintf(str, "CDELT%-3d", l+1);
    FITSREADF(buf, str, wcs->cdelt[l], 1.0);
    sprintf(str, "CRDER%-3d", l+1);
    FITSREADF(buf, str, wcs->crder[l], 0.0);
    sprintf(str, "CSYER%-3d", l+1);
    FITSREADF(buf, str, wcs->csyer[l], 0.0);
    if (fabs(wcs->cdelt[l]) < 1e-30)
      error(EXIT_FAILURE, "*Error*: CDELT parameters out of range in ",
	filename);
    }

  if (fitsfind(buf, "CD?_????")!=RETURN_ERROR)
    {
/*-- If CD keywords exist, use them for the linear mapping terms... */
    for (l=0; l<naxis; l++)
      for (j=0; j<naxis; j++)
        {
        sprintf(str, "CD%d_%d", l+1, j+1);
        FITSREADF(buf, str, wcs->cd[l*naxis+j], l==j?1.0:0.0)
        }
    }
  else if (fitsfind(buf, "PC00?00?")!=RETURN_ERROR)
/*-- ...If PC keywords exist, use them for the linear mapping terms... */
    for (l=0; l<naxis; l++)
      for (j=0; j<naxis; j++)
        {
        sprintf(str, "PC%03d%03d", l+1, j+1);
        FITSREADF(buf, str, wcs->cd[l*naxis+j], l==j?1.0:0.0)
        wcs->cd[l*naxis+j] *= wcs->cdelt[l];
        }
  else
    {
/*-- ...otherwise take the obsolete CROTA2 parameter */
    FITSREADF(buf, "CROTA2  ", drota, 0.0)
    wcs->cd[3] = wcs->cd[0] = cos(drota*DEG);
    wcs->cd[1] = -(wcs->cd[2] = sin(drota*DEG));
    wcs->cd[0] *= wcs->cdelt[0];
    wcs->cd[2] *= wcs->cdelt[0];
    wcs->cd[1] *= wcs->cdelt[1];
    wcs->cd[3] *= wcs->cdelt[1];
    }
  QCALLOC(wcs->wcsprm, struct wcsprm, 1);

/* Test if the WCS is recognized and a celestial pair is found */
  if (!wcsset(wcs->naxis,(const char(*)[9])wcs->ctype, wcs->wcsprm)
	&& wcs->wcsprm->flag<999)
    {
     char	*pstr;
     double	date;
     int	biss, dpar[3];

/*-- Coordinate reference frame */
/*-- Search for an observation date expressed in Julian days */
    FITSREADF(buf, "MJD-OBS ", date, -1.0);
    if (date<0.0)
      FITSREADF(buf, "MJDSTART", date, -1.0);
/*-- Precession date (defined from Ephemerides du Bureau des Longitudes) */
/*-- in Julian years from 2000.0 */
    if (date>0.0)
      wcs->obsdate = 2000.0 - (MJD2000 - date)/365.25;
    else
      {
/*---- Search for an observation date expressed in "civilian" format */
      FITSREADS(buf, "DATE-OBS", str, "");
      if (*str)
        {
/*------ Decode DATE-OBS format: DD/MM/YY or YYYY-MM-DD */
        for (l=0; l<3 && (pstr = strtok_r(l?NULL:str,"/- ", &ptr)); l++)
          dpar[l] = atoi(pstr);
        if (l<3 || !dpar[0] || !dpar[1] || !dpar[2])
          {
/*-------- If DATE-OBS value corrupted or incomplete, assume 2000-1-1 */
          warning("Invalid DATE-OBS value in header: ", str);
          dpar[0] = 2000; dpar[1] = 1; dpar[2] = 1;
          }
        else if (strchr(str, '/') && dpar[0]<32 && dpar[2]<100)
          {
          j = dpar[0];
          dpar[0] = dpar[2]+1900;
          dpar[2] = j;
          }

        biss = (dpar[0]%4)?0:1;
/*------ Convert date to MJD */
        date = -678956 + (365*dpar[0]+dpar[0]/4) - biss
			+ ((dpar[1]>2?((int)((dpar[1]+1)*30.6)-63+biss)
		:((dpar[1]-1)*(63+biss))/2) + dpar[2]);
        wcs->obsdate = 2000.0 - (MJD2000 - date)/365.25;
        }
      else
/*------ Well if really no date is found */
        wcs->obsdate = 0.0;
      }

    FITSREADF(buf, "EPOCH", wcs->epoch, 2000.0);
    FITSREADF(buf, "EQUINOX", wcs->equinox, wcs->epoch);
    if (fitsread(buf, "RADESYS", str, H_STRING,T_STRING) != RETURN_OK)
      FITSREADS(buf, "RADECSYS", str,
	wcs->equinox >= 2000.0? "ICRS" : (wcs->equinox<1984.0? "FK4" : "FK5"));
    if (!strcmp(str, "ICRS"))
      wcs->radecsys = RDSYS_ICRS;
    else if (!strcmp(str, "FK5"))
      wcs->radecsys = RDSYS_FK5;
    else if (!strcmp(str, "FK4"))
      {
      if (wcs->equinox == 2000.0)
        {
        FITSREADF(buf, "EPOCH  ", wcs->equinox, 1950.0);
        FITSREADF(buf, "EQUINOX", wcs->equinox, wcs->equinox);
        }
      wcs->radecsys = RDSYS_FK4;
      warning("FK4 precession formulae not yet implemented:\n",
		"            Astrometry may be slightly inaccurate");
      }
    else if (!strcmp(str, "FK4-NO-E"))
      {
      if (wcs->equinox == 2000.0)
        {
        FITSREADF(buf, "EPOCH", wcs->equinox, 1950.0);
        FITSREADF(buf, "EQUINOX", wcs->equinox, wcs->equinox);
        }
      wcs->radecsys = RDSYS_FK4_NO_E;
      warning("FK4 precession formulae not yet implemented:\n",
		"            Astrometry may be slightly inaccurate");
      }
    else if (!strcmp(str, "GAPPT"))
      {
      wcs->radecsys = RDSYS_GAPPT;
      warning("GAPPT reference frame not yet implemented:\n",
		"            Astrometry may be slightly inaccurate");
      }
    else
      {
      warning("Using ICRS instead of unknown astrometric reference frame: ",
		str);
      wcs->radecsys = RDSYS_ICRS;
      }

/*-- Projection parameters */
    if (!strcmp(wcs->wcsprm->pcode, "TNX"))
      {
/*---- IRAF's TNX projection: decode these #$!?@#!! WAT parameters */
      if (fitsfind(buf, "WAT?????") != RETURN_ERROR)
        {
/*------ First we need to concatenate strings */
        pstr = wstr1;
        sprintf(str, "WAT1_001");
        for (j=2; fitsread(buf,str,pstr,H_STRINGS,T_STRING)==RETURN_OK; j++)
	  {
          sprintf(str, "WAT1_%03d", j);
          pstr += strlen(pstr);
	  }
        pstr = wstr2;
        sprintf(str, "WAT2_001");
        for (j=2; fitsread(buf,str,pstr,H_STRINGS,T_STRING)==RETURN_OK; j++)
	  {
          sprintf(str, "WAT2_%03d", j);
          pstr += strlen(pstr);
	  }
/*------ LONGPOLE defaulted to 180 deg if not found */
        if ((pstr = strstr(wstr1, "longpole"))
		|| (pstr = strstr(wstr2, "longpole")))
          pstr = strpbrk(pstr, "1234567890-+.");
        wcs->longpole = pstr? atof(pstr) : 999.0;
        wcs->latpole = 999.0;
/*------ RO defaulted to 180/PI if not found */
        if ((pstr = strstr(wstr1, "ro"))
		|| (pstr = strstr(wstr2, "ro")))
          pstr = strpbrk(pstr, "1234567890-+.");
        wcs->r0 = pstr? atof(pstr) : 0.0;
/*------ Read the remaining TNX parameters */
        if ((pstr = strstr(wstr1, "lngcor"))
		|| (pstr = strstr(wstr2, "lngcor")))
          wcs->tnx_lngcor = read_tnxaxis(pstr);
        if (!wcs->tnx_lngcor)
          error(EXIT_FAILURE, "*Error*: incorrect TNX parameters in ",
			filename);
        if ((pstr = strstr(wstr1, "latcor"))
		|| (pstr = strstr(wstr2, "latcor")))
          wcs->tnx_latcor = read_tnxaxis(pstr);
        if (!wcs->tnx_latcor)
          error(EXIT_FAILURE, "*Error*: incorrect TNX parameters in ",
			filename);
        }
      }
    else
      {
      FITSREADF(buf, "LONGPOLE", wcs->longpole, 999.0);
      FITSREADF(buf, "LATPOLE ", wcs->latpole, 999.0);
/*---- Old convention */
      if (fitsfind(buf, "PROJP???") != RETURN_ERROR)
        for (j=0; j<10; j++)
          {
          sprintf(str, "PROJP%-3d", j);
          FITSREADF(buf, str, wcs->projp[j], 0.0);
          }
/*---- New convention */
      if (fitsfind(buf, "PV?_????") != RETURN_ERROR)
        for (l=0; l<naxis; l++)
          for (j=0; j<100; j++)
            {
            sprintf(str, "PV%d_%d ", l+1, j);
            FITSREADF(buf, str, wcs->projp[j+l*100], 0.0);
            }
      }
    }

/* Initialize other WCS structures */
  init_wcs(wcs);

/* Find the range of coordinates */
  range_wcs(wcs);
/* Invert projection corrections */
  invert_wcs(wcs);

#undef FITSREADF
#undef FITSREADI
#undef FITSREADS

  return wcs;
  }


/******* write_wcs ***********************************************************
PROTO	void write_wcs(tabstruct *tab, wcsstruct *wcs)
PURPOSE	Write WCS (World Coordinate System) info in the FITS header.
INPUT	tab structure,
	WCS structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	01/09/2010
 ***/
void	write_wcs(tabstruct *tab, wcsstruct *wcs)

  {
   double	mjd;
   char		str[MAXCHARS];
   int		j, l, naxis;

  naxis = wcs->naxis;
  addkeywordto_head(tab, "BITPIX  ", "Bits per pixel");
  fitswrite(tab->headbuf, "BITPIX  ", &tab->bitpix, H_INT, T_LONG);
  addkeywordto_head(tab, "NAXIS   ", "Number of axes");
  fitswrite(tab->headbuf, "NAXIS   ", &wcs->naxis, H_INT, T_LONG);
  for (l=0; l<naxis; l++)
    {
    sprintf(str, "NAXIS%-3d", l+1);
    addkeywordto_head(tab, str, "Number of pixels along this axis");
    fitswrite(tab->headbuf, str, &wcs->naxisn[l], H_INT, T_LONG);
    }
  addkeywordto_head(tab, "EQUINOX ", "Mean equinox");
  fitswrite(tab->headbuf, "EQUINOX ", &wcs->equinox, H_FLOAT, T_DOUBLE);
  if (wcs->obsdate!=0.0)
    {
    mjd = (wcs->obsdate-2000.0)*365.25 + MJD2000;
    addkeywordto_head(tab, "MJD-OBS ", "Modified Julian date at start");
    fitswrite(tab->headbuf, "MJD-OBS ", &mjd, H_EXPO,T_DOUBLE);
    }
  addkeywordto_head(tab, "RADECSYS", "Astrometric system");
  switch(wcs->radecsys)
    {
    case RDSYS_ICRS:
      fitswrite(tab->headbuf, "RADECSYS", "ICRS", H_STRING, T_STRING);
      break;
    case RDSYS_FK5:
      fitswrite(tab->headbuf, "RADECSYS", "FK5", H_STRING, T_STRING);
      break;
    case RDSYS_FK4:
      fitswrite(tab->headbuf, "RADECSYS", "FK4", H_STRING, T_STRING);
      break;
    case RDSYS_FK4_NO_E:
      fitswrite(tab->headbuf, "RADECSYS", "FK4-NO-E", H_STRING, T_STRING);
      break;
    case RDSYS_GAPPT:
      fitswrite(tab->headbuf, "RADECSYS", "GAPPT", H_STRING, T_STRING);
      break;
    default:
      error(EXIT_FAILURE, "*Error*: unknown RADECSYS type in write_wcs()", "");
    }
  for (l=0; l<naxis; l++)
    {
    sprintf(str, "CTYPE%-3d", l+1);
    addkeywordto_head(tab, str, "WCS projection type for this axis");
    fitswrite(tab->headbuf, str, wcs->ctype[l], H_STRING, T_STRING);
    sprintf(str, "CUNIT%-3d", l+1);
    addkeywordto_head(tab, str, "Axis unit");
    fitswrite(tab->headbuf, str, wcs->cunit[l], H_STRING, T_STRING);
    sprintf(str, "CRVAL%-3d", l+1);
    addkeywordto_head(tab, str, "World coordinate on this axis");
    fitswrite(tab->headbuf, str, &wcs->crval[l], H_EXPO, T_DOUBLE);
    sprintf(str, "CRPIX%-3d", l+1);
    addkeywordto_head(tab, str, "Reference pixel on this axis");
    fitswrite(tab->headbuf, str, &wcs->crpix[l], H_EXPO, T_DOUBLE);
    for (j=0; j<naxis; j++)
      {
      sprintf(str, "CD%d_%d", l+1, j+1);
      addkeywordto_head(tab, str, "Linear projection matrix");
      fitswrite(tab->headbuf, str, &wcs->cd[l*naxis+j], H_EXPO, T_DOUBLE);
      }
    for (j=0; j<100; j++)
      if (wcs->projp[j+100*l] != 0.0)
        {
        sprintf(str, "PV%d_%d", l+1, j);
        addkeywordto_head(tab, str, "Projection distortion parameter");
        fitswrite(tab->headbuf, str, &wcs->projp[j+100*l], H_EXPO, T_DOUBLE);
        }
    }

/* Update the tab data */
  readbasic_head(tab);

  return;
  }


/******* end_wcs **************************************************************
PROTO	void end_wcs(wcsstruct *wcs)
PURPOSE	Free WCS (World Coordinate System) infos.
INPUT	WCS structure.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	24/05/2000
 ***/
void	end_wcs(wcsstruct *wcs)

  {
  if (wcs)
    {
    if (wcs->lin)
      {
      free(wcs->lin->cdelt);
      free(wcs->lin->crpix);
      free(wcs->lin->pc);
      free(wcs->lin->piximg);
      free(wcs->lin->imgpix);
      free(wcs->lin);
      }
    free(wcs->cel);
    free(wcs->prj);
    free(wcs->wcsprm);
    free_tnxaxis(wcs->tnx_lngcor);
    free_tnxaxis(wcs->tnx_latcor);
    poly_end(wcs->inv_x);
    poly_end(wcs->inv_y);
    free(wcs->projp);
    free(wcs);
    }

  return;
  }


/******* wcs_supproj *********************************************************
PROTO	int wcs_supproj(char *name)
PURPOSE	Tell if a projection system is supported or not.
INPUT	Proposed projection code name.
OUTPUT	RETURN_OK if projection is supported, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	14/06/2012
 ***/
int	wcs_supproj(char *name)

  {
   char	projcode[28][5] =
	{"TAN", "TPV", "AZP", "SIN", "STG", "ARC", "ZPN", "ZEA", "AIR", "CYP",
	"CAR", "MER", "CEA", "COP", "COD", "COE", "COO", "BON", "PCO", "GLS",
	"PAR", "AIT", "MOL", "CSC", "QSC", "TSC", "TNX", "NONE"};

   int	i;

  for (i=0; i<28; i++)
    if (!strcmp(name, projcode[i]))
      return RETURN_OK;

  return RETURN_ERROR;
  }


/******* invert_wcs ***********************************************************
PROTO	void invert_wcs(wcsstruct *wcs)
PURPOSE	Invert WCS projection mapping (using a polynomial).
INPUT	WCS structure.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	13/07/2012
 ***/
void	invert_wcs(wcsstruct *wcs)

  {
   polystruct		*poly;
   double		pixin[NAXIS],raw[NAXIS],rawmin[NAXIS];
   double		*outpos,*outpost, *lngpos,*lngpost,
			*latpos,*latpost,
			lngstep,latstep, rawsize, epsilon;
   int			group[] = {1,1};
				/* Don't ask, this is needed by poly_init()! */
   int		i,j,lng,lat,deg, tnxflag, maxflag;

/* Check first that inversion is not straightforward */
  lng = wcs->wcsprm->lng;
  lat = wcs->wcsprm->lat;
  if (!strcmp(wcs->wcsprm->pcode, "TNX"))
    tnxflag = 1;
  else if ((!strcmp(wcs->wcsprm->pcode, "TAN")
	|| !strcmp(wcs->wcsprm->pcode, "TPV"))
		&& (wcs->projp[1+lng*100] || wcs->projp[1+lat*100]))
    tnxflag = 0;
  else
    return;

/* We define x as "longitude" and y as "latitude" projections */
/* We assume that PCxx cross-terms with additional dimensions are small */
/* Sample the whole image with a regular grid */
  lngstep = wcs->naxisn[lng]/(WCS_NGRIDPOINTS-1.0);
  latstep = wcs->naxisn[lat]/(WCS_NGRIDPOINTS-1.0);
  QMALLOC(outpos, double, 2*WCS_NGRIDPOINTS2);
  QMALLOC(lngpos, double, WCS_NGRIDPOINTS2);
  QMALLOC(latpos, double, WCS_NGRIDPOINTS2);
  for (i=0; i<wcs->naxis; i++)
    raw[i] = rawmin[i] = 0.5;
  outpost = outpos;
  lngpost = lngpos;
  latpost = latpos;
  for (j=WCS_NGRIDPOINTS; j--; raw[lat]+=latstep)
    {
    raw[lng] = rawmin[lng];
    for (i=WCS_NGRIDPOINTS; i--; raw[lng]+=lngstep)
      {
      if (linrev(raw, wcs->lin, pixin))
        error(EXIT_FAILURE, "*Error*: incorrect linear conversion in ",
		wcs->wcsprm->pcode);
      *(lngpost++) = pixin[lng];
      *(latpost++) = pixin[lat];
      if (tnxflag)
        {
        *(outpost++) = pixin[lng]
			+raw_to_tnxaxis(wcs->tnx_lngcor,pixin[lng],pixin[lat]);
        *(outpost++) = pixin[lat]
			+raw_to_tnxaxis(wcs->tnx_latcor,pixin[lng],pixin[lat]);
        }
      else
        {
        raw_to_pv(wcs->prj,pixin[lng],pixin[lat], outpost, outpost+1);
        outpost += 2;
        }
      }
    }

/* Invert "longitude" */
/* Compute the extent of the pixel in reduced projected coordinates */
  linrev(rawmin, wcs->lin, pixin);
  pixin[lng] += ARCSEC/DEG;
  linfwd(pixin, wcs->lin, raw);
  rawsize = sqrt((raw[lng]-rawmin[lng])*(raw[lng]-rawmin[lng])
		+(raw[lat]-rawmin[lat])*(raw[lat]-rawmin[lat]))*DEG/ARCSEC;
  if (!rawsize)
    error(EXIT_FAILURE, "*Error*: incorrect linear conversion in ",
		wcs->wcsprm->pcode);
  epsilon = WCS_INVACCURACY/rawsize;
/* Find the lowest degree polynom */
  poly = NULL;  /* to avoid gcc -Wall warnings */
  maxflag = 1;
  for (deg=1; deg<=WCS_INVMAXDEG && maxflag; deg++)
    {
    if (deg>1)
      poly_end(poly);
    poly = poly_init(group, 2, &deg, 1);
    poly_fit(poly, outpos, lngpos, NULL, WCS_NGRIDPOINTS2, NULL);
    maxflag = 0;
    outpost = outpos;
    lngpost = lngpos;
    for (i=WCS_NGRIDPOINTS2; i--; outpost+=2)
      if (fabs(poly_func(poly, outpost)-*(lngpost++))>epsilon)
        {
        maxflag = 1;
        break;
        }
    }
  if (maxflag)
    warning("Significant inaccuracy likely to occur in projection","");
/* Now link the created structure */
  wcs->prj->inv_x = wcs->inv_x = poly;

/* Invert "latitude" */
/* Compute the extent of the pixel in reduced projected coordinates */
  linrev(rawmin, wcs->lin, pixin);
  pixin[lat] += ARCSEC/DEG;
  linfwd(pixin, wcs->lin, raw);
  rawsize = sqrt((raw[lng]-rawmin[lng])*(raw[lng]-rawmin[lng])
		+(raw[lat]-rawmin[lat])*(raw[lat]-rawmin[lat]))*DEG/ARCSEC;
  if (!rawsize)
    error(EXIT_FAILURE, "*Error*: incorrect linear conversion in ",
		wcs->wcsprm->pcode);
  epsilon = WCS_INVACCURACY/rawsize;
/* Find the lowest degree polynom */
  maxflag = 1;
  for (deg=1; deg<=WCS_INVMAXDEG && maxflag; deg++)
    {
    if (deg>1)
      poly_end(poly);
    poly = poly_init(group, 2, &deg, 1);
    poly_fit(poly, outpos, latpos, NULL, WCS_NGRIDPOINTS2, NULL);
    maxflag = 0;
    outpost = outpos;
    latpost = latpos;
    for (i=WCS_NGRIDPOINTS2; i--; outpost+=2)
      if (fabs(poly_func(poly, outpost)-*(latpost++))>epsilon)
        {
        maxflag = 1;
        break;
        }
    }
  if (maxflag)
    warning("Significant inaccuracy likely to occur in projection","");
/* Now link the created structure */
  wcs->prj->inv_y = wcs->inv_y = poly;

/* Free memory */
  free(outpos);
  free(lngpos);
  free(latpos);

  return;
  }


/******* range_wcs ***********************************************************
PROTO	void range_wcs(wcsstruct *wcs)
PURPOSE	Find roughly the range of WCS coordinates on all axes,
	and typical pixel scales.
INPUT	WCS structure.
OUTPUT	-.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	24/08/2010
 ***/
void	range_wcs(wcsstruct *wcs)

  {
   double		step[NAXIS], raw[NAXIS], rawmin[NAXIS],
			world[NAXIS], world2[NAXIS];
   double		*worldmin, *worldmax, *scale, *worldc,
			rad, radmax, lc;
   int			linecount[NAXIS];
   int			i,j, naxis, npoints, lng,lat;

  naxis = wcs->naxis;

/* World range */
  npoints = 1;
  worldmin = wcs->wcsmin;
  worldmax = wcs->wcsmax;
/* First, find the center and use it as a reference point for lng */
  lng = wcs->lng;
  lat = wcs->lat;
  for (i=0; i<naxis; i++)
    raw[i] = (wcs->naxisn[i]+1.0)/2.0;
  if (raw_to_wcs(wcs, raw, world))
    {
/*-- Oops no mapping there! So explore the image in an increasingly large */
/*-- domain to find  a better "center" (now we know there must be angular */
/*-- coordinates) */
    for (j=0; j<100; j++)
      {
      for (i=0; i<naxis; i++)
        raw[i] += wcs->naxisn[i]/100.0*(0.5-(double)rand()/RAND_MAX);      
      if (!raw_to_wcs(wcs, raw, world))
        break;
      }
    }

  if (lng!=lat)
    lc = world[lng];
  else
    {
    lc = 0.0;   /* to avoid gcc -Wall warnings */
    lng = -1;
    }

/* Pixel scales at image center */
  scale = wcs->wcsscale;
  for (i=0; i<naxis; i++)
    {
    if ((i==lng || i==lat) && lng!=lat)
      wcs->pixscale = scale[i] = sqrt(wcs_scale(wcs, raw));
    else
      {
      raw[i] += 1.0;
      raw_to_wcs(wcs, raw, world2);
      scale[i] = fabs(world2[i] - world[i]);
      raw[i] -= 1.0;
      if (lng==lat)
        wcs->pixscale = scale[i];
      }
    wcs->wcsscalepos[i] = world[i];
    }


/* Find "World limits" */
  for (i=0; i<naxis; i++)
    {
    raw[i] = rawmin[i] = 0.5;
    step[i] = wcs->naxisn[i]/(WCS_NRANGEPOINTS-1.0);
    npoints *= WCS_NRANGEPOINTS;
    worldmax[i] = -(worldmin[i] = 1e31);
    linecount[i] = 0;
    }

  radmax = 0.0;
  worldc = wcs->wcsscalepos;

  for (j=npoints; j--;)
    {
    raw_to_wcs(wcs, raw, world);
/*-- Compute maximum distance to center */
    if ((rad=wcs_dist(wcs, world, worldc)) > radmax)
      radmax = rad;
    for (i=0; i<naxis; i++)
      {
/*---- Handle longitudes around 0 */
      if (i==lng)
        {
        world[i] -= lc;
        if (world[i]>180.0)
          world[i] -= 360.0;
        else if (world[i] <= -180.0)
          world[i] += 360.0;
        }
      if (world[i]<worldmin[i])
        worldmin[i] = world[i];
      if (world[i]>worldmax[i])
        worldmax[i] = world[i];
      }


    for (i=0; i<naxis; i++)
      {
      raw[i] += step[i];
      if (++linecount[i]<WCS_NRANGEPOINTS)
        break;
      else
        {
        linecount[i] = 0;       /* No need to initialize it to 0! */
        raw[i] = rawmin[i];
        }
      }
    }

  wcs->wcsmaxradius = radmax;

  if (lng!=lat)
    {
    worldmin[lng] = fmod_0_p360(worldmin[lng]+lc);
    worldmax[lng] = fmod_0_p360(worldmax[lng]+lc);
    if (worldmax[lat]<-90.0)
      worldmax[lat] = -90.0;
    if (worldmax[lat]>90.0)
      worldmax[lat] = 90.0;
    }

  return;
  }


/******* frame_wcs ***********************************************************
PROTO	int frame_wcs(wcsstruct *wcsin, wcsstruct *wcsout)
PURPOSE	Find the x and y limits of an input frame in an output image.
INPUT	WCS structure of the input frame,
	WCS structure of the output frame.
OUTPUT	1 if frames overlap, 0 otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/03/2012
 ***/
int	frame_wcs(wcsstruct *wcsin, wcsstruct *wcsout)

  {
   double		rawin[NAXIS], rawout[NAXIS], world[NAXIS];
   int			linecount[NAXIS];
   double		worldc;
   int			*min, *max,
			i,j, naxis, npoints, out, swapflag, overlapflag;

  naxis = wcsin->naxis;

/* World range */
  npoints = 1;
  min = wcsin->outmin;
  max = wcsin->outmax;
  for (i=0; i<naxis; i++)
    {
    rawin[i] = 0.5;	/* Lower pixel limits */
    npoints *= WCS_NRANGEPOINTS;
    max[i] = -(min[i] = 1<<30);
    linecount[i] = 0;
    }

/* Check if lng and lat are swapped between in and out wcs (vicious idea!) */
  swapflag = (((wcsin->lng != wcsout->lng) || (wcsin->lat != wcsout->lat))
	&& (wcsin->lng != wcsin->lat) && (wcsout->lng != wcsout->lat));

  for (j=npoints; j--;)
    {
    if (!raw_to_wcs(wcsin, rawin, world))
      {
      if (swapflag)
        {
        worldc = world[wcsout->lat];
        world[wcsout->lat] = world[wcsin->lat];
        world[wcsin->lat] = worldc;
        }
      if (!wcs_to_raw(wcsout, world, rawout))
        for (i=0; i<naxis; i++)
          {
          if ((out=(int)(rawout[i]+0.499))<min[i])
            min[i] = out;
          if (out>max[i])
            max[i] = out;
          }
      }
    for (i=0; i<naxis; i++)
      {
      rawin[i] = 0.5 + 0.5*wcsin->naxisn[i]
	*(1-cos(PI*(linecount[i]+1.0)/(WCS_NRANGEPOINTS-1)));
      if (++linecount[i]<WCS_NRANGEPOINTS)
        break;
      else
        {
        linecount[i] = 0;       /* No need to initialize it to 0! */
        rawin[i] =  0.5;
        }
      }
    }

/* Add a little margin, just in case... */
  for (i=0; i<naxis; i++)
    {
    if (min[i]>-2147483647)
      min[i] -= 2;
    if (max[i]>2147483647)
      max[i] += 2;
    }

/* Check overlap */
  overlapflag = 1;
  for (i=0; i<naxis; i++)
    if (min[i]>wcsout->naxisn[i] || max[i]<0)
      {
      overlapflag = 0;
      break;
      }

  return overlapflag;
  }


/******* reaxe_wcs ***********************************************************
PROTO	int reaxe_wcs(wcsstruct *wcs, int lng, int lat)
PURPOSE	Reformulate a wcs structure to match lng and lat axis indices
INPUT	WCS structure,
	longitude index,
	latitude index.
OUTPUT	RETURN_OK if successful, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	20/11/2003
 ***/
int	reaxe_wcs(wcsstruct *wcs, int lng, int lat)

  {
   char		strlng[80], strlat[80];
   double	dlng,dlat,dlng2,dlat2;
   int		l, ilng,ilat,olng,olat, naxis;

  olng = wcs->lng;
  olat = wcs->lat;
  if (lng<0 || lat<0 || olng<0 || olat<0)
    return RETURN_ERROR;

  ilng = wcs->naxisn[olng];
  ilat = wcs->naxisn[olat];
  wcs->naxisn[lng] = ilng;
  wcs->naxisn[lat] = ilat;
  strcpy(strlng, wcs->ctype[olng]);
  strcpy(strlat, wcs->ctype[olat]);
  strcpy(wcs->ctype[lng], strlng);
  strcpy(wcs->ctype[lat], strlat);
  dlng = wcs->crval[olng];
  dlat = wcs->crval[olat];
  wcs->crval[lng] = dlng;
  wcs->crval[lat] = dlat;
  naxis = wcs->naxis;
  dlng =  wcs->cd[olng+olng*naxis];
  dlng2 = wcs->cd[olng+olat*naxis];
  dlat =  wcs->cd[olat+olat*naxis];
  dlat2 = wcs->cd[olat+olng*naxis];
  wcs->cd[lng+lng*naxis] = dlng2;
  wcs->cd[lng+lat*naxis] = dlng;
  wcs->cd[lat+lat*naxis] = dlat2;
  wcs->cd[lat+lng*naxis] = dlat;
  for (l=0; l<100; l++)
    {
    dlng = wcs->projp[l+olng*100];
    dlat = wcs->projp[l+olat*100];
    wcs->projp[l+lng*100] = dlng;
    wcs->projp[l+lat*100] = dlat;
    }

/*-- Reinitialize wcs */
    wcsset(wcs->naxis,(const char(*)[9])wcs->ctype, wcs->wcsprm);

/*-- Initialize other WCS structures */
    init_wcs(wcs);
/*-- Find the range of coordinates */
    range_wcs(wcs);

  return RETURN_OK;
  }


/******* celsys_to_eq *********************************************************
PROTO	int celsys_to_eq(wcsstruct *wcs, double *wcspos)
PURPOSE	Convert arbitrary celestial coordinates to equatorial.
INPUT	WCS structure,
	Coordinate vector.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/02/2007
 ***/
int	celsys_to_eq(wcsstruct *wcs, double *wcspos)

  {
   double	*mat,
		a2,d2,sd2,cd2cp,sd,x,y;
   int		lng, lat;

  mat = wcs->celsysmat;
  a2 = wcspos[lng = wcs->wcsprm->lng]*DEG - mat[1];
  d2 = wcspos[lat = wcs->wcsprm->lat]*DEG;
/* A bit of spherical trigonometry... */
/* Compute the latitude... */
  sd2 = sin(d2);
  cd2cp = cos(d2)*mat[2];
  sd = sd2*mat[3]-cd2cp*cos(a2);
/* ...and the longitude */
  y = cd2cp*sin(a2);
  x = sd2 - sd*mat[3];
  wcspos[lng] = fmod((atan2(y,x) + mat[0])/DEG+360.0, 360.0);
  wcspos[lat] = asin(sd)/DEG;

  return RETURN_OK;
  }


/******* eq_to_celsys *********************************************************
PROTO	int eq_to_celsys(wcsstruct *wcs, double *wcspos)
PURPOSE	Convert equatorial to arbitrary celestial coordinates.
INPUT	WCS structure,
	Coordinate vector.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/02/2007
 ***/
int	eq_to_celsys(wcsstruct *wcs, double *wcspos)

  {
   double	*mat,
		a,d,sd2,cdcp,sd,x,y;
   int		lng, lat;

  mat = wcs->celsysmat;
  a = wcspos[lng = wcs->wcsprm->lng]*DEG - mat[0];
  d = wcspos[lat = wcs->wcsprm->lat]*DEG;
/* A bit of spherical trigonometry... */
/* Compute the latitude... */
  sd = sin(d);
  cdcp = cos(d)*mat[2];
  sd2 = sd*mat[3]+cdcp*cos(a);
/* ...and the longitude */
  y = cdcp*sin(a);
  x = sd2*mat[3]-sd;
  wcspos[lng] = fmod((atan2(y,x) + mat[1])/DEG+360.0, 360.0);
  wcspos[lat] = asin(sd2)/DEG;

  return RETURN_OK;
  }


/******* raw_to_wcs ***********************************************************
PROTO	int raw_to_wcs(wcsstruct *, double *, double *)
PURPOSE	Convert raw (pixel) coordinates to WCS (World Coordinate System).
INPUT	WCS structure,
	Pointer to the array of input coordinates,
	Pointer to the array of output coordinates.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/02/2007
 ***/
int	raw_to_wcs(wcsstruct *wcs, double *pixpos, double *wcspos)

  {
   double	imgcrd[NAXIS],
		phi,theta;
   int		i;

  if (wcsrev((const char(*)[9])wcs->ctype, wcs->wcsprm, pixpos,
	wcs->lin,imgcrd, wcs->prj, &phi, &theta, wcs->crval, wcs->cel, wcspos))
    {
    for (i=0; i<wcs->naxis; i++)
      wcspos[i] = WCS_NOCOORD;
    return RETURN_ERROR;
    }

/* If needed, convert from a different coordinate system to equatorial */
  if (wcs->celsysconvflag)
    celsys_to_eq(wcs, wcspos);

  return RETURN_OK;
  }


/******* wcs_to_raw ***********************************************************
PROTO	int wcs_to_raw(wcsstruct *, double *, double *)
PURPOSE	Convert WCS (World Coordinate System) coords to raw (pixel) coords.
INPUT	WCS structure,
	Pointer to the array of input coordinates,
	Pointer to the array of output coordinates.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	08/02/2007
 ***/
int	wcs_to_raw(wcsstruct *wcs, double *wcspos, double *pixpos)

  {
   double	imgcrd[NAXIS],
		phi,theta;
   int		i;

/* If needed, convert to a coordinate system different from equatorial */
  if (wcs->celsysconvflag)
    eq_to_celsys(wcs, wcspos);

  if (wcsfwd((const char(*)[9])wcs->ctype,wcs->wcsprm,wcspos,
	wcs->crval, wcs->cel,&phi,&theta,wcs->prj, imgcrd,wcs->lin,pixpos))
    {
    for (i=0; i<wcs->naxis; i++)
      pixpos[i] = WCS_NOCOORD;
    return RETURN_ERROR;
    }

  return RETURN_OK;
  }


/******* red_to_raw **********************************************************
PROTO	int red_to_raw(wcsstruct *, double *, double *)
PURPOSE	Convert reduced (World Coordinate System) coords to raw (pixel)
	coords.
INPUT	WCS structure,
	Pointer to the array of input (reduced) coordinates,
	Pointer to the array of output (pixel) coordinates.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/10/2003
 ***/
int	red_to_raw(wcsstruct *wcs, double *redpos, double *pixpos)

  {
   struct wcsprm	*wcsprm;
   double		offset;

  wcsprm = wcs->wcsprm;
/* Initialize if required */
  if (wcsprm && wcsprm->flag != WCSSET)
    {
    if (wcsset(wcs->naxis, (const char(*)[9])wcs->ctype, wcsprm))
      return RETURN_ERROR;
    }

  if (wcsprm && wcsprm->flag != 999 && wcsprm->cubeface != -1)
    {
/*-- Separation between faces */
    offset = (wcs->prj->r0 == 0.0 ? 90.0 : wcs->prj->r0*PI/2.0);
/*-- Stack faces in a cube */
    if (redpos[wcs->lat] < -0.5*offset)
      {
      redpos[wcs->lat] += offset;
      redpos[wcsprm->cubeface] = 5.0;
      }
    else if (redpos[wcs->lat] > 0.5*offset)
      {
      redpos[wcs->lat] -= offset;
      redpos[wcsprm->cubeface] = 0.0;
      }
    else if (redpos[wcs->lng] > 2.5*offset)
      {
      redpos[wcs->lng] -= 3.0*offset;
      redpos[wcsprm->cubeface] = 4.0;
      }
    else if (redpos[wcs->lng] > 1.5*offset)
      {
      redpos[wcs->lng] -= 2.0*offset;
      redpos[wcsprm->cubeface] = 3.0;
      }
    else if (redpos[wcs->lng] > 0.5*offset)
      {
      redpos[wcs->lng] -= offset;
      redpos[wcsprm->cubeface] = 2.0;
      }
    else
      redpos[wcsprm->cubeface] = 1.0;
    }

/* Apply forward linear transformation */
  if (linfwd(redpos, wcs->lin, pixpos))
      return RETURN_ERROR;

  return RETURN_OK;
  }


/******* raw_to_red **********************************************************
PROTO	int raw_to_red(wcsstruct *, double *, double *)
PURPOSE	Convert raw (pixel) coordinates to reduced WCS coordinates.
INPUT	WCS structure,
	Pointer to the array of input (pixel) coordinates,
	Pointer to the array of output (reduced) coordinates.
OUTPUT	RETURN_OK if mapping successful, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/10/2003
 ***/
int	raw_to_red(wcsstruct *wcs, double *pixpos, double *redpos)

  {
   struct wcsprm	*wcsprm;
   double		offset;
   int			face;

  wcsprm = wcs->wcsprm;
/* Initialize if required */
  if (wcsprm && wcsprm->flag != WCSSET)
    {
    if (wcsset(wcs->naxis, (const char(*)[9])wcs->ctype, wcsprm))
      return RETURN_ERROR;
    }
   
/* Apply reverse linear transformation */
   if (linrev(pixpos, wcs->lin, redpos))
     return RETURN_ERROR;

  if (wcsprm && wcsprm->flag != 999 && wcsprm->cubeface != -1)
    {
/*-- Do we have a CUBEFACE axis? */
    face = (int)(redpos[wcsprm->cubeface] + 0.5);
    if (fabs(redpos[wcsprm->cubeface]-face) > 1e-10)
      return RETURN_ERROR;

/*-- Separation between faces. */
    offset = (wcs->prj->r0 == 0.0 ? 90.0 : wcs->prj->r0*PI/2.0);
/*-- Lay out faces in a plane. */
    switch (face)
      {
      case 0:
        redpos[wcs->lat] += offset;
        break;
      case 1:
        break;
      case 2:
        redpos[wcs->lng] += offset;
        break;
      case 3:
        redpos[wcs->lng] += offset*2;
        break;
      case 4:
        redpos[wcs->lng] += offset*3;
        break;
      case 5:
        redpos[wcs->lat] -= offset;
        break;
      default:
        return RETURN_ERROR;
      }
    }

  return RETURN_OK;
  }


/******* wcs_dist ***********************************************************
PROTO	double wcs_dist(wcsstruct *wcs, double *wcspos1, double *wcspos2)
PURPOSE	Compute the angular distance between 2 points on the sky.
INPUT	WCS structure,
	Pointer to the first array of world coordinates,
	Pointer to the second array of world coordinates.
OUTPUT	Angular distance (in degrees) between points.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	24/07/2002
 ***/
double	wcs_dist(wcsstruct *wcs, double *wcspos1, double *wcspos2)

  {
  double	d, dp;
  int		i, lng, lat;

  lng = wcs->lng;
  lat = wcs->lat;
  if (lat!=lng)
    {
/*-- We are operating in angular coordinates */
    d = sin(wcspos1[lat]*DEG)*sin(wcspos2[lat]*DEG)
	+ cos(wcspos1[lat]*DEG)*cos(wcspos2[lat]*DEG)
		*cos((wcspos1[lng]-wcspos2[lng])*DEG);
    return d>-1.0? (d<1.0 ? acos(d)/DEG : 0.0) : 180.0;
    }
  else
    {
    d = 0.0;
    for (i=0; i<wcs->naxis; i++)
      {
      dp = wcspos1[i] - wcspos2[i];
      d += dp*dp;
      }
    return sqrt(d);
    }
  }


/******* wcs_scale ***********************************************************
PROTO	double wcs_scale(wcsstruct *wcs, double *pixpos)
PURPOSE	Compute the sky area equivalent to a local pixel.
INPUT	WCS structure,
	Pointer to the array of local raw coordinates,
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2008
 ***/
double	wcs_scale(wcsstruct *wcs, double *pixpos)

  {
   double	wcspos[NAXIS], wcspos1[NAXIS], wcspos2[NAXIS], pixpos2[NAXIS];
   double	dpos1,dpos2;
   int		lng, lat;

  if (raw_to_wcs(wcs, pixpos, wcspos))
    return 0.0;

  lng = wcs->lng;
  lat = wcs->lat;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }

/* Compute pixel scale */
  pixpos2[lng] = pixpos[lng] + 1.0;
  pixpos2[lat] = pixpos[lat];
  if (raw_to_wcs(wcs, pixpos2, wcspos1))
    return 0.0;
  pixpos2[lng] -= 1.0;
  pixpos2[lat] += 1.0;
  if (raw_to_wcs(wcs, pixpos2, wcspos2))
    return 0.0;
  dpos1 = wcspos1[lng]-wcspos[lng];
  dpos2 = wcspos2[lng]-wcspos[lng];
  if (wcs->lng!=wcs->lat)
    {
    if (dpos1>180.0)
      dpos1 -= 360.0;
    else if (dpos1<-180.0)
      dpos1 += 360.0;
    if (dpos2>180.0)
      dpos2 -= 360.0;
    else if (dpos2<-180.0)
      dpos2 += 360.0;
    return fabs((dpos1*(wcspos2[lat]-wcspos[lat])
		-(wcspos1[lat]-wcspos[lat])*dpos2)*cos(wcspos[lat]*DEG));
    }
  else
    return fabs((dpos1*(wcspos2[lat]-wcspos[lat])
		-(wcspos1[lat]-wcspos[lat])*dpos2));
  }


/****** wcs jacobian *********************************************************
PROTO	double wcs_jacobian(wcsstruct *wcs, double *pixpos, double *jacob)
PURPOSE	Compute the local Jacobian matrix of the astrometric deprojection.
INPUT	WCS structure,
	Pointer to the array of local raw coordinates,
	Pointer to the jacobian array (output).
OUTPUT	Determinant over spatial coordinates (=pixel area), or -1.0 if mapping
	was unsuccesful.
NOTES   Memory must have been allocated (naxis*naxis*sizeof(double)) for the
        Jacobian array.
AUTHOR	E. Bertin (IAP)
VERSION	11/10/2007
 ***/
double	wcs_jacobian(wcsstruct *wcs, double *pixpos, double *jacob)
  {
   double	pixpos0[NAXIS], wcspos0[NAXIS], wcspos[NAXIS],
		dpos;
   int		i,j, lng,lat,naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    pixpos0[i] = pixpos[i];
  if (raw_to_wcs(wcs, pixpos0, wcspos0) == RETURN_ERROR)
    return -1.0;
  for (i=0; i<naxis; i++)
    {
    pixpos0[i] += 1.0;
    if (raw_to_wcs(wcs, pixpos0, wcspos) == RETURN_ERROR)
      return -1.0;
    pixpos0[i] -= 1.0;
    for (j=0; j<naxis; j++)
      {
      dpos = wcspos[j]-wcspos0[j];
      if (lng!=lat && j==lng)
        {
        if (dpos>180.0)
          dpos -= 360.0;
        else if (dpos<-180.0)
          dpos += 360.0;
        dpos *= cos(wcspos0[lat]*DEG);
        }
      jacob[j*naxis+i] = dpos;
      }
    }

  if (lng==lat)
    {
    lng = 0;
    lat = 1;
    }

  return fabs(jacob[lng+naxis*lng]*jacob[lat+naxis*lat]
		- jacob[lat+naxis*lng]*jacob[lng+naxis*lat]);
  }


/****** wcs rawtoraw *********************************************************
PROTO	double wcs_rawtoraw(wcsstruct *wcsin, wcsstruct *wcsout,
		double *pixpos, double *jacob)
PURPOSE	Convert raw coordinates from input wcs structure to raw coordinates
	from output wcs structure, and the local Jacobian matrix of the
	reprojection.
INPUT	WCS input structure,
	WCS output structure,
	pointer to the array of local input raw coordinates,
	pointer to the array of local output raw coordinates (output),
	pointer to the jacobian array (output).
OUTPUT	Determinant over spatial coordinates (ratio of pixel areas),
	0.0 if jacob is NULL (Jacobian not computed in that case),
	or -1.0 if mapping was unsuccesful.
NOTES   Memory must have been allocated (naxis*naxis*sizeof(double)) for the
        Jacobian array.
AUTHOR	E. Bertin (IAP)
VERSION	12/06/2012
 ***/
double	wcs_rawtoraw(wcsstruct *wcsin, wcsstruct *wcsout,
		double *pixposin, double *pixposout, double *jacob)
  {
   double	pixpos0[NAXIS], pixpos2[NAXIS], wcspos[NAXIS];
   int		i,j, lng,lat,naxis;

  naxis = wcsin->naxis;
  for (i=0; i<naxis; i++)
    pixpos0[i] = pixposin[i];

/* Coordinate transformation */
  if (raw_to_wcs(wcsin, pixposin, wcspos) == RETURN_ERROR)
    return -1.0;
  if (wcs_to_raw(wcsout, wcspos, pixposout) == RETURN_ERROR)
    return -1.0;
 
/* Jacobian */
  if (!jacob)
    return 0.0;

  lng = wcsin->lng;
  lat = wcsin->lat;
  for (i=0; i<naxis; i++)
    {
    pixpos0[i] += 1.0;
    if (raw_to_wcs(wcsin, pixpos0, wcspos) == RETURN_ERROR)
      return -1.0;
    if (wcs_to_raw(wcsout, wcspos, pixpos2) == RETURN_ERROR)
      return -1.0;
    pixpos0[i] -= 1.0;
    for (j=0; j<naxis; j++)
      jacob[j*naxis+i] = pixpos2[j] - pixposout[j];
    }

  if (lng==lat)
    {
    lng = 0;
    lat = 1;
    }

  return fabs(jacob[lng+naxis*lng]*jacob[lat+naxis*lat]
		- jacob[lat+naxis*lng]*jacob[lng+naxis*lat]);
  }


/******* wcs_chirality *******************************************************
PROTO	int wcs_chirality(wcsstruct *wcs)
PURPOSE	Compute the chirality of a WCS projection.
INPUT	WCS structure.
OUTPUT	+1 if determinant of matrix is positive, -1 if negative, 0 if null.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/09/2006
 ***/
int	wcs_chirality(wcsstruct *wcs)

  {
   double	a;
   int		lng,lat, naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;
  if (lng==lat && naxis>=2)
    {
    lng = 0;
    lat = 1;
    }

  a = wcs->cd[lng*naxis+lng]*wcs->cd[lat*naxis+lat]
	- wcs->cd[lng*naxis+lat]*wcs->cd[lat*naxis+lng];
  return a>TINY? 1 : (a<-TINY? -1 : 0);
  }


/****** precess_wcs **********************************************************
PROTO	void precess_wcs(wcsstruct *wcs, double yearin, double yearout)
PURPOSE	Precess the content of a WCS structure according to the equinox.
INPUT	WCS structure,
	Input year,
	Output year.
OUTPUT	-.
NOTES	Epoch for coordinates should be J2000 (FK5 system).
AUTHOR	E. Bertin (IAP)
VERSION	04/01/2008
 ***/
void	precess_wcs(wcsstruct *wcs, double yearin, double yearout)

  {
   double	crval[NAXIS],a[NAXIS*NAXIS],b[NAXIS*NAXIS],
		*c,*at,
		val, cas, sas, angle, dalpha;
   int		i,j,k, lng,lat, naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  if (lat==lng || yearin==yearout)
    return;
  naxis = wcs->naxis;
/* Precess to year out */
  precess(yearin, wcs->crval[lng], wcs->crval[lat], yearout,
	&crval[lng], &crval[lat]);

  dalpha = (crval[lng] - wcs->crval[lng])*DEG;

/* Compute difference angle with the north axis between start and end */
  angle = (dalpha!=0.0 && (crval[lat] - wcs->crval[lat])*DEG != 0.0) ?
	180.0 - (atan2(sin(dalpha),
	cos(crval[lat]*DEG)*tan(wcs->crval[lat]*DEG)
	- sin(crval[lat]*DEG)*cos(dalpha))
	+ atan2(sin(dalpha),
	cos(wcs->crval[lat]*DEG)*tan(crval[lat]*DEG)
	- sin(wcs->crval[lat]*DEG)*cos(dalpha)))/DEG
	: 0.0;

/* A = C*B */
  c = wcs->cd;
/* The B matrix is made of 2 numbers */

  cas = cos(angle*DEG);
  sas = sin(-angle*DEG);
  for (i=0; i<naxis; i++)
    b[i+i*naxis] = 1.0;
  b[lng+lng*naxis] = cas;
  b[lat+lng*naxis] = -sas;
  b[lng+lat*naxis] = sas;
  b[lat+lat*naxis] = cas;
  at = a;
  for (j=0; j<naxis; j++)
    for (i=0; i<naxis; i++)
      {
      val = 0.0;
      for (k=0; k<naxis; k++)
        val += c[k+j*naxis]*b[i+k*naxis];
      *(at++) = val;
      }

  at = a;

  for (i=0; i<naxis*naxis; i++)
    *(c++) = *(at++);

  wcs->crval[lng] = crval[lng];
  wcs->crval[lat] = crval[lat];
  wcs->equinox = yearout;

/* Initialize other WCS structures */
  init_wcs(wcs);
/* Find the range of coordinates */
  range_wcs(wcs);
/* Invert projection corrections */
  invert_wcs(wcs);

  return;
  }


/********************************* precess ***********************************/
/*
precess equatorial coordinates according to the equinox (from Ephemerides du
Bureau des Longitudes 1992). Epoch for coordinates should be J2000
(FK5 system).
*/
void	precess(double yearin, double alphain, double deltain,
		double yearout, double *alphaout, double *deltaout)

  {
   double	dzeta,theta,z, t1,t1t1, t2,t2t2,t2t2t2,
		cddsadz, cddcadz, cdd, sdd, adz, cdin,sdin,ct,st,caindz;

  alphain *= DEG;
  deltain *= DEG;

  t1 = (yearin - 2000.0)/1000.0;
  t2 = (yearout - yearin)/1000.0;
  t1t1 = t1*t1;
  t2t2t2 = (t2t2 = t2*t2)*t2;
  theta = (97171.735e-06 - 413.691e-06*t1 - 1.052e-06 * t1t1) * t2
	+ (-206.846e-06 - 1.052e-06*t1) * t2t2 - 202.812e-06 * t2t2t2;
  dzeta = (111808.609e-06 + 677.071e-06*t1 - 0.674e-06 * t1t1) * t2
	+ (146.356e-06 - 1.673e-06*t1) * t2t2 + 87.257e-06 * t2t2t2;
  z = (111808.609e-06 +677.071e-06*t1 - 0.674e-06 * t1t1) * t2
	+ (530.716e-06 + 0.320e-06*t1) * t2t2 + 88.251e-06 * t2t2t2;
  cddsadz = (cdin=cos(deltain)) * sin(alphain+dzeta);
  cddcadz = -(sdin=sin(deltain))*(st=sin(theta))
	+cdin*(ct=cos(theta))*(caindz=cos(alphain+dzeta));
  sdd = sdin*ct + cdin*st*caindz;
  cdd = cos(*deltaout = asin(sdd));
  adz = asin(cddsadz/cdd);
  if (cddcadz<0.0)
    adz = PI - adz;
  if (adz<0.0)
    adz += 2.0*PI;
  adz += z;
  *alphaout = adz/DEG;
  *deltaout /= DEG;

  return;
  }


/********************************* b2j ***********************************/
/*
conver equatorial coordinates from equinox and epoch B1950 to equinox and
epoch J2000 for extragalactic sources (from Aoki et al. 1983).
*/
void	b2j(double yearobs, double alphain, double deltain,
		double *alphaout, double *deltaout)
  {
   int		i,j;
   double	a[3] = {-1.62557e-6, -0.31919e-6, -0.13843e-6},
		ap[3] = {1.245e-3, -1.580e-3, -0.659e-3},
		m[6][6] = {
  { 0.9999256782,     -0.0111820611,     -0.0048579477,
    0.00000242395018, -0.00000002710663, -0.00000001177656},
  { 0.0111820610,      0.9999374784,     -0.0000271765,
    0.00000002710663,  0.00000242397878, -0.00000000006587},
  { 0.0048579479,     -0.0000271474,      0.9999881997,
    0.00000001177656, -0.00000000006582,  0.00000242410173},
  {-0.000551,        -0.238565,           0.435739,
    0.99994704,	     -0.01118251,        -0.00485767},
  { 0.238514,        -0.002662,          -0.008541,
    0.01118251,	      0.99995883,        -0.00002718},
  {-0.435623,         0.012254,           0.002117,
    0.00485767,      -0.00002714,         1.00000956}},
 		a1[3], r[3], ro[3], r1[3], r2[3], v1[3], v[3];
   double		cai, sai, cdi, sdi, dotp, rmod, alpha, delta,
			t1 = (yearobs - 1950.0)/100.0;

  alphain *= PI/180.0;
  deltain *= PI/180.0;
  cai = cos(alphain);
  sai = sin(alphain);
  cdi = cos(deltain);
  sdi = sin(deltain);
  ro[0] = cdi*cai;
  ro[1] = cdi*sai;
  ro[2] = sdi;
  dotp = 0.0;
  for (i=0; i<3; i++)
    {
    a1[i] = a[i]+ap[i]*ARCSEC*t1;
    dotp += a1[i]*ro[i];
    }
  for (i=0; i<3; i++)
    {
    r1[i] = ro[i] - a1[i] + dotp*ro[i];
    r[i] = v[i] = v1[i] = 0.0;
    }
  for (j=0; j<6; j++)
    for (i=0; i<6; i++)
      {
      if (j<3)
        r[j] += m[j][i]*(i<3?r1[i]:v1[i-3]);
      else
         v[j-3] += m[j][i]*(i<3?r1[i]:v1[i-3]);
      }
  rmod = 0.0;
  for (i=0; i<3; i++)
    {
    r2[i] = r[i]+v[i]*ARCSEC*(t1-0.5);
    rmod += r2[i]*r2[i];
    }
  rmod = sqrt(rmod);
  delta = asin(r2[2]/rmod);
  alpha = acos(r2[0]/cos(delta)/rmod);
  if (r2[1]<0)
    alpha = 2*PI - alpha;
  *alphaout = alpha*180.0/PI;
  *deltaout = delta*180.0/PI;

  return;			
  }


/*********************************** j2b *************************************/
/*
conver equatorial coordinates from equinox and epoch J2000 to equinox and
epoch B1950 for extragalactic sources (from Aoki et al. 1983, after
inversion of their matrix and some custom arrangements).
*/
void    j2b(double yearobs, double alphain, double deltain,
	double *alphaout, double *deltaout)
  {
   int			i,j;
   double		a[3] = {-1.62557e-6, -0.31919e-6, -0.13843e-6},
                        ap[3] = {1.245e-3, -1.580e-3, -0.659e-3},
                        m[6][6] = {
  { 0.9999256794678425,    0.01118148281196562,   0.004859003848996022,
   -2.423898417033081e-06,-2.710547600126671e-08,-1.177738063266745e-08},
  {-0.01118148272969232,   0.9999374849247641,   -2.717708936468247e-05,
    2.710547578707874e-08,-2.423927042585208e-06, 6.588254898401055e-11},
  {-0.00485900399622881,  -2.715579322970546e-05, 0.999988194643078,
    1.177738102358923e-08, 6.582788892816657e-11,-2.424049920613325e-06},
  {-0.0005508458576414713, 0.2384844384742432,   -0.4356144527773499,
    0.9999043171308133,    0.01118145410120206,   0.004858518651645554},
  {-0.2385354433560954,   -0.002664266996872802,  0.01225282765749546,
   -0.01118145417187502,   0.9999161290795875,   -2.717034576263522e-05},
  { 0.4357269351676567,   -0.008536768476441086,  0.002113420799663768,
   -0.004858518477064975, -2.715994547222661e-05, 0.9999668385070383}},
			a1[3], r[3], ro[3], r1[3], r2[3], v1[3], v[3];
   double		cai, sai, cdi, sdi, dotp, rmod, alpha, delta, t1;

/* Convert Julian years from J2000.0 to tropic centuries from B1950.0 */
  t1 = ((yearobs - 2000.0) + (MJD2000 - MJD1950)/365.25)*JU2TROP/100.0;
  alphain *= DEG;
  deltain *= DEG;
  cai = cos(alphain);
  sai = sin(alphain);
  cdi = cos(deltain);
  sdi = sin(deltain);
  r[0] = cdi*cai;
  r[1] = cdi*sai;
  r[2] = sdi;
  for (i=0; i<3; i++)
    v[i] = r2[i] = v1[i] = 0.0;
  for (j=0; j<6; j++)
    for (i=0; i<6; i++)
      if (j<3)
        r2[j] += m[j][i]*(i<3?r[i]:v[i-3]);
      else
        v1[j-3] += m[j][i]*(i<3?r[i]:v[i-3]);

  for (i=0; i<3; i++)
    r1[i] = r2[i]+v1[i]*ARCSEC*t1;

  dotp = 0.0;
  for (i=0; i<3; i++)
    {
    a1[i] = a[i]+ap[i]*ARCSEC*t1;
    dotp += a1[i]*(r1[i]+a1[i]);
    }
  dotp = 2.0/(sqrt(1+4.0*dotp)+1.0);
  rmod = 0.0;
  for (i=0; i<3; i++)
    {
    ro[i] = dotp*(r1[i]+a1[i]);
    rmod += ro[i]*ro[i];
    }
  rmod = sqrt(rmod);
  delta = asin(ro[2]/rmod);
  alpha = acos(ro[0]/cos(delta)/rmod);
  if (ro[1]<0)
    alpha = 2.0*PI - alpha;
  *alphaout = alpha/DEG;
  *deltaout = delta/DEG;

  return;
  }


/******************************** degtosexal *********************************/
/*
Convert degrees to hh mm ss.xx alpha coordinates.
*/
char    *degtosexal(double alpha, char *str)

  {
   int		hh, mm;
   double	ss;

  if (alpha>=0.0 && alpha <360.0)
    {
    hh = (int)(alpha/15.0);
    mm = (int)(60.0*(alpha/15.0 - hh));
    ss = 60.0*(60.0*(alpha/15.0 - hh) - mm);
    }
  else
    hh = mm = ss = 0.0;
  sprintf(str,"%02d:%02d:%05.2f", hh, mm, ss);

  return str;
  }


/******************************** degtosexde *********************************/
/*
Convert degrees to dd dm ds.x delta coordinates.
*/
char    *degtosexde(double delta, char *str)

  {
   char		sign;
   double	ds;
   int		dd, dm;

  sign = delta<0.0?'-':'+';
  delta = fabs(delta);
  if (delta>=-90.0 && delta <=90.0)
    {
    dd = (int)delta;
    dm = (int)(60.0*(delta - dd));
    ds = 60.0*fabs(60.0*(delta - dd) - dm);
    }
  else
    dd = dm = ds = 0.0;
  sprintf(str,"%c%02d:%02d:%04.1f", sign, dd, dm, ds);
  return str;
  }


/******************************** sextodegal *********************************/
/*
Convert hh mm ss.xxx alpha coordinates to degrees.
*/
double  sextodegal(char *hms)
  {
   double	val;
   char		*ptr;

  val = atof(strtok_r(hms, ": \t", &ptr))*15.0;		/* Hours */
  val += atof(strtok_r(NULL, ": \t", &ptr))/4.0;	/* Minutes */
  val += atof(strtok_r(NULL, ": \t", &ptr))/240.0;	/* Seconds */

  return val;
  }


/******************************** sextodegde *********************************/
/*
Convert dd dm ds.xxx delta coordinates to degrees.
*/
double  sextodegde(char *dms)
  {
   double	val, sgn;
   char		*str, *ptr;

  str = strtok_r(dms, ": \t", &ptr);
  sgn = (strchr(str, '-') ? -1.0:1.0);
  val = atof(dms);					/* Degrees */
  val += atof(strtok_r(NULL, ": \t", &ptr))*sgn/60.0;	/* Minutes */
  val += atof(strtok_r(NULL, ": \t", &ptr))*sgn/3600.0;	/* Seconds */

  return val;
  }


/******************************** fmod_0_p360 *******************************/
/*
Fold input angle in the [0,+360[ domain.
*/
double  fmod_0_p360(double angle)
  {
  return angle>0.0? fmod(angle,360.0) : fmod(angle,360.0)+360.0;
  }


/******************************** fmod_m90_p90 *******************************/
/*
Fold input angle in the [-90,+90[ domain.
*/
double  fmod_m90_p90(double angle)
  {
  return angle>0.0? fmod(angle+90.0,180.0)-90.0 : fmod(angle-90.0,180.0)+90.0;
  }


/********************************* fcmp_0_p360 *******************************/
/*
Compare angles in the [0,+360[ domain: return 1 if anglep>anglem, 0 otherwise.
*/
int  fcmp_0_p360(double anglep, double anglem)
  {
   double dval = anglep - anglem;

  return (int)((dval>0.0 && dval<180.0) || dval<-180.0);
  }


