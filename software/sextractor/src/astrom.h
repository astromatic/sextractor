 /*
 				astrom.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN, IAP & Leiden observatory
*
*	Contents:	Astrometrical stuff.
*
*	Last modify:	13/07/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

#define		DEG	(PI/180.0)	/* 1 deg in radians */
#define		ARCSEC	(DEG/3600.0)	/* 1 arcsec in radians */
#define		MJD2000	51544.50000	/* Modified Julian date for J2000.0 */
#define		MJD1950	33281.92346	/* Modified Julian date for B1950.0 */
#define		JU2TROP	1.0000214	/* 1 Julian century in tropical units*/
#define		NAXIS	3		/* Max number of FITS axes */
#define		MAMA_CORFLEX	3.3e-5	/* MAMA coordinate correction factor */

/*------------------------------- structures --------------------------------*/

typedef struct structastrom
  {
  int		naxis;			/* Number of image axes */

  char		ctype[NAXIS][9];	/* FITS CTYPE strings */
  char		cunit[NAXIS][32];	/* FITS CUNIT strings */
  double	crval[NAXIS];		/* FITS CRVAL parameters */
  double	cdelt[NAXIS];		/* FITS CDELT parameters */
  double	crpix[NAXIS];		/* FITS CRPIX parameters */
  double	projp[100*NAXIS];	/* FITS PROJP parameters */
  double	longpole,latpole;	/* FITS LONGPOLE and LATPOLE */
  double	pc[NAXIS*NAXIS];	/* FITS PC matrix */
  double	linmat[NAXIS*NAXIS];	/* Local linear mapping matrix */
  double	lindet;			/* Determinant of the local matrix */
  double	pixscale;		/* (Local) pixel scale */
  double	ap2000,dp2000;		/* J2000 coordinates of pole */
  double	ap1950,dp1950;		/* B1950 coordinates of pole */
  double	equinox;		/* Equinox of observations */
  enum {RDSYS_ICRS, RDSYS_FK5, RDSYS_FK4, RDSYS_FK4_NO_E, RDSYS_GAPPT}
		radecsys;		/* FITS RADECSYS reference frame */
  int		wcs_flag;		/* WCSLIB: can it be used? */
  int		lat,lng;		/* longitude and latitude axes # */
  double	r0;			/* projection "radius" */
  struct wcsprm	*wcs;			/* WCSLIB's wcsprm structure */
  struct linprm	*lin;			/* WCSLIB's linprm structure */
  struct celprm	*cel;			/* WCSLIB's celprm structure */
  struct prjprm *prj;			/* WCSLIB's prjprm structure */
  struct tnxaxis *tnx_latcor;		/* IRAF's TNX latitude corrections */
  struct tnxaxis *tnx_lngcor;		/* IRAF's TNX longitude corrections */
  }	astromstruct;

/*------------------------------- functions ---------------------------------*/
extern void		astrom_errparam(picstruct *, objstruct *),
			astrom_winerrparam(picstruct *, objstruct *),
			astrom_shapeparam(picstruct *, objstruct *),
			astrom_winshapeparam(picstruct *, objstruct *),
			computeastrom(picstruct *, objstruct *),
			copyastrom(picstruct *infield, picstruct *outfield),
			endastrom(picstruct *),
			initastrom(picstruct *),
			j2b(double, double, double, double *, double *),
			precess(double,double,double,double,double *,double *);

extern double		*compute_wcs(picstruct *, double, double);
