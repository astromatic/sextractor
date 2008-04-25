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
*	Last modify:	24/04/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

/*----------------------------- Internal constants --------------------------*/

#define		DEG	(PI/180.0)	/* 1 deg in radians */
#define		ARCSEC	(DEG/3600.0)	/* 1 arcsec in radians */
#define		MJD2000	51544.50000	/* Modified Julian date for J2000.0 */
#define		MJD1950	33281.92346	/* Modified Julian date for B1950.0 */
#define		JU2TROP	1.0000214	/* 1 Julian century in tropical units*/
#define		MAMA_CORFLEX	3.3e-5	/* MAMA coordinate correction factor */

/*------------------------------- structures --------------------------------*/
/*------------------------------- functions ---------------------------------*/
extern void		astrom_errparam(picstruct *, objstruct *),
			astrom_profshapeparam(picstruct *, objstruct *),
			astrom_shapeparam(picstruct *, objstruct *),
			astrom_winerrparam(picstruct *, objstruct *),
			astrom_winshapeparam(picstruct *, objstruct *),
			computeastrom(picstruct *, objstruct *),
			initastrom(picstruct *),
			j2b(double, double, double, double *, double *),
			precess(double,double,double,double,double *,double *);

extern double		*compute_wcs(picstruct *, double, double);
