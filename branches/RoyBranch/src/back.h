/*
*				back.h
*
* Include file for back.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		09/03/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*----------------------------- Internal constants --------------------------*/
#define	BACK_BUFSIZE		1048576		/* bkgnd buffer */
#define	BACK_MINGOODFRAC	0.5		/* min frac with good weights*/
#define	QUANTIF_NSIGMA		5		/* histogram limits */
#define	QUANTIF_NMAXLEVELS	4096		/* max nb of quantif. levels */
#define	QUANTIF_AMIN		4		/* min nb of "mode pixels" */

#define	BACK_WSCALE		1		/* Activate weight scaling */
#define	BACK_NOWSCALE		0		/* No weight scaling */

/* NOTES:
One must have:		BACK_BUFSIZE >= MAXPICSIZE
			0 < QUANTIF_NSIGMA <= 10
			QUANTIF_AMIN > 0
*/

/*------------------------------- structures --------------------------------*/
/* Background info */
typedef struct structback
  {
  float		mode, mean, sigma;	/* Background mode, mean and sigma */
  LONG		*histo;			/* Pointer to a histogram */
  int		nlevels;		/* Nb of histogram bins */
  float		qzero, qscale;		/* Position of histogram */
  float		lcut, hcut;		/* Histogram cuts */
  int		npix;			/* Number of pixels involved */
  }	backstruct;


/*------------------------------- functions ---------------------------------*/
void		backhisto(backstruct *, backstruct *, PIXTYPE *, PIXTYPE *,
			size_t, int, int, int, PIXTYPE),
		backstat(backstruct *, backstruct *, PIXTYPE *, PIXTYPE *,
			size_t, int, int, int, PIXTYPE),
		backrmsline(picstruct *, int, PIXTYPE *),
		copyback(picstruct *infield, picstruct *outfield),
		endback(picstruct *),
		filterback(picstruct *),
		makeback(picstruct *, picstruct *, int),
		subbackline(picstruct *, int, PIXTYPE *);

float		backguess(backstruct *, float *, float *),
		localback(picstruct *, objstruct *),
		*makebackspline(picstruct *, float *);

extern PIXTYPE	back(picstruct *, int, int);
