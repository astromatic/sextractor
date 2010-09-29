 /*
 				winpos.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN, IAP
*
*	Contents:	Include file for winpos.c.
*
*	Last modify:	16/07/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

#define	WINPOS_NITERMAX	16	/* Maximum number of steps */
#define	WINPOS_NSIG	4	/* Measurement radius */
#define	WINPOS_OVERSAMP	11	/* oversampling in each dimension */
#define	WINPOS_STEPMIN	0.0001	/* Minimum change in position for continueing*/
#define	WINPOS_FAC	2.0	/* Centroid offset factor (2 for a Gaussian) */

/* NOTES:
One must have:
	WINPOS_NITERMAX >= 1
	WINPOS_OVERSAMP >= 1
*/

/*------------------------------- functions ---------------------------------*/
extern void	compute_winpos(picstruct *field, picstruct *wfield,
			       objstruct *obj);
