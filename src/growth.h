 /*
 				growth.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN, IAP & Leiden observatory
*
*	Contents:	Include file for growth.c.
*
*	Last modify:	02/07/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*----------------------------- Internal constants --------------------------*/

#define	GROWTH_NSTEP	64	/* number of growth curve samples */
#define	GROWTH_OVERSAMP	5	/* pixel oversampling in each dimension */
#define	GROWTH_NSIG	3*MARGIN_SCALE	/* MAG_AUTO analysis range (nsigmas) */
#define	GROWTH_MINHLRAD	0.5	/* Minimum internal half-light radius (pixels)*/

/* NOTES:
One must have:	GROWTH_SAMP >= 1
		GROWTH_OVERSAMP >= 1
		GROWTH_OVERSAMPRADIUS >= 0
		GROWTH_NSIG > 0.0
*/

/*------------------------------- functions ---------------------------------*/
extern void	endgrowth(void),
		initgrowth(void),
		makeavergrowth(picstruct *field, picstruct *wfield,
			objstruct *obj);

