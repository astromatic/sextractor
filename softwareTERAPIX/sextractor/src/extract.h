  /*
 				extract.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	functions for extraction of connected pixels from
*			a bitmap.
*
*	Last modify:	29/11/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*------------------------------ definitions --------------------------------*/

#define	NOBJ			256		/* starting number of obj. */
#define	UNKNOWN			-1		/* flag for LUTZ */

/*--------------------------------- typedefs --------------------------------*/

typedef	enum		{COMPLETE, INCOMPLETE, NONOBJECT, OBJECT}
				status;	/* Extraction status */

/*--------------------------------- variables -------------------------------*/
PIXTYPE		*dumscan;

/*------------------------------- structures --------------------------------*/
/* Temporary object parameters during extraction */
typedef struct structinfo
  {
  LONG		pixnb;			/* Number of pixels included */
  LONG		firstpix;		/* Pointer to first pixel of pixlist */
  LONG		lastpix;		/* Pointer to last pixel of pixlist */
  short		flag;			/* Extraction flag */
  }       infostruct;


/*------------------------------- functions ---------------------------------*/
void		lutzalloc(int, int),
		lutzfree(void),
		lutzsort(infostruct *, objliststruct *),
		sortit(picstruct *, picstruct *, picstruct *, picstruct *,
			infostruct *, objliststruct *, PIXTYPE *, PIXTYPE *),
		update(infostruct *, infostruct *, pliststruct *);

int		lutz(objliststruct *, int, objstruct *, objliststruct *); 
