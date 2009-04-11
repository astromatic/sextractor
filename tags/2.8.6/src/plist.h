 /*
 				plist.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN, (IAP)
*
*	Contents:	functions dealing with the handling of pixel lists.
*
*	Last modify:	29/11/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*------------------------------- definitions -------------------------------*/

#define	PLIST(ptr, elem)	(((pbliststruct *)(ptr))->elem)

#define	PLISTEXIST(elem)	(plistexist_##elem)

#define	PLISTPIX(ptr, elem)	(*((PIXTYPE *)((ptr)+plistoff_##elem)))

#define	PLISTFLAG(ptr, elem)	(*((FLAGTYPE *)((ptr)+plistoff_##elem)))

/*------------------------------- structures --------------------------------*/

typedef struct
  {
  int		nextpix;
  int		x, y;
  PIXTYPE       value;
  }	pbliststruct;

/*-------------------------------- globals ----------------------------------*/

int	plistexist_value, plistexist_dvalue, plistexist_cdvalue,
	plistexist_flag, plistexist_wflag, plistexist_dthresh, plistexist_var,
	plistoff_value, plistoff_dvalue, plistoff_cdvalue,
	plistoff_flag[MAXFLAG], plistoff_wflag, plistoff_dthresh, plistoff_var,
	plistsize;

/*------------------------------- functions ---------------------------------*/

void	init_plist(void);

int	createblank(objliststruct *objlist, int n),
	createsubmap(objliststruct *objlist, int n);
