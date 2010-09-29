 /*
 				globals.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	global declarations.
*
*	Last modify:	20/08/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include	"types.h"

/*----------------------- miscellaneous variables ---------------------------*/

sexcatstruct		thecat;
picstruct		thefield1,thefield2, thewfield1,thewfield2;
objstruct		flagobj;
obj2struct		flagobj2;
extern obj2struct	outobj2;
float			ctg[37], stg[37];
char			gstr[MAXCHAR];

/*------------------------------- functions ---------------------------------*/
extern void	alloccatparams(void),
		allocparcelout(void),
		analyse(picstruct *, picstruct *, int, objliststruct *),
		blankit(char *, int),
                endcat(char *error),
                reendcat(void),
		changecatparamarrays(char *keyword, int *axisn, int naxis),
                closecheck(void),
		copydata(picstruct *, int, int),
		dumpparams(void),
		endfield(picstruct *),
		endobject(picstruct *, picstruct *, picstruct *, picstruct *,
			int, objliststruct *),
		examineiso(picstruct *, picstruct *, objstruct *,
			pliststruct *),
		flagcleancrowded(int, objliststruct *),
		freeparcelout(void),
		getnnw(void),
		initcat(void),
		reinitcat(picstruct *),
		initglob(void),
		makeit(void),
		mergeobject(objstruct *, objstruct *),
		neurinit(void),
		neurclose(void),
		neurresp(double *, double *),
		preanalyse(int, objliststruct *, int),
		propagate_covar(double *vi, double *d, double *vo,
				int ni, int no,	double *temp),
		readcatparams(char *),
		readdata(picstruct *, PIXTYPE *, int),
		readidata(picstruct *, FLAGTYPE *, int),
		readimagehead(picstruct *),
		readprefs(char *, char **, char **, int),
		scanimage(picstruct *, picstruct *, picstruct **, int,
			picstruct *, picstruct *),
		sexcircle(PIXTYPE *bmp, int, int, double, double, double,
			PIXTYPE),
		sexdraw(PIXTYPE *bmp, int, int, double, double, PIXTYPE),
		sexellips(PIXTYPE *bmp, int, int, double, double, double,
			double, double, PIXTYPE, int),
		sexmove(double, double),
		updateparamflags(void),
		useprefs(void),
		writecat(int, objliststruct *),
		write_error(char *msg1, char *msg2),
		write_vo_fields(FILE *file);

extern double	counter_seconds(void);

extern float	hmedian(float *, int);

extern int	addobj(int, objliststruct *, objliststruct *),
		belong(int, objliststruct *, int, objliststruct *),
		gatherup(objliststruct *, objliststruct *),
		parcelout(objliststruct *, objliststruct *);

extern void	*loadstrip(picstruct *, picstruct *);

extern char	*readfitshead(FILE *, char *, int *);

extern picstruct	*inheritfield(picstruct *infield, int flags),
			*newfield(char *, int , int);

