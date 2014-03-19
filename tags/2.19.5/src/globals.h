/*
*				globals.h
*
* Global declarations and variables.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 1993-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		12/09/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
		write_vo_fields(FILE *file),
		zerocat(void);

extern double	counter_seconds(void);

extern float	fqmedian(float *, int);

extern int	addobj(int, objliststruct *, objliststruct *),
		belong(int, objliststruct *, int, objliststruct *),
		gatherup(objliststruct *, objliststruct *),
		parcelout(objliststruct *, objliststruct *);

extern void	*loadstrip(picstruct *, picstruct *);

extern char	*readfitshead(FILE *, char *, int *);

extern picstruct	*inheritfield(picstruct *infield, int flags),
			*newfield(char *, int , int);

