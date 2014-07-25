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
*	Last modified:		30/06/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
#include <stdbool.h>

/*----------------------------- Internal constants --------------------------*/
#define	BACK_BUFSIZE		1048576		/* bkgnd buffer */
#define	BACK_MINGOODFRAC	0.5		/* min frac with good weights*/
#define	QUANTIF_NSIGMA		5		/* histogram limits */
#define	QUANTIF_NMAXLEVELS	4096		/* max nb of quantif. levels */
#define	QUANTIF_AMIN		4		/* min nb of "mode pixels" */
#define	BACK_EPS		(1e-4)		/* a small number */

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

  // TODO: put this struct in new file
 typedef struct objmask
  {
    // these are the basic components
    size_t width;   // nx of the mask
    size_t height;  // ny of the mask
    size_t npix;    // nx*ny of the mask
    bool* maskdata; // the central array

    // these components are associated with swapping/mapping
    int  swapflag;          // marks swapping
    size_t size;            // the memory footprint
    char swapname[MAXCHAR]; // the name of the mapped file
    FILE *swapfile;         // file pointer to the mapped file

    // helping variable to sweep through the data
    long int actindex;      // a pointer to the actual position in the data
    bool  *actpos;
  } objmaskstruct;

/*------------------------------- functions ---------------------------------*/
void		back_histo(backstruct *backmesh, backstruct *wbackmesh,
			PIXTYPE *buf, PIXTYPE *wbuf, size_t bufsize,
			int n, int w, int bw, PIXTYPE wthresh),
		back_stat(backstruct *backmesh, backstruct *wbackmesh,
			PIXTYPE *buf, PIXTYPE *wbuf, size_t bufsize,
			int n, int w, int bw, PIXTYPE wthresh),
		back_rmsline(fieldstruct *field, int y, PIXTYPE *line),
		back_copy(fieldstruct *infield, fieldstruct *outfield),
		back_end(fieldstruct *field),
		back_field(fieldstruct *field),
		back_filter(fieldstruct *field),
                back_map(fieldstruct *field, fieldstruct *wfield,
                        int wscale_flag),
                back_map_mask(fieldstruct *field, fieldstruct *wfield,
                        int wscale_flag, objmaskstruct *omask),
		back_subline(fieldstruct *field, int y, int xmin, int width,
			PIXTYPE *line),
		back_printmeshs(const backstruct *backmesh, const int nmeshs, int *tofile),
		back_printmesh(const backstruct *backmesh, FILE * outstream);

float		*back_makespline(fieldstruct *, float *),
		back_guess(backstruct *bkg, float *mean, float *sigma),
		back_local(fieldstruct *field, objstruct *obj, float *sigma);

PIXTYPE		back_interpolate(fieldstruct *field, double x, double y);

// TODO: put everything from here into new file
objmaskstruct *create_objmask(const size_t nx, const size_t ny,const int swapflag, const int index);
void free_objmask(objmaskstruct *omask);
void get_obmask_name(const int index, char *filename);
void set_objmask_value(const size_t xpos, const size_t ypos, const int value, objmaskstruct *omask);
int get_objmask_value(const size_t xpos, const size_t ypos, const objmaskstruct *omask);
void reset_objmask(objmaskstruct *omask);
long int tell_objmask(const objmaskstruct *omask);
void seek_objmask(const long int offset, const int origin, objmaskstruct *omask);
bool *read_objmaskOld(const size_t nobj, objmaskstruct *omask);
void read_objmask(const size_t nobj, bool **boolbuff, objmaskstruct *omask);
void apply_objmask(PIXTYPE *buffer, size_t npix, objmaskstruct *omask);
void print_objmask_elem(const size_t xpos, const size_t ypos, const objmaskstruct *omask);
void print_objmask(const objmaskstruct *omask);
void objmask_info(const objmaskstruct *omask);
void populate_objmask(fieldstruct *dfield, fieldstruct *dwfield, objmaskstruct *omask);
