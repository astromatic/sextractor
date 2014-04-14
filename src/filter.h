/*
*				filter.h
*
* Include file for filter.c.
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
*	Last modified:		11/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*------------------------------- definitions -------------------------------*/

#define	MAXMASK		1024	/* Maximum number of mask elements (=32x32) */
#define POLY_MAXDIM        4    /* Max dimensionality of polynom (from poly.h)*/

/*------------------------------- structures --------------------------------*/

/** @struct varconv
 *  @brief Structure for a 2D variable convolution
 */
typedef struct convvar
  {
    double    pos[POLY_MAXDIM]; /**< member pos[] the independent coords for psf in the right order */
    double    *xpos;            /**< member xpos  points to the correct location of x in pos[] */
    double    *ypos;            /**< member ypos  points to the correct location of y in pos[] */
    psfstruct *psf;             /**< member psf   the psf structure */
    double    *basis;           /**< member basis pointer to store the polynomial coefficients for an image row */
  } varconv;

  /** @struct filterstruct
   *  @brief Structure for the convolution filter
   */
typedef struct structfilter
  {
  float	  *conv;		 /**< member conv pointer to the convolution mask */
  int	  nconv;		 /**< member nconv total number of elements */
  int     convw;	         /**< member convw width of the mask */
  int     convh;                 /**< member convw height of the mask */
  float   varnorm;
  varconv *varpsf;               /**< member varpsf for 2D variable convolution */
  struct  structbpann	*bpann;  /**< member structbpann  neural filtering */
  void    (*filterFunc)(fieldstruct *, PIXTYPE *, int y); /**< member filterFunc */
  }	filterstruct;
filterstruct *thefilter;

/*------------------------------- functions ---------------------------------*/
void  convolve(fieldstruct *, PIXTYPE *, int y);
void  convolve_var(fieldstruct *, PIXTYPE *, int y);
void  convolve_image(fieldstruct *field, float *vig1, float *vig2, int width, int height);
void  filter(fieldstruct *, PIXTYPE *, int y);
void  neurfilter(fieldstruct *, PIXTYPE *, int y);
void  endfilter(void);


extern int  getneurfilter(const char *filename);
extern void getfilter(const char *filename);
extern int  getconv(const char *filename);
extern int  getASCIIconv(const char *filename);
extern int  getConstPSFExConv(const char *filename);
extern int  getPSFExConv(const char *filename);
extern int  optPSFExSize(const psfstruct *psf, int *size);
extern void optPSFExSizeOld(const psfstruct *psf);
extern void getImageStats(const float *pix, const int width, const int height, float stats[]);
