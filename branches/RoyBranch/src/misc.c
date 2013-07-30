/*
*				misc.c
*
* Miscellaneous functions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SExtractor
*
*	Copyright:		(C) 2009-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		11/03/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<stdlib.h>
#include	<time.h>
#include	<sys/time.h>

#include	"define.h"
#include	"globals.h"


/*i**** fqcmp **************************************************************
PROTO	int	fqcmp(const void *p1, const void *p2)
PURPOSE	Sorting function for floats in qsort().
INPUT	Pointer to first element,
	pointer to second element.
OUTPUT	1 if *p1>*p2, 0 if *p1=*p2, and -1 otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	05/10/2010
 ***/
static int	fqcmp(const void *p1, const void *p2)
  {
   double	f1=*((float *)p1),
		f2=*((float *)p2);
  return f1>f2? 1 : (f1<f2? -1 : 0);
  }


/****** fqmedian **************************************************************
PROTO	float   fqmedian(float *ra, int n)
PURPOSE	Compute the median of an array of floats, using qsort().
INPUT	Pointer to the array,
	Number of array elements.
OUTPUT	Median of the array.
NOTES	Warning: the order of input data is modified!.
AUTHOR	E. Bertin (IAP)
VERSION	05/10/2010
 ***/
float	fqmedian(float *ra, int n)

  {
   int dqcmp(const void *p1, const void *p2);

  qsort(ra, n, sizeof(float), fqcmp);
  if (n<2)
    return *ra;
  else
    return n&1? ra[n/2] : (ra[n/2-1]+ra[n/2])/2.0;
  }


/****** propagate_covar ******************************************************
PROTO	void	propagate_covar(double *vi, double *d, double *vo,
				int ni, int no,	double *temp)
PURPOSE	Compute Dt.V.D (propagate covariance matrix errors)
INPUT	Pointer to the original covariance matrix,
	pointer to the matrix of derivatives,
	input number of parameters,
	output number of parameters,
	pointer to a ni*no work array.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/08/2010
 ***/
void propagate_covar(double *vi, double *d, double *vo,
				int ni, int no,	double *temp)
  {
   double	*vit,*dt,*vot,*tempt,
		dval;
   int		i,j,k;

  tempt = temp;
  vit = vi;
  for (j=0; j<ni; j++)
    {
    dt = d;
    for (i=no; i--;)
      {
      vit = vi + j*ni;
      dval = 0.0;
      for (k=ni; k--;)
        dval += *(vit++)**(dt++);
      *(tempt++) = dval;
      }
    }

  vot = vo;
  for (j=0; j<no; j++)
    {
    for (i=0; i<no; i++)
      {
      dt = d + j*ni;
      tempt = temp + i;
      dval = 0.0;
      for (k=ni; k--; tempt+=no)
        dval += *(dt++)**tempt;
      *(vot++) = dval;
      }
    }

  return;
  }


/****** counter_seconds *******************************************************
PROTO	double counter_seconds(void)
PURPOSE	Count the number of seconds (with an arbitrary offset).
INPUT	-.
OUTPUT	Returns a number of seconds.
NOTES	Results are meaningful only for tasks that take one microsec or more.
AUTHOR	E. Bertin (IAP)
VERSION	24/09/2009
 ***/
double	counter_seconds(void)
  {
   struct timeval	tp;
   struct timezone	tzp;
   int			dummy;

  dummy = gettimeofday(&tp,&tzp);
  return (double) tp.tv_sec + (double) tp.tv_usec * 1.0e-6;
  }


