 /*
 				misc.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	miscellaneous functions.
*
*	Last modify:	20/08/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include	<time.h>
#include	<sys/time.h>

#include	"define.h"
#include	"globals.h"


/******************************** hmedian ***********************************/
/*
Median using Heapsort algorithm (for float arrays) (based on Num.Rec algo.).
Warning: changes the order of data!
*/
float	hmedian(float *ra, int n)

  {
   int		l, j, ir, i;
   float	rra;


  if (n<2)
    return *ra;
  ra--;
  for (l = ((ir=n)>>1)+1;;)
    {
    if (l>1)
      rra = ra[--l];
    else
      {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1)
        {
        ra[1] = rra;
        return n&1? ra[n/2+1] : (float)((ra[n/2]+ra[n/2+1])/2.0);
        }
      }
    for (j = (i=l)<<1; j <= ir;)
      {
      if (j < ir && ra[j] < ra[j+1])
        ++j;
      if (rra < ra[j])
        {
        ra[i] = ra[j];
        j += (i=j);
        }
      else
        j = ir + 1;
      }
    ra[i] = rra;
    }

/* (the 'return' is inside the loop!!) */
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


