/*
*				poly.c
*
* Manage polynomials.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic WCS library
*
*	Copyright:		(C) 1998-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		20/12/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
#endif

#include	"poly.h"

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		  qerror("Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		  qerror("Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

/********************************* qerror ************************************/
/*
I hope it will never be used!
*/
void	qerror(char *msg1, char *msg2)
  {
  fprintf(stderr, "\n> %s%s\n\n",msg1,msg2);
  exit(-1);
  }


/****** poly_init ************************************************************
PROTO   polystruct *poly_init(int *group, int ndim, int *degree, int ngroup)
PURPOSE Allocate and initialize a polynom structure.
INPUT   1D array containing the group for each parameter,
        number of dimensions (parameters),
        1D array with the polynomial degree for each group,
        number of groups.
OUTPUT  polystruct pointer.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 30/08/2011
 ***/
polystruct	*poly_init(int *group, int ndim, int *degree, int ngroup)
  {
   void	qerror(char *msg1, char *msg2);
   polystruct	*poly;
   char		str[512];
   int		nd[POLY_MAXDIM];
   int		*groupt,
		d,g,n, num,den, dmax;

  QCALLOC(poly, polystruct, 1);
  if ((poly->ndim=ndim) > POLY_MAXDIM)
    {
    sprintf(str, "The dimensionality of the polynom (%d) exceeds the maximum\n"
		"allowed one (%d)", ndim, POLY_MAXDIM);
    qerror("*Error*: ", str);
    }

  if (ndim)
    QMALLOC(poly->group, int, poly->ndim);
    for (groupt=poly->group, d=ndim; d--;)
      *(groupt++) = *(group++)-1;

  poly->ngroup = ngroup;
  if (ngroup)
    {
    group = poly->group;	/* Forget the original *group */

    QMALLOC(poly->degree, int, poly->ngroup);

/*-- Compute the number of context parameters for each group */
    memset(nd, 0, ngroup*sizeof(int));
    for (d=0; d<ndim; d++)
      {
      if ((g=group[d])>=ngroup)
        qerror("*Error*: polynomial GROUP out of range", "");
      nd[g]++;
      }
    }

/* Compute the total number of coefficients */
  poly->ncoeff = 1;
  for (g=0; g<ngroup; g++)
    {
    if ((dmax=poly->degree[g]=*(degree++))>POLY_MAXDEGREE)
      {
      sprintf(str, "The degree of the polynom (%d) exceeds the maximum\n"
		"allowed one (%d)", poly->degree[g], POLY_MAXDEGREE);
      qerror("*Error*: ", str);
      }

/*-- There are (n+d)!/(n!d!) coeffs per group = Prod_(i<=d)(n+i)/Prod_(i<=d)i */
    n = nd[g];
    d = dmax>n? n: dmax;
    for (num=den=1; d; num*=(n+dmax--), den*=d--);
    poly->ncoeff *= num/den;
    }

  QMALLOC(poly->basis, double, poly->ncoeff);
  QCALLOC(poly->coeff, double, poly->ncoeff);

  return poly;
  }


/****** poly_end *************************************************************
PROTO   void poly_end(polystruct *poly)
PURPOSE Free a polynom structure and everything it contains.
INPUT   polystruct pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 09/04/2000
 ***/
void	poly_end(polystruct *poly)
  {
  if (poly)
    {
    free(poly->coeff);
    free(poly->basis);
    free(poly->degree);
    free(poly->group);
    free(poly);
    }
  }


/****** poly_func ************************************************************
PROTO   double poly_func(polystruct *poly, double *pos)
PURPOSE Evaluate a multidimensional polynom.
INPUT   polystruct pointer,
        pointer to the 1D array of input vector data.
OUTPUT  Polynom value.
NOTES   Values of the basis functions are updated in poly->basis.
AUTHOR  E. Bertin (IAP)
VERSION 03/03/2004
 ***/
double	poly_func(polystruct *poly, double *pos)
  {
   double	xpol[POLY_MAXDIM+1];
   double      	*post, *xpolt, *basis, *coeff, xval;
   long double	val;
   int		expo[POLY_MAXDIM+1], gexpo[POLY_MAXDIM+1];
   int	       	*expot, *degree,*degreet, *group,*groupt, *gexpot,
			d,g,t, ndim;

/* Prepare the vectors and counters */
  ndim = poly->ndim;
  basis = poly->basis;
  coeff = poly->coeff;
  group = poly->group;
  degree = poly->degree;
  if (ndim)
    {
    for (xpolt=xpol, expot=expo, post=pos, d=ndim; --d;)
      {
      *(++xpolt) = 1.0;
      *(++expot) = 0;
      }
    for (gexpot=gexpo, degreet=degree, g=poly->ngroup; g--;)
      *(gexpot++) = *(degreet++);
    if (gexpo[*group])
      gexpo[*group]--;
    }

/* The constant term is handled separately */
  val = *(coeff++);
  *(basis++) = 1.0;
  *expo = 1;
  *xpol = *pos;

/* Compute the rest of the polynom */
  for (t=poly->ncoeff; --t; )
    {
/*-- xpol[0] contains the current product of the x^n's */
    val += (*(basis++)=*xpol)**(coeff++);
/*-- A complex recursion between terms of the polynom speeds up computations */
/*-- Not too good for roundoff errors (prefer Horner's), but much easier for */
/*-- multivariate polynomials: this is why we use a long double accumulator */
    post = pos;
    groupt = group;
    expot = expo;
    xpolt = xpol;
    for (d=0; d<ndim; d++, groupt++)
      if (gexpo[*groupt]--)
        {
        ++*(expot++);
        xval = (*(xpolt--) *= *post);
        while (d--)
          *(xpolt--) = xval;
        break;
        }
      else
        {
        gexpo[*groupt] = *expot;
        *(expot++) = 0;
        *(xpolt++) = 1.0;
        post++;
        }
    }

  return (double)val;
  }


/****** poly_fit *************************************************************
PROTO   double poly_fit(polystruct *poly, double *x, double *y, double *w,
        int ndata, double *extbasis)
PURPOSE Least-Square fit of a multidimensional polynom to weighted data.
INPUT   polystruct pointer,
        pointer to the (pseudo)2D array of inputs to basis functions,
        pointer to the 1D array of data values,
        pointer to the 1D array of data weights,
        number of data points,
        pointer to a (pseudo)2D array of computed basis function values.
OUTPUT  Chi2 of the fit.
NOTES   If different from NULL, extbasis can be provided to store the
        values of the basis functions. If x==NULL and extbasis!=NULL, the
        precomputed basis functions stored in extbasis are used (which saves
        CPU). If w is NULL, all points are given identical weight.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 08/03/2005
 ***/
void	poly_fit(polystruct *poly, double *x, double *y, double *w, int ndata,
		double *extbasis)
  {
   void	qerror(char *msg1, char *msg2);
   double	/*offset[POLY_MAXDIM],*/x2[POLY_MAXDIM],
		*alpha,*alphat, *beta,*betat, *basis,*basis1,*basis2, *coeff,
		*extbasist,*xt,
		val,wval,yval;
   int		ncoeff, ndim, matsize,
		d,i,j,n;

  if (!x && !extbasis)
    qerror("*Internal Error*: One of x or extbasis should be "
	"different from NULL\nin ", "poly_func()");
  ncoeff = poly->ncoeff;
  ndim = poly->ndim;
  matsize = ncoeff*ncoeff;
  basis = poly->basis;
  extbasist = extbasis;
  QCALLOC(alpha, double, matsize);
  QCALLOC(beta, double, ncoeff);

/* Subtract an average offset to maintain precision (droped for now ) */
/*
  if (x)
    {
    for (d=0; d<ndim; d++)
      offset[d] = 0.0;
    xt = x;
    for (n=ndata; n--;)
      for (d=0; d<ndim; d++)
        offset[d] += *(xt++);
    for (d=0; d<ndim; d++)
      offset[d] /= (double)ndata;    
    }
*/ 
/* Build the covariance matrix */
  xt = x;
  for (n=ndata; n--;)
    {
    if (x)
      {
/*---- If x!=NULL, compute the basis functions */
      for (d=0; d<ndim; d++)
        x2[d] = *(xt++)/* - offset[d]*/;     
      poly_func(poly, x2);
/*---- If, in addition, extbasis is provided, then fill it */
      if (extbasis)
        for (basis1=basis,j=ncoeff; j--;)
          *(extbasist++) = *(basis1++);
      }
    else
/*---- If x==NULL, then rely on pre-computed basis functions */
      for (basis1=basis,j=ncoeff; j--;)
        *(basis1++) = *(extbasist++);

    basis1 = basis;
    wval = w? *(w++) : 1.0;
    yval = *(y++);
    betat = beta;
    alphat = alpha;
    for (j=ncoeff; j--;)
      {
      val = *(basis1++)*wval;
      *(betat++) += val*yval;
      for (basis2=basis,i=ncoeff; i--;)
        *(alphat++) += val**(basis2++);
      }
    }

/* Solve the system */
  poly_solve(alpha,beta,ncoeff);

  free(alpha);

/* Now fill the coeff array with the result of the fit */
  betat = beta;
  coeff = poly->coeff;
  for (j=ncoeff; j--;)
    *(coeff++) = *(betat++);
/*
  poly_addcste(poly, offset);
*/
  free(beta);

  return;
  }


/****** poly_addcste *********************************************************
PROTO   void poly_addcste(polystruct *poly, double *cste)
PURPOSE Modify matrix coefficients to mimick the effect of adding a cst to
	the input of a polynomial.
INPUT   Pointer to the polynomial structure,
        Pointer to the vector of cst.
OUTPUT  -.
NOTES   Requires quadruple-precision. **For the time beeing, this function
	returns completely wrong results!!**
AUTHOR  E. Bertin (IAP)
VERSION 05/10/2010
 ***/
void	poly_addcste(polystruct *poly, double *cste)
  {
   long double	*acoeff;
   double	*coeff,*mcoeff,*mcoefft,
		val;
   int		*mpowers,*powers,*powerst,*powerst2,
		i,j,n,p, denum, flag, maxdegree, ncoeff, ndim;

  ncoeff = poly->ncoeff;
  ndim = poly->ndim;
  maxdegree = 0;
  for (j=0; j<poly->ngroup; j++)
    if (maxdegree < poly->degree[j])
      maxdegree = poly->degree[j];
  maxdegree++;		/* Actually we need maxdegree+1 terms */
  QCALLOC(acoeff, long double, ncoeff);
  QCALLOC(mcoeff, double, ndim*maxdegree);
  QCALLOC(mpowers, int, ndim);
  mcoefft = mcoeff;		/* To avoid gcc -Wall warnings */
  powerst = powers = poly_powers(poly);
  coeff = poly->coeff;
  for (i=0; i<ncoeff; i++)
    {
    for (j=0; j<ndim; j++)
      {
      mpowers[j] = n = *(powerst++);
      mcoefft = mcoeff+j*maxdegree+n;
      denum = 1;
      val = 1.0;
      for (p=n+1; p--;)
        {
        *(mcoefft--) = val;
        val *= (cste[j]*(n--))/(denum++);	/* This is C_n^p X^(n-p) */
        }
      }
/*-- Update all valid coefficients */
    powerst2 = powers;
    for (p=0; p<ncoeff; p++)
      {
/*---- Check that this combination of powers is included in the series above */
      flag = 0;
      for (j=0; j<ndim; j++)
        if (mpowers[j] < powerst2[j])
	  {
          flag = 1;
          powerst2 += ndim;
          break;
          }
      if (flag == 1)
        continue;
      val = 1.0;
      mcoefft = mcoeff;
      for (j=ndim; j--; mcoefft += maxdegree)
        val *= mcoefft[*(powerst2++)];
      acoeff[i] += val*coeff[p];
      }
    }

/* Add the new coefficients to the previous ones */

  for (i=0; i<ncoeff; i++)
    coeff[i] = (double)acoeff[i];

  free(acoeff);
  free(mcoeff);
  free(mpowers);
  free(powers);

  return;
  }

/****** poly_solve ************************************************************
PROTO   void poly_solve(double *a, double *b, int n)
PURPOSE Solve a system of linear equations, using Cholesky decomposition.
INPUT   Pointer to the (pseudo 2D) matrix of coefficients,
        pointer to the 1D column vector,
        matrix size.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
VERSION 20/12/2011
 ***/
void	poly_solve(double *a, double *b, int n)
  {
#if defined(HAVE_LAPACKE)
  LAPACKE_dposv(LAPACK_COL_MAJOR, 'L', n, 1, a, n, b, n);
#elif defined(HAVE_ATLAS)
  clapack_dposv(CblasRowMajor, CblasUpper, n, 1, a, n, b, n);
#else
  cholsolve(a,b,n);
#endif

  return;
  }


/****** cholsolve *************************************************************
PROTO	int cholsolve(double *a, double *b, int n)
PURPOSE	Solve a system of linear equations, using Cholesky decomposition.
INPUT	Pointer to the (pseudo 2D) matrix of coefficients,
	pointer to the 1D column vector,
 	matrix size.
OUTPUT	-1 if the matrix is not positive-definite, 0 otherwise.
NOTES	Based on algorithm described in Numerical Recipes, 2nd ed. (Chap 2.9).
	The matrix of coefficients must be symmetric and positive definite.
AUTHOR	E. Bertin (IAP)
VERSION	10/10/2010
 ***/
int	cholsolve(double *a, double *b, int n)
  {
   double	*p, *x, sum;
   int		i,j,k;

/* Allocate memory to store the diagonal elements */
  QMALLOC(p, double, n);

/* Cholesky decomposition */
  for (i=0; i<n; i++)
    for (j=i; j<n; j++)
      {
      sum = a[i*n+j];
      for (k=i; k--;)
        sum -= a[i*n+k]*a[j*n+k];
      if (i==j)
        {
        if (sum <= 0.0)
	  {
          free(p);
          return -1;
          }
        p[i] = sqrt(sum);
        }
      else
        a[j*n+i] = sum/p[i];
      }

/* Solve the system */
  x = b;		/* Just to save memory:  the solution replaces b */
  for (i=0; i<n; i++)
    {
    for (sum=b[i],k=i; k--;)
      sum -= a[i*n+k]*x[k];
    x[i] = sum/p[i];
    }

  for (i=n; i--;)
    {
    sum = x[i];
    for (k=i; ++k<n;)
      sum -= a[k*n+i]*x[k];
    x[i] = sum/p[i];
    }

  free(p);

  return 0;
  }


/****** poly_powers ***********************************************************
PROTO   int *poly_powers(polystruct *poly)
PURPOSE	Return an array of powers of polynom terms
INPUT   polystruct pointer,
OUTPUT  Pointer to an array of polynom powers (int *), (ncoeff*ndim numbers).
NOTES   The returned pointer is mallocated.
AUTHOR  E. Bertin (IAP)
VERSION 23/10/2003
 ***/
int	*poly_powers(polystruct *poly)
  {
   int		expo[POLY_MAXDIM+1], gexpo[POLY_MAXDIM+1];
   int	       	*expot, *degree,*degreet, *group,*groupt, *gexpot,
		*powers, *powerst,
		d,g,t, ndim;

/* Prepare the vectors and counters */
  ndim = poly->ndim;
  group = poly->group;
  degree = poly->degree;
  QMALLOC(powers, int, ndim*poly->ncoeff);
  if (ndim)
    {
    for (expot=expo, d=ndim; --d;)
      *(++expot) = 0;
    for (gexpot=gexpo, degreet=degree, g=poly->ngroup; g--;)
      *(gexpot++) = *(degreet++);
    if (gexpo[*group])
      gexpo[*group]--;
    }

/* The constant term is handled separately */
  powerst = powers;
  for (d=0; d<ndim; d++)
    *(powerst++) = 0;
  *expo = 1;

/* Compute the rest of the polynom */
  for (t=poly->ncoeff; --t; )
    {
    for (d=0; d<ndim; d++)
      *(powerst++) = expo[d];
/*-- A complex recursion between terms of the polynom speeds up computations */
    groupt = group;
    expot = expo;
    for (d=0; d<ndim; d++, groupt++)
      if (gexpo[*groupt]--)
        {
        ++*(expot++);
        break;
        }
      else
        {
        gexpo[*groupt] = *expot;
        *(expot++) = 0;
        }
    }

  return powers;
  }

