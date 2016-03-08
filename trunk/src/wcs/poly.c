/*
*				poly.c
*
* Polynomial functions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1998-2012 IAP/CNRS/UPMC
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
*	Last modified:		20/11/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"poly.h"

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
//#define MATSTORAGE_PACKED 1
#endif

#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		  qerror("Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		  qerror("Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
		  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		    qerror("Not enough memory for ", \
			#ptrout " (" #nel " elements) !"); \
		memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ));};;}

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
AUTHOR  E. Bertin (IAP)
VERSION 04/11/2008
 ***/
void	poly_end(polystruct *poly)
  {
  if (poly)
    {
    free(poly->coeff);
    free(poly->basis);
    free(poly->orthobasis);
    free(poly->degree);
    free(poly->group);
    free(poly->orthomat);
    free(poly->deorthomat);
    free(poly);
    }

  return;
  }


/****** poly_copy *************************************************************
PROTO   polystruct *poly_copy(polystruct *poly)
PURPOSE Copy a polynom structure and everything it contains.
INPUT   polystruct pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 04/11/2008
 ***/
polystruct *poly_copy(polystruct *poly)
  {
   polystruct	*newpoly;

  if (poly)
    {
    QMALLOC(newpoly, polystruct, 1);
    *newpoly = *poly;
    if (poly->ncoeff)
      {
      QMEMCPY(poly->coeff, newpoly->coeff, double, poly->ncoeff);
      QMEMCPY(poly->basis, newpoly->basis, double, poly->ncoeff);
      }
    if (poly->ndim)
      QMEMCPY(poly->group, newpoly->group, int, poly->ndim);
    if (poly->ngroup)
      QMEMCPY(poly->degree, newpoly->degree, int, poly->ngroup);
    if (poly->orthomat)
      {
      QMEMCPY(poly->orthomat, newpoly->orthomat, double,
		poly->ncoeff*poly->ncoeff);
      QMEMCPY(poly->deorthomat, newpoly->deorthomat, double,
		poly->ncoeff*poly->ncoeff);
      QMEMCPY(poly->orthobasis, newpoly->orthobasis, double, poly->ncoeff);
      }

    return newpoly;
    }
  else
    return NULL;
  }


/****** poly_func ************************************************************
PROTO   double poly_func(polystruct *poly, double *pos)
PURPOSE Evaluate a multidimensional polynomial.
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


/****** poly_cfunc ************************************************************
PROTO   double poly_cfunc(polystruct *poly, double *pos)
PURPOSE Evaluate a multidimensional Chebyshev polynomial.
INPUT   polystruct pointer,
        pointer to the 1D array of input vector data.
OUTPUT  Polynom value.
NOTES   Values of the basis functions are updated in poly->basis.
AUTHOR  E. Bertin (IAP)
VERSION 29/01/2013
 ***/
double	poly_cfunc(polystruct *poly, double *pos)
  {
   double	pol[POLY_MAXDIM*(POLY_MAXDEGREE+1)],
	      	*polt, *post, *basis, *coeff, xval;
   long double	val;
   int		expo[POLY_MAXDIM+1], gexpo[POLY_MAXDIM+1];
   int	       	*expot, *degree,*degreet, *group,*groupt, *gexpot,
			d,d2,g,t, ndim;

/* Prepare the vectors and counters */
  ndim = poly->ndim;
  basis = poly->basis;
  coeff = poly->coeff;
  group = poly->group;
  degree = poly->degree;
  if (ndim)
    {
    for (groupt=group, expot=expo, post=pos, d=0; d<ndim; d++)
      {
      *(expot++) = 0;
      polt = pol + d*(POLY_MAXDEGREE+1);
      *(polt++) = 1.0;
      *(polt++) = xval = *(post++);
      for (d2 = degree[*(groupt++)]; --d2 > 0; polt++)
        *polt = 2.0*xval**(polt-1) - *(polt-2);
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

/* Compute the rest of the polynom */
  for (t=poly->ncoeff; --t; )
    {
    polt = pol;
    expot = expo;
/*-- xval contains the current product of the polynomials */
    xval = 1.0;
    for (d=ndim; d--; polt += POLY_MAXDEGREE+1)
      xval *= polt[*(expot++)];
    val += (*(basis++)=xval)**(coeff++);
/*-- A complex recursion between terms of the polynom speeds up computations */
/*-- Not too good for roundoff errors (prefer Horner's), but much easier for */
/*-- multivariate polynomials: this is why we use a long double accumulator */
    expot = expo;
    groupt = group;
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

  return (double)val;
  }


/****** poly_fit *************************************************************
PROTO   int poly_fit(polystruct *poly, double *x, double *y, double *w,
        int ndata, double *extbasis, double regul)
PURPOSE Least-Square fit of a multidimensional polynom to weighted data.
INPUT   polystruct pointer,
        pointer to the (pseudo)2D array of inputs to basis functions,
        pointer to the 1D array of data values,
        pointer to the 1D array of data weights,
        number of data points,
        pointer to a (pseudo)2D array of computed basis function values.
	Tikhonov regularization parameter (0 = no regularization).
OUTPUT  Chi2 of the fit.
NOTES   If different from NULL, extbasis can be provided to store the
        values of the basis functions. If x==NULL and extbasis!=NULL, the
        precomputed basis functions stored in extbasis are used (which saves
        CPU). If w is NULL, all points are given identical weight.
AUTHOR  E. Bertin (IAP)
VERSION 20/11/2012
 ***/
int	poly_fit(polystruct *poly, double *x, double *y, double *w, int ndata,
		double *extbasis, double regul)
  {
   void	qerror(char *msg1, char *msg2);
   double	/*offset[POLY_MAXDIM],*/x2[POLY_MAXDIM],
		*alpha,*alphat, *beta,*betat, *basis,*basis1,*basis2, *coeff,
		*extbasist,*xt,
		val,wval,yval;
   int		ncoeff, ndim, matsize, info,
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

  if (regul>POLY_TINY)
/*-- Simple Tikhonov regularization */
    for (i=0; i<ncoeff; i++)
      alpha[i*(ncoeff+1)] += regul;

/* Solve the system */
  info = poly_solve(alpha,beta,ncoeff);

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

  return info;
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
PROTO   int poly_solve(double *a, double *b, int n)
PURPOSE Solve a system of linear equations, using Cholesky decomposition.
INPUT   Pointer to the (pseudo 2D) matrix of coefficients,
        pointer to the 1D column vector,
        matrix size.
OUTPUT  0 if solution OK, !=0 otherwise (e.g., singular matrix).
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 20/11/2012
 ***/
int	poly_solve(double *a, double *b, int n)
  {
#if defined(HAVE_LAPACKE)
  return LAPACKE_dposv(LAPACK_COL_MAJOR, 'L', n, 1, a, n, b, n);
#elif defined(HAVE_ATLAS)
  return clapack_dposv(CblasRowMajor, CblasUpper, n, 1, a, n, b, n);
#else
  return cholsolve(a,b,n);
#endif
  }


/****** cholsolve *************************************************************
PROTO	void cholsolve(double *a, double *b, int n)
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


/****** poly_initortho ********************************************************
PROTO   void poly_initortho(polystruct *poly, double *data, int ndata)
PURPOSE Compute orthonormalization and de-orthonormalization matrices for a
	polynomial basis on a data set.
INPUT   polystruct pointer,
        pointer to the 1D array of input vector data,
	number of data vectors.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 10/07/2012
 ***/
void	poly_initortho(polystruct *poly, double *data, int ndata)
  {
   double	*basis, *coeff, *invec,*invect0,*invect,*invect02,*invect2,
		*rdiag, *deortho,
		scale,s, dval;
   int		c,i,j,m,n, ndmc, ndim,ncoeff;

/* Prepare the vectors and counters */
  ndim = poly->ndim;
  ncoeff = poly->ncoeff;
  basis = poly->basis;
  coeff = poly->coeff;

/* Allocate memory for orthonormalization matrix and vector */
  QCALLOC(poly->deorthomat, double, ncoeff*ncoeff);
  QMALLOC(poly->orthobasis, double, poly->ncoeff);
  QMALLOC(rdiag, double, ncoeff);

/* Do a QR decomposition of input vector set */
/* Vectors are stored as rows to speed up the Householder transformation */
  n = ncoeff;
  m = ndata;
  invec = data;
  for (c=0; c<ncoeff; c++)
    {
    ndmc = ndata - c;
    scale = 0.0;
    invect = invect0 = data + c*(ndata+1);
    for (i=ndmc; i--; invect++)
      scale = sqrt(scale*scale + *invect**invect);
    if (scale > POLY_TINY)
      {
      if (*invect0 < 0.0)
        scale = -scale;
      invect = invect0;
      for (i=ndmc; i--;)
        *(invect++) /= scale;
      *invect0 += 1.0;
      invect02 = invect0 + ndata;
      for (j=ncoeff-c; --j; invect02+=ndata)
        {
        s = 0.0;
        invect = invect0;
        invect2 = invect02;
        for (i=ndmc; i--;)
          s += *(invect++)**(invect2++);
        s /= -*invect0;
        invect = invect0;
        invect2 = invect02;
        for (i=ndmc; i--;)
          *(invect2++) += s**(invect++);
        }
      }
    rdiag[c] = -scale;
    }

/* Convert to deorthonormalization matrix */
  deortho = poly->deorthomat;
  for (j=0; j<ncoeff; j++)
    for (i=0; i<ncoeff; i++)
      deortho[j*ncoeff+i] = i<j? data[j*ndata+i] : (i==j?rdiag[i] : 0.0);

  free(rdiag);

/* Compute the "unorthonormalization" matrix */
  QMEMCPY(poly->deorthomat, poly->orthomat, double, ncoeff*ncoeff);
#if defined(HAVE_LAPACKE)
  LAPACKE_dtrtri(LAPACK_ROW_MAJOR, 'L', 'N', ncoeff,poly->orthomat,ncoeff);
#elif defined(HAVE_ATLAS)
  clapack_dtrtri(CblasRowMajor, CblasLower, CblasNonUnit, ncoeff,
	poly->orthomat, ncoeff);
#else
  qerror("*Internal Error*: no routine available", " for triangular inverse");
#endif

/* Transpose orthonormalization matrix to speed up later use */
  deortho = poly->deorthomat;
  for (j=0; j<ncoeff; j++)
    for (i=j; i<ncoeff; i++)
      {
      dval = deortho[j*ncoeff+i];
      deortho[j*ncoeff+i] = deortho[i*ncoeff+j];
      deortho[i*ncoeff+j] = dval;
      }

  return;
  }


/****** poly_ortho ************************************************************
PROTO   double *poly_ortho(polystruct *poly, double *datain, double *dataout)
PURPOSE Apply orthonormalization to the poly basis vector ("ket>").
INPUT   polystruct pointer,
	pointer to the input vector,
	pointer to the output vector.
OUTPUT  Pointer to poly->orthobasis, or poly->basis if no ortho. matrix exists.
NOTES   The poly->basis vector must have been updated with poly_func() first.
AUTHOR  E. Bertin (IAP)
VERSION 04/11/2008
 ***/
double	*poly_ortho(polystruct *poly, double *datain, double *dataout)
  {
   double	*omat,*basis,*obasis,
		dval;
   int		i,j, ncoeff;

  if (!poly->orthomat)
    return datain;

  ncoeff = poly->ncoeff;

/* Compute matrix product */
  omat = poly->orthomat;
  obasis = dataout;
  for (j=ncoeff; j--;)
    {
    basis = datain;
    dval = 0.0;
    for (i=ncoeff; i--;)
      dval += *(omat++)**(basis++);
    *(obasis++) = dval;
    }

  return dataout;
  }


/****** poly_deortho **********************************************************
PROTO   void poly_deortho(polystruct *poly, double *datain, double *dataout)
PURPOSE Apply deorthonormalization to the poly basis component vector("<bra|").
INPUT   polystruct pointer,
	pointer to the input vector,
	pointer to the output vector.
OUTPUT  Pointer to poly->basis.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 04/11/2008
 ***/
double	*poly_deortho(polystruct *poly, double *datain, double *dataout)
  {
   double	*omat,*basis,*obasis,
		dval;
   int		i,j, ncoeff;

  if (!poly->deorthomat)
    return datain;

  ncoeff = poly->ncoeff;

/* Compute matrix product */
  omat = poly->deorthomat;
  basis = dataout;
  for (j=ncoeff; j--;)
    {
    obasis = datain;
    dval = 0.0;
    for (i=ncoeff; i--;)
      dval += *(omat++)**(obasis++);
    *(basis++) = dval;
    }

  return dataout;
  }


