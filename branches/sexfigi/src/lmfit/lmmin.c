/*
 * lmfit
 *
 * Solves or minimizes the sum of squares of m nonlinear
 * functions of n variables.
 *
 * From public domain Fortran version
 * of Argonne National Laboratories MINPACK
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more
 * C translation by Steve Moshier
 * Joachim Wuttke converted the source into C++ compatible ANSI style
 * and provided a simplified interface
 */

 
#include <stdlib.h>
#include <math.h>
#include "lmmin.h"
#define _LMDIF

/* *********************** high-level interface **************************** */


void lm_initialize_control( lm_control_type *control )
{
    control->maxcall = 100;
    control->epsilon = 1.e-14;
    control->stepbound = 100.;
    control->ftol = 1.e-14;
    control->xtol = 1.e-14;
    control->gtol = 1.e-14;
}

void lm_minimize( int m_dat, int n_par, double* par,
                  lm_evaluate_ftype *evaluate, lm_print_ftype *printout,
                  void *data, lm_control_type *control )
{

// *** allocate work space.

    double *fvec, *diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
    int *ipvt;

    int n = n_par;
    int m = m_dat;

    if (!(fvec = (double*) malloc(  m*sizeof(double))) ||
        !(diag = (double*) malloc(n*  sizeof(double))) ||
        !(qtf =  (double*) malloc(n*  sizeof(double))) ||
        !(fjac = (double*) malloc(n*m*sizeof(double))) ||
        !(wa1 =  (double*) malloc(n*  sizeof(double))) ||
        !(wa2 =  (double*) malloc(n*  sizeof(double))) ||
        !(wa3 =  (double*) malloc(n*  sizeof(double))) ||
        !(wa4 =  (double*) malloc(  m*sizeof(double))) ||
        !(ipvt = (int*)    malloc(n*  sizeof(int)))) {
            control->info = 9;
            return;
    }

// *** perform fit.
			
    control->info = 0;
    control->nfev = 0;

    // this goes through the modified legacy interface:
    lm_lmdif( m, n, par, fvec, control->ftol, control->xtol, control->gtol,
              control->maxcall*(n+1), control->epsilon, diag, 1,
              control->stepbound, &(control->info),
              &(control->nfev), fjac, ipvt, qtf, wa1, wa2, wa3, wa4,
              evaluate, printout, data );

    (*printout)( n, par, m, fvec, data, -1, 0, control->nfev );
    control->fnorm = lm_enorm(m, fvec);
    if (control->info < 0 ) control->info = 10;

// *** clean up.

    free(fvec);
    free(diag);
    free(qtf); 
    free(fjac);
    free(wa1); 
    free(wa2); 
    free(wa3 );
    free(wa4); 
    free(ipvt);
}


// ***** the following messages are referenced by the variable info.

char *lm_infmsg[] = {
    "improper input parameters",
    "the relative error in the sum of squares is at most tol",
    "the relative error between x and the solution is at most tol",
    "both errors are at most tol",
    "fvec is orthogonal to the columns of the jacobian to machine precision",
    "number of calls to fcn has reached or exceeded 200*(n+1)",
    "ftol is too small. no further reduction in the sum of squares is possible",
    "xtol too small. no further improvement in approximate solution x possible",
    "gtol too small. no further improvement in approximate solution x possible",
    "not enough memory",
    "break requested within function evaluation"
};
 
char *lm_shortmsg[] = {
        "invalid input",
        "success (f)",
        "success (p)",
        "success (f,p)",
        "degenerate",
        "call limit",
        "failed (f)",
        "failed (p)",
        "failed (o)",
        "no memory",
        "user break"
};


/* ************************** implementation ******************************* */


#define BUG 0
#if BUG 
#include <stdio.h>
#endif

// the following values seem good for an x86:
#define LM_MACHEP .555e-16 /* resolution of arithmetic */
#define LM_DWARF  9.9e-324 /* smallest nonzero number */
// the follwoing values should work on any machine:
// #define LM_MACHEP 1.2e-16
// #define LM_DWARF 1.0e-38

// the squares of the following constants shall not under/overflow:
// these values seem good for an x86:
#define LM_SQRT_DWARF 1.e-160
#define LM_SQRT_GIANT 1.e150
// the following values should work on any machine:
// #define LM_SQRT_DWARF 3.834e-20
// #define LM_SQRT_GIANT 1.304e19


void lm_qrfac( int m, int n, double* a, int pivot, int* ipvt,
               double* rdiag, double* acnorm, double* wa);
void lm_qrsolv(int n, double* r, int ldr, int* ipvt, double* diag,
               double* qtb, double* x, double* sdiag, double* wa);
void lm_lmpar( int n, double* r, int ldr, int* ipvt, double* diag, double* qtb,
               double delta, double* par, double* x, double* sdiag,
               double* wa1, double* wa2);

#define MIN(a,b) (((a)<=(b)) ? (a) : (b))
#define MAX(a,b) (((a)>=(b)) ? (a) : (b))
#define SQR(x)   (x)*(x) 


// ***** the low-level legacy interface for full control.

void lm_lmdif( int m, int n, double* x, double* fvec, double ftol, double xtol,
               double gtol, int maxfev, double epsfcn, double* diag, int mode,
               double factor, int *info, int *nfev, 
               double* fjac, int* ipvt, double* qtf,
               double* wa1, double* wa2, double* wa3, double* wa4,
               lm_evaluate_ftype *evaluate, lm_print_ftype *printout,
               void *data )
{
/*
 *   the purpose of lmdif is to minimize the sum of the squares of
 *   m nonlinear functions in n variables by a modification of
 *   the levenberg-marquardt algorithm. the user must provide a
 *   subroutine evaluate which calculates the functions. the jacobian
 *   is then calculated by a forward-difference approximation.
 *
 *   the multi-parameter interface lm_lmdif is for users who want
 *   full control and flexibility. most users will be better off using
 *   the simpler interface lmfit provided above.
 *
 *   the parameters are the same as in the legacy FORTRAN implementation,
 *   with the following exceptions:
 *      the old parameter ldfjac which gave leading dimension of fjac has
 *        been deleted because this C translation makes no use of two-
 *        dimensional arrays;
 *      the old parameter nprint has been deleted; printout is now controlled
 *        by the user-supplied routine *printout;
 *      the parameter field *data and the function parameters *evaluate and
 *        *printout have been added; they help avoiding global variables.
 *
 *   parameters:
 *
 *	m is a positive integer input variable set to the number
 *	  of functions.
 *
 *	n is a positive integer input variable set to the number
 *	  of variables. n must not exceed m.
 *
 *	x is an array of length n. on input x must contain
 *	  an initial estimate of the solution vector. on output x
 *	  contains the final estimate of the solution vector.
 *
 *	fvec is an output array of length m which contains
 *	  the functions evaluated at the output x.
 *
 *	ftol is a nonnegative input variable. termination
 *	  occurs when both the actual and predicted relative
 *	  reductions in the sum of squares are at most ftol.
 *	  therefore, ftol measures the relative error desired
 *	  in the sum of squares.
 *
 *	xtol is a nonnegative input variable. termination
 *	  occurs when the relative error between two consecutive
 *	  iterates is at most xtol. therefore, xtol measures the
 *	  relative error desired in the approximate solution.
 *
 *	gtol is a nonnegative input variable. termination
 *	  occurs when the cosine of the angle between fvec and
 *	  any column of the jacobian is at most gtol in absolute
 *	  value. therefore, gtol measures the orthogonality
 *	  desired between the function vector and the columns
 *	  of the jacobian.
 *
 *	maxfev is a positive integer input variable. termination
 *	  occurs when the number of calls to lm_fcn is at least
 *	  maxfev by the end of an iteration.
 *
 *	epsfcn is an input variable used in determining a suitable
 *	  step length for the forward-difference approximation. this
 *	  approximation assumes that the relative errors in the
 *	  functions are of the order of epsfcn. if epsfcn is less
 *	  than the machine precision, it is assumed that the relative
 *	  errors in the functions are of the order of the machine
 *	  precision.
 *
 *	diag is an array of length n. if mode = 1 (see below), diag is
 *        internally set. if mode = 2, diag must contain positive entries
 *        that serve as multiplicative scale factors for the variables.
 *
 *	mode is an integer input variable. if mode = 1, the
 *	  variables will be scaled internally. if mode = 2,
 *	  the scaling is specified by the input diag. other
 *	  values of mode are equivalent to mode = 1.
 *
 *	factor is a positive input variable used in determining the
 *	  initial step bound. this bound is set to the product of
 *	  factor and the euclidean norm of diag*x if nonzero, or else
 *	  to factor itself. in most cases factor should lie in the
 *	  interval (.1,100.). 100. is a generally recommended value.
 *
 *	info is an integer output variable that indicates the termination
 *        status of lm_lmdif as follows:
 *
 *        info < 0  termination requested by user-supplied routine *evaluate;
 *
 *	  info = 0  improper input parameters;
 *
 *	  info = 1  both actual and predicted relative reductions
 *		    in the sum of squares are at most ftol;
 *
 *	  info = 2  relative error between two consecutive iterates
 *		    is at most xtol;
 *
 *	  info = 3  conditions for info = 1 and info = 2 both hold;
 *
 *	  info = 4  the cosine of the angle between fvec and any
 *		    column of the jacobian is at most gtol in
 *		    absolute value;
 *
 *	  info = 5  number of calls to lm_fcn has reached or
 *		    exceeded maxfev;
 *
 *	  info = 6  ftol is too small. no further reduction in
 *		    the sum of squares is possible;
 *
 *	  info = 7  xtol is too small. no further improvement in
 *		    the approximate solution x is possible;
 *
 *	  info = 8  gtol is too small. fvec is orthogonal to the
 *		    columns of the jacobian to machine precision;
 *
 *	nfev is an output variable set to the number of calls to the
 *        user-supplied routine *evaluate.
 *
 *	fjac is an output m by n array. the upper n by n submatrix
 *	  of fjac contains an upper triangular matrix r with
 *	  diagonal elements of nonincreasing magnitude such that
 *
 *		 t     t	   t
 *		p *(jac *jac)*p = r *r,
 *
 *	  where p is a permutation matrix and jac is the final
 *	  calculated jacobian. column j of p is column ipvt(j)
 *	  (see below) of the identity matrix. the lower trapezoidal
 *	  part of fjac contains information generated during
 *	  the computation of r.
 *
 *	ipvt is an integer output array of length n. ipvt
 *	  defines a permutation matrix p such that jac*p = q*r,
 *	  where jac is the final calculated jacobian, q is
 *	  orthogonal (not stored), and r is upper triangular
 *	  with diagonal elements of nonincreasing magnitude.
 *	  column j of p is column ipvt(j) of the identity matrix.
 *
 *	qtf is an output array of length n which contains
 *	  the first n elements of the vector (q transpose)*fvec.
 *
 *	wa1, wa2, and wa3 are work arrays of length n.
 *
 *	wa4 is a work array of length m.
 *
 *   the following parameters are newly introduced in this C translation:
 *
 *      evaluate is the name of the subroutine which calculates the functions.
 *        a default implementation lm_evaluate_default is provided in lm_eval.c;
 *        alternatively, evaluate can be provided by a user calling program.
 *        it should be written as follows:
 *
 *        void evaluate ( double* par, int m_dat, double* fvec, 
 *                       void *data, int *info )
 *        {
 *           // for ( i=0; i<m_dat; ++i )
 *           //     calculate fvec[i] for given parameters par;
 *           // to stop the minimization, 
 *           //     set *info to a negative integer.
 *        }
 *
 *      printout is the name of the subroutine which nforms about fit progress.
 *        a default implementation lm_print_default is provided in lm_eval.c;
 *        alternatively, printout can be provided by a user calling program.
 *        it should be written as follows:
 *
 *        void printout ( int n_par, double* par, int m_dat, double* fvec, 
 *                       void *data, int iflag, int iter, int nfev )
 *        {
 *           // iflag : 0 (init) 1 (outer loop) 2(inner loop) -1(terminated)
 *           // iter  : outer loop counter
 *           // nfev  : number of calls to *evaluate
 *        }
 *
 *      data is an input pointer to an arbitrary structure that is passed to
 *        evaluate. typically, it contains experimental data to be fitted.
 *
 */
    int i, iter, j;
    double actred, delta, dirder, eps, fnorm, fnorm1, gnorm, par, pnorm,
        prered, ratio, step, sum, temp, temp1, temp2, temp3, xnorm;
    static double p1 = 0.1;
    static double p5 = 0.5;
    static double p25 = 0.25;
    static double p75 = 0.75;
    static double p0001 = 1.0e-4;

    *nfev = 0; // function evaluation counter
    iter = 1;  // outer loop counter
    par = 0;   // levenberg-marquardt parameter 
    delta = 0; // just to prevent a warning (initialization within if-clause)
    xnorm = 0; // dito

    temp = MAX(epsfcn,LM_MACHEP);
    eps = sqrt(temp); // used in calculating the Jacobian by forward differences

// *** check the input parameters for errors.

    if ( (n <= 0) || (m < n) || (ftol < 0.)
	|| (xtol < 0.) || (gtol < 0.) || (maxfev <= 0) || (factor <= 0.) )
    {
        *info = 0; // invalid parameter
        return;
    }
    if ( mode == 2 )  /* scaling by diag[] */
    {
	for ( j=0; j<n; j++ )  /* check for nonpositive elements */
        {
            if ( diag[j] <= 0.0 )
            {
                *info = 0; // invalid parameter
                return;
            }
        }	
    }
#if BUG
    printf( "lmdif\n" );
#endif

// *** evaluate the function at the starting point and calculate its norm.

    *info = 0;
    (*evaluate)( x, m, fvec, data, info );
    (*printout)( n, x, m, fvec, data, 0, 0, ++(*nfev) );
    if ( *info < 0 ) return;
    fnorm = lm_enorm(m,fvec);

// *** the outer loop.

    do { 
#if BUG 
        printf( "lmdif/ outer loop iter=%d nfev=%d fnorm=%.10e\n",
                iter, *nfev, fnorm );
#endif

// O** calculate the jacobian matrix.

        for ( j=0; j<n; j++ )
        {
            temp = x[j];
            step = eps * fabs(temp);
            if (step == 0.) step = eps;
            x[j] = temp + step;
            *info = 0;
            (*evaluate)( x, m, wa4, data, info );
            (*printout)( n, x, m, wa4, data, 1, iter, ++(*nfev) );
            if ( *info < 0 ) return;  // user requested break
            x[j] = temp;
            for ( i=0; i<m; i++ )
                fjac[j*m+i] = (wa4[i] - fvec[i]) / step;
        }
#if BUG>1
        // DEBUG: print the entire matrix
        for ( i=0; i<m; i++ )
        {
            for ( j=0; j<n; j++ )
                printf( "%.5e ", y[j*m+i] );
            printf( "\n" );
        }
#endif

// O** compute the qr factorization of the jacobian.

        lm_qrfac( m, n, fjac, 1, ipvt, wa1, wa2, wa3);

// O** on the first iteration ... 

        if (iter == 1)
        {
            if (mode != 2)
//      ... scale according to the norms of the columns of the initial jacobian.
            {
                for ( j=0; j<n; j++ )
                {
                    diag[j] = wa2[j];
                    if ( wa2[j] == 0. )
                        diag[j] = 1.;
                }
            }

//      ... calculate the norm of the scaled x and 
//          initialize the step bound delta.

            for ( j=0; j<n; j++ )
                wa3[j] = diag[j] * x[j];

            xnorm = lm_enorm( n, wa3 );
            delta = factor*xnorm;
            if (delta == 0.)
                delta = factor;
        }

// O** form (q transpose)*fvec and store the first n components in qtf.

        for ( i=0; i<m; i++ )
            wa4[i] = fvec[i];

        for ( j=0; j<n; j++ )
        {
            temp3 = fjac[j*m+j];
            if (temp3 != 0.)
            {
                sum = 0;
                for ( i=j; i<m; i++ )
                    sum += fjac[j*m+i] * wa4[i];
                temp = -sum / temp3;
                for ( i=j; i<m; i++ )
                    wa4[i] += fjac[j*m+i] * temp;
            }
            fjac[j*m+j] = wa1[j];
            qtf[j] = wa4[j];
        }

// O** compute the norm of the scaled gradient and test for convergence.

        gnorm = 0;
        if ( fnorm != 0 )
        {
            for ( j=0; j<n; j++ )
            {
                if ( wa2[ ipvt[j] ] == 0 ) continue;
                
                sum = 0.;
                for ( i=0; i<=j; i++ )
                    sum += fjac[j*m+i] * qtf[i] / fnorm;
                gnorm = MAX( gnorm, fabs(sum/wa2[ ipvt[j] ]) );
            }
        }

        if ( gnorm <= gtol )
        {
            *info = 4;
            return;
        }

// O** rescale if necessary.

        if ( mode != 2 )
        {
            for ( j=0; j<n; j++ )
                diag[j] = MAX(diag[j],wa2[j]);
        }

// O** the inner loop.

        do {
#if BUG 
            printf( "lmdif/ inner loop iter=%d nfev=%d\n", iter, *nfev );
#endif

// OI* determine the levenberg-marquardt parameter.

            lm_lmpar( n,fjac,m,ipvt,diag,qtf,delta,&par,wa1,wa2,wa3,wa4 );

// OI* store the direction p and x + p. calculate the norm of p.

            for ( j=0; j<n; j++ )
            {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j]*wa1[j];
            }
            pnorm = lm_enorm(n,wa3);

// OI* on the first iteration, adjust the initial step bound.

            if ( *nfev<= 1+n ) // bug corrected by J. Wuttke in 2004
                delta = MIN(delta,pnorm);

// OI* evaluate the function at x + p and calculate its norm.

            *info = 0;
            (*evaluate)( wa2, m, wa4, data, info );
            (*printout)( n, x, m, wa4, data, 2, iter, ++(*nfev) );
            if ( *info < 0 ) return;  // user requested break

            fnorm1 = lm_enorm(m,wa4);
#if BUG 
            printf( "lmdif/ pnorm %.10e  fnorm1 %.10e  fnorm %.10e"
                    " delta=%.10e par=%.10e\n",
                    pnorm, fnorm1, fnorm, delta, par );
#endif

// OI* compute the scaled actual reduction.

            if ( p1*fnorm1 < fnorm )
                actred = 1 - SQR( fnorm1/fnorm );
            else
                actred = -1;

// OI* compute the scaled predicted reduction and 
//     the scaled directional derivative.

            for ( j=0; j<n; j++ )
            {
                wa3[j] = 0;
                for ( i=0; i<=j; i++ )
                    wa3[i] += fjac[j*m+i]*wa1[ ipvt[j] ];
            }
            temp1 = lm_enorm(n,wa3) / fnorm;
            temp2 = sqrt(par) * pnorm / fnorm;
            prered = SQR(temp1) + 2 * SQR(temp2);
            dirder = - ( SQR(temp1) + SQR(temp2) );

// OI* compute the ratio of the actual to the predicted reduction.

            ratio = prered!=0 ? actred/prered : 0;
#if BUG 
            printf( "lmdif/ actred=%.10e prered=%.10e ratio=%.10e"
                    " sq(1)=%.10e sq(2)=%.10e dd=%.10e\n",
                    actred, prered, prered!=0 ? ratio : 0.,
                    SQR(temp1), SQR(temp2), dirder );
#endif

// OI* update the step bound.

            if (ratio <= p25)
            {
                if (actred >= 0.)
                    temp = p5;
                else
                    temp = p5*dirder/(dirder + p5*actred);
                if ( p1*fnorm1 >= fnorm || temp < p1 )
                    temp = p1;
                delta = temp * MIN(delta,pnorm/p1);
                par /= temp;
            }
            else if ( par == 0. || ratio >= p75 )
            {
                delta = pnorm/p5;
                par *= p5;
            }

// OI* test for successful iteration...

            if (ratio >= p0001)
            {

//     ... successful iteration. update x, fvec, and their norms.

                for ( j=0; j<n; j++ )
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j]*x[j];
                }
                for ( i=0; i<m; i++ )
                    fvec[i] = wa4[i];
                xnorm = lm_enorm(n,wa2);
                fnorm = fnorm1;
                iter++;
            }
#if BUG 
            else {
                printf( "ATTN: iteration considered unsuccessful\n" );
            } 
#endif

// OI* tests for convergence ( otherwise *info = 1, 2, or 3 )

            *info = 0; // do not terminate (unless overwritten by nonzero value)
            if ( fabs(actred) <= ftol && prered <= ftol && p5*ratio <= 1 )
                *info = 1;
            if (delta <= xtol*xnorm)
                *info += 2;
            if ( *info != 0)
                return;

// OI* tests for termination and stringent tolerances.

            if ( *nfev >= maxfev)
                *info = 5;
            if ( fabs(actred) <= LM_MACHEP &&
                 prered <= LM_MACHEP && p5*ratio <= 1 )
                *info = 6;
            if (delta <= LM_MACHEP*xnorm)
                *info = 7;
            if (gnorm <= LM_MACHEP)
                *info = 8;
            if ( *info != 0)
                return;

// OI* end of the inner loop. repeat if iteration unsuccessful.

        } while (ratio < p0001);

// O** end of the outer loop.

    } while (1);
	
}



void lm_lmpar(int n, double* r, int ldr, int* ipvt, double* diag, double* qtb,
              double delta, double* par, double* x, double* sdiag,
              double* wa1, double* wa2)
{
/*     given an m by n matrix a, an n by n nonsingular diagonal
 *     matrix d, an m-vector b, and a positive number delta,
 *     the problem is to determine a value for the parameter
 *     par such that if x solves the system
 *
 *	    a*x = b ,	  sqrt(par)*d*x = 0 ,
 *
 *     in the least squares sense, and dxnorm is the euclidean
 *     norm of d*x, then either par is 0. and
 *
 *	    (dxnorm-delta) .le. 0.1*delta ,
 *
 *     or par is positive and
 *
 *	    abs(dxnorm-delta) .le. 0.1*delta .
 *
 *     this subroutine completes the solution of the problem
 *     if it is provided with the necessary information from the
 *     qr factorization, with column pivoting, of a. that is, if
 *     a*p = q*r, where p is a permutation matrix, q has orthogonal
 *     columns, and r is an upper triangular matrix with diagonal
 *     elements of nonincreasing magnitude, then lmpar expects
 *     the full upper triangle of r, the permutation matrix p,
 *     and the first n components of (q transpose)*b. on output
 *     lmpar also provides an upper triangular matrix s such that
 *
 *	     t	 t		     t
 *	    p *(a *a + par*d*d)*p = s *s .
 *
 *     s is employed within lmpar and may be of separate interest.
 *
 *     only a few iterations are generally needed for convergence
 *     of the algorithm. if, however, the limit of 10 iterations
 *     is reached, then the output par will contain the best
 *     value obtained so far.
 *
 *     parameters:
 *
 *	n is a positive integer input variable set to the order of r.
 *
 *	r is an n by n array. on input the full upper triangle
 *	  must contain the full upper triangle of the matrix r.
 *	  on output the full upper triangle is unaltered, and the
 *	  strict lower triangle contains the strict upper triangle
 *	  (transposed) of the upper triangular matrix s.
 *
 *	ldr is a positive integer input variable not less than n
 *	  which specifies the leading dimension of the array r.
 *
 *	ipvt is an integer input array of length n which defines the
 *	  permutation matrix p such that a*p = q*r. column j of p
 *	  is column ipvt(j) of the identity matrix.
 *
 *	diag is an input array of length n which must contain the
 *	  diagonal elements of the matrix d.
 *
 *	qtb is an input array of length n which must contain the first
 *	  n elements of the vector (q transpose)*b.
 *
 *	delta is a positive input variable which specifies an upper
 *	  bound on the euclidean norm of d*x.
 *
 *	par is a nonnegative variable. on input par contains an
 *	  initial estimate of the levenberg-marquardt parameter.
 *	  on output par contains the final estimate.
 *
 *	x is an output array of length n which contains the least
 *	  squares solution of the system a*x = b, sqrt(par)*d*x = 0,
 *	  for the output par.
 *
 *	sdiag is an output array of length n which contains the
 *	  diagonal elements of the upper triangular matrix s.
 *
 *	wa1 and wa2 are work arrays of length n.
 *
 */
    int i, iter, j, nsing;
    double dxnorm, fp, fp_old, gnorm, parc, parl, paru;
    double sum, temp;
    static double p1 = 0.1;
    static double p001 = 0.001;

#if BUG
    printf( "lmpar\n" );
#endif

// *** compute and store in x the gauss-newton direction. if the
//     jacobian is rank-deficient, obtain a least squares solution.

    nsing = n;
    for ( j=0; j<n; j++ )
    {
	wa1[j] = qtb[j];
	if ( r[j*ldr+j] == 0 && nsing == n )
            nsing = j;
	if (nsing < n)
            wa1[j] = 0;
    }
#if BUG
    printf( "nsing %d ", nsing );
#endif
    for ( j=nsing-1; j>=0; j-- )
    {
        wa1[j] = wa1[j]/r[j+ldr*j];
        temp = wa1[j];
        for ( i=0; i<j; i++ )
            wa1[i] -= r[j*ldr+i]*temp;
    }

    for ( j=0; j<n; j++ )
	x[ ipvt[j] ] = wa1[j];

// *** initialize the iteration counter.
//     evaluate the function at the origin, and test
//     for acceptance of the gauss-newton direction.

    iter = 0;
    for ( j=0; j<n; j++ )
	wa2[j] = diag[j]*x[j];
    dxnorm = lm_enorm(n,wa2);
    fp = dxnorm - delta;
    if (fp <= p1*delta)
    {
#if BUG
	printf( "lmpar/ terminate (fp<delta/10\n" );
#endif
        *par = 0;
        return;
    }

// *** if the jacobian is not rank deficient, the newton
//     step provides a lower bound, parl, for the 0. of
//     the function. otherwise set this bound to 0..

    parl = 0;
    if (nsing >= n)
    {
	for ( j=0; j<n; j++ )
            wa1[j] = diag[ ipvt[j] ] * wa2[ ipvt[j] ] / dxnorm;

	for ( j=0; j<n; j++ )
        {
            sum = 0.;
            for ( i=0; i<j; i++ )
                sum += r[j*ldr+i]*wa1[i];
            wa1[j] = (wa1[j] - sum)/r[j+ldr*j];
        }
	temp = lm_enorm(n,wa1);
	parl = fp/delta/temp/temp;
    }

// *** calculate an upper bound, paru, for the 0. of the function.

    for ( j=0; j<n; j++ )
    {
	sum = 0;
	for ( i=0; i<=j; i++ )
            sum += r[j*ldr+i]*qtb[i];
	wa1[j] = sum/diag[ ipvt[j] ];
    }
    gnorm = lm_enorm(n,wa1);
    paru = gnorm/delta;
    if (paru == 0.)
	paru = LM_DWARF/MIN(delta,p1);

// *** if the input par lies outside of the interval (parl,paru),
//     set par to the closer endpoint.

    *par = MAX( *par,parl);
    *par = MIN( *par,paru);
    if ( *par == 0.)
	*par = gnorm/dxnorm;
#if BUG
    printf( "lmpar/ parl %.4e  par %.4e  paru %.4e\n", parl, *par, paru );
#endif

// *** iterate.

    for ( ; ; iter++ ) {

// *** evaluate the function at the current value of par.

        if ( *par == 0.)
            *par = MAX(LM_DWARF,p001*paru);
        temp = sqrt( *par );
        for ( j=0; j<n; j++ )
            wa1[j] = temp*diag[j];
        lm_qrsolv( n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
        for ( j=0; j<n; j++ )
            wa2[j] = diag[j]*x[j];
        dxnorm = lm_enorm(n,wa2);
        fp_old = fp;
        fp = dxnorm - delta;

// ***	 if the function is small enough, accept the current value
//	 of par. also test for the exceptional cases where parl
//	 is 0. or the number of iterations has reached 10.

        if ( fabs(fp) <= p1*delta
             || (parl == 0. && fp <= fp_old && fp_old < 0.)
             || iter == 10 )
            break; // the only exit from this loop

// *** compute the Newton correction.

        for ( j=0; j<n; j++ )
            wa1[j] = diag[ ipvt[j] ] * wa2[ ipvt[j] ] / dxnorm;

        for ( j=0; j<n; j++ )
        {
            wa1[j] = wa1[j]/sdiag[j];
            for ( i=j+1; i<n; i++ )
                wa1[i] -= r[j*ldr+i]*wa1[j];
        }
        temp = lm_enorm( n, wa1);
        parc = fp/delta/temp/temp;

// *** depending on the sign of the function, update parl or paru.

        if (fp > 0)
            parl = MAX(parl, *par);
        else if (fp < 0)
            paru = MIN(paru, *par);
        // the case fp==0 is precluded by the break condition 

// *** compute an improved estimate for par.

        *par = MAX(parl, *par + parc);

    }

}



void lm_qrfac(int m, int n, double* a, int pivot, int* ipvt,
           double* rdiag, double* acnorm, double* wa)
{
/*
 *     this subroutine uses householder transformations with column
 *     pivoting (optional) to compute a qr factorization of the
 *     m by n matrix a. that is, qrfac determines an orthogonal
 *     matrix q, a permutation matrix p, and an upper trapezoidal
 *     matrix r with diagonal elements of nonincreasing magnitude,
 *     such that a*p = q*r. the householder transformation for
 *     column k, k = 1,2,...,min(m,n), is of the form
 *
 *			    t
 *	    i - (1/u(k))*u*u
 *
 *     where u has 0.s in the first k-1 positions. the form of
 *     this transformation and the method of pivoting first
 *     appeared in the corresponding linpack subroutine.
 *
 *     parameters:
 *
 *	m is a positive integer input variable set to the number
 *	  of rows of a.
 *
 *	n is a positive integer input variable set to the number
 *	  of columns of a.
 *
 *	a is an m by n array. on input a contains the matrix for
 *	  which the qr factorization is to be computed. on output
 *	  the strict upper trapezoidal part of a contains the strict
 *	  upper trapezoidal part of r, and the lower trapezoidal
 *	  part of a contains a factored form of q (the non-trivial
 *	  elements of the u vectors described above).
 *
 *	pivot is a logical input variable. if pivot is set true,
 *	  then column pivoting is enforced. if pivot is set false,
 *	  then no column pivoting is done.
 *
 *	ipvt is an integer output array of length lipvt. ipvt
 *	  defines the permutation matrix p such that a*p = q*r.
 *	  column j of p is column ipvt(j) of the identity matrix.
 *	  if pivot is false, ipvt is not referenced.
 *
 *	rdiag is an output array of length n which contains the
 *	  diagonal elements of r.
 *
 *	acnorm is an output array of length n which contains the
 *	  norms of the corresponding columns of the input matrix a.
 *	  if this information is not needed, then acnorm can coincide
 *	  with rdiag.
 *
 *	wa is a work array of length n. if pivot is false, then wa
 *	  can coincide with rdiag.
 *
 */
    int i, j, k, kmax, minmn;
    double ajnorm, sum, temp;
    static double p05 = 0.05;

// *** compute the initial column norms and initialize several arrays.

    for ( j=0; j<n; j++ )
    {
	acnorm[j] = lm_enorm(m, &a[j*m]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	if ( pivot )
            ipvt[j] = j;
    }
#if BUG
    printf( "qrfac\n" );
#endif

// *** reduce a to r with householder transformations.

    minmn = MIN(m,n);
    for ( j=0; j<minmn; j++ )
    {
        if ( !pivot ) goto pivot_ok;

// *** bring the column of largest norm into the pivot position.

        kmax = j;
        for ( k=j+1; k<n; k++ )
            if (rdiag[k] > rdiag[kmax])
		kmax = k;
        if (kmax == j) goto pivot_ok; // bug fixed in rel 2.1

        for ( i=0; i<m; i++ )
	{
            temp        = a[j*m+i];
            a[j*m+i]    = a[kmax*m+i];
            a[kmax*m+i] = temp;
	}
        rdiag[kmax] = rdiag[j];
        wa[kmax] = wa[j];
        k = ipvt[j];
        ipvt[j] = ipvt[kmax];
        ipvt[kmax] = k;

    pivot_ok:

// *** compute the Householder transformation to reduce the
//     j-th column of a to a multiple of the j-th unit vector.

        ajnorm = lm_enorm( m-j, &a[j*m+j] );
        if (ajnorm == 0.)
        {
            rdiag[j] = 0;
            continue;
        }

        if (a[j*m+j] < 0.)
            ajnorm = -ajnorm;
        for ( i=j; i<m; i++ )
            a[j*m+i] /= ajnorm;
        a[j*m+j] += 1;

// *** apply the transformation to the remaining columns
//     and update the norms.

        for ( k=j+1; k<n; k++ )
        {
            sum = 0;

            for ( i=j; i<m; i++ )
                sum += a[j*m+i]*a[k*m+i];

            temp = sum/a[j+m*j];

            for ( i=j; i<m; i++ )
                a[k*m+i] -= temp * a[j*m+i];

            if ( pivot && rdiag[k] != 0. )
            {
                temp = a[m*k+j]/rdiag[k];
                temp = MAX( 0., 1-temp*temp );
                rdiag[k] *= sqrt(temp);
                temp = rdiag[k]/wa[k];
                if ( p05*SQR(temp) <= LM_MACHEP )
                {
                    rdiag[k] = lm_enorm( m-j-1, &a[m*k+j+1]);
                    wa[k] = rdiag[k];
                }
            }
        }

	rdiag[j] = -ajnorm;
    }
}



void lm_qrsolv(int n, double* r, int ldr, int* ipvt, double* diag,
              double* qtb, double* x, double* sdiag, double* wa)
{
/*
 *     given an m by n matrix a, an n by n diagonal matrix d,
 *     and an m-vector b, the problem is to determine an x which
 *     solves the system
 *
 *	    a*x = b ,	  d*x = 0 ,
 *
 *     in the least squares sense.
 *
 *     this subroutine completes the solution of the problem
 *     if it is provided with the necessary information from the
 *     qr factorization, with column pivoting, of a. that is, if
 *     a*p = q*r, where p is a permutation matrix, q has orthogonal
 *     columns, and r is an upper triangular matrix with diagonal
 *     elements of nonincreasing magnitude, then qrsolv expects
 *     the full upper triangle of r, the permutation matrix p,
 *     and the first n components of (q transpose)*b. the system
 *     a*x = b, d*x = 0, is then equivalent to
 *
 *		   t	   t
 *	    r*z = q *b ,  p *d*p*z = 0 ,
 *
 *     where x = p*z. if this system does not have full rank,
 *     then a least squares solution is obtained. on output qrsolv
 *     also provides an upper triangular matrix s such that
 *
 *	     t	 t		 t
 *	    p *(a *a + d*d)*p = s *s .
 *
 *     s is computed within qrsolv and may be of separate interest.
 *
 *     parameters
 *
 *	n is a positive integer input variable set to the order of r.
 *
 *	r is an n by n array. on input the full upper triangle
 *	  must contain the full upper triangle of the matrix r.
 *	  on output the full upper triangle is unaltered, and the
 *	  strict lower triangle contains the strict upper triangle
 *	  (transposed) of the upper triangular matrix s.
 *
 *	ldr is a positive integer input variable not less than n
 *	  which specifies the leading dimension of the array r.
 *
 *	ipvt is an integer input array of length n which defines the
 *	  permutation matrix p such that a*p = q*r. column j of p
 *	  is column ipvt(j) of the identity matrix.
 *
 *	diag is an input array of length n which must contain the
 *	  diagonal elements of the matrix d.
 *
 *	qtb is an input array of length n which must contain the first
 *	  n elements of the vector (q transpose)*b.
 *
 *	x is an output array of length n which contains the least
 *	  squares solution of the system a*x = b, d*x = 0.
 *
 *	sdiag is an output array of length n which contains the
 *	  diagonal elements of the upper triangular matrix s.
 *
 *	wa is a work array of length n.
 *
 */
    int i, kk, j, k, nsing;
    double qtbpj, sum, temp;
    double sin, cos, tan, cotan; // these are local variables, not functions
    static double p25 = 0.25;
    static double p5 = 0.5;

// *** copy r and (q transpose)*b to preserve input and initialize s.
//     in particular, save the diagonal elements of r in x.

    for ( j=0; j<n; j++ )
    {
	for ( i=j; i<n; i++ )
            r[j*ldr+i] = r[i*ldr+j];
	x[j] = r[j*ldr+j];
	wa[j] = qtb[j];
    }
#if BUG
    printf( "qrsolv\n" );
#endif

// *** eliminate the diagonal matrix d using a givens rotation.

    for ( j=0; j<n; j++ )
    {

// ***	 prepare the row of d to be eliminated, locating the
// 	 diagonal element using p from the qr factorization.

        if (diag[ ipvt[j] ] == 0.)
            goto L90;
        for ( k=j; k<n; k++ )
            sdiag[k] = 0.;
        sdiag[j] = diag[ ipvt[j] ];

// ***	 the transformations to eliminate the row of d
//	 modify only a single element of (q transpose)*b
//	 beyond the first n, which is initially 0..

        qtbpj = 0.;
        for ( k=j; k<n; k++ )
	{

//	    determine a givens rotation which eliminates the
//	    appropriate element in the current row of d.

            if (sdiag[k] == 0.)
		continue;
            kk = k + ldr * k; // <! keep this shorthand !>
            if ( fabs(r[kk]) < fabs(sdiag[k]) )
            {
		cotan = r[kk]/sdiag[k];
		sin = p5/sqrt(p25+p25*SQR(cotan));
		cos = sin*cotan;
            }
            else
            {
		tan = sdiag[k]/r[kk];
		cos = p5/sqrt(p25+p25*SQR(tan));
		sin = cos*tan;
            }

// ***	    compute the modified diagonal element of r and
//	    the modified element of ((q transpose)*b,0).

            r[kk] = cos*r[kk] + sin*sdiag[k];
            temp = cos*wa[k] + sin*qtbpj;
            qtbpj = -sin*wa[k] + cos*qtbpj;
            wa[k] = temp;

// *** accumulate the tranformation in the row of s.

            for ( i=k+1; i<n; i++ )
            {
                temp = cos*r[k*ldr+i] + sin*sdiag[i]; 
                sdiag[i] = -sin*r[k*ldr+i] + cos*sdiag[i];
                r[k*ldr+i] = temp;
            }
	}
    L90:

// *** store the diagonal element of s and restore
//     the corresponding diagonal element of r.

	sdiag[j] = r[j*ldr+j];
	r[j*ldr+j] = x[j];
    }

// *** solve the triangular system for z. if the system is
//     singular, then obtain a least squares solution.

    nsing = n;
    for ( j=0; j<n; j++ )
    {
	if ( sdiag[j] == 0. && nsing == n )
            nsing = j;
	if (nsing < n)
            wa[j] = 0;
    }

    for ( j=nsing-1; j>=0; j-- )
    {
	sum = 0;
        for ( i=j+1; i<nsing; i++ )
            sum += r[j*ldr+i]*wa[i];
	wa[j] = (wa[j] - sum)/sdiag[j];
    }

// *** permute the components of z back to components of x.

    for ( j=0; j<n; j++ )
	x[ ipvt[j] ] = wa[j];
}


 
double lm_enorm( int n, double *x )
{
/*     given an n-vector x, this function calculates the
 *     euclidean norm of x.
 *
 *     the euclidean norm is computed by accumulating the sum of
 *     squares in three different sums. the sums of squares for the
 *     small and large components are scaled so that no overflows
 *     occur. non-destructive underflows are permitted. underflows
 *     and overflows do not occur in the computation of the unscaled
 *     sum of squares for the intermediate components.
 *     the definitions of small, intermediate and large components
 *     depend on two constants, LM_SQRT_DWARF and LM_SQRT_GIANT. the main
 *     restrictions on these constants are that LM_SQRT_DWARF**2 not
 *     underflow and LM_SQRT_GIANT**2 not overflow.
 *
 *     parameters
 *
 *	n is a positive integer input variable.
 *
 *	x is an input array of length n.
 */
    int i;
    double agiant, s1, s2, s3, xabs, x1max, x3max, temp;

    s1 = 0;
    s2 = 0;
    s3 = 0;
    x1max = 0;
    x3max = 0;
    agiant = LM_SQRT_GIANT/( (double) n);

    for ( i=0; i<n; i++ )
    {
        xabs = fabs(x[i]);
        if ( xabs > LM_SQRT_DWARF && xabs < agiant )
	{
// **  sum for intermediate components.
            s2 += xabs*xabs;
            continue;
	}

        if ( xabs >  LM_SQRT_DWARF )
	{
// **  sum for large components.
            if (xabs > x1max)
            {
		temp = x1max/xabs;
		s1 = 1 + s1*SQR(temp);
		x1max = xabs;
            }
            else
            {
		temp = xabs/x1max;
		s1 += SQR(temp);
            }
            continue;
	}
// **  sum for small components.
        if (xabs > x3max)
	{
            temp = x3max/xabs;
            s3 = 1 + s3*SQR(temp);
            x3max = xabs;
	}
        else	
	{
            if (xabs != 0.)
            {
		temp = xabs/x3max;
		s3 += SQR(temp);
            }
	}
    }

// *** calculation of norm.

    if (s1 != 0)
	return x1max*sqrt(s1 + (s2/x1max)/x1max);
    if (s2 != 0)
    {
	if (s2 >= x3max)
            return sqrt( s2*(1+(x3max/s2)*(x3max*s3)) );
	else
            return sqrt( x3max*((s2/x3max)+(x3max*s3)) );
    }

    return x3max*sqrt(s3);
}
