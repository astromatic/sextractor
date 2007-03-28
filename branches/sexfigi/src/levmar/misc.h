/////////////////////////////////////////////////////////////////////////////////
// 
//  Levenberg - Marquardt non-linear minimization algorithm
//  Copyright (C) 2004  Manolis Lourakis (lourakis@ics.forth.gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef _MISC_H_
#define _MISC_H_

/* common prefix for BLAS subroutines */
#define LM_BLAS_PREFIX // none
//#define LM_BLAS_PREFIX f2c_   // f2c'd BLAS
//#define LM_BLAS_PREFIX cblas_ // C interface to BLAS

#define LCAT_(a, b)    #a b
#define LCAT(a, b)    LCAT_(a, b) // force substitution
#define RCAT_(a, b)    a #b
#define RCAT(a, b)    RCAT_(a, b) // force substitution

#define __BLOCKSZ__       32 /* block size for cache-friendly matrix-matrix multiply. It should be
                              * such that __BLOCKSZ__^2*sizeof(LM_REAL) is smaller than the CPU (L1)
                              * data cache size. Notice that a value of 32 when LM_REAL=double assumes
                              * an 8Kb L1 data cache (32*32*8=8K). This is a concervative choice since
                              * newer Pentium 4s have a L1 data cache of size 16K, capable of holding
                              * up to 45x45 double blocks.
                              */
#define __BLOCKSZ__SQ    (__BLOCKSZ__)*(__BLOCKSZ__)

#ifdef _MSC_VER
#define inline __inline //MSVC
#elif !defined(__GNUC__)
#define inline //other than MSVC, GCC: define empty
#endif

/* add a prefix in front of a token */
#define LM_CAT__(a, b) a ## b
#define LM_CAT_(a, b) LM_CAT__(a, b) // force substitution
#define LM_ADD_PREFIX(s) LM_CAT_(LM_PREFIX, s)

#ifdef __cplusplus
extern "C" {
#endif

/* blocking-based matrix multiply */
extern void strans_mat_mat_mult(float *a, float *b, int n, int m);
extern void dtrans_mat_mat_mult(double *a, double *b, int n, int m);

/* forward finite differences */
extern void sfdif_forw_jac_approx(void (*func)(float *p, float *hx, int m, int n, void *adata),
					float *p, float *hx, float *hxx, float delta,
					float *jac, int m, int n, void *adata);
extern void dfdif_forw_jac_approx(void (*func)(double *p, double *hx, int m, int n, void *adata),
					double *p, double *hx, double *hxx, double delta,
					double *jac, int m, int n, void *adata);

/* central finite differences */
extern void sfdif_cent_jac_approx(void (*func)(float *p, float *hx, int m, int n, void *adata),
          float *p, float *hxm, float *hxp, float delta,
          float *jac, int m, int n, void *adata);
extern void dfdif_cent_jac_approx(void (*func)(double *p, double *hx, int m, int n, void *adata),
          double *p, double *hxm, double *hxp, double delta,
          double *jac, int m, int n, void *adata);

/* covariance of LS fit */
extern int slevmar_covar(float *JtJ, float *C, float sumsq, int m, int n);
extern int dlevmar_covar(double *JtJ, double *C, double sumsq, int m, int n);

#ifdef __cplusplus
}
#endif

#endif /* _MISC_H */
