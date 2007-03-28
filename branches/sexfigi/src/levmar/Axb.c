/////////////////////////////////////////////////////////////////////////////////
// 
//  Solution of linear systems involved in the Levenberg - Marquardt
//  minimization algorithm
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

/******************************************************************************** 
 * LAPACK-based implementations for various linear system solvers. The same core
 * code is used with appropriate #defines to derive single and double precision
 * solver versions, see also Axb_core.c
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lm.h"
#include "misc.h"

/* double precision definitions */
#define LM_REAL double
#define LM_PREFIX d
#define CNST(x) (x)
#ifndef HAVE_LAPACK
#include <float.h>
#define LM_REAL_EPSILON DBL_EPSILON
#endif

#include "Axb_core.c"

#undef LM_REAL
#undef LM_PREFIX
#undef CNST
#undef LM_REAL_EPSILON

/* single precision (float) definitions */
#define LM_REAL float
#define LM_PREFIX s
#define SUBCNST(x) x##F
#define CNST(x) SUBCNST(x) // force substitution
#ifndef HAVE_LAPACK
#define LM_REAL_EPSILON FLT_EPSILON
#endif

#include "Axb_core.c"

#undef LM_REAL
#undef LM_PREFIX
#undef SUBCNST
#undef CNST
#undef LM_REAL_EPSILON
