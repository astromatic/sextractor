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

/******************************************************************************** 
 * Levenberg-Marquardt nonlinear minimization. The same core code is used with
 * appropriate #defines to derive single and double precision versions, see
 * also lm_core.c
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "lm.h"
#include "misc.h"

#define EPSILON       1E-12
#define ONE_THIRD     0.3333333334 /* 1.0/3.0 */

/* single precision (float) definitions */
#define LM_REAL float
#define LM_PREFIX s

#define LM_REAL_MAX FLT_MAX
#define LM_REAL_MIN FLT_MIN
#define LM_REAL_EPSILON FLT_EPSILON
#define SUBCNST(x) x##F
#define CNST(x) SUBCNST(x) // force substitution

#include "lm_core.c" // read in core code

#undef LM_REAL
#undef LM_PREFIX
#undef LM_REAL_MAX
#undef LM_REAL_EPSILON
#undef LM_REAL_MIN
#undef SUBCNST
#undef CNST

/* double precision definitions */
#define LM_REAL double
#define LM_PREFIX d

#define LM_REAL_MAX DBL_MAX
#define LM_REAL_MIN DBL_MIN
#define LM_REAL_EPSILON DBL_EPSILON
#define CNST(x) (x)

#include "lm_core.c" // read in core code

#undef LM_REAL
#undef LM_PREFIX
#undef LM_REAL_MAX
#undef LM_REAL_EPSILON
#undef LM_REAL_MIN
#undef CNST
