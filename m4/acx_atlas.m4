dnl
dnl				acx_atlas.m4
dnl
dnl Figure out if the ATLAS library and header files are installed.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2003-2016 IAP/CNRS/UPMC
dnl
dnl	License:		GNU General Public License
dnl
dnl	AstrOmatic software is free software: you can redistribute it and/or
dnl	modify it under the terms of the GNU General Public License as
dnl	published by the Free Software Foundation, either version 3 of the
dnl	License, or (at your option) any later version.
dnl	AstrOmatic software is distributed in the hope that it will be useful,
dnl	but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl	GNU General Public License for more details.
dnl	You should have received a copy of the GNU General Public License
dnl	along with AstrOmatic software.
dnl	If not, see <http://www.gnu.org/licenses/>.
dnl
dnl	Last modified:		19/10/2016
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_ATLAS([ATLAS_LIBSDIR, ATLAS_INCDIR, ATLAS_PFLAG,
dnl                     [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$ATLAS_LIBS $LIBS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if BLAS/LAPACK
dnl is found (HAVE_ATLAS is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_ATLAS], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl --------------------

ATLAS_ERROR=""
if test x$2 = x; then
    acx_atlas_incdir=""
    AC_CHECK_HEADERS([cblas.h clapack.h],,
      [
      acx_atlas_incdir="atlas/"
      AC_CHECK_HEADERS([${acx_atlas_incdir}cblas.h ${acx_atlas_incdir}clapack.h],,
        [ATLAS_ERROR="ATLAS header files not found!"])
      ]
    )
  else
    acx_atlas_incdir="$2/"
    AC_CHECK_HEADERS([${acx_atlas_incdir}cblas.h ${acx_atlas_incdir}clapack.h],,
    [
      [acx_atlas_incdir="$2/include/"]
      AC_CHECK_HEADER(
       [${acx_atlas_incdir}cblas.h ${acx_atlas_incdir}clapack.h],,
        [ATLAS_ERROR="ATLAS header files not found in "$2"!"]
      )
    ]
  )
fi

if test x$ATLAS_ERROR = x; then
  AC_DEFINE_UNQUOTED(ATLAS_BLAS_H, "${acx_atlas_incdir}cblas.h", [BLAS header filename.])
  AC_DEFINE_UNQUOTED(ATLAS_LAPACK_H, "${acx_atlas_incdir}clapack.h", [CLAPACK header filename.])

dnl --------------------
dnl Search library files
dnl --------------------

  OLIBS="$LIBS"
  LIBS=""
  if test x$1 = x; then
    if test -d "/usr/lib64/atlas"; then
      acx_atlas_libopt="-L/usr/lib64/atlas"
    elif test -d "/usr/lib/atlas"; then
      acx_atlas_libopt="-L/usr/lib/atlas"
    else
      acx_atlas_libopt=""
    fi
  else
    acx_atlas_libopt="-L$1"
  fi
  if test x$3 == xyes; then
dnl Parallel ATLAS 3.10+
    acx_atlas_newlibs="tatlas"
dnl Older parallel ATLAS
    acx_atlas_oldextralibs="-latlas -lptcblas -lcblas"
  else
dnl Serial ATLAS 3.10+
    acx_atlas_newlibs="satlas"
dnl Older serial ATLAS
    acx_atlas_oldextralibs="-latlas -lcblas"
  fi
  acx_atlas_extralibs=""
  AC_SEARCH_LIBS(
    [clapack_dpotrf], [$acx_atlas_newlibs],,
    [
      unset ac_cv_search_clapack_dpotrf
      acx_atlas_extralibs=$acx_atlas_oldextralibs
      AC_SEARCH_LIBS(
        [clapack_dpotrf], [lapack_atlas lapack],,
        [
          unset ac_cv_search_clapack_dpotrf
          acx_atlas_extralibs=""
          AC_SEARCH_LIBS(
            [clapack_dpotrf], [atlas],
            [ATLAS_WARN="Parallel ATLAS not found, reverting to serial!"],
            [
              unset ac_cv_search_clapack_dpotrf
              acx_atlas_extralibs="-latlas -lcblas"
              AC_SEARCH_LIBS(
                [clapack_dpotrf], [lapack_atlas lapack],
                [ATLAS_WARN="Parallel ATLAS not found, reverting to serial!"],
                [ATLAS_ERROR="ATLAS library files not found!"],
                [$acx_atlas_libopt $acx_atlas_extralibs]
              )
            ],
            $acx_atlas_libopt
          )
        ],
        [$acx_atlas_libopt $acx_atlas_extralibs]
      )
    ],
    $acx_atlas_libopt
  )
  LIBS="$OLIBS"
fi

dnl -------------------------------------------------------------------------
dnl Finally execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND
dnl -------------------------------------------------------------------------

if test "x$ATLAS_ERROR" = "x"; then
  AC_DEFINE(HAVE_ATLAS,1,
	[Define if you have the ATLAS libraries and header files.])
  ATLAS_LIBS="$acx_atlas_libopt $ac_cv_search_clapack_dpotrf"
  AC_SUBST(ATLAS_CFLAGS)
  AC_SUBST(ATLAS_LDFLAGS, "")
  AC_SUBST(ATLAS_LIBS)
  AC_SUBST(ATLAS_WARN)
  $4
else
  AC_SUBST(ATLAS_ERROR)
  $5
fi

])dnl ACX_ATLAS

