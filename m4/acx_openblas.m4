dnl
dnl				acx_openblas.m4
dnl
dnl Figure out if the OpenBLAS library and header files are installed.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2016 IAP/CNRS/UPMC
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
dnl	Last modified:		12/10/2016
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_OPENBLAS([OPENBLAS_LIBSDIR, OPENBLAS_INCDIR, OPENBLAS_PFLAG,
dnl                  ILP64_FLAG, [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$OPENBLAS_LIBS $LIBS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if OPENBLAS
dnl is found (HAVE_OPENBLAS is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_OPENBLAS], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl --------------------

OPENBLAS_ERROR=""
if test x$2 = x; then
  [acx_openblas_incdir="openblas/"]
  AC_CHECK_HEADERS(
    [${acx_openblas_incdir}cblas.h ${acx_openblas_incdir}lapacke.h],,
    [
      [acx_openblas_incdir=""]
      AC_CHECK_HEADER(
        [cblas.h lapacke.h],,
        [OPENBLAS_ERROR="OpenBLAS header files not found!"]
      )
    ]
  )
else
  acx_openblas_incdir="$2/"
  AC_CHECK_HEADER(
    [${acx_openblas_incdir}cblas.h],,
    [
      [acx_openblas_incdir="$2/include/"]
      AC_CHECK_HEADERS(
        [${acx_openblas_incdir}cblas.h ${acx_openblas_incdir}lapacke.h],,
        [OPENBLAS_ERROR="OpenBLAS header files not found in "$2"!"]
    )]
  )
fi

if test "x$OPENBLAS_ERROR" = "x"; then
  AC_DEFINE_UNQUOTED(BLAS_H, "${acx_openblas_incdir}cblas.h", [BLAS header filename.])
  AC_DEFINE_UNQUOTED(LAPACKE_H, "${acx_openblas_incdir}lapacke.h", [LAPACKe header filename.])

dnl ----------------------------
dnl Search OpenBLAS library file
dnl ----------------------------

  OLIBS="$LIBS"
  LIBS=""
  if test x$4 = xyes; then
    acx_openblas_suffix="64"
    OPENBLAS_CFLAGS="-DOPENBLAS_USE64BITINT -DLAPACK_ILP64"
  else
    acx_openblas_suffix=""
    OPENBLAS_CFLAGS=""
  fi
  if test x$1 = x; then
    acx_openblas_libopt=""
  else
    acx_openblas_libopt="-L$1"
  fi
  if test x$3 == xyes; then
    AC_SEARCH_LIBS(
      [LAPACKE_dpotrf], ["openblasp"$acx_openblas_suffix],
      AC_DEFINE(HAVE_OPENBLASP,1,
		[Define if you have the OpenBLAS parallel library and header files.]),
      unset ac_cv_search_LAPACKE_dpotrf
      [AC_SEARCH_LIBS(
        [LAPACKE_dpotrf], ["openblas"$acx_openblas_suffix],
	[OPENBLAS_WARN="parallel OpenBLAS"$acx_openblas_suffix" not found, reverting to scalar OpenBLAS"$acx_openblas_suffix"!"],
        [OPENBLAS_ERROR="OpenBLAS"$acx_openblas_suffix" library file not found!"],
        $acx_openblas_libopt
      )],
      $acx_openblas_libopt
    )
  else
    AC_SEARCH_LIBS(
      [LAPACKE_dpotrf], ["openblas"$acx_openblas_suffix],,
      [OPENBLAS_ERROR="OpenBLAS"$acx_openblas_suffix" library file not found!"],
      $acx_openblas_libopt
    )
  fi
  LIBS="$OLIBS"
fi

dnl -------------------------------------------------------------------------
dnl Finally execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND
dnl -------------------------------------------------------------------------

if test "x$OPENBLAS_ERROR" = "x"; then
  AC_DEFINE(HAVE_OPENBLAS,1, [Define if you have the OpenBLAS library and header files.])
  AC_DEFINE(HAVE_BLAS,1, [Define if you have the BLAS library and header files.])
  AC_DEFINE(HAVE_LAPACKE,1, [Define if you have the LAPACKe library and header files.])
  OPENBLAS_LIBS="$acx_openblas_libopt $ac_cv_search_LAPACKE_dpotrf"
  AC_SUBST(OPENBLAS_CFLAGS)
  AC_SUBST(OPENBLAS_LDFLAGS, "")
  AC_SUBST(OPENBLAS_LIBS)
  AC_SUBST(OPENBLAS_WARN)
  $5
else
  AC_SUBST(OPENBLAS_ERROR)
  $6
fi

])dnl ACX_OPENBLAS

