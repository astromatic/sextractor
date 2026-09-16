dnl
dnl				ax_lapack.m4
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
dnl	Last modified:		2025-07-21
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis AX_LAPACK([LAPACK_LIBSDIR, LAPACK_INCDIR,
dnl                  [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$LAPACK_LIBS $LIBS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if LAPACK
dnl is found (HAVE_LAPACK is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([AX_LAPACK], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl Under openSUSE and Fedora we have /usr/include/openblas/{lapack.h,lapacke.h,lapacke_config.h}
dnl Under Ubuntu we have /usr/include/{lapack.h,lapacke.h,lapacke_64.h}
dnl --------------------

LAPACK_ERROR=""
if test x$2 = x; then
  [acx_lapack_incdir="openblas/"]
  AC_CHECK_HEADERS(
    [${acx_lapack_incdir}lapacke.h],,
    [
      [acx_lapack_incdir=""]
      AC_CHECK_HEADER(
        [lapacke.h],,
        [LAPACK_ERROR="lapacke header files not found!"]
      )
    ]
  )
else
  acx_lapack_incdir="$2/"
  AC_CHECK_HEADER(
    [${acx_lapack_incdir}lapacke.h],,
    [
      [acx_lapack_incdir="$2/include/"]
      AC_CHECK_HEADERS(
        [${acx_lapack_incdir}lapacke.h],,
        [LAPACK_ERROR="lapacke header files not found in "$2"!"]
    )]
  )
fi

if test "x$LAPACK_ERROR" = "x"; then
  AC_DEFINE_UNQUOTED(LAPACKE_H, "${acx_lapack_incdir}lapacke.h", [LAPACKe header filename.])
  LAPACK_CFLAGS="-I ${acx_lapack_incdir}"

dnl ----------------------------
dnl Search lapack library file
dnl On Fedora we have (dnf install lapack-devel) /usr/lib64/{liblapack64*,liblapack* ,liblapacke*}
dnl On Ubuntu we have /usr/lib/x86_64-linux-gnu/liblapacke* 
dnl    and /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so
dnl On openSUSE we have /usr/lib64/lapack/liblapack.so.3, /usr/lib64/liblapack{e}.so
dnl ----------------------------

  OLIBS="$LIBS"
  LIBS=""
  if test x$1 = x; then
    acx_lapack_libopt=""
  else
    acx_lapack_libopt="-L$1"
  fi

  AC_SEARCH_LIBS(
    [LAPACKE_dpotrf], ["lapacke" "lapack" "openblas"],,
    [LAPACK_ERROR="lapacke library file not found!"],
    $acx_lapack_libopt
  )
  LIBS="$OLIBS"
fi

dnl -------------------------------------------------------------------------
dnl Finally execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND
dnl -------------------------------------------------------------------------

if test "x$LAPACK_ERROR" = "x"; then
  AC_DEFINE(HAVE_LAPACKE,1, [Define if you have the lapacke library and header files.])
  LAPACK_LIBS="$ac_cv_search_LAPACKE_dpotrf"
  AC_SUBST(LAPACK_CFLAGS)
  AC_SUBST(LAPACK_LDFLAGS, $acx_lapack_libopt)
  AC_SUBST(LAPACK_LIBS)
  AC_SUBST(LAPACK_WARN)
  $3
else
  AC_SUBST(LAPACK_ERROR)
  $4
fi

])dnl AX_LAPACK

