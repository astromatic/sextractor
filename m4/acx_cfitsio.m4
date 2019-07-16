dnl
dnl				acx_cfitsio.m4
dnl
dnl Figure out if the CFITSIO library and header files are installed.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2019 IAP/CNRS/UPMC
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
dnl	Last modified:		16/07/2019
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_CFITSIO([CFITSIO_LIBSDIR, CFITSIO_INCDIR, CFITSIO_PFLAG,
dnl                  ILP64_FLAG, [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$CFITSIO_LIBS $LIBS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if CFITSIO
dnl is found (HAVE_CFITSIO is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_CFITSIO], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl --------------------

CFITSIO_ERROR=""
if test x$2 = x; then
  [acx_cfitsio_incdir="/"]
  AC_CHECK_HEADERS(
    [${acx_cfitsio_incdir}fitsio.h],,
    [
      [acx_cfitsio_incdir=""]
      AC_CHECK_HEADER(
        [fitsio.h],,
        [CFITSIO_ERROR="CFITSIO header files not found!"]
      )
    ]
  )
else
  acx_cfitsio_incdir="$2/"
  AC_CHECK_HEADER(
    [${acx_cfitsio_incdir}fitsio.h],,
    [
      [acx_cfitsio_incdir="$2/include/"]
      AC_CHECK_HEADERS(
        [${acx_cfitsio_incdir}fitsio.h],,
        [CFITSIO_ERROR="CFITSIO header files not found in "$2"!"]
    )]
  )
fi

if test "x$CFITSIO_ERROR" = "x"; then
  AC_DEFINE_UNQUOTED(FITSIO_H, "${acx_cfitsio_incdir}fitsio.h", [CFITSIO header filename.])

dnl ----------------------------
dnl Search CFITSIO library file
dnl ----------------------------

  OLIBS="$LIBS"
  LIBS=""
  if test x$4 = xyes; then
    acx_cfitsio_suffix="64"
    CFITSIO_CFLAGS="-DCFITSIO_USE64BITINT -DLAPACK_ILP64"
  else
    acx_cfitsio_suffix=""
    CFITSIO_CFLAGS=""
  fi
  if test x$1 = x; then
    acx_cfitsio_libopt=""
  else
    acx_cfitsio_libopt="-L$1"
  fi
  AC_SEARCH_LIBS(
      [ffmahd], ["cfitsio"$acx_cfitsio_suffix],,
      [CFITSIO_ERROR="CFITSIO"$acx_cfitsio_suffix" library file not found!"],
      $acx_cfitsio_libopt
    )
  LIBS="$OLIBS"
fi

dnl -------------------------------------------------------------------------
dnl Finally execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND
dnl -------------------------------------------------------------------------

if test "x$CFITSIO_ERROR" = "x"; then
  AC_DEFINE(HAVE_CFITSIO,1, [Define if you have the CFITSIO library and header files.])
  CFITSIO_LIBS="$acx_cfitsio_libopt $ac_cv_search_ffmahd"
  AC_SUBST(CFITSIO_CFLAGS)
  AC_SUBST(CFITSIO_LDFLAGS, "")
  AC_SUBST(CFITSIO_LIBS)
  AC_SUBST(CFITSIO_WARN)
  $5
else
  AC_SUBST(CFITSIO_ERROR)
  $6
fi

])dnl ACX_CFITSIO

