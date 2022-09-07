dnl
dnl				acx_cfitsio.m4
dnl
dnl Figure out if the CFITSIO library and header files are installed.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
dnl	Last modified:		27/02/2013
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_CFITSIO([CFITSIO_DIR, CFITSIO_INCDIR,
dnl                    [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the CFITSIO libraries and header
dnl files are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$CFITSIO_LIBS $LIBS"
dnl
dnl You may wish to use these variables in your default CFLAGS:
dnl
dnl        CFLAGS="$CFLAGS $CFITSIO_CFLAGS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if CFITSIO
dnl is found (HAVE_CFITSIO are defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_CFITSIO], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl --------------------

acx_cfitsio_ok=no
if test x$2 = x; then
  AC_CHECK_HEADER(fitsio.h,[acx_cfitsio_ok=yes])
  if test x$acx_cfitsio_ok = xyes; then
    AC_DEFINE(CFITSIO_H, "fitsio.h", [CFITSIO header filename.])
  else
    AC_CHECK_HEADER(cfitsio/fitsio.h,[acx_cfitsio_ok=yes])
    if test x$acx_cfitsio_ok = xyes; then
      AC_DEFINE(CFITSIO_H, "cfitsio/fitsio.h", [CFITSIO header filename.])
    else
      CFITSIO_ERROR="CFITSIO include files not found at default location!"
    fi
  fi
else
  AC_CHECK_HEADER($2/fitsio.h,[acx_cfitsio_ok=yes])
  if test x$acx_cfitsio_ok = xyes; then
    AC_DEFINE_UNQUOTED(CFITSIO_H, "$2/fitsio.h", [CFITSIO header filename.])
  else
    CFITSIO_ERROR="CFITSIO include files not found in $2!"
  fi
fi

dnl --------------------
dnl Search library files
dnl --------------------

CFITSIO_LIBS=""
OLIBS="$LIBS"
LIBS=""

if test x$acx_cfitsio_ok = xyes; then
  if test x$1 = x; then
    AC_CHECK_LIB(cfitsio, ffopen, [acx_cfitsio_ok=yes], [acx_cfitsio_ok=no])
    if test x$acx_cfitsio_ok = xyes; then
      AC_DEFINE(HAVE_CFITSIO,1, [Define if you have the CFITSIO libraries and header files.])
      CFITSIO_LIBS="-lcfitsio"
    else
      CFITSIO_ERROR="CFITSIO library files not found at usual locations!"
    fi
  else
dnl -------------------------
dnl Specific libdir specified
dnl -------------------------
    AC_CHECK_LIB(cfitsio, ffopen, [acx_cfitsio_ok=yes], [acx_cfitsio_ok=no], [-L$1])
    if test x$acx_cfitsio_ok = xyes; then
      AC_DEFINE(HAVE_CFITSIO,1, [Define if you have the CFITSIO libraries and header files.])
        CFITSIO_LIBS="-L$1 -lcfitsio"
    else
      CFITSIO_ERROR="CFITSIO library files not found in $1!"
    fi
  fi
fi

LIBS="$OLIBS"
if test x$acx_cfitsio_ok = xyes; then
  AC_SUBST(CFITSIO_LIBS)
  $3
else
  AC_SUBST(CFITSIO_ERROR)
  $4
fi

])dnl ACX_CFITSIO
