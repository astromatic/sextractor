dnl
dnl				acx_mkl.m4
dnl
dnl Set up options for using the INTEL MKL library.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2003-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
dnl	Last modified:		09/07/2012
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_MKL()
dnl
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$MKL_LIBS $LIBS"
dnl
dnl You may wish to use these variables in your default CFLAGS:
dnl
dnl        CFLAGS="$CFLAGS $MKL_CFLAGS"
dnl

AC_DEFUN([ACX_MKL], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl ----------------------
dnl Set architecture flags
dnl ----------------------

dnl Try to find INTEL architecture (Intel 64 or ia32)
if icc -V 2>&1 | grep -i "Intel(R) 64" > /dev/null 2>&1; then
  AC_SUBST(MKL_CFLAGS, "-DMKL_ILP64")
  AC_SUBST(MKL_LIBS, "-mkl")
elif icc -V 2>&1 | grep -i "Intel(R)" > /dev/null 2>&1; then
  AC_SUBST(MKL_CFLAGS, "")
  AC_SUBST(MKL_LIBS, "-mkl")
else
  AC_SUBST(MKL_CFLAGS, "")
  AC_SUBST(MKL_LIBS, "")
  MKL_WARN="INTEL compiler not detected"
  AC_SUBST(MKL_WARN)
fi

dnl --------------------
dnl Set internal flags
dnl --------------------

AC_DEFINE(HAVE_FFTW,1, [Define if you have the FFTW libraries.])
AC_DEFINE(HAVE_LAPACK,1, [Define if you have the LAPACK libraries.])
AC_DEFINE(HAVE_LAPACKE,1, [Define if you have the LAPACKe libraries.])

dnl --------------------
dnl Set include files
dnl --------------------

AC_DEFINE_UNQUOTED(FFTW_H, "fftw/fftw3_mkl.h", [FFTW header filename.])
AC_DEFINE_UNQUOTED(LAPACK_H, "mkl_lapack.h", [LAPACK header filename.])
AC_DEFINE_UNQUOTED(LAPACKE_H, "mkl_lapacke.h", [LAPACKe header filename.])


])dnl ACX_MKL
