dnl
dnl				acx_mkl.m4
dnl
dnl Set up options for using the INTEL MKL library.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2003-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
dnl	Last modified:		17/04/2013
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_MKL([MKL_DIR, ILP64_FLAG, STATIC_FLAG, CONV_LIBS])
dnl
dnl This macro sets the MKL_CFLAGS, MKL_LDFLAGS and MKL_LIBS variables to
dnl for compiling and linking with INTEL's MKL. A coma-separated list of
dnl convenience libraries may be included in the linked group for static linking.
dnl You may wish to use these variables in your default CFLAGS:
dnl
dnl        CFLAGS="$CFLAGS $MKL_CFLAGS"
dnl
dnl You may wish to use these variables in your default LDFLAGS:
dnl
dnl        LDFLAGS="$LDFLAGS $MKL_LDLAGS"
dnl
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$LIBS $MKL_LIBS"
dnl

AC_DEFUN([ACX_MKL], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl ------------------------
dnl Set MKL's root directory
dnl ------------------------

if test x$1 = x; then
  mklroot=${MKLROOT}
else
  mklroot=$1
fi

dnl -----------------------------
dnl Include convenience libraries
dnl -----------------------------

if test x$4 = x; then
  startgroup="-Wl,--start-group"
else
  startgroup="-Wl,--start-group,$4"
fi

dnl ----------------------
dnl Set architecture flags
dnl ----------------------

dnl check if INTEL compiler is present
icc -V 2>&1 | grep -i "Intel" > /dev/null 2>&1 && flagicc=yes
dnl check if INTEL compiler uses x86_64 architecture
icc -V 2>&1 | grep -i "Intel(R) 64" > /dev/null 2>&1 && flag64=yes
dnl check if the platform is OSX
icc -dumpmachine 2>&1 | grep -i "darwin" > /dev/null 2>&1 && flagosx=yes

dnl ----------------------
dnl Exit if INTEL compiler is not found
dnl ----------------------
if test x$flagicc = x; then
  AC_SUBST(MKL_CFLAGS, "")
  AC_SUBST(MKL_LDFLAGS, "")
  AC_SUBST(MKL_LIBS, "")
  MKL_WARN="INTEL compiler not detected"
  AC_SUBST(MKL_WARN)
  exit
fi

if test x$flagosx = xyes; then
dnl MacOSX
  if test x$flag64 = xyes; then
dnl INTEL compiler uses Intel64 architecture
    if test x$2 = xyes; then
dnl 64 bit pointers
      AC_SUBST(MKL_CFLAGS, "-openmp -DMKL_ILP64 -I$mklroot/include")
      if test x$3 = xyes; then
dnl Static linking
        AC_SUBST(MKL_LIBS, ["$mklroot/lib/libmkl_intel_ilp64.a \
		$mklroot/lib/libmkl_intel_thread.a \
		$mklroot/lib/libmkl_core.a -lpthread -lm"])
      else
dnl Dynamic linking
        AC_SUBST(MKL_LIBS, "-L$mklroot/lib  -lmkl_intel_ilp64 \
	-lmkl_intel_thread -lmkl_core -lpthread -lm")
      fi
    else
dnl 32 bit pointers
      AC_SUBST(MKL_CFLAGS, "-openmp -I$mklroot/include")
      if test x$3 = xyes; then
dnl Static linking
        AC_SUBST(MKL_LIBS, ["$mklroot/lib/libmkl_intel_lp64.a \
		$mklroot/lib/libmkl_intel_thread.a \
		$mklroot/lib/libmkl_core.a -lpthread -lm"])
      else
dnl Dynamic linking
        AC_SUBST(MKL_LIBS, "-L$mklroot/lib -lmkl_intel_lp64 \
		-lmkl_intel_thread -lmkl_core -lpthread -lm")
      fi
    fi
  else
dnl INTEL compiler uses IA32 architecture
    AC_SUBST(MKL_CFLAGS, "-openmp -I$mklroot/include")
    if test x$3 = xyes; then
dnl Static linking
    AC_SUBST(MKL_LIBS, ["$mklroot/lib/libmkl_intel.a \
	$mklroot/lib/libmkl_intel_thread.a \
	$mklroot/lib/libmkl_core.a -lpthread -lm"])
    else
dnl Dynamic linking
      AC_SUBST(MKL_LIBS, "-L$mklroot/lib -lmkl_intel -lmkl_intel_thread \
	-lmkl_core -lpthread -lm")
    fi
  fi
else
dnl Linux
  if test x$flag64 = xyes; then
dnl INTEL compiler uses Intel64 architecture
    if test x$2 = xyes; then
dnl 64 bit pointers
      AC_SUBST(MKL_CFLAGS, "-openmp -DMKL_ILP64 -I$mklroot/include")
      if test x$3 = xyes; then
dnl Static linking
      AC_SUBST(MKL_LIBS,
	["$startgroup,$mklroot/lib/intel64/libmkl_intel_ilp64.a,\
$mklroot/lib/intel64/libmkl_intel_thread.a,\
$mklroot/lib/intel64/libmkl_core.a,-end-group -lpthread -lm"])
      else
dnl Dynamic linking
        AC_SUBST(MKL_LIBS, "-L$mklroot/lib/intel64  -lmkl_intel_ilp64 \
		-lmkl_intel_thread -lmkl_core -lpthread -lm")
      fi
    else
dnl 32 bit pointers
      AC_SUBST(MKL_CFLAGS, "-openmp -I$mklroot/include")
      if test x$3 = xyes; then
dnl Static linking
        AC_SUBST(MKL_LIBS,
		["$startgroup,$mklroot/lib/intel64/libmkl_intel_lp64.a,\
$mklroot/lib/intel64/libmkl_intel_thread.a,\
$mklroot/lib/intel64/libmkl_core.a,--end-group -lpthread -lm"])
      else
dnl Dynamic linking
        AC_SUBST(MKL_LIBS, "-L$mklroot/lib/intel64 -lmkl_intel_lp64 \
	-lmkl_intel_thread -lmkl_core -lpthread -lm")
      fi
    fi
  else
dnl INTEL compiler uses IA32 architecture
    AC_SUBST(MKL_CFLAGS, "-openmp -I$mklroot/include")
    if test x$3 = xyes; then
dnl Static linking
      AC_SUBST(MKL_LIBS, ["$startgroup,$mklroot/lib/ia32/libmkl_intel.a,\
$mklroot/lib/ia32/libmkl_intel_thread.a,\
$mklroot/lib/ia32/libmkl_core.a,--end-group -lpthread -lm"])
    else
dnl Dynamic linking
      AC_SUBST(MKL_LIBS, "-L$mklroot/lib/ia32 -lmkl_intel -lmkl_intel_thread \
	-lmkl_core -lpthread -lm")
    fi
  fi
fi

AC_SUBST(MKL_LDFLAGS, "")

dnl --------------------
dnl Set internal flags
dnl --------------------

AC_DEFINE(HAVE_MKL,1, [Define if you have the MKL libraries.])
AC_DEFINE(HAVE_FFTW,1, [Define if you have the FFTW libraries.])
AC_DEFINE(HAVE_LAPACK,1, [Define if you have the LAPACK libraries.])
AC_DEFINE(HAVE_LAPACKE,1, [Define if you have the LAPACKe libraries.])

dnl --------------------
dnl Set include files
dnl --------------------

AC_DEFINE_UNQUOTED(MKL_H, "mkl.h", [MKL header filename.])
AC_DEFINE_UNQUOTED(FFTW_H, "fftw/fftw3_mkl.h", [FFTW header filename.])
AC_DEFINE_UNQUOTED(LAPACK_H, "mkl_lapack.h", [LAPACK header filename.])
AC_DEFINE_UNQUOTED(LAPACKE_H, "mkl_lapacke.h", [LAPACKe header filename.])

])dnl ACX_MKL
