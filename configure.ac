#
#				configure.ac
#
# Process this file with autoconf to generate a configure script.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SExtractor
#
#	Copyright:		(C) 2002-2018 IAP/CNRS/UPMC
#
#	License:		GNU General Public License
#
#	SExtractor is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SExtractor is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SExtractor. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		08/02/2018
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT(SExtractor, 2.25.0,[astromatic@astromatic.net],
	sextractor, [http://astromatic.net/software/sextractor])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_SRCDIR(src/makeit.c)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE
date=`date +%Y-%m-%d`
date1=`date -R`
date2=`date +"%a %b %d %Y"`
date3=`date +"%B %Y"`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])
AC_SUBST(PACKAGER, "AstrOmatic")
AC_SUBST(DATE1, "$date1")
AC_SUBST(DATE2, "$date2")
AC_SUBST(DATE3, "$date3")

# Include macros
sinclude(acx_atlas.m4)
sinclude(acx_fftw.m4)
sinclude(acx_mkl.m4)
sinclude(acx_prog_cc_optim.m4)
sinclude(acx_pthread.m4)
sinclude(acx_urbi_resolve_dir.m4)

# Display pakage and version number
AC_MSG_RESULT([******** Configuring:  $PACKAGE_NAME $VERSION ($date) ********])

# Provide special option for choosing automatically the compilation flags.
AC_MSG_CHECKING([if compilation flags are set automatically])
AC_ARG_ENABLE(auto-flags,
	[AS_HELP_STRING([--enable-auto-flags],
	[Let the configure script choose the compilation flags (default = no)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special option for the Linux Intel C compiler
AC_MSG_CHECKING([whether the INTEL compiler is enabled])
AC_ARG_ENABLE(icc,
	[AS_HELP_STRING([--enable-icc],
	[Use the Intel compiler (default = no)])],
        CC="icc"
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special options for INTEL MKL
# We force the use of icc
AC_MSG_CHECKING([whether INTEL's MKL is enabled])
AC_ARG_ENABLE(mkl,
	[AS_HELP_STRING([--enable-mkl],
	[Use INTEL's MKL for solvers and FFTs (default = no)])],
        enable_icc="yes"
        CC="icc"
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Checks for programs.
AC_LANG(C)

AC_SEARCH_LIBS([strerror],[cposix])
if test x$enable_auto_flags = xyes; then
  CFLAGS=""
  LDFLAGS=""
  ACX_PROG_CC_OPTIM
fi
LT_INIT
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB(m, sin)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stdlib.h string.h sys/mman.h \
		sys/types.h unistd.h])
# Checks for INTEL math header files.
if test "$enable_icc" = "yes"; then
  AC_CHECK_HEADERS(mathimf.h)
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_STRUCT_TM
AC_TYPE_UID_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MMAP
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([atexit getenv gettimeofday isinf isnan logf memcpy memmove \
	memset mkdir munmap posix_memalign setlinebuf sincosf strstr sysconf])

# Check support for large files
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Provide special options for ATLAS
AC_ARG_WITH(atlas-libdir,
	[AS_HELP_STRING([--with-atlas-libdir=<ATLAS library path>],
	[Provide an alternative path to the ATLAS library])])
AC_ARG_WITH(atlas-incdir,
	[AS_HELP_STRING([--with-atlas-incdir=<ATLAS include dir>],
	[Provide an alternative path to the ATLAS include directory])])

# Provide special options for FFTW
AC_ARG_WITH(fftw-libdir,
	[AS_HELP_STRING([--with-fftw-libdir=<FFTW library path>],
	[Provide an alternative path to the FFTW library])])
AC_ARG_WITH(fftw-incdir,
	[AS_HELP_STRING([--with-fftw-incdir=<FFTW include dir>],
	[Provide an alternative path to the FFTW include directory])])

# Provide special options for the MKL library
AC_ARG_WITH(mkl-dir,
	[AS_HELP_STRING([--with-mkl-dir=<MKL path>],
	[Provide an alternative path to the MKL library])])

# Provide special options for OpenBLAS
AC_MSG_CHECKING([whether OpenBLAS is enabled])
AC_ARG_ENABLE(openblas,
	[AS_HELP_STRING([--enable-openblas],
	[Use the OpenBLAS library instead of ATLAS (default = no)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))
AC_ARG_WITH(openblas-libdir,
	[AS_HELP_STRING([--with-openblas-libdir=<OpenBLAS library path>],
	[Provide an alternative path to the OpenBLAS library])])
AC_ARG_WITH(openblas-incdir,
	[AS_HELP_STRING([--with-openblas-incdir=<OpenBLAS header dir>],
	[Provide an alternative path to the OpenBLAS header directory])])

# Provide a special option for the default XSLT URL
with_xsl_url="file://"$(URBI_RESOLVE_DIR([$datadir]))"/$PACKAGE_TARNAME/$PACKAGE_TARNAME.xsl"
AC_ARG_WITH(xsl_url,
	[AS_HELP_STRING([--with-xsl_url=<default URL for XSLT filter>],
	[Provide an alternative default URL of the XSLT filter])])

AC_DEFINE_UNQUOTED([XSL_URL], "$with_xsl_url",[Default URL of the XSLT filter])

# Provide special option to disable model-fitting (enabled by default)
AC_MSG_CHECKING([if model-fitting should be disabled (default=enabled)])
AC_ARG_ENABLE([model-fitting],
	[AS_HELP_STRING([--disable-model-fitting],
	[Disable model-fitting and library dependencies])],
	if test "$enable_model_fitting" = "no"; then
	  AC_MSG_RESULT([yes])
	else
	  AC_MSG_RESULT([no])
	fi,
	AC_MSG_RESULT([no]))

# Set flags for multithreading
n_pthreads=1024
AC_ARG_ENABLE(threads,
	[AS_HELP_STRING([--enable-threads@<:@=<max_number_of_threads>@:>@],
	[Enable multhreading (on with up to 1024 threads by default)])],
    if test "$enableval" = "no"; then
      use_pthreads="no"
    else
      use_pthreads="yes"
      if test "$enableval" != "yes"; then
        n_pthreads=$enableval
      fi
    fi,
    use_pthreads=yes
    )
# Deactivate multithreading for now
use_pthreads="no"

# Provide special option for profiling
AC_MSG_CHECKING([for profiler mode])
AC_ARG_ENABLE(profiling,
	[AS_HELP_STRING([--enable-profiling],
	[Enable special mode for profiling (default = no)])]
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Enable linking options for making the executable as portable as possible.
AC_MSG_CHECKING([best linking option])
AC_ARG_ENABLE(best-link,
	[AS_HELP_STRING([--enable-best-link],
	[Choose the right combination of static and dynamic linking to make \
the executable as portable as possible (default = no)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))


################# Actions to complete in case of multhreading ################
AC_DEFINE_UNQUOTED(THREADS_NMAX, $n_pthreads,[Maximum number of POSIX threads])
if test "$use_pthreads" = "yes"; then
  AC_MSG_CHECKING([for multithreading])
  AC_MSG_RESULT([maximum of $n_pthreads thread(s)])
  AC_DEFINE(USE_THREADS, 1, [Triggers multhreading])
# CC, CFLAGS and LIBS are system and compiler-dependent
  ACX_PTHREAD
  CC="$PTHREAD_CC"
  [AM_CFLAGS="$AM_CFLAGS $PTHREAD_CFLAGS -D_REENTRANT"]
  LIBS="$LIBS $PTHREAD_LIBS"
fi
AM_CONDITIONAL(USE_THREADS, test $use_pthreads = "yes")

############# Actions to complete if model-fitting is activated ##############
if test "$enable_model_fitting" != "no"; then
  AC_DEFINE(USE_MODEL, 1, [Triggers model-fitting])
############ handle the INTEL MKL library (FFTW + LAPACKe) ###########
  if test "$enable_mkl" = "yes"; then
    convlibs="${srcdir}/wcs/libwcs_c.a,${srcdir}/levmar/liblevmar.a"
    ACX_MKL($with_mkl_dir,,$enable_best_link,$convlibs)
    if test "$MKL_WARN" == ""; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      AC_MSG_WARN([$MKL_WARN])
    fi
    AM_CFLAGS="$AM_CFLAGS $MKL_CFLAGS "
    AM_LDFLAGS="$AM_LDFLAGS $MKL_LDFLAGS "
    LIBS="$LIBS $MKL_LIBS"
  else
################ handle the FFTW library (Fourier transforms) ################
    ACX_FFTW($with_fftw_libdir, $with_fftw_incdir, $use_pthreads,yes,
      [
        [LIBS="$FFTW_LIBS $LIBS"]
        if test "$FFTW_WARN" != ""; then
          AC_MSG_WARN([$FFTW_WARN])
        fi
      ],
      AC_MSG_ERROR([$FFTW_ERROR Exiting.])
    )

    if test "x$enable_openblas" = "xyes"; then
######## Handle the OpenBLAS library (linear algebra: BLAS + LAPACKe) ########
      ACX_OPENBLAS($with_openblas_libdir, $with_openblas_incdir, $use_pthreads, no,
        [
          AM_CFLAGS="$AM_CFLAGS $OPENBLAS_CFLAGS "
          AM_LDFLAGS="$AM_LDFLAGS $OPENBLAS_LDFLAGS "
          LIBS="$OPENBLAS_LIBS $LIBS"
          if test "$OPENBLAS_WARN" != ""; then
            AC_MSG_WARN([$OPENBLAS_WARN])
          fi
        ],
        AC_MSG_ERROR([$OPENBLAS_ERROR Exiting.])
      )
    else
######### handle the ATLAS library (linear algebra: BLAS + cLAPACK) ##########
      ACX_ATLAS($with_atlas_libdir, $with_atlas_incdir, $use_pthreads,
        [
          [LIBS="$ATLAS_LIBS $LIBS"]
          if test "$ATLAS_WARN" != ""; then
            AC_MSG_WARN([$ATLAS_WARN])
	  fi
        ],
        AC_MSG_ERROR([$ATLAS_ERROR Exiting.])
      )
    fi
  fi
fi

AM_CONDITIONAL(USE_MODEL, [test "$enable_model_fitting" != "no"])

# Compile with profiling option
if test "$enable_profiling" = "yes"; then
  if test "$enable_icc" = "yes"; then
    AM_CFLAGS="$AM_CFLAGS -pq"
  else
    AM_CFLAGS="$AM_CFLAGS -pg"
  fi
  enable_best_link="no"
fi

# "Best" linking option
if test "$enable_best_link" = "yes"; then
  if test "$enable_icc" = "yes"; then
    AM_LDFLAGS="-static-intel -qopenmp-link static -shared-libgcc \
	-static-libtool-libs -no-intel-extensions -avoid-version $AM_LDFLAGS"
  else
    AM_LDFLAGS="-shared-libgcc -static-libtool-libs -avoid-version $AM_LDFLAGS"
  fi
fi

AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)

# Display compiler and linker flags
AC_MSG_RESULT([***************************************************************])
AC_MSG_RESULT([Compile cmdline:  $CC $AM_CPPFLAGS $CPPFLAGS $AM_CFLAGS $CFLAGS])
AC_MSG_RESULT([Link    cmdline:  $CC $AM_LDFLAGS $LDFLAGS $LIBS])
AC_MSG_RESULT([Default XSLT URL: $xsl_url])
AC_MSG_RESULT([***************************************************************])

AC_CONFIG_FILES([
Makefile
sextractor.spec
debian/Makefile
debian/changelog
doc/Makefile
doc/src/conf.py
man/Makefile
man/sex.1
src/Makefile
src/fits/Makefile
src/levmar/Makefile
src/wcs/Makefile
tests/Makefile
])
AC_OUTPUT
