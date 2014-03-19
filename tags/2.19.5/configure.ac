#
#				configure.ac
#
# Process this file with autoconf to generate a configure script.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SExtractor
#
#	Copyright:		(C) 2002-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
#	Last modified:		14/12/2013
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT(sextractor, 2.19.5, [bertin@iap.fr])
AC_CONFIG_SRCDIR(src/makeit.c)
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE
date=`date +%Y-%m-%d`
date2=`date +"%a %b %d %Y"`
date3=`date +"%B %Y"`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])
AC_SUBST(PACKAGER, "Emmanuel Bertin")
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
AC_MSG_RESULT([*********** Configuring:  $PACKAGE $VERSION ($date) **********])

# Provide special option for choosing automatically the compilation flags.
AC_MSG_CHECKING([if compilation flags should be set automatically])
AC_ARG_ENABLE(auto-flags,
	[AS_HELP_STRING([--enable-auto-flags],
	[Let the configure script choose the compilation flags \
(off by default)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special option for the Linux Intel C compiler
AC_MSG_CHECKING([whether we should use the INTEL compiler])
AC_ARG_ENABLE(icc,
	[AS_HELP_STRING([--enable-icc],
	[Enable special mode for compilation with the Intel compiler \
(off by default)])],
        CC="icc"
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special options for INTEL MKL
# We force the use of icc
AC_MSG_CHECKING([whether we should use INTEL's MKL])
AC_ARG_ENABLE(mkl,
	[AS_HELP_STRING([--enable-mkl],
	[Use INTEL's MKL for solvers and FFTs (off by default)])],
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
if test x$enable_icc = xyes; then
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

# Provide a special option for the default XSLT URL
with_xsl_url="file://"$(URBI_RESOLVE_DIR([$datadir]))"/$PACKAGE_NAME/$PACKAGE_NAME.xsl"
AC_ARG_WITH(xsl_url,
	[AS_HELP_STRING([--with-xsl_url=<default URL for XSLT filter>],
	[Provide an alternative default URL of the XSLT filter])])

AC_DEFINE_UNQUOTED([XSL_URL], "$with_xsl_url",[Default URL of the XSLT filter])

# Provide special option to disable model-fitting (enabled by default)
AC_MSG_CHECKING([if model-fitting should be disabled])
AC_ARG_ENABLE([model-fitting],
	[AS_HELP_STRING([--disable-model-fitting],
	[Disable model-fitting and library dependencies])],
	if test x$enable_model_fitting = xno; then
	  AC_MSG_RESULT([yes])
	else
	  AC_MSG_RESULT([no])
	fi,
	AC_MSG_RESULT([no]))

# Set flags for multithreading
n_pthreads=16
AC_ARG_ENABLE(threads,
        [AS_HELP_STRING([--enable-threads@<:@=<max_number_of_threads>@:>@],
	[Enable multhreading (on with up to 16 threads by default)])],
    if test x$enableval = xno; then
      use_pthreads="no"
    else
      use_pthreads="yes"
      if test x$enableval != xyes; then
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
	[Enable special mode for profiling (off by default)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Enable linking options for making the executable as portable as possible.
AC_MSG_CHECKING([best linking option])
AC_ARG_ENABLE(best-link,
	[AS_HELP_STRING([--enable-best-link],
	[Choose the right combination of static and dynamic linking to make \
the executable as portable as possible (off by default)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))


################# Actions to complete in case of multhreading ################
AC_DEFINE_UNQUOTED(THREADS_NMAX, $n_pthreads,[Maximum number of POSIX threads])
if test x$use_pthreads = xyes; then
  AC_MSG_CHECKING([for multithreading])
  AC_MSG_RESULT([maximum of $n_pthreads thread(s)])
  AC_DEFINE(USE_THREADS, 1, [Triggers multhreading])
# CC, CFLAGS and LIBS are system and compiler-dependent
  ACX_PTHREAD
  CC="$PTHREAD_CC"
  [AM_CFLAGS="$AM_CFLAGS $PTHREAD_CFLAGS -D_REENTRANT"]
  LIBS="$PTHREAD_LIBS $LIBS"
fi
AM_CONDITIONAL(USE_THREADS, [test x$use_pthreads = xyes])

############# Actions to complete if model-fitting is activated ##############
if test x$enable_model_fitting != xno; then
  AC_DEFINE(USE_MODEL, 1, [Triggers model-fitting])
############ handle the INTEL MKL library (FFTW + LAPACK) ###########
  if test x$enable_mkl = xyes; then
    convlibs="${srcdir}/wcs/libwcs_c.a,${srcdir}/levmar/liblevmar.a"
    ACX_MKL($with_mkl_dir,,$enable_best_link,$convlibs)
    if test x$MKL_WARN != x; then
      AC_MSG_WARN([$MKL_WARN])
    fi
    AM_CFLAGS="$AM_CFLAGS $MKL_CFLAGS "
    AM_LDFLAGS="$AM_LDFLAGS $MKL_LDFLAGS"
    LIBS="$LIBS $MKL_LIBS"
  else
################ handle the FFTW library (Fourier transforms) ################
    ACX_FFTW($with_fftw_libdir,$with_fftw_incdir,$use_pthreads,yes,
	LIBS="$FFTW_LIBS $LIBS"
	if test x$FFTW_WARN != x; then
	  AC_MSG_WARN([$FFTW_WARN])
	fi,
	AC_MSG_ERROR([$FFTW_ERROR Exiting.]))
################## handle the ATLAS library(linear algebra) ##################
    ACX_ATLAS($with_atlas_libdir,$with_atlas_incdir,$use_pthreads,
	LIBS="$ATLAS_LIBS $LIBS"
	if test x$ATLAS_WARN != x; then
	  AC_MSG_WARN([$ATLAS_WARN])
	fi,
	AC_MSG_ERROR([$ATLAS_ERROR Exiting.]))
  fi
fi
AM_CONDITIONAL(USE_MODEL, [test x$enable_model_fitting != xno])

# Compile with profiling option
if test x$enable_profiling = xyes; then
  if test x$enable_icc = xyes; then
    AM_CFLAGS="$AM_CFLAGS -pq"
  else
    AM_CFLAGS="$AM_CFLAGS -pg"
  fi
  enable_best_link="no"
fi

# "Best" linking option
if test x$enable_best_link = xyes; then
  if test x$enable_icc = xyes; then
    AM_LDFLAGS="-static-intel -openmp-link static -shared-libgcc \
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
src/Makefile
src/fits/Makefile
src/levmar/Makefile
src/wcs/Makefile
man/Makefile
tests/Makefile
sextractor.spec
man/sex.1])
AC_OUTPUT
