dnl @synopsis ACX_FFTW([FFTW_DIR, FFTW_INCDIR, FFTW_PFLAG, FFTW_FFLAG, [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the FFTW3 libraries and header
dnl files are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$FFTW_LIBS $LIBS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if FFTW
dnl is found (HAVE_FFTWx are defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.
dnl
dnl @version $Id: acx_fftw.m4,v 1.0 2007/10/19 21:30:17 bertin Exp $
dnl @author Emmanuel Bertin <bertin@iap.fr>

AC_DEFUN([ACX_FFTW], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl --------------------

acx_fftw_ok=no
if test x$2 = x; then
  if test x$1 = x; then
    AC_CHECK_HEADER(fftw3.h,[acx_fftw_ok=yes])
    if test x$acx_fftw_ok = xyes; then
      AC_DEFINE(FFTW_H, "fftw.h", [FFTW header filename.])
    else
      AC_CHECK_HEADER(fftw/fftw3.h,[acx_fftw_ok=yes])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(FFTW_H, "fftw/fftw3.h", [FFTW header filename.])
      else
        FFTW_ERROR="FFTW include files not found in default location!"
      fi
    fi
  else
    AC_CHECK_HEADER($1/include/fftw3.h,[acx_fftw_ok=yes])
    if test x$acx_fftw_ok = xyes; then
      AC_DEFINE(FFTW_H, "$1/include/fftw3.h", [FFTW header filename.])
    else
      AC_CHECK_HEADER(fftw3.h,[acx_fftw_ok=yes])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(FFTW_H, "fftw.h", [FFTW header filename.])
      else
        FFTW_ERROR="FFTW include files not found in $1/include!"
      fi
    fi
  fi
else
  AC_CHECK_HEADER($2/fftw3.h,[acx_fftw_ok=yes])
  if test x$acx_fftw_ok = xyes; then
    AC_DEFINE(FFTW_H, "$2/fftw3.h", [FFTW header filename.])
  else
    FFTW_ERROR="FFTW include files not found in $2!"
  fi
fi

dnl --------------------
dnl Search library files
dnl --------------------

FFTW_LIBS=""
OLIBS="$LIBS"
LIBS=""

if test x$acx_fftw_ok = xyes; then
  if test x$1 = x; then
    if test x$4 = xyes; then
      AC_CHECK_LIB(fftw3f, fftwf_execute, [acx_fftw_ok=yes],
		[acx_fftw_ok=no], [-lm])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(HAVE_FFTWF,1,
    [Define if you have the FFTW single precision libraries and header files.])
        FFTW_LIBS="-lfftw3f"
      else
        FFTW_ERROR="FFTW single precision library files not found at usual locations!"
      fi
    else
      AC_CHECK_LIB(fftw3, fftw_execute, [acx_fftw_ok=yes],
		[acx_fftw_ok=no], [-lm])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(HAVE_FFTW,1,
    [Define if you have the FFTW double precision libraries and header files.])
        FFTW_LIBS="-lfftw3"
      else
        FFTW_ERROR="FFTW double precision library files not found at usual locations!"
      fi
    fi
    if test x$acx_fftw_ok = xyes && test x$3 = xyes; then
      if test x$4 = xyes; then
        AC_CHECK_LIB(fftw3f_threads, fftwf_init_threads,
		[acx_fftw_ok=yes], [acx_fftw_ok=no], [-lfftw3f -lm])
        if test x$acx_fftw_ok = xyes; then
          AC_DEFINE(HAVE_FFTWFT,1,
    [Define if you have the FFTW single precision multithreaded libraries and header files.])
          FFTW_LIBS="-lfftw3f_threads -lfftw3f"
        else
          FFTW_ERROR="FFTW single precision library was compiled without multithreading support!"
        fi
      else
        AC_CHECK_LIB(fftw3_threads, fftw_init_threads,
		[acx_fftw_ok=yes], [acx_fftw_ok=no], [-lfftw3 -lm])
        if test x$acx_fftw_ok = xyes; then
          AC_DEFINE(HAVE_FFTWT,1,
    [Define if you have the FFTW double precision multithreaded libraries and header files.])
          FFTW_LIBS="-lfftw3_threads -lfftw3"
        else
          FFTW_ERROR="FFTW double precision library was compiled without multithreading support!"
        fi
      fi
    fi
  else
dnl -------------------------
dnl Specific libdir specified
dnl -------------------------
    if test x$4 = xyes; then
      AC_CHECK_LIB(fftw3f, fftwf_execute, [acx_fftwf_ok=yes],
		[acx_fftw_ok=no], [-lm])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(HAVE_FFTWF,1,
    [Define if you have the FFTW single precision libraries and header files.])
        FFTW_LIBS="-L$1 -lfftw3f"
      else
        FFTW_ERROR="FFTW single precision library files not found in $1!"
      fi
    else
      AC_CHECK_LIB(fftw3, fftw_execute, [acx_fftw_ok=yes],
		[acx_fftw_ok=no], [-lm])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(HAVE_FFTW,1,
    [Define if you have the FFTW double precision libraries and header files.])
        FFTW_LIBS="-L$1 -lfftw3"
      else
        FFTW_ERROR="FFTW double precision library files not found in $1!"
      fi
    fi
    if test x$acx_fftw_ok = xyes && test x$3 = xyes; then
      if test x$4 = xyes; then
        AC_CHECK_LIB(fftw3f_threads, fftwf_init_threads,
		[acx_fftw_ok=yes], [acx_fftw_ok=no], [-lfftw3f -lm])
        if test x$acx_fftwft_ok = xyes; then
          AC_DEFINE(HAVE_FFTWFT,1,
    [Define if you have the FFTW single precision multithreaded libraries and header files.])
          FFTW_LIBS="-L$1 -lfftw3f_threads -lfftw3f"
        else
          FFTW_ERROR="FFTW single precision library in $1 was compiled without multithreading support!"
        fi
      else
        AC_CHECK_LIB(fftw3_threads, fftw_init_threads,
		[acx_fftw_ok=yes], [acx_fftw_ok=no], [-lfftw3 -lm])
        if test x$acx_fftw_ok = xyes; then
          AC_DEFINE(HAVE_FFTWT,1,
    [Define if you have the FFTW double precision multithreaded libraries and header files.])
          FFTW_LIBS="-L$1 -lfftw3_threads -lfftw3"
        else
          FFTW_ERROR="FFTW double precision library in $1 was compiled without multithreading support!"
        fi
      fi
    fi
  fi
fi

LIBS="$OLIBS"
if test x$acx_fftw_ok = xyes; then
  AC_SUBST(FFTW_LIBS)
  $5
else
  AC_SUBST(FFTW_ERROR)
  $6
fi

])dnl ACX_FFTW
