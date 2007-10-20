dnl @synopsis ACX_FFTW([FFTW_DIR, FFTW_INCDIR, FFTW_PFLAG, FFTW_FFLAG, [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the FFTW3 libraries and header
dnl files are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$FFTWF_LIBS $LIBS"
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
      AC_DEFINE(ATLAS_FFTW_H, "fftw.h", [FFTW header filename.])
    else
      AC_CHECK_HEADER(fftw/fftw3.h,[acx_fftw_ok=yes])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(ATLAS_FFTW_H, "fftw/fftw3.h", [FFTW header filename.])
      else
        FFTW_ERROR="FFTW include files not found in default location!"
      fi
    fi
  else
    AC_CHECK_HEADER($1/include/fftw3.h,[acx_fftw_ok=yes])
    if test x$acx_fftw_ok = xyes; then
      AC_DEFINE(ATLAS_FFTW_H, "$1/include/fftw3.h", [FFTW header filename.])
    else
      AC_CHECK_HEADER(fftw3.h,[acx_fftw_ok=yes])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(ATLAS_FFTW_H, "fftw.h", [FFTW header filename.])
      else
        FFTW_ERROR="FFTW include files not found in $1/include!"
      fi
    fi
  fi
else
  AC_CHECK_HEADER($2/fftw3.h,[acx_fftw_ok=yes])
  if test x$acx_fftw_ok = xyes; then
    AC_DEFINE(ATLAS_FFTW_H, "$2/fftw3.h", [FFTW header filename.])
  else
    FFTW_ERROR="FFTW include files not found in $2!"
  fi
fi

dnl --------------------
dnl Search library files
dnl --------------------

FFTW_LIBS=""
if test x$acx_fftw_ok = xyes; then
  OLIBS="$LIBS"
  LIBS=""
  if test x$1 = x; then
    if test x$4 = xyes; then
      AC_CHECK_LIB(fftw3f, fftwf_execute, acx_fftwf_ok=yes FFTW_LIBS="-lfftw3f",
		[acx_fftwf_ok=no], [-lm])
      if test x$acx_fftwf_ok = xyes; then
        AC_DEFINE(HAVE_FFTWF,1,
    [Define if you have the FFTW single precision libraries and header files.])
        FFTW_LIBPATH=""
        FFTW_LIBS="-lfftw3f"
        AC_SUBST(FFTW_LIBPATH)
        AC_SUBST(FFTW_LIBS)
        $5
      else
        ATLAS_ERROR="FFTW single precision library files not found at usual locations!"
        $6
      fi
    else
      AC_CHECK_LIB(fftw3, fftw_execute, acx_fftwd_ok=yes FFTW_LIBS="-lfftw3",
		[acx_fftwd_ok=no], [-lm])
      if test x$acx_fftwd_ok = xyes; then
        AC_DEFINE(HAVE_FFTW,1,
    [Define if you have the FFTW double precision libraries and header files.])
        FFTW_LIBPATH=""
        FFTW_LIBS="-lfftw3"
        AC_SUBST(FFTW_LIBPATH)
        AC_SUBST(FFTW_LIBS)
        $5
      else
        ATLAS_ERROR="FFTW double precision library files not found at usual locations!"
        $6
      fi
    fi
    if test x$3 = xyes; then
      if test x$4 = xyes; then
        AC_CHECK_LIB(fftw3f_threads, fftwf_init_threads,
		acx_fftwft_ok=yes FFTW_LIBS="-lfftw3f_threads -lfftw3f",
		[acx_fftwft_ok=no], [-lm])
        if test x$acx_fftwft_ok = xyes; then
          AC_DEFINE(HAVE_FFTWFT,1,
    [Define if you have the FFTW single precision multithreaded libraries and header files.])
          FFTW_LIBPATH=""
          FFTW_LIBS="-lfftw3f_threads -lfftw3f"
          AC_SUBST(FFTW_LIBPATH)
          AC_SUBST(FFTW_LIBS)
          $5
        else
          ATLAS_ERROR="FFTW single precision library was compiled without multithreading support!"
          $6
        fi
      else
        AC_CHECK_LIB(fftw3_threads, fftw_init_threads,
		acx_fftwdt_ok=yes FFTW_LIBS="-lfftw3_threads -lfftw3",
		[acx_fftwdt_ok=no], [-lm])
        if test x$acx_fftwdt_ok = xyes; then
          AC_DEFINE(HAVE_FFTWT,1,
    [Define if you have the FFTW double precision multithreaded libraries and header files.])
          FFTW_LIBPATH=""
          FFTW_LIBS="-lfftw3_threads -lfftw3"
          AC_SUBST(FFTW_LIBPATH)
          AC_SUBST(FFTW_LIBS)
          $5
        else
          ATLAS_ERROR="FFTW double precision library was compiled without multithreading support!"
          $6
        fi
      fi
    fi
  else
dnl -------------------------
dnl Specific libdir specified
dnl -------------------------
    if test x$4 = xyes; then
      AC_CHECK_LIB(fftw3f, fftwf_execute, acx_fftwf_ok=yes FFTW_LIBS="-L$1 -lfftw3f",
		[acx_fftwf_ok=no], [-lm])
      if test x$acx_fftwf_ok = xyes; then
        AC_DEFINE(HAVE_FFTWF,1,
    [Define if you have the FFTW single precision libraries and header files.])
        FFTW_LIBPATH=""
        FFTW_LIBS="-L$1 -lfftw3f"
      else
        ATLAS_ERROR="FFTW single precision library files not found in $1!"
        $6
      fi
    else
      AC_CHECK_LIB(fftw3, fftw_execute, acx_fftwd_ok=yes FFTW_LIBS="-L$1 -lfftw3",
		[acx_fftwd_ok=no], [-lm])
      if test x$acx_fftwd_ok = xyes; then
        AC_DEFINE(HAVE_FFTW,1,
    [Define if you have the FFTW double precision libraries and header files.])
        FFTW_LIBPATH=""
        FFTW_LIBS="-L$1 -lfftw3"
      else
        ATLAS_ERROR="FFTW double precision library files not found in $1!"
        $6
      fi
    fi
    if test x$3 = xyes; then
      if test x$4 = xyes; then
        AC_CHECK_LIB(fftw3f_threads, fftwf_init_threads,
		acx_fftwft_ok=yes FFTW_LIBS="-L$1 -lfftw3f_threads -lfftw3f",
		[acx_fftwft_ok=no], [-lm])
        if test x$acx_fftwft_ok = xyes; then
          AC_DEFINE(HAVE_FFTWFT,1,
    [Define if you have the FFTW single precision multithreaded libraries and header files.])
          FFTW_LIBPATH=""
          FFTW_LIBS="-L$1 -lfftw3f_threads -lfftw3f"
        else
          ATLAS_ERROR="FFTW single precision library in $1 was compiled without multithreading support!"
          $6
        fi
      else
        AC_CHECK_LIB(fftw3_threads, fftw_init_threads,
		acx_fftwdt_ok=yes FFTW_LIBS="-L$1 -lfftw3_threads -lfftw3",
		[acx_fftwdt_ok=no], [-lm])
        if test x$acx_fftwdt_ok = xyes; then
          AC_DEFINE(HAVE_FFTWT,1,
    [Define if you have the FFTW double precision multithreaded libraries and header files.])
          FFTW_LIBPATH=""
          FFTW_LIBS="-L$1 -lfftw3_threads -lfftw3"
        else
          ATLAS_ERROR="FFTW double precision library in $1 was compiled without multithreading support!"
          $6
        fi
      fi
    fi
  fi
  LIBS="$OLIBS"
  AC_SUBST(FFTW_LIBPATH)
  AC_SUBST(FFTW_LIBS)
  $5
else
  $6
fi


])dnl ACX_FFTW
