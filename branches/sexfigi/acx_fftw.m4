dnl @synopsis ACX_FFTW([FFTW_DIR,[ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the FFTW3 libraries and header
dnl files are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$FFTWF_LIBS $LIBS"
dnl        CFLAGS="$CFLAGS $FFTW_CFLAGS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if FFTW
dnl is found (HAVE_FFTWx are defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.
dnl
dnl @version $Id: acx_fftw.m4,v 1.0 2004/06/02 21:30:17 bertin Exp $
dnl @author Emmanuel Bertin <bertin@iap.fr>

AC_DEFUN([ACX_FFTW], [
AC_REQUIRE([AC_CANONICAL_HOST])

acx_fftw_ok=yes
AC_CHECK_HEADER(fftw3.h,,[acx_fftw_ok=no])

if test x$acx_fftw_ok = xyes; then
    AC_CHECK_LIB(fftw3f, fftwf_execute, 
		acx_fftwf_ok=yes
		FFTWF_LIBS="-lfftw3f",
		[acx_fftwf_ok=no],
		[-lm])
    AC_CHECK_LIB(fftw3, fftw_execute,
		acx_fftwd_ok=yes
		FFTW_LIBS="-lfftw3",
		[acx_fftwd_ok=no],
		[-lm])
    AC_CHECK_LIB(fftw3f_threads, fftwf_init_threads,
		acx_fftwft_ok=yes
		FFTWFT_LIBS="-lfftw3f_threads -lfftw3f",
		[acx_fftwft_ok=no],
		[-lfftw3f -lm])
    AC_CHECK_LIB(fftw3_threads, fftw_init_threads,
		acx_fftwdt_ok=yes
		FFTWT_LIBS="-lfftw3_threads -lfftw3",
		[acx_fftwdt_ok=no],
		[-lfftw3 -lm])
fi

AC_SUBST(FFTW_LIBS)
AC_SUBST(FFTW_CFLAGS)
AC_SUBST(acx_fftwft_ok)
AC_SUBST(acx_fftwdt_ok)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_fftwf_ok" = xyes; then
  AC_DEFINE(HAVE_FFTWF,1,
    [Define if you have the FFTW single precision libraries and header files.])
  if test x"$acx_fftwft_ok" = xyes; then
    AC_DEFINE(HAVE_FFTWFT,1,
    [Define if you have the FFTW threaded single precision libraries and header files.])
  fi
  $2
elif test x"$acx_fftwd_ok" = xyes; then
  AC_DEFINE(HAVE_FFTW,1,
    [Define if you have the FFTW double precision libraries and header files.])
  if test x"$acx_fftwdt_ok" = xyes; then
    AC_DEFINE(HAVE_FFTWT,1,
    [Define if you have the FFTW threaded double precision libraries and header files.])
  fi
  $2
else
  $3
fi

])dnl ACX_FFTW
