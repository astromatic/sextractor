dnl @synopsis ACX_ATLAS([ATLAS_LIBDIR,[ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the ATLAS library and header files
dnl are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$ATLAS_LIBPATH -llapack -lcblas -latlas  $LIBS"
dnl        CFLAGS="$CFLAGS $ATLAS_CFLAGS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if BLAS/LAPACK
dnl is found (HAVE_ATLAS is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.
dnl
dnl @version $Id: acx_atlas.m4,v 1.0 2004/06/01 21:30:17 bertin Exp $
dnl @author Emmanuel Bertin <bertin@iap.fr>

AC_DEFUN([ACX_ATLAS], [
AC_REQUIRE([AC_CANONICAL_HOST])

acx_atlas_ok=no
AC_CHECK_HEADERS([cblas.h clapack.h],[acx_atlas_ok=yes])
if test x$acx_atlas_ok = xno; then
    AC_CHECK_HEADERS([atlas/cblas.h atlas/clapack.h],[acx_atlas_ok=yes])
    if test x$acx_atlas_ok = xyes; then
      AC_DEFINE(ATLAS_BLAS_H, "atlas/cblas.h",
        [BLAS header filename.])
      AC_DEFINE(ATLAS_LAPACK_H, "atlas/clapack.h",
        [LAPACK header filename.])
    else
      AC_DEFINE_UNQUOTED(ATLAS_BLAS_H, "$1/include/cblas.h",
        [BLAS header filename.])
      AC_DEFINE_UNQUOTED(ATLAS_LAPACK_H, "$1/include/clapack.h",
        [LAPACK header filename.])
      AC_CHECK_HEADERS([$1/include/cblas.h $1/include/clapack.h],[acx_atlas_ok=yes])
    fi
else
    AC_DEFINE(ATLAS_BLAS_H, "cblas.h",
        [BLAS header filename.])
    AC_DEFINE(ATLAS_LAPACK_H, "clapack.h",
        [LAPACK header filename.])
fi
if test x$acx_atlas_ok = xyes; then
    OLIBS="$LIBS"
    LIBS=""
    if test x$1 = x; then
      AC_CHECK_LIB(lapack, [clapack_dpotrf],, [acx_atlas_ok=no],
		[-lcblas -latlas -lm])
      AC_CHECK_LIB(cblas, cblas_sgemm,, [acx_atlas_ok=no],
		[-latlas -lm])
      ATLAS_LIBPATH=""
    else
      AC_CHECK_LIB(lapack, [clapack_dpotrf],, [acx_atlas_ok=no],
		[-L$1 -lcblas -latlas -lm])
      AC_CHECK_LIB(cblas, cblas_sgemm,, [acx_atlas_ok=no],
		[-L$1 -latlas -lm])
      if test x$acx_atlas_ok = xyes; then
        ATLAS_LIBPATH="-L$1"
      else
        acx_atlas_ok=yes
        unset ac_cv_lib_lapack_clapack_dpotrf
        AC_CHECK_LIB(lapack, [clapack_dpotrf],, [acx_atlas_ok=no],
		[-L$1/lib -lcblas -latlas -lm])
        unset ac_cv_lib_cblas_cblas_sgemm
        AC_CHECK_LIB(cblas, cblas_sgemm,, [acx_atlas_ok=no],
		[-L$1/lib -latlas -lm])
        ATLAS_LIBPATH="-L$1/lib"
      fi
    fi
    LIBS="$OLIBS"
fi

AC_SUBST(ATLAS_LIBPATH)
AC_SUBST(ATLAS_CFLAGS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_atlas_ok" = xyes; then
        AC_DEFINE(HAVE_ATLAS,1,
        [Define if you have the ATLAS libraries and header files.])
        $2
else
        $3
fi

])dnl ACX_ATLAS
