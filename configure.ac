#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([rational], [2.13.4], [heiko@rangun.de])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign dist-xz no-dist-gzip tar-pax no-installinfo no-texinfo.tex])
AM_MAINTAINER_MODE
LT_INIT

CFLAGS="${CXXFLAGS} -std=c99"

# Checks for programs.
AC_PROG_CXX
PKG_PROG_PKG_CONFIG([0.22])
AX_CXX_COMPILE_STDCXX_11([noext], [optional])
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

# Checks for libraries.
AM_CONDITIONAL([WITH_CLN], [false])
AC_ARG_ENABLE([cln], [AS_HELP_STRING([--enable-cln], [enable CLN extension])], [
if test x$enable_cln = xyes; then
  PKG_CHECK_MODULES([CLN], [cln >= 1.3.2], [
  AM_CONDITIONAL([WITH_CLN], [true])
  AC_DEFINE([HAVE_CLN_RATIONAL], [1], [define to 1 if CLN is available])
fi
])])

PKG_CHECK_MODULES([CPPUNIT], [cppunit >= 1.14])

# Checks for header files.
AC_LANG_PUSH([C++])
AC_HEADER_STDBOOL
AM_CONDITIONAL([WITH_GMPXX], [false])
AM_CONDITIONAL([HAVE_MPREAL_H], [false])
AC_ARG_ENABLE([gmp], [AS_HELP_STRING([--enable-gmp], [enable GMP extension])], [
if test x$enable_gmp = xyes; then
  if test -n "${GMP_CFLAGS}" ; then
    ac_save_cppflags="${CPPFLAGS}"
    CPPFLAGS="${CPPFLAGS} ${GMP_CFLAGS}"
  fi
  AC_CHECK_HEADER([gmpxx.h], [
    AM_CONDITIONAL([WITH_GMPXX], [true])
    AC_DEFINE([HAVE_GMPXX_H], [1], [define to 1 if gmpxx.h is available])])
  if test -n "${GMP_CFLAGS}" ; then
    CPPFLAGS=$ac_save_cppflags
  fi
  AC_ARG_WITH([mpreal], [AS_HELP_STRING([--with-mpreal],
    [if gmp is enabled, use mpreal.h and mpfr for floating points])], [
       if test x$with_mpreal = xyes; then
       ac_save_cppflags="${CPPFLAGS}"
       CPPFLAGS="${CPPFLAGS} -std=c++0x"
       AC_CHECK_HEADER([mpreal.h],
         [ AC_LANG_PUSH([C])
           AC_CHECK_LIB([mpfr], [mpfr_set_str], [
           AM_CONDITIONAL([HAVE_MPREAL_H], [true])
	       GMP_CFLAGS="${GMP_CFLAGS} -DHAVE_MPREAL_H=1"
           GMP_LIBS="${GMP_LIBS} -lmpfr"
           AC_MSG_NOTICE([found mpreal.h, will use mpfr and mpfrc++ for floating points])
        ])
       AC_LANG_POP([C])
      ])
    CPPFLAGS=$ac_save_cppflags
  fi])
  if test -z "${GMP_LIBS}" -o "${GMP_LIBS}" = " -lmpfr"; then
    GMP_LIBS="${GMP_LIBS}"
    AC_LANG_PUSH([C])
    AC_CHECK_LIB([gmpxx], main, [GMP_LIBS="${GMP_LIBS} -lgmpxx"],
      [AC_MSG_ERROR([GNU MP++ not found, see https://gmplib.org/])])
    AC_CHECK_LIB([gmp], [__gmpz_gcd], [GMP_LIBS="${GMP_LIBS} -lgmp"],
      [AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])
    AC_LANG_POP([C])
  fi
fi
])
AM_CONDITIONAL([WITH_INFINT], [false])
AC_ARG_ENABLE([infint], [AS_HELP_STRING([--enable-infint],
  [enable InfInt extension (deprecated)])], [
if test x$enable_infint = xyes; then
  if test -n "${INFINT_CFLAGS}" ; then
    ac_save_cppflags="${CPPFLAGS}"
    CPPFLAGS="${CPPFLAGS} ${INFINT_CFLAGS}"
  fi
  AC_CHECK_HEADER([InfInt.h], [
    AM_CONDITIONAL([WITH_INFINT], [true])
    AC_DEFINE([HAVE_INFINT_H], [1], [define to 1 if InfInt.h is available])])
  if test -n "${INFINT_CFLAGS}" ; then
    CPPFLAGS=$ac_save_cppflags
  fi
fi
])
AM_CONDITIONAL([WITH_UNLOGICAL], [false])
AC_ARG_ENABLE([unlogical], [AS_HELP_STRING([--enable-unlogical],
  [enable unlogical extension (experimental)])], [
if test x$enable_unlogical = xyes; then
  if test -n "${UNLOGICAL_CFLAGS}" ; then
    ac_save_cppflags="${CPPFLAGS}"
    CPPFLAGS="${CPPFLAGS} ${UNLOGICAL_CFLAGS}"
  fi
  AC_CHECK_HEADER([number.h], [
    AM_CONDITIONAL([WITH_UNLOGICAL], [true])
    AC_DEFINE([HAVE_NUMBER_H], [1], [define to 1 if number.h is available])])
  if test -n "${UNLOGICAL_CFLAGS}" ; then
    CPPFLAGS=$ac_save_cppflags
  fi
fi
])
AC_LANG_POP([C++])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_LANG_PUSH([C++])
AC_CHECK_FUNCS([modf])
AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([sqrt])
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([strndup])
AC_LANG_POP([C++])

AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([docs/Doxyfile])])

AC_ARG_VAR([GMP_CFLAGS], [C compiler flags for GMP, useful if installed in a non-standard location])
AC_ARG_VAR([GMP_LIBS], [linker flags for GMP, useful if installed in a non-standard location])
AC_ARG_VAR([INFINT_CFLAGS], [C compiler flags for InfInt, useful if installed in a non-standard location])
AC_ARG_VAR([UNLOGICAL_CFLAGS], [C compiler flags for unlogical, useful if installed in a non-standard location])
AC_ARG_VAR([UNLOGICAL_LIBS], [linker flags for unlogical, useful if installed in a non-standard location])

AC_OUTPUT([
	Makefile
	docs/Makefile
	src/Makefile
	src/rational/Makefile
	src/test/Makefile
])
