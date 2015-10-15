#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([rational], [1.5], [heiko@rangun.de])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign dist-xz no-dist-gzip tar-pax no-installinfo no-texinfo.tex])
AM_MAINTAINER_MODE
LT_INIT

# Checks for programs.
AC_PROG_CXX
AM_PATH_CPPUNIT([1.12])
PKG_PROG_PKG_CONFIG([0.22])

# Checks for libraries.
AM_CONDITIONAL([WITH_CLN], [false])
AC_ARG_ENABLE([cln], [AS_HELP_STRING([--enable-cln], [enable CLN extension])], [
if test x$enable_cln = xyes; then
  PKG_CHECK_MODULES([CLN], [cln >= 1.3.2], [
  AM_CONDITIONAL([WITH_CLN], [true])
  AC_DEFINE([HAVE_CLN_RATIONAL], [1], [define to 1 if CLN is available])
fi
])])

# Checks for header files.
AC_LANG_PUSH([C++])
AC_HEADER_STDBOOL
AM_CONDITIONAL([WITH_GMPXX], [false])
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
  if test -z "${GMP_LIBS}" ; then
    GMP_LIBS="-lgmpxx -lgmp"
  fi
fi
])
AC_LANG_POP([C++])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
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
AC_CHECK_FUNCS([ffs])
AC_LANG_POP([C++])

AC_ARG_VAR([GMP_CFLAGS], [C compiler flags for GMP, useful if installed in a non-standard location])
AC_ARG_VAR([GMP_LIBS], [linker flags for GMP, useful if installed in a non-standard location])

AC_OUTPUT([
	Makefile
	src/Makefile
	src/librational/Makefile
	src/test/Makefile
])
