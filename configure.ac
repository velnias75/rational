#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([rational], [0.0], [heiko@rangun.de])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([foreign dist-xz no-dist-gzip tar-pax no-installinfo no-texinfo.tex])
AM_MAINTAINER_MODE
LT_INIT

# Checks for programs.
AC_PROG_CXX
AM_PATH_CPPUNIT(1.13)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.
AC_LANG_PUSH([C++])
AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([sqrt])
AC_LANG_POP([C++])

AC_OUTPUT([
	Makefile
	src/Makefile
	src/librational/Makefile
	src/test/Makefile
])
