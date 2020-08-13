AC_INIT([jellyfishsim], [0.0.1], [gmarcais@cs.cmu.edu])
AC_CANONICAL_HOST
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects foreign parallel-tests color-tests])
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([.])
AC_CONFIG_HEADERS([config.h])

# Change default compilation flags
AC_LANG(C++)
AC_PROG_CXX

PKG_CHECK_MODULES([JELLYFISH], [jellyfish-2.0])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
