#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE
# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lconfuse':
AC_CHECK_LIB([confuse], [cfg_parse])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [pow])

AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
#AC_CHECK_LIB([xml2],[xmlNodeListGetString])
PKG_CHECK_MODULES([libxml2], [libxml-2.0])
AC_CHECK_LIB([z],[deflate])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

dnl Check if we have enable debug support.
AC_MSG_CHECKING(whether to enable debugging)
debug_default="yes"
AC_ARG_ENABLE(debug, [ --enable-debug=[no/yes] turn on debugging
[default=$debug_default]],, enable_debug=$debug_default)
dnl Yes, shell scripts can be used
if test "x$enable_debug" = "xyes"; then
CFLAGS="-g -O2 -DDEBUG"
AC_MSG_RESULT(yes)
else
CFLAGS="-O3 -ffast-math"
AC_MSG_RESULT(no)
fi

#detect commit version and date
#-                   ["Trisurf $VERSION on $host, compiled by $cc_string, `expr $ac_cv_sizeof_void_p \* 8`-bit"],
#                   ["Trisurf (commit `cd $srcdir && git log -1 --format=format:%h`) on $host, compiled by $cc_string, `expr $ac_cv_sizeof_void_p \* 8`-bit"],
AC_DEFINE_UNQUOTED(TS_VERSION,
                    ["Trisurf (commit `git describe  --abbrev=6 --dirty --always --tags`) on `hostname`."],
                    [A string containing the version number])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt strndup])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 test/Makefile])
AC_OUTPUT
