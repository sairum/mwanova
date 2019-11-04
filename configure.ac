#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(mwanova, 2.1.3, amsantos@fc.up.pt)
AC_CONFIG_SRCDIR(src/main.cpp)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [pow])

# Checks for header files.
AC_STDC_HEADERS

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_HEADER_STDC
AC_FUNC_MEMCMP
AC_CHECK_FUNCS([memchr memset pow sqrt strchr strstr log])

dnl this macro is used to get the arguments supplied
dnl to the configure script (./configure --enable-cgi)
dnl Check if we have enable cgi compilation.
cgi_default="no"
AC_MSG_CHECKING(whether to enable cgi compilation)
AC_ARG_ENABLE(cgi, [  --enable-cgi=[no/yes] turn on cgi compilation [default=no]],,enable_cgi=$cgi_default)
dnl Yes, shell scripts can be used
if test "x$enable_cgi" = "xyes"; then
         CPPFLAGS="$CPPFLAGS -DCGI "
AC_MSG_RESULT(yes)
else
CPPFLAGS="$CPPFLAGS"
AC_MSG_RESULT(no)
fi
AC_SUBST(CPPFLAGS)

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
