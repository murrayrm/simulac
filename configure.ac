dnl configure.ac - autoconf configuration for Simulac
dnl RMM, 16 Mar 2010

AC_INIT([simulac], [1.0b])
AM_INIT_AUTOMAKE			dnl initialize automake
AC_PROG_MAKE_SET			dnl set SET_MAKE in toplevel Makefile
AC_PROG_CC				dnl configure C compiler
AC_PROG_INSTALL				dnl configuration install
AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile examples/lambda/Makefile])
AC_CHECK_LIB([m], [cos])		dnl include math library

dnl Check if we have enable debug support.
AC_MSG_CHECKING(whether to enable debugging)
debug_default="no"
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                       [default=$debug_default]],, enable_debug=$debug_default)
if test "x$enable_debug" = "xyes"; then
         CFLAGS="-g"
AC_MSG_RESULT(yes)
else
CFLAGS="-O2"
AC_MSG_RESULT(no)
fi

AC_OUTPUT
