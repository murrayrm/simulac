dnl configure.ac - autoconf configuration for Simulac
dnl RMM, 16 Mar 2010

AC_INIT([simulac], [1.0b])
AM_INIT_AUTOMAKE			dnl initialize automake
AC_PROG_MAKE_SET			dnl set SET_MAKE in toplevel Makefile
AC_PROG_CC				dnl configure C compiler
AC_PROG_INSTALL				dnl configuration install
AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile examples/lambda/Makefile])
AC_CHECK_LIB([m], [cos])		dnl include math library
AC_OUTPUT
