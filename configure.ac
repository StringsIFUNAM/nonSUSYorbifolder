AC_INIT([Orbifolder], [1.0], [orbifolder@projects.hepforge.org ], [orbifolder], [http://projects.hepforge.org/orbifolder/])

m4_include([m4/m4_ax_boost_base.m4])

AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall])
LT_INIT

AC_CONFIG_SRCDIR([src/io/io.cc])

AC_SUBST([foo1dir], ['$(PWD)'])
AC_SUBST([foo2dir], ['$(PWD)/../..'])

AC_PROG_CXX
AC_PROG_LIBTOOL

AX_PATH_GSL(1.9, [], [AC_MSG_ERROR(GSL - version >= 1.9 not found!)])

AX_BOOST_BASE([1.0], [], [AC_MSG_ERROR(BOOST - version >= 1.0 not found!)])

AC_CONFIG_HEADERS([Config/config.h])
AC_CONFIG_FILES([Makefile src/io/Makefile src/linalg/Makefile src/groupTheory/Makefile src/orbifolder/Makefile])

AC_OUTPUT
