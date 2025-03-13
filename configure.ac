AC_INIT([Orbifolder],[1.0],[enriquescalante@gmail.com],[orbifolder],[https://github.com/StringsIFUNAM/nonSUSYorbifolder])

m4_include([m4/m4_ax_boost_base.m4])

AC_CONFIG_AUX_DIR([Config])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall])
LT_INIT

AC_CONFIG_SRCDIR([src/io/io.cc])

AC_SUBST([foo1dir], ['$(PWD)'])
AC_SUBST([foo2dir], ['$(PWD)/../..'])

AC_PROG_CXX
LT_INIT

AX_PATH_GSL(1.9, [], [AC_MSG_ERROR(GSL - version >= 1.9 not found!)])
AX_BOOST_BASE([1.0], [], [AC_MSG_ERROR(BOOST - version >= 1.0 not found!)])

# Check for readline library
AC_CHECK_LIB([readline], [readline], [], [AC_MSG_ERROR([readline library not found!])])
AC_CHECK_HEADER([readline/readline.h], [], [AC_MSG_ERROR([readline header not found!])])

AC_CONFIG_HEADERS([Config/config.h])
AC_CONFIG_FILES([Makefile src/io/Makefile src/linalg/Makefile src/groupTheory/Makefile src/orbifolder/Makefile])

# Add the compiler version
AC_SUBST(CXXFLAGS, "-std=c++11 -g -O2")



AC_OUTPUT
