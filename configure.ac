AC_INIT([hapfuse], [1.6.2], [wkretzsch@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

# Use -O3 as default optimization level
: ${CXXFLAGS="-O3"}

AC_CONFIG_MACRO_DIR([m4])

AM_PROG_AR

AC_MSG_CHECKING([for C++ compiler])
AC_PROG_CXX
AC_LANG(C++)

# Check for C++11
AC_MSG_CHECKING([for C++11 support in compiler])
AX_CXX_COMPILE_STDCXX_11(,mandatory)

AC_PROG_LIBTOOL

# Check for Boost >= 1.34.1 and boost_iostreams
AC_MSG_CHECKING([for Boost])
BOOST_REQUIRE([1.34.1])
BOOST_IOSTREAMS

# Checks for zlib and adds -lz to LIBS and defined HAVE_LIBZ
AC_ARG_VAR([ZLIB_CFLAGS],[C compiler flags for ZLIB])
AC_ARG_VAR([ZLIB_LIBS],[linker flags for ZLIB])
AC_MSG_CHECKING([for zlib])
AC_CHECK_LIB([z], [zlibVersion], [], [AC_MSG_FAILURE([zlib is required but check for zlibVersion function failed! (is ZLIB_LIBS set correctly?)])], [${ZLIB_LIBS}])

# Check for htslib (which requires zlib)
AC_ARG_VAR([HTSLIB_CFLAGS],[C compiler flags for HTSLIB])
AC_ARG_VAR([HTSLIB_LIBS],[linker flags for HTSLIB])
AC_MSG_CHECKING([for htslib])
AC_CHECK_LIB([hts], [hts_version], [], [AC_MSG_FAILURE([htslib is required but check for hts_version function failed! (is HTSLIB_LIBS set correctly?)])], [${LIBS} ${HTSLIB_LIBS} ${ZLIB_LIBS} ${LTLIBMULTITHREAD}])

AC_MSG_CHECKING([for Perl])
AC_PATH_PROG(PERL, perl)
if test -z "$PERL"; then
   AC_MSG_ERROR([perl not found])
fi
if test -n "`$PERL -v | fgrep 'version 5.001' 2> /dev/null`"; then
   AC_MSG_ERROR([perl 5.001 has bug which causes automake to fail])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT


