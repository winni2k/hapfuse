AC_INIT([hapfuse], [0.12], [wkretzsch@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Use -O3 as default optimization level
: ${CXXFLAGS="-O3"}

AC_CONFIG_MACRO_DIR([m4])

echo "Testing for a c++ complier"
AC_PROG_CXX

AC_LANG(C++)
AC_PROG_CXX

# Check for C++11
AX_CXX_COMPILE_STDCXX_11()

AC_PROG_LIBTOOL

#echo "Testing for boost"
#BOOST_REQUIRE([1.34.1])
#BOOST_IOSTREAMS

AC_PATH_PROG(PERL, perl)
if test -z "$PERL"; then
   AC_MSG_ERROR([perl not found])
fi
if test -n "`$PERL -v | fgrep 'version 5.001' 2> /dev/null`"; then
   AC_MSG_ERROR([perl 5.001 has bug which causes automake to fail])
fi


#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

#AC_CHECK_HEADERS(iostream)
#AC_CHECK_LIB(m, cos)

