AC_INIT([hapfuse], [1.3], [wkretzsch@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Use -O3 as default optimization level
: ${CXXFLAGS="-O3"}

AC_CONFIG_MACRO_DIR([m4])

echo "Testing for a c++ complier"
AC_PROG_CXX

echo "Setting language to C++"
AC_LANG(C++)
AC_PROG_CXX

# Check for C++11
echo "Check for C++11 support in compiler"
AX_CXX_COMPILE_STDCXX_11(,mandatory)

AC_PROG_LIBTOOL

# Checks for libraries
echo "Check for other libraries"
AC_SEARCH_LIBS([hts_version],[hts],[],[
        AC_MSG_ERROR([unable to find the hts_version() function from libhts])
        ])

echo "Check for Boost"
BOOST_REQUIRE([1.34.1])
BOOST_IOSTREAMS

echo "Check for Perl"
AC_PATH_PROG(PERL, perl)
if test -z "$PERL"; then
   AC_MSG_ERROR([perl not found])
fi
if test -n "`$PERL -v | fgrep 'version 5.001' 2> /dev/null`"; then
   AC_MSG_ERROR([perl 5.001 has bug which causes automake to fail])
fi

AC_CONFIG_FILES([Makefile])
AC_OUTPUT


