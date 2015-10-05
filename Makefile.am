dist_doc_DATA = README.md
ACLOCAL_AMFLAGS = -I m4

testsubdir=results
TESTS = src/testing/tests/05-exceptions.pl src/testing/tests/10-testPhaseAPP.pl src/testing/tests/105-linearPrior.pl src/testing/tests/11-testPhaseAPP.pl src/testing/tests/15-simplePhase.pl src/testing/tests/20-bigPhaseAPP.pl src/testing/tests/30-bigPhaseGP.pl src/testing/tests/40-bigPhaseWTCCC.pl gtest/hapSampT

TEST_EXTENSIONS = .pl
PL_LOG_COMPILER = $(PERL)

distclean-local:
	-rm -rf $(testsubdir) gtest


hfheaders=src/hapfuse.hpp src/hapSamp.hpp src/site.hpp src/writer.hpp src/hfHelper.hpp 
hfcpp=src/hapfuse.cpp  src/hapSamp.cpp src/utils.cpp src/writer.cpp src/hfHelper.cpp

lib_LIBRARIES=libhapfuse.a
libhapfuse_a_SOURCES=$(hfheaders) $(hfcpp)

bin_PROGRAMS = hapfuse
hapfuse_SOURCES = $(hfheaders) $(hfcpp) src/main.cpp

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(ZLIB_CFLAGS) $(HTSLIB_CFLAGS)
AM_LDFLAGS = $(BOOST_IOSTREAMS_LDFLAGS)

LIBS += $(BOOST_IOSTREAMS_LIBS) $(ZLIB_LIBS) $(HTSLIB_LIBS)


# gtest
checkflags=-I../gtest/include -I../src -pthread
gtestdir=gtest
check_PROGRAMS = gtest/hapSampT
gtest_hapSampT_LDADD=../gtest/make/gtest_main.a libhapfuse.a
gtest_hapSampT_SOURCES = src/testing/gtest/hapSampT.cpp
gtest_hapSampT_CXXFLAGS=$(checkflags)

