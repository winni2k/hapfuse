dist_doc_DATA = README.md
ACLOCAL_AMFLAGS = -I m4

testsubdir=results
TESTS = src/testing/tests/05-exceptions.pl src/testing/tests/10-testPhaseAPP.pl src/testing/tests/105-linearPrior.pl src/testing/tests/11-testPhaseAPP.pl src/testing/tests/15-simplePhase.pl src/testing/tests/20-bigPhaseAPP.pl src/testing/tests/30-bigPhaseGP.pl src/testing/tests/40-bigPhaseWTCCC.pl

TEST_EXTENSIONS = .pl
PL_LOG_COMPILER = $(PERL)

distclean-local:
	-rm -rf $(testsubdir)

bin_PROGRAMS = hapfuse
hapfuse_SOURCES = src/hapfuse.cpp src/hapfuse.hpp src/hapSamp.cpp src/hapSamp.hpp src/main.cpp src/site.hpp src/utils.cpp src/writer.cpp src/writer.hpp src/hfHelper.hpp src/hfHelper.cpp

AM_CPPFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_IOSTREAMS_LDFLAGS)
LIBS += $(BOOST_IOSTREAMS_LIBS)
