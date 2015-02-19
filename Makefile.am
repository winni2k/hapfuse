dist_doc_DATA = README.md
ACLOCAL_AMFLAGS = -I m4

testsubdir=results
TESTS = src/testing/tests/10-testPhaseAPP.pl src/testing/tests/105-linearPrior.pl src/testing/tests/11-testPhaseAPP.pl src/testing/tests/15-simplePhase.pl src/testing/tests/20-bigPhaseAPP.pl src/testing/tests/30-bigPhaseGP.pl
#src/testing/runtests.pl
TEST_EXTENSIONS = .pl
PL_LOG_COMPILER = $(PERL)

distclean-local:
	-rm -rf $(testsubdir)

bin_PROGRAMS = hapfuse
hapfuse_SOURCES = src/hapfuse.cpp src/hapfuse.hpp src/hapSamp.cpp src/hapSamp.hpp src/main.cpp src/site.hpp src/utils.cpp

hapfuse_LDFLAGS = -lboost_iostreams -lgsl -lgslcblas -lbz2 -lm -lstdc++ -lz -lhts

AM_CXXFLAGS = -std=c++11
