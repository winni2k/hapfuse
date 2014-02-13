export CXX= g++-4.7+

debug: DEBUGFLAGS= -ggdb
debug: CXXFLAGS += -O0

#OBJ=objects
BIN=../bin
#SRC=src

export CXXFLAGS= -Wall -Wextra -std=c++11 -std=gnu++11 $(DEBUGFLAGS)
BOOSTLIB = -lboost_iostreams-mt #alternatively: /opt/boost/boost-1.54.0-intel/lib/libboost_iostreams.a
export LIBS += -lgsl -lgslcblas -lbz2 $(BOOSTLIB) -lm -lstdc++ -lz

# path variables
#export CPATH += 
export CPLUS_INCLUDE_PATH += :src/olivier
#export LIBRARY_PATH += :u-wkretzsch-snptools/samtools-0.1.16:u-wkretzsch-snptools/tabix-0.2.5
#export LD_LIBRARY_PATH += $(LIBRARY_PATH)

all: $(BIN)/hapfuse

$(BIN)/hapfuse: CXXFLAGS += -O3

oxford: $(BIN)/hapfuse
	cp $< ~/bin

debug:
	$(CXX) -E -x c++ - -v < /dev/null

test: $(BIN)/hapfuse
	make -C ../testing
#test : insti.o impute.o emcchain.o olivier/utils.o haplotype.o kNN.o relationship.hpp.gch
#	$(MAKE) -C testing

#testDebug:insti.o impute.o emcchain.o olivier/utils.o haplotype.o kNN.o relationship.hpp.gch
#	$(MAKE) -C testing debug

clean:
	rm -f *.o *.h.gch

.PHONY: clean

.DELETE_ON_ERROR:

### normal compilation

## precompile headers
%.hpp.gch: %.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $< -o $@

utls/utils.hpp.gch: utils/utils.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $< -o $@

vcfParser/vcf_parser.hpp.gch: vcfParser/vcf_parser.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $< -o $@

#hapfuse.hpp.gch: hainsti.hpp impute.hpp.gch emcchain.hpp.gch olivier/utils.hpp.gch relationship.hpp.gch version.hpp.gch enums.hpp.gch snp.hpp.gch hapPanel.hpp.gch
#	$(CXX) $(CXXFLAGS) -x c++-header $<

## object files
#%.o: %.cpp %.hpp.gch
#	$(CXX) $(CXXFLAGS) -c $< $(LIBS) -o $@

utils/utils.o: utils/utils.cpp utils/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -Wno-sign-compare -c $< -o $@ $(LIBS)

#main.o: main.cpp version.hpp.gch impute.hpp.gch insti.hpp.gch
#	$(CXX) $(CXXFLAGS) -c $< $(LIBS)


OBJECT_FILES = utils/utils.o
HEADER_FILES = utils/utils.hpp.gch vcfParser/vcf_parser.hpp.gch version.hpp.gch

$(BIN)/hapfuse: hapfuse.cpp $(OBJECT_FILES) $(HEADER_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJECT_FILES) $(LIBS)

#instiO0: $(OBJECT_FILES)
#	$(CXX) $(CXXFLAGS) -o $@ $(OBJECT_FILES) $(LIBS)
