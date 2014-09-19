export CXX= g++-4.7+

debug: DEBUGFLAGS= -g
CXXFLAGS +=-Wall -Wextra -std=c++11 -std=gnu++11 $(DEBUGFLAGS)

debug: CXXFLAGS += -Og
all: CXXFLAGS += -O3

VERSION := $(shell grep VERSION_MAJOR version.hpp | cut -f3 -d' ').$(shell grep VERSION_MINOR version.hpp | cut -f3 -d' ')

#OBJ=objects
BIN=../bin
#SRC=src

BOOSTLIB = -lboost_iostreams #alternatively: /opt/boost/boost-1.54.0-intel/lib/libboost_iostreams.a
export LIBS += -lgsl -lgslcblas -lbz2 $(BOOSTLIB) -lm -lstdc++ -lz -lhts

# path variables
#export CPATH += 
export CPLUS_INCLUDE_PATH += :utils
#export LIBRARY_PATH += :u-wkretzsch-snptools/samtools-0.1.16:u-wkretzsch-snptools/tabix-0.2.5
#export LD_LIBRARY_PATH += $(LIBRARY_PATH)

EFILE := $(BIN)/hapfuse.$(VERSION)
DEFILE := $(EFILE).debug

all: $(EFILE)

oxford: $(EFILE)
	cd $(BIN) && cp hapfuse $@ ~/bin/

debug: $(DEFILE)

test: $(BIN)/hapfuse.$(VERSION)
	make -C ../testing

clean:
	rm -f $(OBJECT_FILES) $(HEADER_FILES)

.PHONY: clean

.DELETE_ON_ERROR:

### normal compilation

## precompile headers
%.hpp.gch: %.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $< -o $@

utils/utils.hpp.gch: utils/utils.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $< -o $@

utils/utils.o: utils/utils.cpp utils/utils.hpp.gch
	$(CXX) $(CXXFLAGS) -Wno-sign-compare -c $< -o $@ $(LIBS)


OBJECT_FILES = utils/utils.o
HEADER_FILES = utils/utils.hpp.gch


$(EFILE): hapfuse.cpp $(OBJECT_FILES) $(HEADER_FILES)
	mkdir -p $(BIN)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJECT_FILES) $(LIBS)
	cd $(BIN) && rm -f hapfuse && ln -s hapfuse.$(VERSION) hapfuse

$(DEFILE): hapfuse.cpp $(OBJECT_FILES) $(HEADER_FILES)
	mkdir -p $(BIN)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJECT_FILES) $(LIBS)
