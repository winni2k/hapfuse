export CXX= g++-4.7+

debug: DEBUGFLAGS= -ggdb
debug: CXXFLAGS += -O0

VERSION := $(shell grep VERSION_MAJOR version.hpp | cut -f3 -d' ').$(shell grep VERSION_MINOR version.hpp | cut -f3 -d' ')

#OBJ=objects
BIN=../bin
#SRC=src

export CXXFLAGS= -Wall -Wextra -std=c++11 -std=gnu++11 $(DEBUGFLAGS)
BOOSTLIB = -lboost_iostreams-mt #alternatively: /opt/boost/boost-1.54.0-intel/lib/libboost_iostreams.a
export LIBS += -lgsl -lgslcblas -lbz2 $(BOOSTLIB) -lm -lstdc++ -lz

# path variables
#export CPATH += 
export CPLUS_INCLUDE_PATH += :utils
#export LIBRARY_PATH += :u-wkretzsch-snptools/samtools-0.1.16:u-wkretzsch-snptools/tabix-0.2.5
#export LD_LIBRARY_PATH += $(LIBRARY_PATH)

all: $(BIN)/hapfuse.$(VERSION)

$(BIN)/hapfuse: CXXFLAGS += -O3

oxford: $(BIN)/hapfuse.$(VERSION)
	cd $(BIN) && cp hapfuse.$(VERSION) ~/bin/
	cd ~/bin && rm -f hapfuse && ln -s hapfuse.$(VERSION) hapfuse

debug:
	$(CXX) -E -x c++ - -v < /dev/null

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


$(BIN)/hapfuse.$(VERSION): hapfuse.cpp $(OBJECT_FILES) $(HEADER_FILES)
	mkdir -p $(BIN)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJECT_FILES) $(LIBS)
	cd $(BIN) && rm -f hapfuse && ln -s hapfuse.$(VERSION) hapfuse
