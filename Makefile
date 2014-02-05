export CXX= g++-4.7+
export DEBUGFLAGS= -ggdb
export CXXFLAGS= -Wall -Wextra -std=c++11 -std=gnu++11 -mpopcnt $(DEBUGFLAGS)
BOOSTLIB = -lboost_program_options-mt #alternatively: /opt/boost/boost-1.54.0-intel/lib/libboost_iostreams.a
export LIBS += -lgsl -lgslcblas -lbz2 $(BOOSTLIB) -lm -lstdc++ -lz

# path variables
#export CPATH += :u-wkretzsch-snptools/samtools-0.1.16:u-wkretzsch-snptools/tabix-0.2.5
#export CPLUS_INCLUDE_PATH += :olivier
#export LIBRARY_PATH += :u-wkretzsch-snptools/samtools-0.1.16:u-wkretzsch-snptools/tabix-0.2.5
export LD_LIBRARY_PATH += $(LIBRARY_PATH)

all: hapFuse

hapFuse: CXXFLAGS += -O3

hapFuseO0: CXXFLAGS += -O0

debug:
	$(CXX) -E -x c++ - -v < /dev/null

clean:
	rm -f *.o *.hpp.gch hapFuse hapFuseO0

.PHONY: clean

.DELETE_ON_ERROR:

### normal compilation

## precompile headers
%.hpp.gch: %.hpp
	$(CXX) $(CXXFLAGS) -x c++-header $<

#hapFuse.hpp.gch: hapFuse.hpp region.hpp.gch
#	$(CXX) $(CXXFLAGS) -x c++-header $<

## object files
%.o: %.cpp %.hpp.gch
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)

main.o: main.cpp hapFuse.hpp.gch
	$(CXX) $(CXXFLAGS) -c $< $(LIBS)


OBJECT_FILES = main.o hapFuse.o

hapFuse: $(OBJECT_FILES)
	$(CXX) $(CXXFLAGS) $(OBJECT_FILES)  -o $@ $(LIBS)

hapFuseO0: $(OBJECT_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECT_FILES) $(LIBS)


