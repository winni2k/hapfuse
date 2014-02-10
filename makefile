export CXX= g++-4.7+
export DEBUGFLAGS= -ggdb
export CXXFLAGS= -Wall -Wextra -std=c++11 -std=gnu++11


all : hapfuse

src/leo/vcf_parser.hpp.gch: src/leo/vcf_parser.hpp
	$(CXX)

hapfuse=g++ hapfuse.cpp -o hapfuse -std=c++11 -lz -O3

debug : hapfuse.cpp
	$(hapfuse) -ggdb

hapfuse : hapfuse.cpp
	$(hapfuse)

oxford : hapfuse
	cp hapfuse ~/bin/.
