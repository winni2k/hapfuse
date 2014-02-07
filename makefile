

all : hapfuse

hapfuse=g++ hapfuse.cpp -o hapfuse -std=c++11 -lz -O3

debug : hapfuse.cpp
	$(hapfuse) -ggdb

hapfuse : hapfuse.cpp
	$(hapfuse)

oxford : hapfuse
	cp hapfuse ~/bin/.
