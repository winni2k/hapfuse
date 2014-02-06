

all : hapfuse

hapfuse : hapfuse.cpp
	g++ hapfuse.cpp -o hapfuse -std=c++11 -lz -O3

oxford : hapfuse
	cp hapfuse ~/bin/.
