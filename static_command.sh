

/bin/bash ./libtool  --tag=CXX   --mode=link g++  -O3 -std=gnu++11  -o hapfuse  src/hapfuse.o src/hapSamp.o src/utils.o src/writer.o src/hfHelper.o src/main.o /usr/lib/x86_64-linux-gnu/libhts.a /usr/lib/x86_64-linux-gnu/libz.a  /usr/lib/x86_64-linux-gnu/libboost_iostreams.a /usr/lib/x86_64-linux-gnu/libbz2.a -lpthread
