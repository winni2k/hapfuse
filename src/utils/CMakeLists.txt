<<<<<<< HEAD
add_library(utils utils.cpp)
target_link_libraries(utils boost_iostreams)
=======
add_library(utils STATIC utils.cpp)
target_link_libraries(utils ${ZLIB_LIBRARIES} ${IOSTREAMS} -lgsl -lgslcblas -lbz2 -lm -lstdc++ -lz)
>>>>>>> refs/remotes/origin/synced/master
