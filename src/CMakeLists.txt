cmake_minimum_required (VERSION 2.6)
project (hapfuse)

# find and include boost!
#set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
#set(Boost_DEBUG ON)
FIND_PACKAGE(Boost)
#find_package(Boost 1.49 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARIES})

#message(status "** Boost Include: ${Boost_INCLUDE_DIR}")
#message(status "** Boost Libraries: ${Boost_LIBRARY_DIRS}")
#message(status "** Boost Libraries: ${Boost_LIBRARIES}")
#message(status "** Boost Libraries: ${Boost_LIBRARY_DIR}")


# The version number.
set (hapfuse_VERSION_MAJOR 0)
set (hapfuse_VERSION_MINOR 1)

SET(CMAKE_CXX_FLAGS "-std=c++11")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/hapfuseConfig.h.in"
  "${PROJECT_BINARY_DIR}/hapfuseConfig.h"
  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")
 

# add vcfParser include
include_directories ("${PROJECT_SOURCE_DIR}/vcfParser")
add_subdirectory (vcfParser)

# add utils include
include_directories ("${PROJECT_SOURCE_DIR}/utils")
add_subdirectory (utils)

# add executable
add_executable( hapfuse hapfuse.cpp)
target_link_libraries (hapfuse utils)
target_link_libraries (hapfuse ${Boost_LIBRARIES})


# add the install targets
#install (TARGETS hapfuse DESTINATION bin)
#install (FILES "${PROJECT_BINARY_DIR}/hapfuseConfig.h"        
#         DESTINATION include)

#set(HEADER_FILES utils.hpp vcf_parser.hpp)

#set(SOURCES file.cpp file2.cpp ${YOUR_DIRECTORY}/file1.h ${YOUR_DIRECTORY}/file2.h)
#add_executable(test ${SOURCES})


# testing
enable_testing()
add_test (hapfuseRuns hapfuse)