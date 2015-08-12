# CMakeLists.txt for CAFFE/src/Parallel
# Author: Adam O'Brien

set(HEADERS Parallel.h)

set(SOURCES Parallel.cpp)

add_library(Parallel ${HEADERS} ${SOURCES})
target_link_libraries(Parallel mpi)
