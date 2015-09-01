# CMakeLists.txt for CAFFE/src/Parallel
# Author: Adam O'Brien

set(HEADERS Parallel.h
            Vector3DBuffer.h)

set(SOURCES Parallel.cpp
            Vector3DBuffer.cpp)

add_library(Parallel ${HEADERS} ${SOURCES})
target_link_libraries(Parallel mpi Math)
