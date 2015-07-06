# CMakeLists.txt for CAFFE/src/Parallel
# Author: Adam O'Brien

set(HEADERS Parallel.h)

set(SOURCES Parallel.cc)

add_library(Parallel ${HEADERS} ${SOURCES})
