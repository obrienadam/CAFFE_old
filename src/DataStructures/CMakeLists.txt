# CMakeLists.txt for CAFFE/src/DataStructures
# Author: Adam O'Brien

set(HEADERS DataStructures.h
            Array3D.h)

set(SOURCE  DataStructures.cpp
            Array3D.tpp)

add_library(DataStructures ${HEADERS} ${SOURCE})
target_link_libraries(DataStructures Math)
