# CMakeLists.txt for CAFFE/src/Schemes
# Author: Adam O'Brien

set(HEADERS FvScheme.h
            FvScalarScheme.h
            FvVectorScheme.h)

set(SOURCE FvSchemeI.h
           FvScalarScheme.cc
           FvVectorScheme.cc)

add_library(Schemes ${HEADERS} ${SOURCE})
target_link_libraries(Schemes HexaFvmMesh)
