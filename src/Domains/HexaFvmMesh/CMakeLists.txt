# CMakeLists.txt for CAFFE/src/Domains/HexaFvmMesh
# Author: Adam O'Brien

set(HEADERS StructuredMesh.h
	    HexaFvmMesh.h
            ParallelHexaFvmMesh.h
            HexaMeshGen.h
            IndexMap.h)

set(SOURCE StructuredMesh.cpp
           HexaFvmMesh.cpp
           ParallelHexaFvmMesh.cpp
           HexaMeshGen.cpp
           IndexMap.cpp)

add_library(HexaFvmMesh ${HEADERS} ${SOURCE})

target_link_libraries(HexaFvmMesh RunControl
                                  Math)
