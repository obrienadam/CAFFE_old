# CMakeLists.txt for CAFFE/src/Domains/HexaFvmMesh
# Author: Adam O'Brien

set(HEADERS StructuredMesh.h
	    HexaFvmMesh.h
            MultiBlockHexaFvmMesh.h
            HexaMeshGen.h
            MultiBlockHexaMeshGen.h
            IndexMap.h)

set(SOURCE StructuredMesh.cc
           HexaFvmMesh.cc
           MultiBlockHexaFvmMesh.cc
           HexaMeshGen.cc
           MultiBlockHexaMeshGen.cc
           IndexMap.cc)

add_library(HexaFvmMesh ${HEADERS} ${SOURCE})

target_link_libraries(HexaFvmMesh RunControl
                                  Math)
