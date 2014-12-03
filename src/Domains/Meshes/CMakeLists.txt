# CMakeLists.txt for CAFFE/src/Domains/Meshes
# Author: Adam O'Brien

set(HEADERS StructuredMesh.h
	    StructuredMeshI.h
	    HexaFdmMesh.h
            HexaMeshGen.h)

set(SOURCE HexaMeshGen.cc)

add_library(HexaMesh ${HEADERS} ${SOURCE})

target_link_libraries(HexaMesh RunControl
                               Math)
