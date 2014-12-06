# CMakeLists.txt for CAFFE/src/Domains/Meshes
# Author: Adam O'Brien

set(HEADERS DomainIncludes.h
	    ../DomainInterface.h
            StructuredMesh.h
	    HexaFvmMesh.h
            HexaMeshGen.h
	    Field.h)

set(SOURCE HexaMeshGen.cc
	   StructuredMesh.cc
           HexaFvmMesh.cc)

add_library(HexaMesh ${HEADERS} ${SOURCE})

target_link_libraries(HexaMesh RunControl
                               Math)
