# CMakeLists.txt for CAFFE/src/Domains/Meshes
# Author: Adam O'Brien

set(HEADERS DomainIncludes.h
	    ../DomainInterface.h
            StructuredMesh.h
	    HexaFvmMesh.h
            MultiBlockHexaFvmMesh.h
            HexaMeshGen.h
            Field.h
            InitialConditions.h
            IndexMap.h)

set(SOURCE HexaMeshGen.cc
	   StructuredMesh.cc
           HexaFvmMesh.cc
           MultiBlockHexaFvmMesh.cc
           FieldI.h
           InitialConditions.cc
           InitialConditionsI.h
           IndexMap.cc)

add_library(HexaMesh ${HEADERS} ${SOURCE})

target_link_libraries(HexaMesh RunControl
                               Math)
