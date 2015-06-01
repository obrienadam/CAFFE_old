# CMakeLists.txt for CAFFE/src/Schemes
# Author: Adam O'Brien

set(HEADERS FvScheme.h
            FvSchemeI.h
	    Diffusion.h
	    LinearAdvection.h
            Simple.h
            Piso.h
            IbSimple.h
            MultiphaseSimple.h)

set(SOURCE FvScheme.cc
	   Diffusion.cc
	   LinearAdvection.cc
           Simple.cc
           Piso.cc
           IbSimple.cc
           MultiphaseSimple.cc)

add_library(Schemes ${HEADERS} ${SOURCE})
target_link_libraries(Schemes HexaMesh)
