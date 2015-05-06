# CMakeLists.txt for CAFFE/src/Schemes
# Author: Adam O'Brien

set(HEADERS FvScheme.h
	    Diffusion.h
	    LinearAdvection.h
	    Simple.h)

set(SOURCE FvScheme.cc
		FvSchemeI.h
	   Diffusion.cc
	   LinearAdvection.cc
	   Simple.cc)

add_library(Schemes ${HEADERS} ${SOURCE})
