# CMakeLists.txt for CAFFE/src/Schemes
# Author: Adam O'Brien

set(HEADERS FvScheme.h
	    Diffusion.h
	    LinearAdvection.h)

set(SOURCE FvScheme.cc
	   Diffusion.cc
	   LinearAdvection.cc)

add_library(Schemes ${HEADERS} ${SOURCE})
