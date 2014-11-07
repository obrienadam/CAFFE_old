# CMakeLists.txt for CAFFE/src/Math
# Author: Adam O'Brien

set(HEADERS Vector3D.h
            Tensor3D.h
	    FiniteDifference.h
	    FiniteDifferenceI.h)

set(SOURCES Vector3D.cc
	    Tensor3D.cc)

add_library(Math ${HEADERS} ${SOURCES})
