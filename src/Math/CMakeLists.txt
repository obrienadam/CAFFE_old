# CMakeLists.txt for CAFFE/src/Math
# Author: Adam O'Brien

set(HEADERS MathConstants.h
            Vector3D.h
            Point3D.h
            Tensor3D.h
	    Geometry.h
            Interpolation.h
	    FiniteDifference.h
	    FiniteDifferenceI.h)

set(SOURCES Vector3D.cc
            Tensor3D.cc
	    Geometry.cc
            Interpolation.cc)

add_library(Math ${HEADERS} ${SOURCES})
