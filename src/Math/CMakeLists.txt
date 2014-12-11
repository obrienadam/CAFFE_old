# CMakeLists.txt for CAFFE/src/Math
# Author: Adam O'Brien

set(HEADERS MathConstants.h
            Vector3D.h
            Point3D.h
            Tensor3D.h
	    Geometry.h
	    Matrix.h
            Interpolation.h
	    FiniteDifference.h
	    FiniteDifferenceI.h)

set(SOURCES Vector3D.cc
            Tensor3D.cc
	    Geometry.cc
	    Matrix.cc
            Interpolation.cc)

add_library(Math ${HEADERS} ${SOURCES})
target_link_libraries(Math lapacke)
