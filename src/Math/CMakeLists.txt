# CMakeLists.txt for CAFFE/src/Math
# Author: Adam O'Brien

set(SOURCES Vector3D.cc
	    Tensor3D.cc)

add_library(Math ${SOURCES})
