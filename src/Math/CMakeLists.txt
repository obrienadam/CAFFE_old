# CMakeLists.txt for CAFFE/src/Math
# Author: Adam O'Brien

set(HEADERS Vector3D.h
            Point3D.h
            Tensor3D.h
            Geometry.h
            Sphere.h
            Matrix.h
            SparseMatrix.h
            SparseVector.h
            Interpolation.h
            Kernel.h)

set(SOURCES Vector3D.cpp
            Tensor3D.cpp
            Geometry.cpp
            Sphere.cpp
            Matrix.cpp
            SparseMatrix.cpp
            SparseVector.cpp
            Interpolation.cpp
            Kernel.cpp)

add_library(Math ${HEADERS} ${SOURCES})
target_link_libraries(Math cblas lapack atlas petsc)
