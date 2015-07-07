# CMakeLists.txt for CAFFE/test
# Author: Adam O'Brien

# Build geometry test unit

add_executable(geometryTest GeometryTest.cc)
target_link_libraries(geometryTest Math)

# Build matrix test unit

add_executable(matrixTest MatrixTest.cc)
target_link_libraries(matrixTest Math RunControl)

# Build sparse matrix test unit

add_executable(sparseMatrixTest SparseMatrixTest.cc)
target_link_libraries(sparseMatrixTest Math RunControl)

# Build parallel test unit

add_executable(parallelTest ParallelTest.cc)
target_link_libraries(parallelTest Parallel Math)
