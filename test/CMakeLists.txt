# CMakeLists.txt for CAFFE/test
# Author: Adam O'Brien

# Build geometry test unit

add_executable(geometryTest GeometryTest.cpp)
target_link_libraries(geometryTest Math)

# Build matrix test unit
add_executable(matrixTest MatrixTest.cpp)
target_link_libraries(matrixTest Math RunControl)

# Build sparse matrix test unit
add_executable(sparseMatrixTest SparseMatrixTest.cpp)
target_link_libraries(sparseMatrixTest Math RunControl)

# Build parallel test unit
add_executable(parallelTest ParallelTest.cpp)
target_link_libraries(parallelTest Parallel Math HexaFvmMesh)

# Build multi-block hexa mesh test unit
add_executable(multiBlockHexaMeshTest multiBlockHexaMeshTest.cpp)
target_link_libraries(multiBlockHexaMeshTest Parallel Math HexaFvmMesh)
