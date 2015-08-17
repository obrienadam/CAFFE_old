# CMakeLists.txt for CAFFE/test
# Author: Adam O'Brien

# Build Vector3D test unit

add_executable(Vector3DTest Vector3DTest.cpp)
target_link_libraries(Vector3DTest Math)

# Build Matrix test unit
add_executable(MatrixTest MatrixTest.cpp)
target_link_libraries(MatrixTest Math RunControl)

# Build geometry test unit

add_executable(geometryTest GeometryTest.cpp)
target_link_libraries(geometryTest Math)

# Build field test unit

add_executable(FieldTest FieldTest.cpp)
target_link_libraries(FieldTest Field HexaFvmMesh)

# Build sparse matrix test unit
add_executable(SparseMatrixTest SparseMatrixTest.cpp)
target_link_libraries(SparseMatrixTest Math RunControl)

# Build parallel test unit
add_executable(parallelTest ParallelTest.cpp)
target_link_libraries(parallelTest Parallel Math HexaFvmMesh)

# Build multi-block hexa mesh test unit
add_executable(MultiBlockHexaMeshTest MultiBlockHexaMeshTest.cpp)
target_link_libraries(MultiBlockHexaMeshTest Parallel Math HexaFvmMesh)
