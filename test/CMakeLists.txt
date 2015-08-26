# CMakeLists.txt for CAFFE/test
# Author: Adam O'Brien

# Build Vector3D test unit

add_executable(Vector3DTest Vector3DTest.cpp)
target_link_libraries(Vector3DTest Math)

# Build Matrix test unit
add_executable(MatrixTest MatrixTest.cpp)
target_link_libraries(MatrixTest Math RunControl)

# Build Array3D test unit
add_executable(Array3DTest Array3DTest.cpp)

# Build field test unit
add_executable(FieldTest FieldTest.cpp)
target_link_libraries(FieldTest Field HexaFvmMesh)

# Build index map test unit
add_executable(IndexMapTest IndexMapTest.cpp)
target_link_libraries(IndexMapTest Parallel HexaFvmMesh)

# Build ParallelHexaFvmMesh test unit
add_executable(ParallelHexaFvmMeshTest ParallelHexaFvmMeshTest.cpp)
target_link_libraries(ParallelHexaFvmMeshTest Parallel HexaFvmMesh)

# Build geometry test unit
add_executable(geometryTest GeometryTest.cpp)
target_link_libraries(geometryTest Math)

# Build sparse matrix test unit
add_executable(SparseMatrixTest SparseMatrixTest.cpp)
target_link_libraries(SparseMatrixTest Math RunControl)

# Build parallel test unit
add_executable(ParallelTest ParallelTest.cpp)
target_link_libraries(ParallelTest Parallel RunControl)

# Build multi-block hexa mesh test unit
add_executable(MultiBlockHexaMeshTest MultiBlockHexaMeshTest.cpp)
target_link_libraries(MultiBlockHexaMeshTest Parallel Math HexaFvmMesh)
