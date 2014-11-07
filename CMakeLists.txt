# CMakeLists.txt for CAFFE
# Author: Adam O'Brien

cmake_minimum_required (VERSION 2.6.2)

project(CAFFE)

# Includes

set(INCLUDE_DIR src/DataStructures
		src/Math
		src/RunControl
		src/Domains
		src/Domains/Meshes
		src/Solvers
		src/Schemes
		src/State)

include_directories(${INCLUDE_DIR})

# Output paths

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Libraries

add_subdirectory(src)

# Modules

add_subdirectory(Modules)
