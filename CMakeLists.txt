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

# Libraries

add_subdirectory(src)

# Modules

add_subdirectory(Modules)
