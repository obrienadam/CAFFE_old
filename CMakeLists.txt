# CMakeLists.txt for CAFFE
# Author: Adam O'Brien

cmake_minimum_required (VERSION 2.8 FATAL_ERROR)

project(CAFFE)

# Find boost package

find_package(Boost REQUIRED COMPONENTS program_options system date_time)

if(Boost_FOUND)

    include_directories(${Boost_INCLUDE_DIRS})

endif()

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

# Compiler flags

set(CMAKE_CC_FLAGS "-g -Wall -O3")
set(CMAKE_CXX_FLAGS "-g -Wall -Wno-reorder -Wno-sign-compare -std=c++11 -O3")

# Output paths

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Libraries

add_subdirectory(src)

# Modules

add_subdirectory(Modules)

# Testing

add_subdirectory(test)
