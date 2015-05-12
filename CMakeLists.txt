# CMakeLists.txt for CAFFE
# Author: Adam O'Brien

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(CAFFE)

# Version number

set(CAFFE_VERSION_MAJOR 1)
set(CAFFE_VERSION_MINOR 0)

# Compiler

set(CMAKE_CXX_COMPILER "g++")

# Compiler flags

set(CMAKE_CXX_FLAGS "-g -Wall -Wno-reorder -Wno-sign-compare -Wno-switch -std=c++11 -O3")

# Output paths

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# Find boost package

find_package(Boost REQUIRED COMPONENTS program_options system date_time filesystem)

if(Boost_FOUND)

    include_directories(${Boost_INCLUDE_DIRS})

endif()

# Find lapack package

find_package(LAPACK)

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

# Testing

add_subdirectory(test)
